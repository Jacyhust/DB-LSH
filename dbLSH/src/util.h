#pragma once

#include <memory>
#include <iostream>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>

#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <cstdint> 
#include <vector>
#include <random>
#include <sstream>


//#include "def.h"
//#include "pri_queue.h"

#if defined(unix) || defined(__unix__)
//under POSIX system, we use clock_gettime()
//remember we have to use linker option "-lrt"
#include <sys/time.h>
#include <unistd.h>
#else
#include <time.h>
#include <windows.h>
#include <io.h>
#include <direct.h>

//using access = _access;
#define F_OK 0 //值为0，判断文件是否存在
#define X_OK 1 //值为1，判断对文件是可执行权限
#define W_OK 2 //值为2，判断对文件是否有写权限
#define R_OK 4 //值为4，判断对文件是否有读权限
//int access(const char* pathname, int mode);
#endif

#define Scalar float

inline int get_num_bits8(uint8_t x)				//get the number of 1 in the binary representation of u
{
	x = (x&0x55) + ((x>>1)&0x55);
	x = (x&0x33) + ((x>>2)&0x33);
	x = (x&0x0f) + ((x>>4)&0x0f);
	return x;
}

inline int get_num_bits64(uint64_t x)			////get the number of 1 in the binary representation of x
{
	x = x - ((x >> 1) & 0x5555555555555555);
    x = (x & 0x3333333333333333) +
        ((x >> 2) & 0x3333333333333333);
    x = ((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F);
    return (x*(0x0101010101010101))>>56;
}

int get_num_prefix(uint16_t u);
int get_num_prefix(uint32_t u);
int get_num_prefix(uint64_t u);

inline uint64_t hash_combine(uint64_t h0, uint64_t h1)
{
    return h0 ^ (0x9e3779b9 + (h0<<6) + (h0>>2) + h1);
}

#if defined(unix) || defined(__unix__)
//under POSIX system, we use clock_gettime()
//remember we have to use linker option "-lrt"
inline int log2i(unsigned x) {
    return sizeof(int) * 8 - __builtin_clz(x) - 1;
}
inline int log2ll(unsigned long long x) {
    return sizeof(int) * 8 - __builtin_clzll(x) - 1;
}
#else
#endif


// -----------------------------------------------------------------------------
template<class ScalarType>
ScalarType sqr(ScalarType x)
{
	return x*x;
}

//FProd :: ScalarType -> ScalarType -> ScalarType
template<class ScalarType, class FProd, class FSum> 
inline ScalarType fast_reduce(int dim, const ScalarType* x, const ScalarType* y, const FProd& fp, const FSum& fs)
{
    unsigned d = dim & ~unsigned(7);
    const ScalarType *aa = x, *end_a = aa + d;
    const ScalarType *bb = y, *end_b = bb + d;
#ifdef __GNUC__
    __builtin_prefetch(aa, 0, 3);
    __builtin_prefetch(bb, 0, 0);
#endif
    ScalarType r = 0.0;
    ScalarType r0, r1, r2, r3, r4, r5, r6, r7;

    const ScalarType *a = end_a, *b = end_b;

    r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0.0;

    switch (dim & 7) {
        case 7: r6 = fp(a[6], b[6]);
  		// fall through
        case 6: r5 = fp(a[5], b[5]);
  		// fall through
        case 5: r4 = fp(a[4], b[4]);
  		// fall through
        case 4: r3 = fp(a[3], b[3]);
  		// fall through
        case 3: r2 = fp(a[2], b[2]);
  		// fall through
        case 2: r1 = fp(a[1], b[1]);
  		// fall through
        case 1: r0 = fp(a[0], b[0]);
    }

    a = aa; b = bb;
	const auto fsum8 = [&](){
		auto r01 = fs(r0, r1);
		auto r23 = fs(r2, r3);
		auto r45 = fs(r4, r5);
		auto r67 = fs(r6, r7);
		auto r0123 = fs(r01, r23);
		auto r4567 = fs(r45, r67);
		return fs(r0123, r4567);
	};

    for (; a < end_a; a += 8, b += 8) {
#ifdef __GNUC__
        __builtin_prefetch(a + 32, 0, 3);
        __builtin_prefetch(b + 32, 0, 0);
#endif
		r = fs(r, fsum8() );
		r0 = fp(a[0], b[0]);
		r1 = fp(a[1], b[1]);
		r2 = fp(a[2], b[2]);
		r3 = fp(a[3], b[3]);
		r4 = fp(a[4], b[4]);
		r5 = fp(a[5], b[5]);
		r6 = fp(a[6], b[6]);
		r7 = fp(a[7], b[7]);
    }

    r = fs(r, fsum8() );
    return r;
}


template<class ScalarType> 
inline ScalarType calc_l2_sqr(int dim, const ScalarType* x, const ScalarType* y)
{
	const auto fProd = [](ScalarType a, ScalarType b){
		return sqr(a-b);
	};
	const auto fSum = [](ScalarType a, ScalarType b){
		return a+b;
	};
	return fast_reduce(dim, x, y, fProd, fSum);
}


template<class ScalarType> 
inline ScalarType calc_l1_dist(int dim, const ScalarType* x, const ScalarType* y)
{
	const auto fProd = [](ScalarType a, ScalarType b){
		return abs(a-b);
	};
	const auto fSum = [](ScalarType a, ScalarType b){
		return a+b;
	};
	return fast_reduce(dim, x, y, fProd, fSum);
}

template<class ScalarType> 
inline ScalarType calc_lp_dist_p(int dim, const ScalarType p, const ScalarType* x, const ScalarType* y)
{
	const auto fProd = [p](ScalarType a, ScalarType b){
		return pow(abs(a-b), p);
	};
	const auto fSum = [](ScalarType a, ScalarType b){
		return a+b;
	};
	return pow(fast_reduce(dim, x, y, fProd, fSum), 1./p);
}

template<class ScalarType> 
inline ScalarType calc_inner_product(int dim, const ScalarType* x, const ScalarType* y)
{
	const auto fProd = [](ScalarType a, ScalarType b){
		return a*b;
	};
	const auto fSum = [](ScalarType a, ScalarType b){
		return a+b;
	};
	return fast_reduce(dim, x, y, fProd, fSum);
}

// -----------------------------------------------------------------------------
template<class ScalarType>
ScalarType calc_l2_dist(					// calc L2 distance
	int   dim,							// dimension
	const ScalarType *p1,					// 1st point
	const ScalarType *p2)					// 2nd point
{
	return sqrt(calc_l2_sqr(dim, p1, p2));
}


template<class Iter>
void printVec(const Iter &its, int dim)
{
	Iter it = its;
	for(int i=0;i<dim;i++){
		std::cout << *it;
		++it;
		if(i==dim-1){
			std::cout << std::endl;
		} else{
			std::cout << ", ";
		}
	}
}

template<class Iter>
void printVec(const Iter& its, const Iter& ite, const std::string &header=""){
//	if(verbosity){
	std::cout << header;
	for(Iter it=its; ; ){
		std::cout << *it;
		++it;
		if(it!=ite){
			std::cout << ", ";
		} else{
			std::cout << std::endl;
			break;
		}
	}
//	}
}
template<class Iter, class OStream>
void printVec(const Iter& its, const Iter& ite, const std::string &header="", OStream &os=std::cout){
//	if(verbosity){
	os << header;
	for(Iter it=its; ; ){
		os << *it;
		++it;
		if(it!=ite){
			os << ", ";
		} else{
			os << std::endl;
			break;
		}
	}
//	}
}

//return idx such that xs[idx[i]] is sorted for i in range(beign, end)
//require begin and end are randomly accessable!!
template<typename T, typename F>
std::vector<int> argsort(const T& begin, const T& end, const F& cmp)
{
	size_t len = distance(begin, end);
	std::vector<int> idx(len);
	for(int i=0;i<idx.size();i++){
		idx[i] = i;
	}
	std::sort(idx.begin(), idx.end(), [&](int a, int b){
		return cmp(*(begin+a), *(begin+b));
	});
	return idx;
}
template<typename T>
std::vector<int> argsort(const T& begin, const T& end)
{
	size_t len = distance(begin, end);
	std::vector<int> idx(len);
	for(int i=0;i<idx.size();i++){
		idx[i] = i;
	}
	std::sort(idx.begin(), idx.end(), [&](int a, int b){
		return *(begin+a)<*(begin+b);
	});
	return idx;
}

template<typename Iter>
void scatter(const Iter& begin, std::vector<int> &idx)
{
    using T = typename std::iterator_traits<Iter>::value_type;
	std::vector<T> tmpxs(begin, begin+idx.size());
	for(int i=0;i<idx.size();i++){
		*(begin+i) = tmpxs[idx[i]];
	}
}

template<typename T> 
std::vector<int> getExtent(std::vector<T>& toSplit, std::vector<T>& refs)
{
	std::vector<int> ret;
	ret.reserve(refs.size()+1);

	int curLoc = 0;
	for(T& r:refs){
		while(curLoc<toSplit.size() && r>toSplit[curLoc]){
			curLoc++;
		}
		ret.push_back(curLoc);
	}
	ret.push_back(toSplit.size());
	return ret;
}

template<class uintt>
struct CountMarkerU
{
	std::vector<uintt> markCount;
	uintt curCnt;
	CountMarkerU(int sz=0):markCount(sz), curCnt(1){}


	void resize(int n){
		markCount.resize(n);
		fill(markCount.begin(), markCount.end(), 0);
	}

	void mark(int n){
		markCount[n] = curCnt;
	}
	bool isMarked(int n){
		return markCount[n] >= curCnt;
	}
	void clear(){
		if(curCnt==~uintt(0)){
			curCnt=1;
			markCount.clear();
		} else{
			curCnt++;
		}
	}
};

using CountMarker = CountMarkerU<unsigned>;

//array version of CountMatker
template<class T> 
struct CountArray
{
	CountMarker marker;
	CountArray(int sz):count(sz), marker(sz){}

	std::vector<T> count;

	void resize(int n){
		marker.resize(n);
		count.resize(n);
		marker.clear();
	}

	T& operator[](int i){
		if(!marker.isMarked(i)){
			marker.mark(i);
			count[i] = 0;
			return count[i];
		}
		return count[i];
	}

	void clear(){
		marker.clear();
	}
};



