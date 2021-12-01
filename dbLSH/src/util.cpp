#include "util.h"


using namespace std;

//#if defined(unix) || defined(__unix__)
////under POSIX system, we use clock_gettime()
////remember we have to use linker option "-lrt"
//#else
//int access(const char* pathname, int mode)
//{
//	return _access(pathname, mode);
//}
//#endif

// -----------------------------------------------------------------------------
float calc_angle_normalized(				// calc angle
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	return acos(calc_inner_product(dim, p1, p2));
}

// -----------------------------------------------------------------------------
float calc_cosangle(				// calc cos(angle)
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	double ret = 0.0f;
	double norm0 = 0., norm1 = 0.;
	for (int i = 0; i < dim; ++i) {
		ret += (p1[i] * p2[i]);
		norm0 += p1[i] * p1[i];
		norm1 += p2[i] * p2[i];
	}
	if(norm0==0 || norm1==0){
		return 0;
	}
	return ret/sqrt(norm0*norm1);
}

// // -----------------------------------------------------------------------------
// float calc_l2_sqr(					// calc L2 square distance
// 	int   dim,							// dimension
// 	const float *p1,					// 1st point
// 	const float *p2)					// 2nd point
// {
// 	float diff = 0.0f;
// 	float ret  = 0.0f;
// 	for (int i = 0; i < dim; ++i) {
// 		diff = p1[i] - p2[i];
// 		ret += diff * diff;
// 	}
// 	return ret;
// }

const int PrefixTableSize = 1 << 16;
//std::array<uint8_t, PrefixTableSize> _prefix_table;

#if defined(unix) || defined(__unix__)
	std::array<uint8_t, PrefixTableSize> _prefix_table;
#else
#include <vector>
std::vector<uint8_t> _prefix_table(PrefixTableSize);
#endif




bool init_prefix_table()
{
	for (int i = 0; i < PrefixTableSize; i++) {
		//calculate the prefix-1 of i, since it will be run only once, implement using stupid way
		_prefix_table[i] = 0;
		for (int j = 0; j < 16; j++) {
			int mask = 1 << (15 - j);
			if (i&mask) {
				_prefix_table[i]++;
			}
			else {
				break;
			}
		}
	}
	return true;
}

int get_num_prefix(uint16_t u)
{
	static bool initialized = init_prefix_table();
	return _prefix_table[u];
}

int get_num_prefix(uint32_t u)
{
	int a = get_num_prefix(uint16_t(u >> 16));
	if (a != 16) {
		return a;
	}
	int b = get_num_prefix(uint16_t(u & 0xffff));
	return a + b;
}

int get_num_prefix(uint64_t u)
{
	int a = get_num_prefix(uint16_t(u >> 48));
	if (a != 16) {
		return a;
	}
	int b = get_num_prefix(uint16_t((u >> 32) & 0xffff));
	if (b != 16) {
		return a + b;
	}
	int c = get_num_prefix(uint16_t((u >> 16) & 0xffff));
	if (c != 16) {
		return a + b + c;
	}
	int d = get_num_prefix(uint16_t(u & 0xffff));
	return a + b + c + d;
}
