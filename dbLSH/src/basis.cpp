#include "basis.h"
#include <cmath>
#include <vector>
#include <random>
#include <set>
#include <assert.h>
#include <omp.h>

float cal_inner_product(float* v1, float* v2, int dim)
{
	float res = 0.0;
//#pragma omp parallel for
	for (int i = 0; i < dim; ++i) {
		res += v1[i] * v2[i];
	}
	return res;
}

float cal_dist(float* v1, float* v2, int dim)
{
	float res = 0.0;
	int i = 0;
#pragma omp parallel for reduction(+:res) default(shared) num_threads(4)
	for (i = 0; i < dim; ++i) {
		auto& e1 = v1[i]; auto& e2 = v2[i];
		res += (e1 - e2) * (e1 - e2);
		//res += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	}
//#pragma omp critical
	return sqrt(res);
}

float cal_l2_dist(std::vector<float>& a, std::vector<float>& b)
{
	//assert(a.size() == b.size());
	if (a.size() != b.size()) {
		std::cout << "the size don't match in Cal_Euclidean_distance" << std::endl;
		std::cout << a.size() << "," << b.size() << "\n";
		system("pause");
	}
	float sum = 0;
	for (int i = 0; i < a.size(); ++i) {
		sum += (a[i] - b[i]) * (a[i] - b[i]);
	}
	return sqrt(sum);
}

float cal_l2_dist(const std::vector<float>& a, const std::vector<float>& b)
{
	assert(a.size() == b.size());
	if (a.size() != b.size()) { 
		std::cout << "the size don't match in Cal_Euclidean_distance" << std::endl;
		std::cout << a.size() << "," << b.size() << "\n";
	}
	float sum = 0;
	for (int i = 0; i < a.size(); ++i) {
		sum += (a[i] - b[i]) * (a[i] - b[i]);
	}
	return sqrt(sum);
}

std::vector<std::vector<float>> Random_pivot(float** raw_data_, int pivot_num_, int N_, int dim, int repeat_num)
{
	std::default_random_engine rng(0);
	std::uniform_int_distribution<int> u_random_int(0, N_ - 1);
	std::vector<std::vector<std::vector<float>>> N_pivot_label(repeat_num);
	for (int i = 0; i < repeat_num; ++i) {
		std::set<int> repeat_flag_set;
		std::vector<std::vector<float>> pivot_label(pivot_num_);
		for (int j = 0; j < pivot_num_; ++j) {
			int label = u_random_int(rng);
			if (repeat_flag_set.count(label) == 0)
			{
				pivot_label[j].resize(dim);
				for (int i = 0; i < dim; ++i) {
					pivot_label[j][i] = raw_data_[label][i];
				}
					
				repeat_flag_set.insert(label);
			}
			else
				--j;
		}
		N_pivot_label.push_back(pivot_label);
	}

	/*maximum distantce*/
	std::pair<float, int> res_label(INT16_MIN, -1);
	for (int i = 0; i < N_pivot_label.size(); ++i) {
		float dist = distance_Between_Piovt_Vector(N_pivot_label[i]);
		if (dist > res_label.first)
			res_label.second = i;
	}
	return N_pivot_label[res_label.second];
}

float distance_Between_Piovt_Vector(std::vector<std::vector<float>>& pivot_vector_)
{
	float all_dist = 0;
	for (int i = 0; i < pivot_vector_.size(); ++i) {
		for (int j = i + 1; j < pivot_vector_.size(); ++j) {
			all_dist += cal_l2_dist(pivot_vector_[i], pivot_vector_[j]);
		}
	}
	return all_dist;
}

void set_rmin(std::string& datatsetName, float& R_min)
{
	if (datatsetName == "audio") {
		R_min = 0.3f;
	}
	else if (datatsetName == "mnist") {
		R_min = 500.0f;
	}
	else if (datatsetName == "cifar") {
		R_min = 300.0f;
	}
	else if (datatsetName == "deep1m") {
		R_min = 0.2f;
	}
	else if (datatsetName == "NUS") {
		R_min = 6.50f;
	}
	else if (datatsetName == "Trevi") {
		R_min = 700.0f;
	}
	else if (datatsetName == "gist") {
		R_min = 0.3f;
	}
}

#if defined(unix) || defined(__unix__)
void showMemoryInfo(){}
#else
	#include <iostream>
	#include <windows.h>
	#include <psapi.h>
	#pragma comment(lib,"psapi.lib")
	using namespace std;
	void showMemoryInfo()
	{
		HANDLE handle = GetCurrentProcess();
		PROCESS_MEMORY_COUNTERS pmc;
		GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
		cout << "Memory Usage: " << pmc.WorkingSetSize / (1024 * 1024) << "M/" << pmc.PeakWorkingSetSize / (1024 * 1024) << "M + " << pmc.PagefileUsage / (1024 * 1024) << "M/" << pmc.PeakPagefileUsage / (1024 * 1024) << "M." << endl;
	}
#endif
