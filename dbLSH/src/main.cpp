// PM_LSH.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <fstream>
#include "Preprocess.h"
#include "dblsh.h"
#include "Query.h"
#include <time.h>
#include "basis.h"

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

void lshknn(float c, int k, Hash& myslsh, Preprocess& prep, float beta, std::string& datasetName, std::string& data_fold);
void expe_k(float c, Hash& myslsh, Preprocess& prep, float beta, std::string& datasetName, std::string& data_fold);

int main(int argc, char const* argv[])
{
	float c = 1.5;
	unsigned k = 50;
	unsigned L = 8, K = 10;//NUS
	L = 5, K = 10;
	float beta = 0.1;
	unsigned Qnum = 100;
	float R_min = 1.0f;

	std::string datasetName;
	if (argc == 2) {
		datasetName = argv[1];
		set_rmin(datasetName, R_min);
	}
	else if (argc > 2) {
		if ((argc < 6 || argc > 8)) {
			std::cerr << "Usage: ./dblsh datasetName approx_ratio k L K beta R_min(optinal)\n\n";
			exit(-1);
		}
		datasetName = argv[1];
		c = std::atof(argv[2]);
		k = std::atoi(argv[3]);
		L = std::atoi(argv[4]);
		K = std::atoi(argv[5]);
		beta = std::atof(argv[6]);
		if (argc == 8) {
			R_min = std::atof(argv[7]);
		}
		else {
			set_rmin(datasetName, R_min);
		}
	}
	else//only for debug, not advised for user
	{
		/*std::cout << BOLDYELLOW << "Usage:" << YELLOW << " ./dblsh datasetName approx_ratio k L K beta R_min(optinal)\n\n" << RESET;
		exit(-1);*/

		const std::string datas[] = { "audio","mnist","cifar","NUS","Trevi","gist","deep1m" };
		datasetName = datas[0];
		datasetName = "audio";

		set_rmin(datasetName, R_min);
		std::cout << "Using the default configuration!\n\n";
	}

	/// <summary>
	/// Show the configuration
	/// </summary>
	/// <param name="argc"></param>
	/// <param name="argv"></param>
	/// <returns></returns>
	std::cout << "Using DB-LSH for " << datasetName << " ..." << std::endl;
	std::cout << "c=        " << c << std::endl;
	std::cout << "k=        " << k << std::endl;
	std::cout << "L=        " << L << std::endl;
	std::cout << "K=        " << K << std::endl;
	std::cout << "beta=     " << beta << std::endl;
	//std::cout << "R_min=    " << R_min << std::endl << std::endl;

	#if defined(unix) || defined(__unix__)
		std::string data_fold = "./../dataset/", index_fold = "";
	#else
		std::string data_fold = "E:/Dataset_for_c/", index_fold = "";
	#endif

	Preprocess prep(data_fold + datasetName + ".data", data_fold + "ANN/" + datasetName + ".bench");

	showMemoryInfo();

	Parameter param(prep, L, K, R_min);
	Hash myslsh(prep, param, index_fold.append(datasetName));

	if (beta > 0) {
		lshknn(c, k, myslsh, prep, beta, datasetName, data_fold);
	}
	else if (beta <= 0 && beta > -10) {
		std::vector<float> Betas = { 0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.2 };
		
		for (auto& x : Betas) {
			if (datasetName == "sift10M") {
				x = x / 4;
			}
			lshknn(c, k, myslsh, prep, x, datasetName, data_fold);
		}
	}
	else {
		beta = 0.1;
		expe_k(c, myslsh, prep, beta, datasetName, data_fold);
	}

	

	return 0;
}

void lshknn(float c, int k, Hash& myslsh, Preprocess& prep, float beta, std::string& datasetName, std::string& data_fold) {
	lsh::timer timer;
	std::cout << std::endl << "RUNNING QUERY ..." << std::endl;
	int Qnum = 100;
	lsh::progress_display pd(Qnum);
	Performance perform;
	for (unsigned j = 0; j < Qnum; j++)
	{
		Query query(j, c, k, myslsh, prep, beta);
		perform.update(query, prep);
		++pd;
	}

	showMemoryInfo();

	float mean_time = (float)perform.time_total / perform.num;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl << std::endl;
	//std::cout << "SORT TIME:         " << ((float)perform.time_sift) / (perform.num) << std::endl;
	//std::cout << "AVG QUERY TIME:    " << (float)perform.time_verify / perform.num * 1000 << "ms." << std::endl << std::endl;
	std::cout << "AVG RECALL:        " << ((float)perform.NN_num) / (perform.num * k) << std::endl;
	std::cout << "AVG RATIO:         " << ((float)perform.ratio) / (perform.res_num) << std::endl;
	std::cout << "AVG COST:          " << ((float)perform.cost) / ((float)perform.num * prep.data.N) << std::endl;
	std::cout << "AVG ROUNDS:        " << ((float)perform.rounds) / (perform.num) << std::endl;
	//std::cout << "AVG ACCESSES:      " << ((float)perform.num_access_in_RTree) / (perform.num) << std::endl;
	std::cout << "\nQUERY FINISH... \n\n\n";

	time_t now = std::time(0);
	time_t zero_point = 1635153971;//Let me set the time at 2021.10.25. 17:27 as the zero point
	float date = ((float)(now - zero_point)) / 86400;

	std::string fpath = data_fold + "ANN/";

	if (!GenericTool::CheckPathExistence(fpath.c_str())) {
		GenericTool::EnsurePathExistence(fpath.c_str());
		std::cout << BOLDGREEN << "WARNING:\n" << GREEN << "Could not find the path of result file. Have created it. \n"
			<< "The query result will be stored in: " << fpath.c_str() << RESET;
	}
	std::ofstream os(fpath + "DB-LSH_result.csv", std::ios_base::app);
	if (os) {
		os.seekp(0, std::ios_base::end); // move to the end of file
		int tmp = (int)os.tellp();
		if (tmp == 0) {
			os << "Dataset,c,k,L,K,R_min,RATIO,RECALL,AVG_TIME,COST,DATE" << std::endl;
		}
		std::string dataset = datasetName;
		os << dataset << ',' << c << ',' << k << ',' << myslsh.L << ',' << myslsh.K << ',' << myslsh.R_min << ','
			<< ((float)perform.ratio) / (perform.res_num) << ','
			<< ((float)perform.NN_num) / (perform.num * k) << ','
			<< mean_time * 1000 << ','
			<< ((float)perform.cost) / (perform.num * prep.data.N) << ','
			<< date << ','
			<< std::endl;
		os.close();
	}
}

void expe_k(float c, Hash& myslsh, Preprocess& prep, float beta, std::string& datasetName, std::string& data_fold) {
	lsh::timer timer;
	std::cout << std::endl << "RUNNING QUERY ..." << std::endl;
	int Qnum = 100;
	lsh::progress_display pd(Qnum);
	

	int k = 100;
	Performance performs[11];
	int Ks[] = { 1,10,20,30,40,50,60,70,80,90,100 };
	for (unsigned j = 0; j < Qnum; j++)
	{
		Query query(j, c, k, myslsh, prep, beta);
		for (int m = 0; m < 11; ++m) {
			query.k = Ks[m];
			performs[m].update(query, prep);
		}
		++pd;
	}

	showMemoryInfo();

	for (int m = 0; m < 11; ++m) {
		Performance perform = performs[m];
		k = Ks[m];
		float mean_time = (float)perform.time_total / perform.num;
		std::cout << "k=                 " << k << std::endl << std::endl;
		std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl << std::endl;
		//std::cout << "SORT TIME:         " << ((float)perform.time_sift) / (perform.num) << std::endl;
		//std::cout << "AVG QUERY TIME:    " << (float)perform.time_verify / perform.num * 1000 << "ms." << std::endl << std::endl;
		std::cout << "AVG RECALL:        " << ((float)perform.NN_num) / (perform.num * k) << std::endl;
		std::cout << "AVG RATIO:         " << ((float)perform.ratio) / (perform.res_num) << std::endl;
		std::cout << "AVG COST:          " << ((float)perform.cost) / ((float)perform.num * prep.data.N) << std::endl;
		std::cout << "AVG ROUNDS:        " << ((float)perform.rounds) / (perform.num) << std::endl;
		std::cout << "\nFINISH... \n\n\n";

		time_t now = std::time(0);
		time_t zero_point = 1635153971;//Let me set the time at 2021.10.25. 17:27 as the zero point
		float date = ((float)(now - zero_point)) / 86400;
		std::ofstream os(data_fold + "ANN/DB-LSH_result.csv", std::ios_base::app);
		os.seekp(0, std::ios_base::end); // move to the end of file
		int tmp = (int)os.tellp();
		if (tmp == 0) {
			os << "Dataset,c,k,L,K,R_min,RATIO,RECALL,AVG_TIME,COST,DATE" << std::endl;
		}
		std::string dataset = datasetName;
		os << dataset << ',' << c << ',' << k << ',' << myslsh.L << ',' << myslsh.K << ',' << myslsh.R_min << ','
			<< ((float)perform.ratio) / (perform.res_num) << ','
			<< ((float)perform.NN_num) / (perform.num * k) << ','
			<< mean_time * 1000 << ','
			<< ((float)perform.cost) / (perform.num * prep.data.N) << ','
			<< date << ','
			<< std::endl;
		os.close();
	}
	
}
