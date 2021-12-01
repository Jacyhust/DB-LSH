#pragma once
#include "StructType.h"
#include "Preprocess.h"
#include "RStarTree.h"
#include <cmath>
#include <assert.h>
#include <unordered_map>
#include <vector>
#include <queue>

class Hash
{
private:
	std::string index_file;
	std::vector<TreeDataP<float>**> hs;
public:
	int N = 0;
	int dim = 0;
	// Number of hash functions
	int S = 0;
	int L = 0;
	int K = 0;
	int MaxSize = 0;
	float R_min = 1.0f;
	
	HashParam hashpar;
	RStarTree<TreeDataP, float>** myIndexes = NULL;

public:
	Hash(Preprocess& prep_, Parameter& param_, const std::string& file);
	void SetHash();
	void GetHash(Preprocess& prep);
	void GetTables(Preprocess& prep);
	bool IsBuilt(const std::string& file);
	~Hash();
};

