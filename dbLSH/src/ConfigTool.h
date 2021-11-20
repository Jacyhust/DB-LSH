//This file includes the interface for handling config input for a program
//by GS
#ifndef _CONFIG_TOOL_H_
#define _CONFIG_TOOL_H_

#include "GenericTool.h"

#include <cstdlib>
#include <cstring>

#include <cstring>
#include <string>
#include <fstream>
#include <map>

using namespace std;

const int defaultBufferMax=4096;

class ConfigTool
{
public:
	ConfigTool();
	ConfigTool(int custom_line_max);
	~ConfigTool();

	static void TrimSpace(char *str);
	static bool IfNumChar(char ch);
	char *CombinedKey(const char *config_key, const char *prefix);

	bool IfExist(const char *config_key, const char *prefix=NULL);
	bool AddConfigFromFile(const char *filename, const char *prefix=NULL);
	void AddConfigFromCommandLine(int argc, char **argv, const char *prefix=NULL);
	void ModifyConfig(const char *config_key, const char *config_value, const char *prefix=NULL);
	const char *StringForm(char seprator='\t');
	void ListConfig();

	int GetConfigInt(const char* config_key, const char *prefix=NULL);
	float GetConfigFloat(const char* config_key, const char *prefix=NULL);
	double GetConfigDouble(const char *config_key, const char *prefix=NULL);
	const char *GetConfigStr(const char* config_key, const char *prefix=NULL);

	int GetConfigIntArray(const char* config_key, int *data_array, const char *prefix=NULL);
	int GetConfigFloatArray(const char* config_key, float *data_array, const char *prefix=NULL);
	int GetConfigDoubleArray(const char* config_key, double *data_array, const char *prefix=NULL);
protected:
	map<string, string> the_map;
	int line_max;
	bool map_str_match;
	char *map_str;
	char *cfg_key;
};

#endif