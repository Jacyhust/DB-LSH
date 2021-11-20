#include "ConfigTool.h"

ConfigTool::ConfigTool()
{
	line_max=defaultBufferMax;
	map_str_match=true;
	map_str=new char[1];
	map_str[0]='\0';
	cfg_key=new char[line_max];
}

ConfigTool::ConfigTool(int custom_line_max)
{
	line_max=custom_line_max;
	map_str_match=true;
	map_str=new char[1];
	map_str[0]='\0';
	cfg_key=new char[line_max];
}

ConfigTool::~ConfigTool()
{
	if(map_str!=NULL)
	{
		delete[] map_str;
		map_str=NULL;
	}
	if(cfg_key!=NULL)
	{
		delete[] cfg_key;
		cfg_key=NULL;
	}
}

void ConfigTool::TrimSpace(char *str)
{
	if(str==NULL) return;

	int len=(int)strlen(str);
	int pos=0;
	for(int i=0;i<len;i++)
	{
		switch(str[i])
		{
		case '\t':
		case '\n':
		case '\f':
		case '\r':
		case ' ':
			continue;
		default:
			str[pos++]=str[i];
		}
	}
	str[pos]='\0';
}

bool ConfigTool::IfNumChar(char ch)
{
	if(((ch>='0')&&(ch<='9'))||(ch=='.')||(ch=='+')||(ch=='-')) return true;
	else return false;
}

char *ConfigTool::CombinedKey(const char *config_key, const char *prefix)
{
	if(cfg_key==NULL) cfg_key=new char[line_max];

	if(prefix==NULL)
	{
		memcpy(cfg_key, config_key, strlen(config_key)+1);
		return cfg_key;
	}

	int cfg_len=strlen(config_key);
	int cur=strlen(prefix);

	if((cur+cfg_len+2)>line_max)
	{
		printf("Combined key too long!\n");
		exit(-1);
	}

	memcpy(cfg_key, prefix, cur);
	cfg_key[cur++]='.';
	memcpy(cfg_key+cur, config_key, cfg_len+1); //we have to add one because we want to copy '\0'

	return cfg_key;
}

bool ConfigTool::IfExist(const char *config_key, const char *prefix)
{
	return (the_map.count(CombinedKey(config_key, prefix))>0);
}

bool ConfigTool::AddConfigFromFile(const char *filename, const char *prefix)
{
	char *line=new char[line_max];
	char *key=new char[line_max];
	char *value=new char[line_max];

	ifstream config_file(filename);
	if(!config_file.is_open())
	{
		char *fn_guess=new char[line_max];
		int ret=c99_snprintf(fn_guess, line_max, "../%s", filename);
		if(ret>=line_max) //error handling
		{
			printf("c99_snprintf overflow in ConfigTool::AddConfigFromFile\n");
			exit(-1);
		}
		config_file.clear();
		config_file.open(fn_guess);
		if(!config_file.is_open())
		{
			printf("Cannot find \"%s\" nor \"%s\"\n", filename, fn_guess);
			config_file.close();
			delete[] value;
			delete[] key;
			delete[] line;
			return false;
		}
		delete[] fn_guess;
	}

	while(config_file.getline(line, line_max))
	{
		if(strstr(line, "//")!=NULL) continue; //ignore comments
		char *sep_pos=strchr(line, '=');
		if(sep_pos!=NULL)
		{
			int pos=int(sep_pos-line)/sizeof(char);
			int key_len=pos;
			int value_len=(int)strlen(line)-1-key_len;
			memcpy(key, line, sizeof(char)*key_len);
			key[key_len]='\0';
			memcpy(value, sep_pos+1, sizeof(char)*value_len);
			value[value_len]='\0';
			TrimSpace(key);
			TrimSpace(value);
			the_map[CombinedKey(key, prefix)]=value;
		}
	}

	config_file.close();
	map_str_match=false; //invalidate string representation buffer

	delete[] value;
	delete[] key;
	delete[] line;
	return true;
}

void ConfigTool::AddConfigFromCommandLine(int argc, char **argv, const char *prefix)
{
	int i=0;
	while(i<argc)
	{
		while((i<argc)&&(argv[i][0]!='-')) i++;
		if(i+1<argc)
		{
			char *key=argv[i]+1;
			char *value=argv[i+1];
			TrimSpace(key);
			TrimSpace(value);
			the_map[CombinedKey(key, prefix)]=value;
			i+=2;
		}
		else return;
	}
	map_str_match=false; //invalidate string representation buffer
}

void ConfigTool::ModifyConfig(const char *config_key, const char *config_value, const char *prefix)
{
	the_map[CombinedKey(config_key, prefix)]=config_value;

	map_str_match=false; //dirty flag
}

const char *ConfigTool::StringForm(char seprator)
{
	int len=0;
	int cur_len=0;

	if(!map_str_match)
	{
		//use c99 printf to calculate total output length
		//cannot simply replace c99_snprintf due to cross-platform compatibility
		map<string, string>::iterator iter=the_map.begin();
		len+=c99_snprintf(map_str, 1, "%s=%s", iter->first.c_str(), iter->second.c_str());
		++iter;
		for(;iter!=the_map.end();++iter) len+=c99_snprintf(map_str, 1, "%c%s=%s", seprator, iter->first.c_str(), iter->second.c_str());
		len++;

		delete[] map_str;
		map_str=new char[len];
		iter=the_map.begin();
		cur_len+=c99_snprintf(map_str, len, "%s=%s", iter->first.c_str(), iter->second.c_str());
		++iter;
		for(;iter!=the_map.end();++iter) cur_len+=c99_snprintf(map_str+cur_len, len-cur_len, "%c%s=%s", seprator, iter->first.c_str(), iter->second.c_str());

		map_str_match=true;
	}

	return (const char *)map_str;
}

void ConfigTool::ListConfig()
{
	printf("%s\n", StringForm());
}

int ConfigTool::GetConfigInt(const char* config_key, const char *prefix)
{
	int value=0;

	CombinedKey(config_key, prefix);

	if(the_map.count(cfg_key)>0)
	{
		const char *val=the_map[cfg_key].c_str();
		if(IfNumChar(val[0])) value=atoi(the_map[cfg_key].c_str());
		else value=GetConfigInt(val); //else we search other gobal key without current prefix
	}
	else
	{
		printf("Config key \"%s\" missing!\n", cfg_key);
		exit(-1);
	}

	return value;
}

float ConfigTool::GetConfigFloat(const char* config_key, const char *prefix)
{
	float value=0;

	CombinedKey(config_key, prefix);

	if(the_map.count(cfg_key)>0)
	{
		const char *val=the_map[cfg_key].c_str();
		if(IfNumChar(val[0])) value=(float)atof(the_map[cfg_key].c_str());
		else value=GetConfigFloat(val); //else we search other gobal key without current prefix
	}
	else
	{
		printf("Config key \"%s\" missing!\n", cfg_key);
		exit(-1);
	}

	return value;
}

double ConfigTool::GetConfigDouble(const char* config_key, const char *prefix)
{
	double value=0;

	CombinedKey(config_key, prefix);

	if(the_map.count(cfg_key)>0)
	{
		const char *val=the_map[cfg_key].c_str();
		if(IfNumChar(val[0])) value=atof(the_map[cfg_key].c_str());
		else value=GetConfigDouble(val); //else we search other gobal key without current prefix
	}
	else
	{
		printf("Config key \"%s\" missing!\n", cfg_key);
		exit(-1);
	}

	return value;
}

const char *ConfigTool::GetConfigStr(const char* config_key, const char *prefix)
{
	const char *value=NULL;

	CombinedKey(config_key, prefix);

	if(the_map.count(cfg_key)>0) value=the_map[cfg_key].c_str();
	else
	{
		printf("Config key \"%s\" missing!\n", cfg_key);
		exit(-1);
	}

	return value;	
}

int ConfigTool::GetConfigIntArray(const char* config_key, int *data_array, const char *prefix)
{
	char *str=new char[line_max];

	int count=0;

	CombinedKey(config_key, prefix);

	if(the_map.count(cfg_key)>0)
	{
		memcpy(str, the_map[cfg_key].c_str(), (the_map[cfg_key].length()+1)*sizeof(char));
	}
	else
	{
		printf("Config key \"%s\" missing!\n", cfg_key);
		delete str;
		exit(-1);
	}

	char *p=str, *q;
	while((q=strchr(p, ','))!=NULL)
	{	   
		*q='\0';
		data_array[count++]=atoi(p);
		p=q+1;
	}
	data_array[count++]=atoi(p);

	delete[] str;

	return count;
}

int ConfigTool::GetConfigFloatArray(const char* config_key, float *data_array, const char *prefix)
{
	char *str=new char[line_max];

	int count=0;

	CombinedKey(config_key, prefix);

	if(the_map.count(cfg_key)>0)
	{
		memcpy(str, the_map[cfg_key].c_str(), (the_map[cfg_key].length()+1)*sizeof(char));
	}
	else
	{
		printf("Config key \"%s\" missing!\n", cfg_key);
		delete str;
		exit(-1);
	}

	char *p=str, *q;
	while((q=strchr(p, ','))!=NULL)
	{	   
		*q='\0';
		data_array[count++]=(float)atof(p);
		p=q+1;
	}
	data_array[count++]=(float)atof(p);

	delete[] str;

	return count;
}

int ConfigTool::GetConfigDoubleArray(const char* config_key, double *data_array, const char *prefix)
{
	char *str=new char[line_max];

	int count=0;

	CombinedKey(config_key, prefix);

	if(the_map.count(cfg_key)>0)
	{
		memcpy(str, the_map[cfg_key].c_str(), (the_map[cfg_key].length()+1)*sizeof(char));
	}
	else
	{
		printf("Config key \"%s\" missing!\n", cfg_key);
		delete str;
		exit(-1);
	}

	char *p=str, *q;
	while((q=strchr(p, ','))!=NULL)
	{	   
		*q='\0';
		data_array[count++]=atof(p);
		p=q+1;
	}
	data_array[count++]=atof(p);

	delete[] str;

	return count;
}