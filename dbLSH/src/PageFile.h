#ifndef _PAGE_FILE_H_
#define _PAGE_FILE_H_

#include <cstring>
#include <list>
#include <vector>
#include <fstream>
using namespace std;

#include "GenericTool.h"

//Provide functionality for cached page file
//we use one page to put header for maintennance
//and the header of objects using this file as storage
class PageFile
{
public:
	fstream *f; // lower level file object
	char *fname; // stored filed name
	bool new_file; // is this a new file?

	int page_len; // the page length
	int page_num; // total number of pages

	//cache
	enum c_policy {C_FULLMEM=0, C_LRU=1, C_MRU=2, C_NO_CACHE=3};

	c_policy policy; // current policy
	int cache_size; // the size of cache
	char *buf_page; // one page buffer used to store page when no cache is used
	int num_cached; // the number of cached pages
	vector<char *> *cache; // the cache map from page id to buffer pointer
	vector<bool> *c_dirty; // the cache map from page id to if it is dirty (been modified)
	list<int> *use_list; // stores the page id list fron (most recent) -> end (least recent)
	vector<list<int>::iterator> *list_pos; // the cache map from page id to buffer pointer

	//statistics
	int page_access;

	// p_len is effective only for new files
	// pol is cache policy of {C_FULLMEM, C_LRU, C_MRU, C_NO_CACHE}
	// c_size is size of cache
	PageFile(const char *name, int p_len, const char *path=NULL, c_policy pol=C_FULLMEM, int c_size=0, bool force_new=false);
	~PageFile();

	void Construct(const char *name, int p_len, const char *path, c_policy pol, int c_size, bool force_new);
	void EmptyValue();
	void CleanValue();

	const char *read_header(); // read header block without blockfile internal header
	void set_header(const char *header); // set block fileheader without blockfile internal header
	int get_header_size() { return page_len-2*sizeof(int); } // get actual header length

	//the index is the number except header
	//i.e, header is -1, and the actual page after header is 0
	const char *read_page(int index); // reads page
	bool write_page(const char *p, int index); // write page
	int append_page(const char *p); // returns page number
	bool truncate_pages(int num); // deletes the last num pages

	//write all page back to disk
	void flush();

protected:
	//avoid copy
	PageFile(const PageFile &pf);
};

#endif