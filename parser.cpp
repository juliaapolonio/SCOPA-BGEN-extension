#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
// #include <zlib.h>
#include <cstring>

#include "helper.h"

using std::string;
using namespace std;

enum BGEN_RETURN_CODE
{
	E_OK = 0,
	E_WRITE_ERROR = -1
};

struct file_info{
	file_info(){}
	file_info(uint32_t header_offset, uint32_t size_of_header, uint32_t num_variants, uint32_t num_samples, string magic, int freedata_size, string freedata, uint32_t flags, uint32_t empty_space){
		this->header_offset=header_offset;
		this->size_of_header=size_of_header;
		this->num_variants=num_variants;
		this->num_samples=num_samples;
		this->magic=magic;
		this->freedata_size=freedata_size;
		this->freedata=freedata;
		this->flags=flags;
		this->empty_space=empty_space;

		this->compression=(flags & 3);
		this->bgen_version=((flags & 60)>>2)-1;
	}
	uint32_t header_offset;
	uint32_t size_of_header;
	uint32_t num_variants;
	uint32_t num_samples;
	string magic;
	int freedata_size;
	string freedata;
	uint32_t flags;
	unsigned short int compression;
	bool bgen_version;
	bool sample_info_included;
	uint32_t empty_space;
};

file_info parse_header(ifstream & bgenfile, bool VERBOSE){
	uint32_t header_offset;
	uint32_t size_of_header;
	uint32_t num_variants;
	uint32_t num_samples;
	bgenfile.read(reinterpret_cast<char *>(&header_offset), sizeof(header_offset));
	bgenfile.read(reinterpret_cast<char *>(&size_of_header), sizeof(size_of_header));
	bgenfile.read(reinterpret_cast<char *>(&num_variants), sizeof(num_variants));
	bgenfile.read(reinterpret_cast<char *>(&num_samples), sizeof(num_samples));
	char magic[4];
	bgenfile.read(magic, 4);
	if(VERBOSE){info("Offset is: ",header_offset);}
	if(VERBOSE){info("Size of header is: ", size_of_header);}
	if(VERBOSE){info("Number of variants: ", num_variants);}
	if(VERBOSE){info("Number of samples: ", num_samples);}
	if(VERBOSE){info("Magic number: ", magic);}
	int freedata_size=0;
	string freedata_s;
	if(size_of_header-20>0){
		freedata_size=size_of_header-20;
		char freedata[freedata_size];
		bgenfile.read(freedata, freedata_size);
		if(VERBOSE){info("Free data detected: ", freedata);}
		freedata_s=std::string(freedata);
	}
	uint32_t flags;
	bgenfile.read(reinterpret_cast<char *>(&flags), sizeof(flags));
	switch(flags & 3){
		case 0 : 
		if(VERBOSE){info("SNP Block probability data is not compressed.");}
		break;
		case 1 :
		if(VERBOSE){info("SNP Block probability data compressed using zlib.");}
		break;
		case 2 :
		if(VERBOSE){info("SNP Block probability data compressed using zstandard.");}
		break;
	}
	switch((flags & 60)>>2){
		case 1:
		if(VERBOSE){info("Bgen v1.1 detected.");}
		break;
		case 2:
		if(VERBOSE){info("Bgen v1.2 detected.");}
		break;
	}

	uint32_t skip=header_offset-size_of_header;
	if(VERBOSE){info("Skipping ", skip, "bytes.");}
	bgenfile.ignore(skip);
	file_info bgen_file_info(header_offset, size_of_header, num_variants, num_samples, std::string(magic), freedata_size, freedata_s, flags, skip);

	return(bgen_file_info);
}

BGEN_RETURN_CODE write_header(ofstream & bgenfile, file_info f_info){

	// TBD: the header data needs to be calculated at the end as it contains numbers of variants and samples
	// TBD: 

	// Write basic header
	try{
	bgenfile.write((char *) &(f_info.header_offset), 4);
	bgenfile.write((char *) &(f_info.size_of_header), 4);
	bgenfile.write((char *) &(f_info.num_variants), 4);
	bgenfile.write((char *) &(f_info.num_samples), 4);
	bgenfile.write((char *) &(f_info.magic), 4);
	bgenfile.write((char *) &(f_info.freedata), f_info.freedata_size);
	bgenfile.write((char *) &(f_info.flags), 4);
	}
	catch (...) {
		return(E_WRITE_ERROR);
	}
	// If sample info is present, write sampe identifier block


	bgenfile.flush();
	return(E_OK);
}

int main(int argc, char* argv[])
{
	bool VERBOSE=0;
	ifstream bgenfile("C:/Users/jiaq8/OneDrive/FYP/example.10bits.bgen", ios::in | ios::binary);
	file_info bgen_finfo;
	bgen_finfo=parse_header(bgenfile, 0);
    cout << bgen_finfo.header_offset << endl;
    cout << bgen_finfo.size_of_header << endl;
    cout << bgen_finfo.num_variants << endl;
    cout << bgen_finfo.num_samples << endl;
    cout << bgen_finfo.flags << endl;
    cout << bgen_finfo.compression << endl;
	ofstream out("C:/Users/jiaq8/OneDrive/FYP", ios :: out | ios :: binary);
	write_header(out, bgen_finfo);
	bgenfile.close();
	exit(0);
	unsigned short id_len;
	string id;
	string rsid;
	string chr;
	uint32_t pos;
}
