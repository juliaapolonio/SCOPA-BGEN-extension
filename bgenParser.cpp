#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
// #include <zlib.h>

using std::string;
using namespace std;

int main(int argc, char* argv[])
{
    uint32_t header_offset;
	uint32_t size_of_header;
	uint32_t num_variants;
	uint32_t num_samples;
    char magic_number[4];

	ifstream bgenfile("C:/Users/jiaq8/OneDrive/FYP/example.10bits.bgen", ios::in | ios::binary);
	bgenfile.read(reinterpret_cast<char *>(&header_offset), sizeof(header_offset));
	bgenfile.read(reinterpret_cast<char *>(&size_of_header), sizeof(size_of_header));
	bgenfile.read(reinterpret_cast<char *>(&num_variants), sizeof(num_variants));
	bgenfile.read(reinterpret_cast<char *>(&num_samples), sizeof(num_samples));
    bgenfile.read(reinterpret_cast<char *>(&magic_number), 32);

    uint32_t free_data_size;
    free_data_size = size_of_header - 20;
    char free_data[free_data_size];
    bgenfile.read(reinterpret_cast<char *>(&free_data), free_data_size);
    
    // flags
    uint32_t flags;
    bgenfile.read(reinterpret_cast<char *>(&flags), sizeof(flags));
    switch(flags & 3){
        case 0:
            cout << "Not compressed" << endl;
        break;
        case 1:
            cout << "Compressed with zlib" << endl;
        break;
        case 2:
            cout << "Compressed with zstandard" << endl;
        break;
    }

    switch((flags << 2) >> 30){
		case 1:
            cout << "v1.1" << endl;
		break;

		case 2:
            cout << "v1.2" << endl;
		break;
	}

    switch(flags >> 30){
		case 1:
            cout << "Sample Identifiers unavailable" << endl;
		break;

		case 2:
            cout << "Sample Identifiers available" << endl;
		break;
	}

	bgenfile.close();

    cout << "Offset = " + to_string(header_offset) << endl;
    cout << "Size of header = " + to_string(size_of_header) << endl;
    cout << "Number of variants = " + to_string(num_variants) << endl;
    cout << "Number of samples = " + to_string(num_samples) << endl;
    // cout << magic_number << endl;
    // cout << free_data_size << endl;
    // cout << free_data << endl;
    // cout << flags << endl;

	exit(0);
}
