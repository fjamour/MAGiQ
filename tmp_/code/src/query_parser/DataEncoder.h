/* 
 * File:   DataEncoder.h
 * Author: abdelai
 *
 * Created on November 13, 2017, 3:08 PM
 */

#ifndef DATAENCODER_H
#define	DATAENCODER_H

#include "partitioner_store.h"
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
using namespace std;
int main_d_encoder(int argc, const char** argv) {
	if(argc < 4){
		print_to_screen("USAGE: RDFEncoder input_dir output_file_path remove_double_edges");
		exit(-1);
	}
	string input_dir = string(argv[1]);
	string output_file_name = string(argv[2]);
        string param = string(argv[3]);
	partitioner_store * store;

        bool remove_double_edges;
        if (param == "1") 
            remove_double_edges = true;
        else
            remove_double_edges = false;
        	
        print_to_screen("remove_double_edges: ");
        cout<<remove_double_edges<<"\n";
        
	store = new partitioner_store();
	store->load_encode_rdf_data(input_dir, output_file_name, remove_double_edges);
	//store->dump_encoded_data(file_name);

        
        
	delete store;
	print_to_screen(part_string);
	print_to_screen("Run Finished Successfully");

}
#endif	/* DATAENCODER_H */

