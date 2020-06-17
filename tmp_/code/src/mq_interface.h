
/*
 File Header TODO have a nice header 
 */


#include "mq_query_exec.h"
#include "DataEncoder.h"
#include "QueryEncoder.h"

#include "mq_matlab_compiler.h"

int interface_main(int argc, char**argv) {
    //./dist/Release/GNU-Linux-x86/magic 1 data/lubm_160/ data/lubm_160/encoded.nt
    //./dist/Release/GNU-Linux-x86/magic 2 data/lubm_160/encoded.ntpredicate_map.dic data/lubm_160/encoded.ntverts_map.dic  data/lubm_160/queries/
    //./dist/Release/GNU-Linux-x86/magic 3 data/lubm_160/encoded.nt data/lubm_160/queries/encoded/ data/lubm_160/queries/magiq_res/
    
    if(argc < 2){
	cout<<"USAGE: MAGIC <operator> <param1> <param2> ..."<<endl;
        cout<<"operator values: (1 > encode data, 2 > encode queries, 3 > evaluate queries)"<<endl;
	exit(-1);
    }
    
    string oper = string(argv[1]);
    if(oper == "1")
    {
        if(argc < 4){
            cout<<"USAGE for encoding data: MAGIC <1> <input_dir> <output_file_name> <remove_double_edges>"<<endl;
            cout<<"Example: ./dist/Release/GNU-Linux-x86/magic 1 data/lubm_160/ data/lubm_160/encoded.nt 0"<<endl;
            exit(-1);
        }
        
        string input_dir = string(argv[2]);    
        string output_file_name = string(argv[3]);    
        string remove_double_edges = string(argv[4]);    
                 
        string tool = "DataEncoder";
        const char* data_arg[4];
        data_arg[0] = tool.c_str();
        data_arg[1] = input_dir.c_str();
        data_arg[2] = output_file_name.c_str();
        data_arg[3] = remove_double_edges.c_str();
        main_d_encoder(4, data_arg);
    
    }
    
    if(oper == "2")
    {
        if(argc < 5){
            cout<<"USAGE for encoding queries: MAGIC <2> <pred_map> <vert_map> <queries_folder>"<<endl;
             cout<<"Example: ./dist/Release/GNU-Linux-x86/magic 2 data/lubm_160/encoded.ntpredicate_map.dic data/lubm_160/encoded.ntverts_map.dic  data/lubm_160/queries/"<<endl;
            exit(-1);
        }
        
        string pred_map = string(argv[2]);    
        string vert_map = string(argv[3]);    
        string queries = string(argv[4]);    
        string tool = "QueryLoadEncoder";
        const char* query_arg[4];
        query_arg[0] = tool.c_str();
        query_arg[1] = pred_map.c_str();
        query_arg[2] = vert_map.c_str();
        query_arg[3] = queries.c_str();
        main_q_encoder(4, query_arg);
    }

    if(oper == "SS_run")
    {  
        if(argc < 5){
            cout<<"USAGE for running queries: MAGIC SS_run <graph_path> <queries_folder> <out_res_folder> "<<endl;
            cout<<"Example: ./dist/Release/GNU-Linux-x86/magic 3 data/lubm_160/encoded.nt data/lubm_160/queries/encoded/ data/lubm_160/queries/magiq_res/"<<endl;
            exit(-1);
        }
          
        string graph_path = string(argv[2]);    
        string queries_folder = string(argv[3]);  
        string out_res_folder = string(argv[4]);   
        
        if(queries_folder.back() != '/') queries_folder += "/";
        if(out_res_folder.back() != '/') out_res_folder += "/";
        
        vector<string> query_path_vec;
        
        vector<string> files = vector<string>();
        getdir(queries_folder,files);
	cout<<"Will evaluate "<<(files.size())<<" query files."<<endl;
        sort(files.begin(), files.end()); 
        for(unsigned q = 0 ; q < files.size(); q++){
            query_path_vec.push_back(queries_folder + files[q]);
            cout<<"Adding query: "<<queries_folder + files[q]<<endl;
	}
        mq_exec_query_list(graph_path, query_path_vec, out_res_folder);
          
    }
    
    // makes a script for each query and a script to load the dataset once
    // and execute all the queries on CPU and GPU
    if(oper == "qrun2matlab")
    {  
        if(argc < 5) {
            cout << "USAGE: MAGIC qrun2matlab <graph> <queries> <output>" << endl;
            exit(-1);
        }
        
        string graph_path     = string(argv[2]);    
        string queries_folder = string(argv[3]);  
        string out_path       = string(argv[4]);   
        
        if(queries_folder.back() != '/') queries_folder += "/";
        if(out_path.back() != '/') out_path += "/";
        
        vector<string> query_path_vec;
        vector<string> files = vector<string>();
        getdir(queries_folder, files);
	cout<<"Will compile "<<(files.size())<<" query files."<<endl;
        sort(files.begin(), files.end()); 
        for(unsigned q = 0 ; q < files.size(); q++){
            query_path_vec.push_back(queries_folder + files[q]);
            cout<<"Adding query: "<< queries_folder + files[q]<< endl;
	}
        
        bool gpu = false;
        q_to_matlab_dataset(graph_path, query_path_vec, out_path, gpu);
        q_to_matlab_dataset(graph_path, query_path_vec, out_path, !gpu);
    }
    
    if(oper == "q2matlab")
    {  
        if(argc < 4) {
            cout << "USAGE: MAGIC q2matlab <queries> <output>" << endl;
            exit(-1);
        }
            
        string queries_folder = string(argv[2]);  
        string out_path       = string(argv[3]);   
        
        if(queries_folder.back() != '/') queries_folder += "/";
        if(out_path.back() != '/') out_path += "/";
        
        vector<string> query_path_vec;
        vector<string> files = vector<string>();
        getdir(queries_folder, files);
	cout<<"Will compile "<<(files.size())<<" query files."<<endl;
        sort(files.begin(), files.end()); 
        for(unsigned q = 0 ; q < files.size(); q++){
            query_path_vec.push_back(queries_folder + files[q]);
            cout<<"Adding query: "<< queries_folder + files[q]<< endl;
	}
        
        q_to_matlab_list(query_path_vec, out_path);
    }
    


    return 0;
}
