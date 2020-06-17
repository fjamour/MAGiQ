/*
 File Header TODO have a nice header 
 */

#include <iostream>
#include "mq_simple_logger.h"

#include "mq_query_exec.h"
#include "mq_interface.h"
#include "DataEncoder.h"
#include "QueryEncoder.h"

#include "mq_matlab_compiler.h"
extern "C" 
{
    #include "GraphBLAS.h"
}

void run_data_encoder(){
    string tool = "DataEncoder";
    string data_dir_path = "data/lubm_160/";
    string output = "data/lubm_160/encoded.nt"; //output file name
    string remove_double_edges = "0";
    
    const char* data_arg[4];
    data_arg[0] = tool.c_str();
    data_arg[1] = data_dir_path.c_str();
    data_arg[2] = output.c_str();
    data_arg[3] = remove_double_edges.c_str();
    
    main_d_encoder(4, data_arg);
}


void run_query_encoder(){
    string tool = "QueryLoadEncoder";
    string pred_map = "data/lubm_160/encoded.ntpredicate_map.dic";
    string vert_map = "data/lubm_160/encoded.ntverts_map.dic";
    string queries = "data/lubm_160/queries/";
    
    const char* query_arg[4];
    query_arg[0] = tool.c_str();
    query_arg[1] = pred_map.c_str();
    query_arg[2] = vert_map.c_str();
    query_arg[3] = queries.c_str();
    main_q_encoder(4, query_arg);
}

void run_magiq(){
    string graph_path;
    vector<string> query_path_vec;
    
//    graph_path = "./data/fig1.nt";
//    query_path_vec.push_back(string("./queries/dev_qg2.q"));
//    query_path_vec.push_back(string("./queries/dev_qg3.q"));
    
//    graph_path = "./data/lubm320/encoded.nt";
//    query_path_vec.push_back(string("./data/lubm320/queries/l1.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm320/queries/l2.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm320/queries/l3.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm320/queries/l4.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm320/queries/l5.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm320/queries/l6.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm320/queries/l7.q_encoded"));

//    graph_path = "./data/lubm2560/encoded.nt";
//    query_path_vec.push_back(string("./data/lubm2560/queries/l1.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm2560/queries/l2.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm2560/queries/l3.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm2560/queries/l4.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm2560/queries/l5.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm2560/queries/l6.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm2560/queries/l7.q_encoded"));
    
//    graph_path = "./data/watdiv100M/encoded.nt";
//    query_path_vec.push_back(string("./data/watdiv100M/queries/S1.in_encoded"));
//    query_path_vec.push_back(string("./data/watdiv100M/queries/SF1.in_encoded"));
//    query_path_vec.push_back(string("./data/watdiv100M/queries/L1.in_encoded"));
//    query_path_vec.push_back(string("./data/watdiv100M/queries/C1.in_encoded")); 
    
//    graph_path = "./data/yago/encoded.nt";
//    query_path_vec.push_back(string("./data/yago/queries/1.q_encoded"));
//    query_path_vec.push_back(string("./data/yago/queries/2.q_encoded"));
//    query_path_vec.push_back(string("./data/yago/queries/3.q_encoded"));
//    query_path_vec.push_back(string("./data/yago/queries/4.q_encoded"));
    
//    graph_path = "./data/watdiv1B/encoded.nt";
//    query_path_vec.push_back(string("./data/watdiv1B/queries/S1.in_encoded"));
//    query_path_vec.push_back(string("./data/watdiv1B/queries/SF1.in_encoded"));
//    query_path_vec.push_back(string("./data/watdiv1B/queries/C1.in_encoded"));
//    query_path_vec.push_back(string("./data/watdiv1B/queries/C6.in_encoded"));
    
//    graph_path = "./data/lubm10240/encoded.nt";
//    query_path_vec.push_back(string("./data/lubm10240/queries/encoded/l1.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm10240/queries/encoded/l2.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm10240/queries/encoded/l3.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm10240/queries/encoded/l4.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm10240/queries/encoded/l5.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm10240/queries/encoded/l6.q_encoded"));
//    query_path_vec.push_back(string("./data/lubm10240/queries/encoded/l7.q_encoded"));
    
//    graph_path = "./data/bio2rdf/encoded.nt";
//    query_path_vec.push_back(string("./data/bio2rdf/queries/b1.q_encoded"));
//    query_path_vec.push_back(string("./data/bio2rdf/queries/b2.q_encoded"));
//    query_path_vec.push_back(string("./data/bio2rdf/queries/b3.q_encoded"));
//    query_path_vec.push_back(string("./data/bio2rdf/queries/b4.q_encoded"));
//    query_path_vec.push_back(string("./data/bio2rdf/queries/b5.q_encoded"));

//    string qres_dir_path = "./scratch/";
//    mq_exec_query_list(graph_path, query_path_vec, qres_dir_path);
    
//    string matlabq_path = "./scratch/";
//    q_to_matlab_list(query_path_vec, matlabq_path);
    
//    string out_path = "./matlab/query_scripts/bio2rdf/";
//    
//    bool gpu = false;
//    q_to_matlab_dataset(graph_path, query_path_vec, out_path, gpu);
//    q_to_matlab_dataset(graph_path, query_path_vec, out_path, !gpu);
}




#include "TurtleParser.hpp"
void adhoc_yago_fix()
{
    string in_path = "/home/fuad/scratch/yago_fix/yago.nt";
    ifstream fin(in_path.c_str());
    TurtleParser parser(fin);
   
    string subject,predicate,object,objectSubType;
    while (true) {
        string subject,predicate,object,objectSubType;
        Type::ID objectType;
        try {
            if(!parser.parse(subject,predicate,object,objectType,objectSubType)) {
                break;
            }
        } catch (const TurtleParser::Exception& e) {
            cerr << e.message << endl;
            while (fin.get()!='\n') ;
            continue;
        }
        if(objectType == 1) {
            object.erase(remove(object.begin(), object.end(), '\"' ),  object.end());
            object.erase(remove(object.begin(), object.end(), ' ' ),  object.end());
            object.erase(remove(object.begin(), object.end(), '\t' ),  object.end());
        }
        subject.erase(remove(subject.begin(), subject.end(), ' ' ),  subject.end());
        subject.erase(remove(subject.begin(), subject.end(), '\t' ),  subject.end());
        predicate.erase(remove(predicate.begin(), predicate.end(), ' ' ),  predicate.end());
        predicate.erase(remove(predicate.begin(), predicate.end(), '\t' ),  predicate.end());
        if(objectType == 1){
            printf("<%s> <%s> \"%s\" .\n", subject.c_str(), predicate.c_str(), object.c_str());
        } else {
            printf("<%s> <%s> <%s> .\n", subject.c_str(), predicate.c_str(), object.c_str());
        }
    }
    
}





int main(int argc, char**argv) {
    setbuf(stdout, NULL); // make printf display immediately
    
    //adhoc_yago_fix();
    interface_main(argc, argv);
    
//    run_data_encoder();
//    run_query_encoder();
//    
//    run_magiq();

    return 0;
}
