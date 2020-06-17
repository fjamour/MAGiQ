/* 
 * File:   mq_matlab_compiler.h
 * Author: fuad
 *
 * Created on February 1, 2018, 4:08 PM
 */

#ifndef MQ_MATLAB_COMPILER_H
#define	MQ_MATLAB_COMPILER_H

#include "mq_qplan.h"
#include "mq_simple_timer.h"
#include "mq_util.h"
#include "mq_matrix.h"
#include "mq_global_params.h"
#include <unordered_map>
#include <unordered_set>
#include <queue>

using namespace std;

/*
 *
 */
string
q_to_matlab
(
    const   Query& query,
    bool    load_layers = true
)
{
    // 1. construct predicate graphs needed by query (or make sure preloaded)
    // 2. make the matrix vector multiplication to get variable bindings
    // 3. make binding matrices (sorted to be ready for joins)
    // 4. construct the results by doing R_join and R_filter (if needed)
    // XXX make sure the cycle breaking filters appear as early as possible
    string  out = "";
    char    buff[2048];
    QPlan  qplan;
    QGraph qg;
    qplan.make_sqplan(query);
    qg = qplan.get_qgraph();
    
    auto walk     = qg.extract_dfs_dwalk_klbe();
    auto qpreds   = qg.get_predicates();
    auto eval_seq = qplan.get_eval_seq();
    
    
    
    // construct predicate graphs
    if(load_layers) {
        out += "t_layers = tic;\n";
        for(auto p : qpreds) {
            sprintf(
                buff,
                "[hI%d, hJ%d] = make_predicate_graph(I, J, V, %d);\n",
                p, p, p);
            out += buff;
        }
        out += "fprintf('  Time (query_layers): \%f\\n', toc(t_layers));\n";
        out += "\n";
    } else { //layers are global, and are loaded and defined as global by caller
        out += "global";
        for(auto p : qpreds) {
            sprintf( buff, " hI%d hJ%d", p, p);
            out += buff;
        }
        out += "\n\n";
    }
    // copy predicate graphs to GPU
    out += "t_all_wcopy = tic;\n";
    out += "use_gpu = ismember('use_gpu', varargin);\n";
    out += "if use_gpu\n";
    if(load_layers) out += "    reset(GPU)\n";
    out += "    t_gpu_copy = tic;\n";
    for(auto p: qpreds) {
        sprintf(
            buff,
            "    I%d = gpuArray(hI%d); J%d = gpuArray(hJ%d);\n",
            p, p, p, p);
        out += buff;
    }
    out += "    fprintf('  Time (gpu_layer_copy): \%f\\n', toc(t_gpu_copy));\n";
    out += "else\n";
    out += "    t_copy = tic;\n";
    for(auto p: qpreds) {
        sprintf(buff, "    I%d = hI%d; J%d = hJ%d;\n", p, p, p, p);
        out += buff;
    }
    out += "    fprintf('  Time (local_layer_copy): \%f\\n', toc(t_copy));\n";
    out += "end\n";
    out += "\n";
    // find edge types: forward and backward
    // true: forwards, false: backward
    vector<bool> etype_vec(walk.size(), true);
    for(int i = 1; i < walk.size(); ++i) {
        auto end = walk.begin() + i + 1;
        auto e   = walk[i];
        auto e_r = make_pair(e.second, e.first);
        if(find(walk.begin(), end, e_r) != end)
            etype_vec[i] = false;
    }
    // do first node
    out += "t_all = tic;\n";
    out += "t_mxv = tic;\n";
    bool v1_is_c, v2_is_c, e_dir_match;
    dnode_t c1, c2;
    edge_label_t edge_scalar;
    qedge_t e;
    e = walk[0];
    qg.get_node_labels(v1_is_c, v2_is_c, c1, c2, e);
    qg.get_edge_label(edge_scalar, e_dir_match, e);
    if(v1_is_c) {
        sprintf(buff, "fv%d = %d;\n", e.first, c1+1);
    } else {
        if(e_dir_match) {
            sprintf(buff, "fv%d = I%d;\n", e.first, edge_scalar);
        } else {
            sprintf(buff, "fv%d = J%d;\n", e.first, edge_scalar);
        }
    }
    out += buff;
    // do all mat vec mul
    auto eqp_vec = qg.get_eq_var_pairs();
    map<qnode_t, qnode_t> sn_map;
    for(auto p : eqp_vec) {sn_map[p.first] = p.second;}
    for(int i = 0; i < walk.size(); ++i) {
        e = walk[i];
        qg.get_node_labels(v1_is_c, v2_is_c, c1, c2, e);
        qg.get_edge_label(edge_scalar, e_dir_match, e);
        string trans = "";
        if(e_dir_match) trans = ", 'trans'";
        if(etype_vec[i]) { // forward edge
            if(v2_is_c) {
                sprintf(buff, "fv%d = %d;\n", e.second, c2+1);
            } else {
                sprintf(
                    buff,
                    "fv%d = mxv(I%d, J%d, fv%d%s);\n",
                    e.second, edge_scalar, edge_scalar, e.first, trans.c_str());
            }
            out += buff;
            if(sn_map.find(e.second) != sn_map.end()) { // cycle filter
                sprintf(buff, "fv%d = fv%d;\n", sn_map[e.second], e.second);
                out += buff;
            }
        } else { // backward edge
            sprintf(buff,
                "fv%d = vandv(fv%d, mxv(I%d, J%d, fv%d%s));\n",
                e.second, e.second,
                edge_scalar, edge_scalar,
                e.first, trans.c_str());
            out += buff;
        }
    }
    out += "fprintf('  Time (var_bindings): \%f\\n', toc(t_mxv));\n";
    out += "\n";
    // construct sorted binding matrices
    out += "t_bm = tic;\n";
    for(auto e: eval_seq) {
        qg.get_edge_label(edge_scalar, e_dir_match, e);
        string trans = "";
        if(!e_dir_match) trans = ", 'trans'";
        sprintf(
            buff,
            "sB%d%d = make_bm(I%d, J%d, fv%d, fv%d%s);\n",
            e.first, e.second, edge_scalar, edge_scalar,
            e.first, e.second, trans.c_str());
        out += buff;
    }
    out += "fprintf('  Time (binding_mats): \%f\\n', toc(t_bm));\n";
    out += "\n";
    // construct result set R with joins and filters
    out += "t_R = tic;\n";
    map<int, int> v2c_map;
    int last_c = 1;
    e = eval_seq[0];
    sprintf(buff, "R = sB%d%d;\n", e.first, e.second);
    v2c_map[e.first] = 1;
    v2c_map[e.second] = 2;
    last_c = 2;
    out += buff;
    for(int i = 1; i < eval_seq.size(); ++i) {
        e = eval_seq[i];
        bool is_filter = sn_map.find(e.second) != sn_map.end();
        if(is_filter) {
            sprintf(
                buff,
                "R = R_filter(R, sB%d%d, %d, %d);\n",
                e.first, e.second,
                v2c_map[e.first], v2c_map[sn_map[e.second]]);
        } else {
            sprintf(
                buff,
                "R = R_join(R, sB%d%d, %d);\n",
                e.first, e.second, v2c_map[e.first]);
            v2c_map[e.second] = ++last_c;
        }
        out += buff;
    }
    out += "fprintf('  Time (join_res): \%f\\n', toc(t_R));\n";
    out += "fprintf('  Time (total_query): \%f\\n', toc(t_all));\n";
    out += "fprintf('  Time (total_query_wcopy): \%f\\n', toc(t_all_wcopy));\n";
    out += "fprintf('  Number of results: \%d\\n', size(R, 1));\n";
    out += "\n";
    return out;
}

/*
 *
 */
void
q_to_matlab_list
(
    vector<string>  query_path_vec,
    string          output_path
)
{
    vector<Query>  raw_q_vec;
    vector<string> raw_qstr_vec;
    for(auto qpath : query_path_vec) {
        Query q;
        string qstr;
        read_query_file(q, qstr, qpath);
        raw_q_vec.push_back(q);
        raw_qstr_vec.push_back(qstr);
        //printf("%s", qstr.c_str());
    }
    for(int i = 0; i < raw_q_vec.size(); ++i) {
        Query query = raw_q_vec[i];
        string qstr = raw_qstr_vec[i];
        string qpath = query_path_vec[i];
        string qname = get_qname(qpath);
        printf("#######################################################\n");
        printf("Query: [%s]\n", qpath.c_str());
        QPlan qplan;
        qplan.make_sqplan(query);
        //printf("========\n");
        //qplan.print_qgraph();
        //printf("========\n");
        //qplan.print_qwalk();
        //printf("=========================================================\n");
        //qplan.print_qplan();
        qplan.make_sqplan_csc_friendly(query);
        //printf("=========================================================\n");
        //qplan.print_qplan();
        //printf("=========================================================\n");
        string mfile = q_to_matlab(query);
        char mfun_str[100000];
        sprintf(
            mfun_str,
            "function R = mq_%s(I, J, V, N, GPU, varargin)\n%send\n",
            qname.c_str(), mfile.c_str());
        string out_path = output_path + "mq_" + qname + string(".m");
        FILE* fout = fopen(out_path.c_str(), "w");
        fprintf(fout, "%s", mfun_str);
        fclose(fout);
        //printf("Matlab program:\n%s", mfun_str);
        //printf("=========================================================\n");
    }
}

/*
 * Given a dataset and its queries, produce matlab function for each query,
 * and a matlab script to load graph run queries.
 */
void
q_to_matlab_dataset
(
    string          graph_path,
    vector<string>  query_path_vec,
    string          out_path,
    bool            gpu_backend = true
)
{
    char buff[100000];
    vector<Query>  raw_q_vec;
    vector<string> raw_qstr_vec;
    for(auto qpath : query_path_vec) {
        Query q;
        string qstr;
        read_query_file(q, qstr, qpath);
        raw_q_vec.push_back(q);
        raw_qstr_vec.push_back(qstr);
        //printf("%s", qstr.c_str());
    }
    vector<string> qfun_name_vec;
    // get dataset prefix
    string gname = get_dataset_prefix(graph_path);
    for(int i = 0; i < raw_q_vec.size(); ++i) {
        Query query = raw_q_vec[i];
        string qstr = raw_qstr_vec[i];
        string qpath = query_path_vec[i];
        string qname = get_qname(qpath);
        string out_qname  = gname + "_" + qname;
        string out_qpath  = out_path + out_qname + ".m";
        
        bool load_layers = false;
        string mfile = q_to_matlab(query, load_layers);
        qfun_name_vec.push_back(out_qname);
        char mfun_str[100000];
        sprintf(
            mfun_str,
            "function R = %s(varargin)\n%send\n",
            out_qname.c_str(), mfile.c_str());
        FILE* fout = fopen(out_qpath.c_str(), "w");
        fprintf(fout, "%s", mfun_str);
        fclose(fout);
        //printf("#######################################################\n");
        //printf("Query: [%s]\n", qpath.c_str());
        //printf("Matlab program:\n%s", mfun_str);
        //printf("=========================================================\n");
    }
    // make the run script
    // 1. get all the needed predicates
    set<edge_label_t> all_p_set;
    for(int i = 0; i < raw_q_vec.size(); ++i) {
        Query query = raw_q_vec[i];
        QPlan qplan;
        qplan.make_sqplan(query);
        QGraph qg = qplan.get_qgraph();
        auto pred_vec = qg.get_predicates();
        all_p_set.insert(pred_vec.begin(), pred_vec.end());
    }
    string run_script = "";
    run_script += "global";
    for(auto p : all_p_set) {
        sprintf(buff, " hI%d hJ%d", p, p);
        run_script += buff;
    }
    run_script += "\n\n";
    // 2. read the graph
    run_script += "if ~graph_is_loaded\n";
    run_script += "    [I J V N] = read_rdf_graph(PATH_TO_GRAPH, 'uint32', 'sort_I');\n\n";
    // 3. make the predicate graphs
    run_script += "    t_layers = tic;\n";
    for(auto p : all_p_set) {
        sprintf(
            buff,
            "    [hI%d, hJ%d] = make_predicate_graph(I, J, V, %d);\n",
            p, p, p);
            run_script += buff;
    }
    run_script += "    fprintf('Time (all_layers): \%f\\n', toc(t_layers));\n";
    run_script += "    graph_is_loaded = true;\n";
    run_script += "end\n\n";
    // 4. run the queries
    if(gpu_backend) {
        run_script += "GPU = gpuDevice;\n";
        run_script += "reset(GPU);\n\n";
    }
    for(auto q : qfun_name_vec) {
        string use_gpu = "";
        if(gpu_backend) use_gpu = "'use_gpu'";
        sprintf(buff, "fprintf('Query: %s\\n');\n", q.c_str());
        run_script += buff;
        sprintf(buff, "R = %s(%s);\n", q.c_str(), use_gpu.c_str());
        run_script += buff;
        run_script += "fprintf('=========================\\n');\n";
    }
    run_script += "\n";
    
    string backend_label;
    if(gpu_backend) backend_label = "_gpu_run.m";
    else            backend_label = "_cpu_run.m";
    string out_script_path = out_path + gname + backend_label;
    FILE* fout = fopen(out_script_path.c_str(), "w");
    fprintf(fout, "%s", run_script.c_str());
    fclose(fout);
}




#endif	/* MQ_MATLAB_COMPILER_H */

