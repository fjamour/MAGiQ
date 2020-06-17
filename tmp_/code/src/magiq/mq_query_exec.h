/* 
 * File:   mq_query_exec.h
 * Author: fuad
 *
 * Created on December 2, 2017, 5:59 PM
 */

#ifndef MQ_QUERY_EXEC_H
#define	MQ_QUERY_EXEC_H

#include "mq_qplan.h"
#include "mq_simple_timer.h"
#include "mq_util.h"
#include "mq_matrix.h"
#include "mq_global_params.h"
#include <unordered_map>
#include <unordered_set>
#include <queue>


#include <omp.h>



void
mq_exec_query
(
    vector<vector<int> >&   res_vec,
    QPlan                   qplan,
    MQGraph                 G
);


/*
 *
 */
void
write_qresults
(
    vector<vector<int> >&   qresults,
    QGraph                  qgraph,
    string                  out_path
)
{
    FILE* fout = fopen(out_path.c_str(), "w");
    if(!fout) {
        LOGGER.log(string("Can't open file: ") + out_path, true, false);
        exit(1);
    }
    auto proj_vars = qgraph.get_proj_vars();
    for(auto& rec : qresults) {
        for(auto vp : proj_vars) {
            fprintf(fout, "<%d> ", rec[vp.second]);
        }
        fprintf(fout, "\n");
   }
    fclose(fout);
}

/*
 *
 */
void
mq_exec_query_list
(
    string          graph_path,
    vector<string>  query_path_vec,
    string          qres_dir_path
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
    timer tm;
    MQGraph G;
    tm.start();
    printf("#######################################################\n");
    printf("Reading graph [%s]...\n", graph_path.c_str());
    G.load_graph(graph_path);
    tm.stop();
    printf("    Graph load (total):   %f\n", tm.interval());
    printf("#######################################################\n\n");
    
    for(int i = 0; i < raw_q_vec.size(); ++i) {
        Query query = raw_q_vec[i];
        string qstr = raw_qstr_vec[i];
        string qpath = query_path_vec[i];
        printf("#######################################################\n");
        printf("Query: [%s]\n", qpath.c_str());
        QPlan qplan;
        qplan.make_sqplan(query);
        printf("========\n");
        qplan.print_qgraph();
        printf("========\n");
        qplan.print_qwalk();
        printf("=========================================================\n");
        qplan.print_qplan();
        qplan.make_sqplan_csc_friendly(query);
        printf("=========================================================\n");
        qplan.print_qplan();
        printf("=========================================================\n");
        vector<vector<int> > res_vec;
        tm.start();
        mq_exec_query(res_vec, qplan, G);
        tm.stop();
        printf("Query time (total):     %f\n", tm.interval());
        printf("#######################################################\n\n");
        
        // Write results to disk
        if(WRITE_RES_TO_DISK) {
            string qres_path = qres_dir_path + get_qname(qpath) + string(".mq_res");
            write_qresults(res_vec, qplan.get_qgraph(), qres_path);
        }
        
    }
    G.destroy();
}

/*
 *
 */
void
mq_enum_res_col_major
(
    vector<vector<int> >&  qresults,
    mat_store_t            mat_store,
    vector<qedge_t>        eval_seq,
    QPlan                  qplan
)
{
    //======================================================================
    //=experimental dag stuff:==============================================
//    for(auto p : eval_seq) {
//        printf("%d, %d \n", p.first, p.second);
//    }
//    // reorder eval_seq
//    vector<qedge_t> eval_tree;
//    eval_tree.push_back(eval_seq[0]);
//    for(int i = 0; i < eval_tree.size(); ++i) {
//        qedge_t nd = eval_tree[i];
//        for(int j = i; j < eval_seq.size(); ++j) {
//            if(nd.first == eval_seq[j].second)
//                eval_tree.push_back(eval_seq[j]);
//        }
//    }
    //======================================================================
   
    //======================================================================
    //=experimental cycle cleaning stuff:===================================
    typedef vector<pair<qnode_t, qnode_t> > eq_pairs_t ;
    eq_pairs_t eq_pairs;
    eq_pairs = qplan.get_eq_var_pairs();
    auto valid_eqp = [](vector<int>& rec, eq_pairs_t& eq_pairs)
    {
        for(auto& p:eq_pairs) {
            if(rec[p.first] == -1 || rec[p.second] == -1)
                return true;
            if(rec[p.first] != rec[p.second])
                return false;
        }
        return true;
    };
    //======================================================================
    
    
    // Prepare qresults
    int num_vars = qplan.get_num_var();
    qresults.clear();
    //qresults.reserve(10000);
    
    // Intermediate arrays
    size_t num_nodes = mat_store[eval_seq[0]].get_nnodes(); // Max possible
    GrB_Index   *I = (GrB_Index*)malloc(num_nodes * sizeof(GrB_Index));
    GrB_Index   *J = (GrB_Index*)malloc(num_nodes * sizeof(GrB_Index));
    int8_t      *X = (int8_t*)malloc(num_nodes * sizeof(int8_t));   
    
    
    // Fill records from the first matrix
    timer ret_tm, vec_tm;
    qedge_t     e = eval_seq[0];
    MQMatrix    m = mat_store[e];
    size_t      nnz;
    int ii, jj;
    ii = e.first;
    jj = e.second;
    ret_tm.start();
    nnz = m.get_all_nz(I, J, X);
    ret_tm.stop();
    vec_tm.start();
    for(int i = 0; i < nnz; ++i) {
        vector<int> rec(num_vars, -1);
        rec[ii] = I[i];
        rec[jj] = J[i];
        qresults.push_back(rec);
    }
    vec_tm.stop();
    printf(
        "Initial fill from m(%d,%d): num_nnz[%lu] get_nz[%f] vec_fill[%f]\n",
        ii, jj, nnz,
        ret_tm.interval(), vec_tm.interval());
    // Fill next var from next mat
    timer cp_tm;
    timer fill_tm;
    for(int e_i = 1; e_i < eval_seq.size(); ++e_i) {
        e = eval_seq[e_i];
        m = mat_store[e];
        int var_src_i = e.second;
        int var_dst_i = e.first;
        vector<vector<int> > updated_qresults;
        updated_qresults.reserve(1000000);
        fill_tm.start();
        //======================================================================
        // EXP parallel
        int num_thread = 20;
        int rec_per_thread = qresults.size()/num_thread;
        omp_set_num_threads(num_thread);
        vector<vector<int> > new_qresults;
#pragma omp parallel for
        for(int block = 0; block < num_thread; ++block) {
            int t_start, t_end;
            t_start = rec_per_thread*block;
            t_end   = t_start + rec_per_thread;
            if(block == num_thread-1)
                t_end = qresults.size();
            GrB_Index   *I = (GrB_Index*)malloc(num_nodes * sizeof(GrB_Index));
            int8_t      *X = (int8_t*)malloc(num_nodes * sizeof(int8_t)); 
            vector<vector<int> > updated_qresults;
            updated_qresults.reserve(1000000);
            for(int rec_i = t_start; rec_i < t_end; ++rec_i) {
                vector<int> rec = qresults[rec_i];
                size_t nnz = m.get_col(I, X, rec[var_src_i]);
                for(int i = 0; i < nnz; ++i) {
                    rec[var_dst_i] = I[i];
                    if(valid_eqp(rec, eq_pairs)) {
                        updated_qresults.push_back(rec);
                    }
                }
            }
            free(I);
            free(X);
            #pragma omp critical
            new_qresults.insert(
                    new_qresults.end(),
                    updated_qresults.begin(),
                    updated_qresults.end());
        }
        //======================================================================
        
//        for(auto& rec : qresults) {
//            nnz = m.get_col(I, X, rec[var_src_i]);
//            for(int i = 0; i < nnz; ++i) {
//                rec[var_dst_i] = I[i];
//                if(valid_eqp(rec, eq_pairs))
//                    updated_qresults.push_back(rec);
//            }
//        }
        fill_tm.stop();
        cp_tm.start();
//        qresults = updated_qresults;
        qresults = new_qresults;
        cp_tm.stop();
        printf( "  Fill from m(%d,%d): num_rec[%10lu] get_nz[%f] vec_cp[%f]\n",
                e.first, e.second,
                qresults.size(),
                fill_tm.interval(), cp_tm.interval());
    }
//    for(auto p : qresults) {
//        for(int i = 0; i < num_vars; ++i)
//            printf("%d ", p[i]);
//        printf("\n");
//    }    
    free(I);
    free(J);
    free(X);
}

/*
 *
 */
void
EXP_mq_enum_res_col_major
(
    vector<vector<int> >&  qresults,
    mat_store_t            mat_store,
    vector<qedge_t>        eval_seq,
    QPlan                  qplan
)
{
    //======================================================================
    //=experimental dag stuff:==============================================
//    for(auto p : eval_seq) {
//        printf("%d, %d \n", p.first, p.second);
//    }
//    // reorder eval_seq
//    vector<qedge_t> eval_tree;
//    eval_tree.push_back(eval_seq[0]);
//    for(int i = 0; i < eval_tree.size(); ++i) {
//        qedge_t nd = eval_tree[i];
//        for(int j = i; j < eval_seq.size(); ++j) {
//            if(nd.first == eval_seq[j].second)
//                eval_tree.push_back(eval_seq[j]);
//        }
//    }
    //======================================================================
   
    
    // rename equal vars to actual names (in eval_seq and mat_store)
    typedef vector<pair<qnode_t, qnode_t> > eq_pairs_t ;
    eq_pairs_t eq_pairs;
    eq_pairs = qplan.get_eq_var_pairs();
    map<qnode_t, qnode_t> eq_var_map;
    for(auto ep : eq_pairs) {
        eq_var_map.insert(ep);
    }
    for(auto& e : eval_seq) {
        if(eq_var_map.find(e.first) != eq_var_map.end()) {
            auto old_e = e;
            e = make_pair(e.second, eq_var_map[e.first]);
            MQMatrix m;
            m.deep_copy(mat_store[old_e], true);
            mat_store[e] = m;
            mat_store.erase(old_e);
        }
    }
    // reorder eval_seq
    vector<qedge_t> eval_tree;
    eval_tree.push_back(eval_seq[0]);
    eval_seq.erase(eval_seq.begin());
    for(int i = 0; i < eval_tree.size(); ++i) {
        qedge_t nd = eval_tree[i];
        for(int j = 0; j < eval_seq.size(); ++j) {
            if(nd.first == eval_seq[j].second) {
                eval_tree.push_back(eval_seq[j]);
                eval_seq.erase(eval_seq.begin() + j);
                --j;
            }
        }
    }
    eval_seq = eval_tree;
    for(auto p : eval_seq) {
        printf("%d, %d \n", p.first, p.second);
    }
    eval_seq.clear();
    eval_seq.push_back(make_pair(4, 0));
    eval_seq.push_back(make_pair(5, 4));
    eval_seq.push_back(make_pair(3, 4));
    eval_seq.push_back(make_pair(5, 3));
    eval_seq.push_back(make_pair(2, 3));
    eval_seq.push_back(make_pair(1, 5));
    
 
    
    // Prepare qresults
    int num_vars = qplan.get_num_actual_vars();
    qresults.clear();
    //qresults.reserve(10000);
    
    // Intermediate arrays
    size_t num_nodes = mat_store[eval_seq[0]].get_nnodes(); // Max possible
    GrB_Index   *I = (GrB_Index*)malloc(num_nodes * sizeof(GrB_Index));
    GrB_Index   *J = (GrB_Index*)malloc(num_nodes * sizeof(GrB_Index));
    int8_t      *X = (int8_t*)malloc(num_nodes * sizeof(int8_t));   
    
        printf("Mat store content:\n");
        for(auto& e : eval_seq) {
            printf
            (
                "  m(%d, %d): %9lu [%9lu, %9lu]\n",
                e.first, e.second,
                mat_store[e].get_nnz(),
                mat_store[e].get_nnz(MQ_FROM_RROWS),
                mat_store[e].get_nnz(MQ_FROM_RCOLS)
            );
        }
    // Fill records from the first matrix
    timer ret_tm, vec_tm;
    qedge_t     e = eval_seq[0];
    MQMatrix    m = mat_store[e];
    size_t      nnz;
    int ii, jj;
    ii = e.first;
    jj = e.second;
    ret_tm.start();
    nnz = m.get_all_nz(I, J, X);
    ret_tm.stop();
    vec_tm.start();
    for(int i = 0; i < nnz; ++i) {
        vector<int> rec(num_vars, -1);
        rec[ii] = I[i];
        rec[jj] = J[i];
        qresults.push_back(rec);
    }
    vec_tm.stop();
    printf(
        "Initial fill from m(%d,%d): num_nnz[%lu] get_nz[%f] vec_fill[%f]\n",
        ii, jj, nnz,
        ret_tm.interval(), vec_tm.interval());
    // Fill next var from next mat
    timer cp_tm;
    timer fill_tm;
    for(int e_i = 1; e_i < eval_seq.size(); ++e_i) {
        e = eval_seq[e_i];
        m = mat_store[e];
        int var_src_i = e.second;
        int var_dst_i = e.first;
        vector<vector<int> > updated_qresults;
        updated_qresults.reserve(1000000);
        fill_tm.start();
        if
        (   e_i+1 < eval_seq.size()
         && qresults.size() > 0
         && var_dst_i == eval_seq[e_i+1].second
         && qresults[0][eval_seq[e_i+1].first] != -1
        )
        {
            for(auto& rec : qresults) {
                nnz = m.get_col(I, X, rec[var_src_i]);
                for(int i = 0; i < nnz; ++i) {
                    qedge_t ee;
                    ee.first = qresults[0][eval_seq[e_i+1].first];
                    ee.second = I[i];
                    if(mat_store[eval_seq[e_i+1]].get_element(ee) != 0) {
                        rec[var_dst_i] = I[i];
                        updated_qresults.push_back(rec);
                    }
                }
            }
            e_i++;
        } else {
            for(auto& rec : qresults) {
                nnz = m.get_col(I, X, rec[var_src_i]);
                for(int i = 0; i < nnz; ++i) {
                    rec[var_dst_i] = I[i];
                    updated_qresults.push_back(rec);
                }
            }
        }
        fill_tm.stop();
        cp_tm.start();
        qresults = updated_qresults;
        cp_tm.stop();
        printf( "  Fill from m(%d,%d): num_rec[%10lu] get_nz[%f] vec_cp[%f]\n",
                e.first, e.second,
                qresults.size(),
                fill_tm.interval(), cp_tm.interval());
    }
//    for(auto p : qresults) {
//        for(int i = 0; i < num_vars; ++i)
//            printf("%d ", p[i]);
//        printf("\n");
//    }    
    free(I);
    free(J);
    free(X);
}

/*
 *
 */
void
mq_exec_query
(
    vector<vector<int> >&   qresults,
    QPlan                   qplan,
    MQGraph                 G
)
{
    // NOTES
    // 1. ZERO aware i.e., any ZERO matrix results in immediate termination.
    // 2. 'Matrix Algebra' functionality fully abstracted
    // 3. Graph is a Matrix
    // 4. leave a comment to remember that many matrices are created and destroyed
    //    (probably unnecessarily )
    GrB_init (GrB_NONBLOCKING);
    
    mat_store_t mat_store;
    
    // Fill mat_store with empty matrices
    auto op_vec = qplan.get_op_vec();
    for(auto op : op_vec) {
        if(op.dst_desc == MQ_CONSTRUCT_MAT) {
            mat_store[op.dst_mat_e] = MQMatrix();
            mat_store[op.dst_mat_e].init_size(G.get_V());
        }
    }
    
    // Execute the query plan
    timer op_tm, grp_tm, ste_tm, rdc_tm;
    timer m1_tm, m2_tm, mul_tm, mxm_tm, trn_tm, zrc_tm;
    timer muls_tm;
    muls_tm.start();
    char str_[1024];
    bool empty_res = false;
    for(auto op : op_vec) {
        if(empty_res) {
            break;
        }
        str_[0] = 0;
        op_tm.start();
        qedge_t e = op.dst_mat_e;
        switch(op.type) {
            case MQ_TWOCB:
            {
                bool T = op.G_transpose;
                dedge_t ee = make_pair(op.c1, op.c2);
                bool G_has_edge;
                grp_tm.start();
                G_has_edge =
                        (G.get_mat().get_element(ee, T) == op.edge_label);
                grp_tm.stop();
                ste_tm.start();
                if(G_has_edge) {
                    mat_store[e].set_element(ee, 1);
                }
                ste_tm.stop();
                sprintf(str_, "G_access[%f] M_set[%f]",
                        grp_tm.interval(),
                        ste_tm.interval());
            }
            break;
            case MQ_TWOCM:
            {
                bool valid_rc;
                rdc_tm.start();
                if(   op.m1.mat_transform == MQ_TRANSPOSE
                   || op.m1.mat_transform == MQ_TRANS_DIAG)
                {
                    valid_rc = mat_store[op.m1.mat_e].reduce_col_lor(op.mc);
                } else {
                    valid_rc = mat_store[op.m1.mat_e].reduce_row_lor(op.mc);
                }
                rdc_tm.stop();
                grp_tm.start();
                if(valid_rc) {
                    bool T = op.G_transpose;
                    dedge_t ee = make_pair(op.c1, op.c2);
                    bool G_has_edge;
                    G_has_edge =
                            (G.get_mat().get_element(ee, T) == op.edge_label);
                    if(G_has_edge) {
                        mat_store[e].set_element(ee, 1);
                    }
                }
                grp_tm.stop();
                sprintf(str_, "Reduction[%f] G_access+M_set[%f]",
                        rdc_tm.interval(),
                        grp_tm.interval());
            }
            break;
            case MQ_MUL:
            {
                MQMatrix m1, m2;
                m1.meta = op.m1;
                m2.meta = op.m2;
                m1_tm.start();
                if(op.m1.mat_content == MQ_GRAPH)
                {
                    bool T = op.m1.mat_transform == MQ_TRANSPOSE;
                    m1 = G.get_mat(T);
                    
                    // layered
//                    m1 = G.get_mat(T, op.m2.scalar);
//                    op.m2.scalar = 1;
                } else if(   op.m1.mat_content == MQ_STORED
                          && op.m1.mat_transform == MQ_NOTRANSFORM)
                {
                    m1.deep_copy(mat_store[op.m1.mat_e]);
                } else
                {
                    m1.construct_from_meta(op.m1, mat_store, G.get_V());
                }
                m1_tm.stop();
                m2_tm.start();
                if(op.m2.mat_content == MQ_GRAPH)
                {
                    bool T = op.m2.mat_transform == MQ_TRANSPOSE;
                    m2 = G.get_mat(T);
                } else if(   op.m2.mat_content == MQ_STORED
                          && op.m2.mat_transform == MQ_NOTRANSFORM)
                {
                    m2.deep_copy(mat_store[op.m2.mat_e]);
                } else
                {
                    m2.construct_from_meta(op.m2, mat_store, G.get_V());
                }
                m2_tm.stop();
                mul_tm.start();
                switch(op.dst_desc) {
                    case MQ_CONSTRUCT_MAT:
                        mat_store[e].construct_from_mul
                        (mxm_tm, trn_tm, zrc_tm, m1, m2, op.trans_res);
                        break;
                    case MQ_UPDATE_MAT:
                        mat_store[e].update_from_mul
                        (mxm_tm, trn_tm, zrc_tm, m1, m2, op.trans_res);
                        break;
                }
                if(mat_store[e].get_nnz() == 0) {
                    empty_res = true;
                }
                mul_tm.stop();                
                if(op.m1.mat_content != MQ_GRAPH) m1.destroy();
                if(op.m2.mat_content != MQ_GRAPH) m2.destroy();
                
                sprintf(str_,
                        "m1[%f] m2[%f] mul[%f] SS_mxm[%f] T[%f] zero_clean[%f]",
                        m1_tm.interval(), m2_tm.interval(),
                        mul_tm.interval(),
                        mxm_tm.interval(),
                        trn_tm.interval(),
                        zrc_tm.interval());
            }
            break;
        }
        op_tm.stop();
        printf("[%f] %s\n", op_tm.interval(), op.stringify().c_str());
        if(ENABLE_PROFILE) {
            printf
            (
                "    MUL(%d, %d): %lu [%lu, %lu]\n",
                op.dst_mat_e.first, op.dst_mat_e.second,
                mat_store[e].get_nnz(),
                mat_store[e].get_nnz(MQ_FROM_RROWS),
                mat_store[e].get_nnz(MQ_FROM_RCOLS)
            );
        }
        printf("    %s\n", str_);
        printf("--------------------------------------------------------\n");
        
    
    }
    printf("=========================================================\n");
    GrB_wait();
    muls_tm.stop();
    

    
    // Enumerate the results
    timer enum_tm;
    if(!empty_res) {
        vector<qedge_t> eval_seq = qplan.get_eval_seq();
        if(qplan.is_csc_friend()) {
            // reverse eval_seq edges
            for(auto& p : eval_seq) p = make_pair(p.second, p.first);
            if(ENABLE_PROFILE) {
                printf("Mat store content:\n");
                for(auto& e : eval_seq) {
                    printf
                    (
                        "  m(%d, %d): %9lu [%9lu, %9lu]\n",
                        e.first, e.second,
                        mat_store[e].get_nnz(),
                        mat_store[e].get_nnz(MQ_FROM_RROWS),
                        mat_store[e].get_nnz(MQ_FROM_RCOLS)
                    );
                }
            }
            printf("--------------------------------------------------------\n");
            enum_tm.start();
            mq_enum_res_col_major(qresults, mat_store, eval_seq, qplan);
            //EXP_mq_enum_res_col_major(qresults, mat_store, eval_seq, qplan);
            enum_tm.stop();
        } else {
            //enum_res_row_major(mat_store, eval_seq, qresults, qplan.get_num_var(), G.get_V());
        }
    } else {
        enum_tm.start();
        enum_tm.stop();
    }
    printf("--------------------------------------------------------\n");
    printf("Num recs: %lu\n", qresults.size());
    printf("--------------------------------------------------------\n");
    printf("Query muls time:        %f\n", muls_tm.interval());
    printf("Query result enum time: %f\n", enum_tm.interval());
    
    // Print results
    if(PRINT_RESULTS) {
        auto proj_vars = qplan.get_proj_vars();
        printf("--------------------------------------------------------\n");
        printf("Results:\n  ");
        for(auto vp : proj_vars) {
            printf("%s ", vp.first.c_str());
        }
        printf("\n");
        for(auto& rec : qresults) {
            printf("  ");
            for(auto vp : proj_vars) {
                printf("<%d> ", rec[vp.second]);
            }
            printf("\n");
        }
        printf("--------------------------------------------------------\n");
    }
    
    for(auto p : mat_store) {
        p.second.destroy();
    }
    GrB_finalize();
}


#endif	/* MQ_QUERY_EXEC_H */

