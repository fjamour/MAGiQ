/* 
 * File:   query_exec.h
 * Author: fuad
 *
 * Created on November 20, 2017, 1:14 PM
 */

#ifndef QUERY_EXEC_H
#define	QUERY_EXEC_H

#include "mq_qplan.h"

#include "mq_simple_logger.h"
#include "mq_simple_timer.h"

#include <cstdio>

extern "C" 
{
    #include "GraphBLAS.h"
}


typedef pair<int, int> edge_t;
typedef map<edge_t, GrB_Matrix> old_mat_store_t;

void
print_mat_SSGBLAS
(const GrB_Matrix mat)
{
    GrB_Index   nnz;
    GrB_Matrix_nvals(&nnz, mat);
    GrB_Index   *I = (GrB_Index*)malloc(nnz * sizeof(GrB_Index));
    GrB_Index   *J = (GrB_Index*)malloc(nnz * sizeof(GrB_Index));
    int8_t      *X = (int8_t*)malloc(nnz * sizeof(int8_t));
    GrB_Matrix_extractTuples_INT8(I, J, X, &nnz, mat);
    printf("Number of edges: %d\n", (int)nnz);
    for(int i = 0; i < nnz; ++i) {
        printf("(%ld, %ld): %d\n", I[i], J[i], X[i]);
    }
    free(I);
    free(J);
    free(X);
}

/*
 *
 */
int
num_rows_wnnz_SSGBLAS
(const GrB_Matrix mat)
{
    GrB_Vector vec;
    GrB_Index num_nodes;
    GrB_Matrix_nrows(&num_nodes, mat);
    GrB_Vector_new(&vec, GrB_INT8, num_nodes);
    GrB_Descriptor desc;
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
    //GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
    GrB_Matrix_reduce_Monoid(vec, NULL, NULL, GrB_LOR_BOOL_MONOID, mat, desc);
    GrB_Index vec_nval;
    GrB_Vector_nvals(&vec_nval, vec);
    GrB_Vector_free(&vec);
    GrB_Descriptor_free(&desc);
    return (int)vec_nval;
}
/*
 *
 */
int
num_cols_wnnz_SSGBLAS
(const GrB_Matrix mat)
{
    GrB_Vector vec;
    GrB_Index num_nodes;
    GrB_Matrix_nrows(&num_nodes, mat);
    GrB_Vector_new(&vec, GrB_INT8, num_nodes);
    GrB_Descriptor desc;
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
    GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
    GrB_Matrix_reduce_Monoid(vec, NULL, NULL, GrB_LOR_BOOL_MONOID, mat, desc);
    GrB_Index vec_nval;
    GrB_Vector_nvals(&vec_nval, vec);
    GrB_Vector_free(&vec);
    GrB_Descriptor_free(&desc);
    return (int)vec_nval;
}

class mat_SSGBLAS_t
{
    GrB_Matrix m;
public:
    /*
     * returns a *shallow copy* of m
     */
    GrB_Matrix
    get_mat()
    {
        return m;
    }
    /*
     *
     */
    bool
    reduce_row
    (
        mat_meta_t& mat_meta,
        old_mat_store_t& mat_store,
        int num_nodes,
        int row_index
    )
    {
        GrB_Matrix m = mat_store[mat_meta.mat_e];
        GrB_Vector vec;
        GrB_Vector_new(&vec, GrB_INT8, num_nodes);
        GrB_Descriptor desc;
        GrB_Descriptor_new(&desc);
        GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
        if(!(mat_meta.mat_transform == MQ_TRANSPOSE || mat_meta.mat_transform == MQ_TRANS_DIAG))
            GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
        GrB_Col_extract(vec, NULL, NULL, m, GrB_ALL, num_nodes, row_index, desc);
        int8_t reduction = 0;
        GrB_Vector_reduce_INT8(&reduction, NULL, GrB_LOR_BOOL_MONOID, vec, NULL);
        bool res = (reduction != 0);
        GrB_Descriptor_free(&desc);
        GrB_Vector_free(&vec);
        return res;
    }
    /*
     *
     */
    void
    construct_mat
    (
        mat_meta_t& mat_meta,
        old_mat_store_t& mat_store,
        int num_nodes
    )
    {
        GrB_Matrix_new(&m, GrB_INT8, num_nodes, num_nodes);
        switch(mat_meta.mat_content) {
            case MQ_ZERO: // ZERO  means there is a one in one place
            {
                int v1 = mat_meta.one_location.first;
                int v2 = mat_meta.one_location.second;
                GrB_Matrix_setElement_INT8(m, mat_meta.scalar, v1, v2);
            }
                break;
            case MQ_IDENTITY:
            {
                for(int i = 0; i < num_nodes; ++i) {
                    GrB_Matrix_setElement_INT8(m, mat_meta.scalar, i, i);
                }
            }
                break;
            case MQ_STORED:
            {
                // if we reach this point, there has to be DIAGONALIZE
                // , may be with TRASNPOSE may be not
                GrB_Vector vec;
                GrB_Vector_new(&vec, GrB_INT8, num_nodes);
                GrB_Descriptor desc;
                GrB_Descriptor_new(&desc);
                GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
                if(mat_meta.mat_transform == MQ_TRANS_DIAG)
                    GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
                GrB_Matrix stored_m = mat_store[mat_meta.mat_e];
                GrB_Matrix_reduce_Monoid(vec, NULL, NULL, GrB_LOR_BOOL_MONOID, stored_m, desc);
                GrB_Index vec_nval;
                GrB_Vector_nvals(&vec_nval, vec);
                GrB_Index   *I = (GrB_Index*)malloc(vec_nval * sizeof(GrB_Index));
                int8_t      *X = (int8_t*)malloc(vec_nval * sizeof(int8_t));
                GrB_Vector_extractTuples_INT8(I, X, &vec_nval, vec);
                for(int i = 0; i < vec_nval; ++i) {
                    GrB_Index   ind = I[i];
                    int8_t      x = X[i];
                    GrB_Matrix_setElement_INT8(m, mat_meta.scalar, ind, ind);
                }
                free(I);
                free(X);
                GrB_Vector_free(&vec);
                GrB_Descriptor_free(&desc);
            }
                break;
        }
    }
    /*
     *
     */
    void
    init
    (int num_nodes)
    {
        GrB_Matrix_new(&m, GrB_INT8, num_nodes, num_nodes);
    }
    /*
     *
     */
    void
    print()
    {
        print_mat_SSGBLAS(m);
    }
    /*
     *
     */
    mat_SSGBLAS_t()
    {
        m = NULL;
    }
    /*
     *
     */
    ~mat_SSGBLAS_t()
    {
        if(m) GrB_Matrix_free(&m);
    }
};

class data_graph_SSGBLAS_t
{
    GrB_Matrix dgraph;
    GrB_Matrix dgraph_T;
    int num_nodes;
    int num_edges;
public:
    /*
     *
     */
    int
    get_edge_label
    (
        bool from_transpose,
        edge_t e
    )
    {
        int res = 0;
        int8_t x = -1;
        GrB_Index i, j;
        if(from_transpose)  {i = e.second; j = e.first;}
        else                {i = e.first; j = e.second;}
        GrB_Matrix_extractElement_INT8(&x, dgraph, i, j);
        if(x != -1)
            res = x;
        return res;
    }
    /*
     *
     */
    int
    get_V()
    {
        return num_nodes;
    }
    /*
     * returns a *shallow copy* of dgraph
     */
    GrB_Matrix
    get_mat(bool transpose = false)
    {
        if(transpose)
            return dgraph_T;
        else
            return dgraph;
    }
    /*
     *
     */
    void
    load_graph
    (string encoded_graph_path)
    {
        // 1. get the number of nodes from the input file
        // 2. construct the matrix
        // 3. fill the matrix
        FILE *fin = fopen(encoded_graph_path.c_str(), "r");
        if(!fin) {
            LOGGER.log("Can't open graph file!\n", true, false);
            exit(1);
        }
        int ne_rd = fscanf(fin, "%d\n", &num_edges);
        int nn_rd = fscanf(fin, "%d\n", &num_nodes);
        GrB_Matrix_new(&dgraph, GrB_INT8, num_nodes, num_nodes);
        for(int i = 0; i < num_edges; ++i) {
            int s_i;
            int p_i;
            int o_i;
            int triple_rd = fscanf(fin, "<%d> <%d> <%d> .\n", &s_i, &p_i, &o_i);
            GrB_Matrix_setElement_INT8 (dgraph, p_i, s_i, o_i);
        }
        fclose(fin);
        //to make sure mat is actually filled (NON_BLOCKING mode needs barriers)
        GrB_Index nval;
        GrB_Matrix_nvals(&nval, dgraph);
        
        GrB_Matrix_new(&dgraph_T, GrB_INT8, num_nodes, num_nodes);
        GrB_transpose(dgraph_T, NULL, NULL, dgraph, NULL);
    }
    /*
     *
     */
    void
    print()
    {
        print_mat_SSGBLAS(dgraph);
    }
    /*
     *
     */
    data_graph_SSGBLAS_t()
    {
        dgraph = NULL;
    }
    /*
     *
     */
    ~data_graph_SSGBLAS_t()
    {
        if(dgraph) GrB_Matrix_free(&dgraph);
    }

};


/*
 *
 */
void
enum_res_row_major
(
    old_mat_store_t& mat_store,
    vector<edge_t> eval_seq,
    vector<vector<int> >& qresults,
    int num_vars,
    int num_nodes
)
{
    // 1. use the first matrix (all nnz) to make the first bunch of records
    // 2. use one var a time to fill one more slot in each record
    edge_t e;
    GrB_Matrix m;
    e = eval_seq[0];
    m = mat_store[e];
    GrB_Index   nnz;
    GrB_Matrix_nvals(&nnz, m);
    GrB_Index   *I = (GrB_Index*)malloc(nnz * sizeof(GrB_Index));
    GrB_Index   *J = (GrB_Index*)malloc(nnz * sizeof(GrB_Index));
    int8_t      *X = (int8_t*)malloc(nnz * sizeof(int8_t));
    GrB_Matrix_extractTuples_INT8(I, J, X, &nnz, m);
    
    int ii, jj;
    ii = e.first;
    jj = e.second;
    for(int i = 0; i < nnz; ++i) {
        vector<int> rec(num_vars, -1);
        rec[ii] = I[i];
        rec[jj] = J[i];
        qresults.push_back(rec);
    }
    free(I);
    free(J);
    free(X);
    
    int res_start;
    int res_stop = 0;
    for(int e_i = 1; e_i < eval_seq.size(); ++e_i) {
        edge_t e = eval_seq[e_i];
        m = mat_store[e];
        ii = e.first;
        jj = e.second;
//        res_start = res_stop;
//        res_stop = qresults.size();
        res_start = 0;
        res_stop = qresults.size();
        for(int rr = res_start; rr < res_stop; ++rr) {
            //vector<int> rec = qresults[rr];
            int row_index = qresults[rr][ii];
            if(row_index == -1) // already filled
                continue;
            //extracl row row_index from m
            GrB_Vector vec;
            GrB_Vector_new(&vec, GrB_INT8, num_nodes);
            GrB_Descriptor desc;
            GrB_Descriptor_new(&desc);
            GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
            GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
            GrB_Col_extract(vec, NULL, NULL, m, GrB_ALL, num_nodes, row_index, desc);
            GrB_Vector_nvals(&nnz, vec);
            GrB_Index   *I = (GrB_Index*)malloc(nnz * sizeof(GrB_Index));
            int8_t      *X = (int8_t*)malloc(nnz * sizeof(int8_t));
            GrB_Vector_extractTuples_INT8(I, X, &nnz, vec);
            for(int i = 0; i < nnz; ++i) {
                if(qresults[rr][jj] == -1) {
                    qresults[rr][jj] = I[i];
                } else {
                    vector<int> new_rec = qresults[rr];
                    new_rec[jj] = I[i];
                    qresults.push_back(new_rec);
                }
            }
            GrB_Descriptor_free(&desc);
            GrB_Vector_free(&vec);
            free(I);
            free(X);
        }
    }
}


/*
 *
 */
void
enum_res_col_major
(
    old_mat_store_t& mat_store,
    vector<edge_t> eval_seq,
    vector<vector<int> >& qresults,
    int num_vars,
    int num_nodes
)
{
    // 1. use the first matrix (all nnz) to make the first bunch of records
    // 2. use one var a time to fill one more slot in each record
    
    //reverse eval_seq edges
    for(auto& p : eval_seq) {
        p = make_pair(p.second, p.first);
    }
    
    edge_t e;
    GrB_Matrix m;
    e = eval_seq[0];
    m = mat_store[e];
    GrB_Index   nnz;
    
//    GrB_Matrix_nvals(&nnz, m);
//    GrB_Index   *I = (GrB_Index*)malloc(nnz * sizeof(GrB_Index));
//    GrB_Index   *J = (GrB_Index*)malloc(nnz * sizeof(GrB_Index));
//    int8_t      *X = (int8_t*)malloc(nnz * sizeof(int8_t));
//    GrB_Matrix_extractTuples_INT8(I, J, X, &nnz, m);
    
    
    GrB_Index   *I = (GrB_Index*)malloc(num_nodes * sizeof(GrB_Index));
    GrB_Index   *J = (GrB_Index*)malloc(num_nodes * sizeof(GrB_Index));
    int8_t      *X = (int8_t*)malloc(num_nodes * sizeof(int8_t));
    GrB_Matrix_nvals(&nnz, m);
    GrB_Matrix_extractTuples_INT8(I, J, X, &nnz, m);    
    
    int ii, jj;
    ii = e.first;
    jj = e.second;
    for(int i = 0; i < nnz; ++i) {
        vector<int> rec(num_vars, -1);
        rec[ii] = I[i];
        rec[jj] = J[i];
        qresults.push_back(rec);
    }
//    free(I);
//    free(J);
//    free(X);

    GrB_Vector vec;
    GrB_Vector_new(&vec, GrB_INT8, num_nodes);
    int res_start;
    int res_stop = 0;
    for(int e_i = 1; e_i < eval_seq.size(); ++e_i) {
        edge_t e = eval_seq[e_i];
        m = mat_store[e];
        ii = e.first;
        jj = e.second;
//        res_start = res_stop;
//        res_stop = qresults.size();
        res_start = 0;
        res_stop = qresults.size();
        for(int rr = res_start; rr < res_stop; ++rr) {
            //vector<int> rec = qresults[rr];
            int col_index = qresults[rr][jj];
            if(col_index == -1) // already filled
                continue;
            //extracl row row_index from m

            GrB_Descriptor desc;
            GrB_Descriptor_new(&desc);
            GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
            //GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
            GrB_Col_extract(vec, NULL, NULL, m, GrB_ALL, num_nodes, col_index, desc);
            GrB_Vector_nvals(&nnz, vec);
//            GrB_Index   *I = (GrB_Index*)malloc(nnz * sizeof(GrB_Index));
//            int8_t      *X = (int8_t*)malloc(nnz * sizeof(int8_t));
            GrB_Vector_extractTuples_INT8(I, X, &nnz, vec);
            for(int i = 0; i < nnz; ++i) {
                if(qresults[rr][ii] == -1) {
                    qresults[rr][ii] = I[i];
                } else {
                    vector<int> new_rec = qresults[rr];
                    new_rec[ii] = I[i];
                    qresults.push_back(new_rec);
                }
            }
            GrB_Descriptor_free(&desc);
//            free(I);
//            free(X);
        }
    }
    GrB_Vector_free(&vec);
    free(I);
    free(J);
    free(X);
}

/*
 * 
 */
void
execute_query_plan
(
    vector<vector<int> >& res_vec,
    QPlan qplan,
    data_graph_SSGBLAS_t& G
)
{
    timer tm;
    GrB_init (GrB_NONBLOCKING);
    old_mat_store_t mat_store;
    tm.start();
    auto op_vec = qplan.get_op_vec();
    for(auto op : op_vec) {
        edge_t e = op.dst_mat_e;
        switch(op.type) {
            case MQ_TWOCB:
            {
                GrB_Matrix store_m;
                // init and construct store_m
                GrB_Matrix_new(&store_m, GrB_INT8, G.get_V(), G.get_V());
                if(G.get_edge_label(op.G_transpose, make_pair(op.c1, op.c2)) == op.edge_label) {
                    GrB_Matrix_setElement_INT8(store_m, 1, op.c1, op.c2);
                }
                mat_store[e] = store_m;
//                printf("###m_(%d, %d)\n", e.first, e.second);
//                print_mat_SSGBLAS(store_m);
            }
            break;
            case MQ_TWOCM:
            {
                GrB_Matrix store_m;
                GrB_Matrix_new(&store_m, GrB_INT8, G.get_V(), G.get_V());
                mat_SSGBLAS_t m;
                //int c = ( qplan.csc_friendly ? op.c2 : op.c1 );
                bool valid_row = m.reduce_row(op.m1, mat_store, G.get_V(), op.mc);
                if(valid_row) {
                    bool G_got_edge =
                    (G.get_edge_label(op.G_transpose, make_pair(op.c1, op.c2)) == op.edge_label);
                    if(G_got_edge) {
                        GrB_Matrix_setElement_INT8(store_m, 1, op.c1, op.c2);
                    }
                }
                mat_store[e] = store_m;
//                printf("###m_(%d, %d)\n", e.first, e.second);
//                print_mat_SSGBLAS(store_m);
            }
            break;
            case MQ_MUL:
            {
                GrB_Matrix store_m;
                GrB_Matrix m1, m2;
                mat_SSGBLAS_t M1, M2;
                GrB_Descriptor desc;
                GrB_Descriptor_new(&desc);
                if(op.m1.mat_content == MQ_GRAPH) {
                    m1 = G.get_mat();
                    if(op.m1.mat_transform == MQ_TRANSPOSE) {
                        //GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
                        m1 = G.get_mat(true);
                    }
                } else if(op.m1.mat_content == MQ_STORED && op.m1.mat_transform == MQ_NOTRANSFORM) {
                    m1 = mat_store[op.m1.mat_e];
                } else {
                    M1.construct_mat(op.m1, mat_store, G.get_V());
                    m1 = M1.get_mat();
                }
                if(op.m2.mat_content == MQ_GRAPH) {
                    m2 = G.get_mat();
                    if(op.m2.mat_transform == MQ_TRANSPOSE) {
                        //GrB_Descriptor_set(desc, GrB_INP1, GrB_TRAN);
                        m2 = G.get_mat(true);
                    }
                } else if(op.m2.mat_content == MQ_STORED && op.m2.mat_transform == MQ_NOTRANSFORM) {
                    m2 = mat_store[op.m2.mat_e];
                } else {
                    M2.construct_mat(op.m2, mat_store, G.get_V());
                    m2 = M2.get_mat();
                }
                
                
                GrB_Matrix tmp_m;
                GrB_Matrix_new(&tmp_m, GrB_INT8, G.get_V(), G.get_V());
                GrB_Matrix_new(&store_m, GrB_INT8, G.get_V(), G.get_V());
                //GrB_mxm(tmp_m, NULL, NULL, GrB_MAX_ISEQ_INT8, m1, m2, desc);
                GrB_mxm(tmp_m, NULL, NULL, GrB_LOR_EQ_INT8, m1, m2, desc);
                
                if(op.trans_res) {
                    GrB_Descriptor d;
                    GrB_Descriptor_new(&d);
                    GrB_Descriptor_set (d, GrB_OUTP, GrB_REPLACE);
                    GrB_transpose(tmp_m, NULL, NULL, tmp_m, d);
                    GrB_Descriptor_free(&d);
                }
                
                //clean from explicit zeros
                GrB_Descriptor desc_clean;
                GrB_Descriptor_new (&desc_clean);
                GrB_Descriptor_set (desc_clean, GrB_INP0, GrB_TRAN);
                GrB_transpose(store_m, tmp_m, NULL, tmp_m, desc_clean);
                
//                printf("###%s\n", print_mat(op.m1).c_str());
//                print_mat_SSGBLAS(m1);
//                printf("###%s\n", print_mat(op.m2).c_str());
//                print_mat_SSGBLAS(m2);
//                printf("###m_(%d, %d)\n", e.first, e.second);
//                print_mat_SSGBLAS(store_m);

                GrB_Index nnz_stored;
                GrB_Index nnz_new;
                if(op.dst_desc == MQ_UPDATE_MAT) {
                    GrB_Matrix_nvals(&nnz_stored, mat_store[e]);
                    GrB_Matrix_nvals(&nnz_new, store_m);
                    int br = num_rows_wnnz_SSGBLAS(mat_store[e]);
                    int bc = num_cols_wnnz_SSGBLAS(mat_store[e]);
                    int ar = num_rows_wnnz_SSGBLAS(store_m);
                    int ac = num_cols_wnnz_SSGBLAS(store_m);
                    printf("---MUL(%d, %d): %d [%d, %d] to %d [%d, %d]\n",
                            op.dst_mat_e.first, op.dst_mat_e.second,
                            (int)nnz_stored, br, bc,
                            (int)nnz_new, ar, ac);
                    
                    GrB_Matrix mm = mat_store[e];
                    GrB_Matrix_free(&mm);
                } else {
                    GrB_Matrix_nvals(&nnz_new, store_m);
                    int ar = num_rows_wnnz_SSGBLAS(store_m);
                    int ac = num_cols_wnnz_SSGBLAS(store_m);
                    printf("---MUL(%d, %d): %d [%d, %d]\n",
                            op.dst_mat_e.first, op.dst_mat_e.second,
                            (int)nnz_new, ar, ac);
                }
                

                mat_store[e] = store_m;
                GrB_Matrix_free(&tmp_m);
                GrB_Descriptor_free(&desc);
                GrB_Descriptor_free(&desc_clean);
            }
            break;
        }
    }
    
    //make sure all ops are done
    GrB_wait();
    tm.stop();
    printf("Time to do query operations: %f\n", tm.interval());

//    printf("*******************************************************\n");
//    for(auto p : mat_store) {
//        GrB_Matrix m = p.second;
//        printf("\n# m_(%d, %d):\n", p.first.first, p.first.second);
//        print_mat_SSGBLAS(m);
//    }
//    printf("*******************************************************\n");
    tm.start();
    vector<qedge_t> eval_seq = qplan.get_eval_seq();
    vector<vector<int> > qresults;
    if(qplan.is_csc_friend()) {
        enum_res_col_major(mat_store, eval_seq, qresults, qplan.get_num_var(), G.get_V());
    } else {
        enum_res_row_major(mat_store, eval_seq, qresults, qplan.get_num_var(), G.get_V());
    }
    tm.stop();
    printf("Time to enumerate results: %f\n", tm.interval());
    // clean up results
    tm.start();
    vector<vector<int> > actual_res;
    int last_var = eval_seq.back().second;
    for(int i = 0; i < qresults.size(); ++i) {
        vector<int> rec = qresults[i];
        if(rec[last_var] == -1)
            continue;
        vector<pair<int, int> > eq_pairs;
        eq_pairs = qplan.get_eq_var_pairs();
        bool valid_rec = true;
        for(auto p : eq_pairs) {
            int i1 = p.first;
            int i2 = p.second;
            if(rec[i1] != rec[i2]) {
                valid_rec = false;
                break;
            }
        }
        if(valid_rec)
            actual_res.push_back(rec);
    }
    tm.stop();
    printf("Time to clean results: %f\n", tm.interval());
    // print results
    tm.start();
    int first_var_index = -1;
    printf("*******************************************************:\n");
    printf("Results:\n");
    for(int i = 0; i < qplan.get_num_actual_vars(); ++i) {
        string s = qplan.get_var_name(i);
        if(s[0] == '?') {
            if(first_var_index == -1) first_var_index = i;
            printf("%s ", s.c_str());
        }
    }
    printf("\n");
//    if(actual_res.size() < 200) {
//        for(auto& rec : actual_res) {
//            string s = "[";
//            for(int i = first_var_index; i < qplan.qgraph.get_num_actual_vars(); ++i) {
//                if(i != first_var_index) s += ", ";
//                s += to_string(rec[i]);
//            }
//            s += "]";
//            printf("%s\n", s.c_str());
//        }
//    }
    printf("*******************************************************:\n");
    printf("Number of recs: %d\n", (int)actual_res.size());
    tm.stop();
    printf("Time to print results: %f\n", tm.interval());
    // clean up
    for(auto p : mat_store) {
        GrB_Matrix m = p.second;
        GrB_Matrix_free(&m);
    }
    
    GrB_finalize();
}




#endif	/* QUERY_EXEC_H */

