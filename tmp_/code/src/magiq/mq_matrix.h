/* 
 * File:   mq_matrix.h
 * Author: fuad
 *
 * Created on December 2, 2017, 5:58 PM
 */

/*
 * TODO study edge type in the matrix store, automatic 
 *      for now assumed int8
 */

#ifndef MQ_MATRIX_H
#define	MQ_MATRIX_H

#include "mq_qplan.h"
extern "C"
{
    #include "GraphBLAS.h"
#include "mq_simple_timer.h"
}

class MQMatrix;
class MQGraph;

typedef map<qedge_t, MQMatrix> mat_store_t;


    bool
    f
    (
        const GrB_Index i,
        const GrB_Index j,
        const GrB_Index nrows,
        const GrB_Index ncols,
        const void *x,
        const void *k
    )
    {
        int row = *((int*)k);
        return i == row;
    }

/*
 *
 */
enum nnz_desc_t
{
    MQ_FROM_MAT,
    MQ_FROM_RROWS,
    MQ_FROM_RCOLS
};

/*
 *
 */
class MQMatrix
{
    friend class MQGraph;
    GrB_Matrix  im;
    size_t      num_nodes;
//    bool        zero_clean;
    
public:
    mat_meta_t  meta;
    
    /*
     *
     */
    MQMatrix()
    {
        im = NULL;
        num_nodes = 0;
//        zero_clean = false;
    }
    
    /*
     *
     */
    void
    destroy()
    {
        GrB_Matrix_free(&im);
    }
    
    
    /*
     *
     */
    void
    deep_copy
    (
        const MQMatrix& m,
        bool            T = false
    )
    {
        num_nodes = m.num_nodes;
        init_size(num_nodes);
        GrB_Descriptor desc;
        GrB_Descriptor_new(&desc);
        GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
        if(!T)
            GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
        GrB_transpose(im, NULL, NULL, m.im, desc);
        GrB_Descriptor_free(&desc);
        
    }

    
    /*
     *
     */
    void
    load_from_disk
    (const string& path)
    {
        FILE *fin = fopen(path.c_str(), "r");
        if(!fin) {
            LOGGER.log(string("Can't open file: ") + path, true, false);
            exit(1);
        }
        size_t num_edges = 0;
        int ne_rd = fscanf(fin, "%lu\n", &num_edges);
        int nn_rd = fscanf(fin, "%lu\n", &num_nodes);
        init_size(num_nodes);
        //=================
        // DBG
//        set<pair<int, int> > dbg_set;
        //=================
        for(int i = 0; i < num_edges; ++i) {
            int s_i;
            int p_i;
            int o_i;
            int triple_rd = fscanf(fin, "<%d> <%d> <%d> .\n", &s_i, &p_i, &o_i);
            // DBG ////////////////
//            if(dbg_set.find(make_pair(s_i, o_i)) != dbg_set.end()) {
//                printf("DUPLICATE EDGE (%d, %d)\n", s_i, o_i);
//                exit(1);
//            } else {
//                dbg_set.insert(make_pair(s_i, o_i));
//            }
            ///////////////////////
            GrB_Matrix_setElement_INT8 (im, p_i, s_i, o_i);
        }
        fclose(fin);
        //to make sure mat is actually filled (NON_BLOCKING mode needs barriers)
        GrB_Index nval;
        GrB_Matrix_nvals(&nval, im);
    }
    
    /*
     *
     */
    void
    print()
    {
        GrB_Index   nnz;
        GrB_Matrix_nvals(&nnz, im);
        GrB_Index   *I = (GrB_Index*)malloc(nnz * sizeof(GrB_Index));
        GrB_Index   *J = (GrB_Index*)malloc(nnz * sizeof(GrB_Index));
        int8_t      *X = (int8_t*)malloc(nnz * sizeof(int8_t));
        GrB_Matrix_extractTuples_INT8(I, J, X, &nnz, im);
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
    void
    dbg_write_to_disk(string path)
    {
        FILE* fout = fopen(path.c_str(), "w");
        GrB_Index   nnz;
        GrB_Matrix_nvals(&nnz, im);
        GrB_Index   *I = (GrB_Index*)malloc(nnz * sizeof(GrB_Index));
        GrB_Index   *J = (GrB_Index*)malloc(nnz * sizeof(GrB_Index));
        int8_t      *X = (int8_t*)malloc(nnz * sizeof(int8_t));
        GrB_Matrix_extractTuples_INT8(I, J, X, &nnz, im);
        fprintf(fout, "Number of edges: %d\n", (int)nnz);
        for(int i = 0; i < nnz; ++i) {
            fprintf(fout, "(%ld, %ld): %d\n", I[i], J[i], X[i]);
        }
        free(I);
        free(J);
        free(X);
        fclose(fout);
    }
    
    /*
     *
     */
    size_t
    get_nnz
    (nnz_desc_t nnz_desc = MQ_FROM_MAT)
    {
        if(nnz_desc == MQ_FROM_MAT) {
            GrB_Index nnz;
            GrB_Matrix_nvals(&nnz, im);
            return nnz;
        } else {
            GrB_Vector vec;
            GrB_Vector_new(&vec, GrB_INT8, num_nodes);
            GrB_Descriptor desc;
            GrB_Descriptor_new(&desc);
            GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
            if(nnz_desc == MQ_FROM_RCOLS) { // SSGBLAS does reduction by row
                GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
            }
            GrB_Matrix_reduce_Monoid
                    (vec, NULL, NULL, GxB_LOR_BOOL_MONOID, im, desc);
            GrB_Index vec_nval;
            GrB_Vector_nvals(&vec_nval, vec);
            GrB_Vector_free(&vec);
            GrB_Descriptor_free(&desc);
            return vec_nval;
        }
    }
    
    /*
     *
     */
    size_t
    get_nnodes()
    {
        return num_nodes;
    }
    
    /*
     * Assumes arrays are preallocated
     */
    size_t
    get_row
    (
        GrB_Index*  J,
        int8_t*     X,
        dnode_t     row
    )
    {
        GrB_Vector vec;
        GrB_Vector_new(&vec, GrB_INT8, num_nodes);
        GrB_Descriptor desc;
        GrB_Descriptor_new(&desc);
        GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
        GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
        GrB_Col_extract
            (vec, NULL, NULL, im, GrB_ALL, num_nodes, row, desc);
        GrB_Index nnz;
        GrB_Vector_nvals(&nnz, vec);
        GrB_Vector_extractTuples_INT8(J, X, &nnz, vec);
        GrB_Vector_free(&vec);
        GrB_Descriptor_free(&desc);
        
        // clean from zeros (in case cleaning was not done before)
//        if(!zero_clean) {
//            int i, j;
//            i = 0; j = 0;
//            while(i < nnz) {
//                if(X[i] != 0) {X[j]=X[i]; J[j]=J[i]; j++;}
//                i++;
//            }
//            nnz = j;
//        }
        return nnz;
    }
    
    /*
     * Assumes arrays are preallocated
     */
    size_t
    get_col
    (
        GrB_Index*  I,
        int8_t*     X,
        dnode_t     col
    )
    {
        GrB_Vector vec;
        GrB_Vector_new(&vec, GrB_INT8, num_nodes);
        GrB_Descriptor desc;
        GrB_Descriptor_new(&desc);
        GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
        //GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
        GrB_Col_extract
            (vec, NULL, NULL, im, GrB_ALL, num_nodes, col, desc);
        GrB_Index nnz;
        GrB_Vector_nvals(&nnz, vec);
        GrB_Vector_extractTuples_INT8(I, X, &nnz, vec);
        GrB_Vector_free(&vec);
        GrB_Descriptor_free(&desc);
        
        // clean from zeros (in case cleaning was not done before)
//        if(!zero_clean) {
//            int i, j;
//            i = 0; j = 0;
//            while(i < nnz) {
//                if(X[i] != 0) {X[j]=X[i]; I[j]=I[i]; j++;}
//                i++;
//            }
//            nnz = j;
//        }
        return nnz;
    }
    
    /*
     * Assumes arrays are preallocated
     */
    size_t
    get_all_nz
    (
        GrB_Index*  I,
        GrB_Index*  J,
        int8_t*     X
    )
    {
        GrB_Index nnz;
        GrB_Matrix_nvals(&nnz, im);
        GrB_Matrix_extractTuples_INT8(I, J, X, &nnz, im);
        
        // clean from zeros (in case cleaning was not done before)
//        if(!zero_clean) {
//            int i, j;
//            i = 0; j = 0;
//            while(i < nnz) {
//                if(X[i] != 0) {X[j]=X[i]; I[j]=I[i]; J[j]=J[i]; j++;}
//                i++;
//            }
//            nnz = j;
//        }
        return nnz;
    }
    
    /*
     *
     */
    edge_label_t
    get_element
    (
        const dedge_t&  e,
        bool            from_trans = false
    )
    {
        edge_label_t res = 0;
        int8_t x = -1;
        GrB_Index i, j;
        if(from_trans)  {i = e.second; j = e.first;}
        else            {i = e.first; j = e.second;}
        GrB_Matrix_extractElement_INT8(&x, im, i, j);
        if(x != -1 && x != 0)
            res = x;
        return res;
    }
    
    /*
     *
     */
    void
    set_element
    (
        const dedge_t&  e,
        edge_label_t    val
    )
    {
        GrB_Matrix_setElement_INT8(im, val, e.first, e.second);
    }
    
    /*
     *
     */
    void
    init_size
    (const size_t& num_nodes)
    {
        this->num_nodes = num_nodes;
        GrB_Matrix_new(&im, GrB_INT8, num_nodes, num_nodes);
    }
    
    /*
     *
     */
    bool
    reduce_col_lor
    (const dnode_t& col)
    {
        GrB_Vector vec;
        GrB_Vector_new(&vec, GrB_INT8, num_nodes);
        GrB_Descriptor desc;
        GrB_Descriptor_new(&desc);
        GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
        //GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
        GrB_Col_extract(vec, NULL, NULL, im, GrB_ALL, num_nodes, col, desc);
        int8_t reduction = 0;
        GrB_Vector_reduce_INT8(&reduction, NULL, GxB_LOR_BOOL_MONOID, vec, NULL);
        bool res = (reduction != 0);
        GrB_Vector_free(&vec);
        GrB_Descriptor_free(&desc);
        return res;
    }
    
    /*
     *
     */
    bool
    reduce_row_lor
    (const dnode_t& row)
    {
        GrB_Vector vec;
        GrB_Vector_new(&vec, GrB_INT8, num_nodes);
        GrB_Descriptor desc;
        GrB_Descriptor_new(&desc);
        GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
        GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
        GrB_Col_extract(vec, NULL, NULL, im, GrB_ALL, num_nodes, row, desc);
        int8_t reduction = 0;
        GrB_Vector_reduce_INT8(&reduction, NULL, GxB_LOR_BOOL_MONOID, vec, NULL);
        bool res = (reduction != 0);
        GrB_Vector_free(&vec);
        GrB_Descriptor_free(&desc);
        return res;
    }
    
    /*
     *
     */
    void
    construct_from_meta
    (
        const mat_meta_t&   mat_meta,
        const mat_store_t&  mat_store,
        size_t              nnodes
    )
    {
        init_size(nnodes);
        switch(mat_meta.mat_content) {
            case MQ_ZERO: // ZERO  means there is a one in one place
            {
                dnode_t v1 = mat_meta.one_location.first;
                dnode_t v2 = mat_meta.one_location.second;
                GrB_Matrix_setElement_INT8(im, mat_meta.scalar, v1, v2);
            }
                break;
            case MQ_IDENTITY: // places mat_meta.scalar on diagonal
            {
                for(dnode_t i = 0; i < num_nodes; ++i) {
                    GrB_Matrix_setElement_INT8(im, mat_meta.scalar, i, i);
                }
            }
                break;
            case MQ_STORED: // Means DIAGONALIZE existing matrix
            {               // may be with TRASNPOSE may be not
                GrB_Vector vec;
                GrB_Vector_new(&vec, GrB_INT8, num_nodes);
                GrB_Descriptor desc;
                GrB_Descriptor_new(&desc);
                GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);
                if(mat_meta.mat_transform == MQ_TRANS_DIAG)
                    GrB_Descriptor_set(desc, GrB_INP0, GrB_TRAN);
                GrB_Matrix stored_m =
                        mat_store.find(mat_meta.mat_e)->second.im;
                GrB_Matrix_reduce_Monoid
                        (vec, NULL, NULL, GxB_LOR_BOOL_MONOID, stored_m, desc);
                GrB_Index vec_nval;
                GrB_Vector_nvals(&vec_nval, vec);
                GrB_Index   *I = (GrB_Index*)malloc(vec_nval * sizeof(GrB_Index));
                int8_t      *X = (int8_t*)malloc(vec_nval * sizeof(int8_t));
                GrB_Vector_extractTuples_INT8(I, X, &vec_nval, vec);
                for(int i = 0; i < vec_nval; ++i) {
                    GrB_Index   ind = I[i];
                    int8_t      x = X[i];
                    if(x != 0) // just in case there are zeros
                        GrB_Matrix_setElement_INT8(im, mat_meta.scalar, ind, ind);
                }
                free(I);
                free(X);
                GrB_Vector_free(&vec);
                GrB_Descriptor_free(&desc);
                
                //==================================
                // DBG remove
//                GrB_Matrix dm;
//                GrB_Matrix_new(&dm, GrB_BOOL, nnodes, nnodes);
//                GrB_Descriptor dd;
//                GrB_Descriptor_new(&dd);
//                GrB_Descriptor_set(desc, GrB_INP1, GrB_TRAN);
//                GrB_Matrix MQ_I;
//                GrB_Matrix_new(&MQ_I, GrB_BOOL, nnodes, nnodes);
//                for(dnode_t i = 0; i < num_nodes; ++i) {
//                    GrB_Matrix_setElement_BOOL(MQ_I, 1, i, i);
//                }
//                GrB_wait();
//                timer dbg_tm;
//                dbg_tm.start();
//                GrB_mxm(dm, MQ_I, NULL, GxB_LOR_LAND_BOOL, stored_m, stored_m, dd);
//                dbg_tm.stop();
//                printf("###### TIME TO DIAG: %f\n", dbg_tm.interval());
//                GrB_Matrix_free(&dm);
//                GrB_Matrix_free(&MQ_I);
                //==================================
            }
                break;
        }
    }
    
    /*
     *
     */
    void
    construct_from_mul
    (
        timer&          mxm_tm,
        timer&          trn_tm,
        timer&          zrc_tm,
        const MQMatrix& m1,
        const MQMatrix& m2,
        bool            trans_res = false
    )
    {
        GrB_Matrix tmp_m;
        GrB_Matrix_new(&tmp_m, GrB_INT8, num_nodes, num_nodes);
        mxm_tm.start();
        GrB_mxm(tmp_m, NULL, NULL, GxB_LOR_EQ_INT8, m1.im, m2.im, NULL);    
//        GrB_mxm(tmp_m, NULL, NULL, GxB_PLUS_TIMES_INT8, m1.im, m2.im, NULL);
        mxm_tm.stop();
        zrc_tm.start();
        //clean from explicit zeros
//        GrB_Descriptor desc_clean;
//        GrB_Descriptor_new (&desc_clean);
//        GrB_Descriptor_set (desc_clean, GrB_INP0, GrB_TRAN);
//        GrB_transpose(im, tmp_m, NULL, tmp_m, desc_clean);
//        GrB_Matrix_free(&tmp_m);
        GxB_Matrix_select(im, NULL, NULL, GxB_NONZERO, tmp_m, NULL, NULL);
        zrc_tm.stop();
        trn_tm.start();
        if(trans_res) {
            GrB_Descriptor d;
            GrB_Descriptor_new(&d);
            GrB_Descriptor_set (d, GrB_OUTP, GrB_REPLACE);
            GrB_transpose(im, NULL, NULL, im, d);
            GrB_Descriptor_free(&d);
        }
        trn_tm.stop();
        GrB_Matrix_free(&tmp_m);
    }
    
    /*
     *
     */
    void
    update_from_mul
    (
        timer&          mxm_tm,
        timer&          trn_tm,
        timer&          zrc_tm,
        const MQMatrix& m1,
        const MQMatrix& m2,
        bool            trans_res = false
    )
    {
        // TODO does replacing this with a mxm with replace make things faster?
        GrB_Matrix_free(&im);
        init_size(num_nodes);
        construct_from_mul(mxm_tm, trn_tm, zrc_tm, m1, m2, trans_res);
    }
    
    
    /*
     *
     */
    void
    EXP_construct_from_mul
    (
        timer&          mxm_tm,
        timer&          trn_tm,
        timer&          zrc_tm,
        const MQMatrix& m1,
        const MQMatrix& m2,
        bool            trans_res = false
    )
    {
    }
    
    //==========================================================================
    // EXP

    /*
     *
     */
    void
    EXP_update_from_mul
    (
        timer&          mxm_tm,
        timer&          trn_tm,
        timer&          zrc_tm,
        const MQMatrix& m1,
        MQMatrix& m2,
        bool            trans_res = false
    )
    {
        if(m1.meta.mat_content == MQ_ZERO) { // simple row selection
            dnode_t     row = m1.meta.one_location.first;
            GrB_Index   *J = (GrB_Index*)malloc(num_nodes * sizeof(GrB_Index));
            int8_t      *X = (int8_t*)malloc(num_nodes * sizeof(int8_t));
            size_t      nnz = m2.get_row(J, X, row);
            GrB_Index I[1];
            I[0] = row;
            GrB_Matrix_free(&im);
            init_size(num_nodes);
            GxB_Matrix_subassign_INT8(im, NULL, NULL, 1, I, 1, J, nnz, NULL);
            
            free(J);
            free(X);
            
            
//            dnode_t     row = m1.meta.one_location.first;
//            GrB_Index   I[1];
//            I[0] = row;
//            GrB_Matrix_free(&im);
//            init_size(num_nodes);
//            GrB_Matrix buffer;
//            GrB_Matrix_new(&buffer, GrB_INT8, 1, num_nodes);
//            GrB_Matrix_extract(buffer, NULL, NULL, m2.im, I, 1, GrB_ALL, 0, NULL);
//            GxB_Matrix_subassign(im, NULL, NULL, buffer, I, 1, GrB_ALL, 0, NULL);
//            GrB_Matrix_free(&buffer);
            
//            GxB_SelectOp op;
//            bool (*fp)(const GrB_Index,
//                    const GrB_Index, 
//                    const GrB_Index, 
//                    const GrB_Index, 
//                    const void*, 
//                    const void*);
//            fp = f;
//            GxB_SelectOp_new(&op, (bool*)fp, NULL);
//            GxB_Matrix_select(im, NULL, NULL, op, m2.im, (const dnode_t*)&row, NULL);
//            GxB_SelectOp_free(&op);
        }
//        else if(m1.meta.mat_transform == MQ_DIAGONALIZE || m1.meta.mat_transform == MQ_TRANS_DIAG) {
//            GrB_Vector vec1, vec2;
//            GrB_Vector_new(&vec1, GrB_INT8, num_nodes);
//            GrB_Vector_new(&vec2, GrB_INT8, num_nodes);
//            GrB_Matrix_reduce_Monoid
//                (vec1, NULL, NULL, GxB_LOR_BOOL_MONOID, m1.im, NULL);
//            GxB_Vector_select(vec2, NULL, NULL, GxB_NONZERO, vec1, NULL, NULL);
//            GrB_Index vec_nval;
//            GrB_Vector_nvals(&vec_nval, vec2);
//            GrB_Index   *I = (GrB_Index*)malloc(vec_nval * sizeof(GrB_Index));
//            int8_t      *X = (int8_t*)malloc(vec_nval * sizeof(int8_t));
//            GrB_Vector_extractTuples_INT8(I, X, &vec_nval, vec2);                
//            
//            GrB_Matrix_free(&im);
//            init_size(num_nodes);
//            GrB_Matrix buffer;
//            GrB_Matrix_new(&buffer, GrB_INT8, vec_nval, num_nodes);
//            GrB_Matrix_extract(buffer, NULL, NULL, m2.im, I, vec_nval, GrB_ALL, 0, NULL);
//            GxB_Matrix_subassign(im, NULL, NULL, buffer, I, vec_nval, GrB_ALL, 0, NULL);
//            
//            free(I);
//            free(X);
//            GrB_Vector_free(&vec1);
//            GrB_Vector_free(&vec2);
//            GrB_Matrix_free(&buffer);          
//        }
        else {
            GrB_Matrix_free(&im);
            init_size(num_nodes);
            construct_from_mul(mxm_tm, trn_tm, zrc_tm, m1, m2, trans_res);
        }
    }
};


class MQGraph
{
    MQMatrix dgraph;
    MQMatrix dgraph_T;
    size_t   num_nodes;
    size_t   num_edges;
    
    //==================
    map<edge_label_t, pair<MQMatrix, MQMatrix> > dgraph_pmap;
    //==================
    
public:
    /*
     *
     */
    void
    destroy()
    {
        dgraph.destroy();
        dgraph_T.destroy();
        
        for(auto& p : dgraph_pmap) {
            p.second.first.destroy();
            p.second.second.destroy();
        }
    }
    
    /*
     *
     */
    void
    load_graph
    (const string& graph_path)
    {
//        timer tm;
//        tm.start();
//        dgraph.load_from_disk(graph_path);
//        tm.stop();
//        printf("    Graph load time:%s%f\n", string(6, ' ').c_str(), tm.interval());
//        tm.start();
//        dgraph_T.deep_copy(dgraph, true);
//        tm.stop();
//        printf("    Graph transpose time:%s%f\n", string(1, ' ').c_str(), tm.interval());
//        num_nodes = dgraph.get_nnodes();
//        num_edges = dgraph.get_nnz();
        
        //==============
        // Fill individual predicate graph layers
        timer tm;
        FILE *fin = fopen(graph_path.c_str(), "r");
        if(!fin) {
            LOGGER.log(string("Can't open file: ") + graph_path, true, false);
            exit(1);
        }
        size_t num_edges = 0;
        int ne_rd = fscanf(fin, "%lu\n", &num_edges);
        int nn_rd = fscanf(fin, "%lu\n", &num_nodes);
        printf("    Will read %lu lines...\n", num_edges);
        dgraph.init_size(num_nodes);
        dgraph_T.init_size(num_nodes);
        tm.start();
 ////////////////////////////////////////////////////////       
        // DBG
        FILE* fout = fopen("dbg_read_graph", "w");
/////////////////////////////////////////////////////////
        for(int i = 0; i < num_edges; ++i) {
            int s_i;
            int p_i;
            int o_i;
            int triple_rd = fscanf(fin, "<%d> <%d> <%d> .\n", &s_i, &p_i, &o_i);
            GrB_Matrix_setElement_INT8 (dgraph.im, p_i, s_i, o_i);
            GrB_Matrix_setElement_INT8 (dgraph_T.im, p_i, o_i, s_i);
            fprintf(fout, "%d %d %d\n", s_i, o_i, p_i);
//            if(dgraph_pmap.find(p_i) == dgraph_pmap.end()) {
//                // new layer
//                dgraph_pmap[p_i] = make_pair(MQMatrix(), MQMatrix());
//                dgraph_pmap[p_i].first.init_size(num_nodes);
//                dgraph_pmap[p_i].second.init_size(num_nodes);
//            }
//            dgraph_pmap[p_i].first.set_element(make_pair(s_i, o_i), 1);
//            dgraph_pmap[p_i].second.set_element(make_pair(o_i, s_i), 1);
        }
        
        ///////////////////////
        fclose(fout);
        ///////////////////////
        
        fclose(fin);
        tm.stop();
        printf("    File read time:%s%f\n", string(6, ' ').c_str(), tm.interval());
        //to make sure mat is actually filled (NON_BLOCKING mode needs barriers)
        num_nodes = dgraph.get_nnodes();
        
        tm.start();
        num_edges = dgraph.get_nnz();
        tm.stop();
        printf("    Graph construction time:%s%f\n", string(1, ' ').c_str(), tm.interval());
        size_t dum;
        tm.start();
        num_edges = dgraph_T.get_nnz();
        tm.stop();
        printf("    Transpose construction time:%s%f\n", string(1, ' ').c_str(), tm.interval());
        printf("    SS:GraphBLAS: num edges: %lu num nodes: %lu\n", num_edges, num_nodes);
//        tm.start();
//        for(auto& p : dgraph_pmap) {
//            dum = p.second.first.get_nnz();
//            dum = p.second.second.get_nnz();
//            printf("        P(%3d) has [%10lu] edges\n", p.first, dum);
//        }
//        tm.stop();
//        printf("    Layers construction time:%s%f\n", string(1, ' ').c_str(), tm.interval());
        //==============
    }
    
    /*
     *
     */
    size_t
    get_V()
    {
        return num_nodes;
    }
    
    /*
     *
     */
    size_t
    get_E()
    {
        return num_edges;
    }
    
    /*
     *
     */
    MQMatrix
    get_mat
    (
        bool         T = false,
        edge_label_t p = -1
    )
    {
        if(p == -1) {
            return ( T ? dgraph_T : dgraph );
        } else {
            return ( T ? dgraph_pmap[p].second : dgraph_pmap[p].first ); 
        }
    }
    
};



#endif	/* MQ_MATRIX_H */

