/* 
 * File:   mq_qplan.h
 * Author: fuad
 *
 * Created on December 2, 2017, 10:33 AM
 */

#ifndef MQ_QPLAN_H
#define	MQ_QPLAN_H


#include "mq_qgraph.h"

using namespace std;

/*
 * 
 */
enum mat_content_t {
    MQ_GRAPH,
    MQ_ZERO,
    MQ_IDENTITY,
    MQ_STORED,
    MQ_NONE
};

/*
 * 
 */
enum mat_transform_t {
    MQ_TRANSPOSE,
    MQ_DIAGONALIZE,
    MQ_TRANS_DIAG,
    MQ_NOTRANSFORM
};

/*
 * Defines one operation (mxm). Collectively with other query ops defines
 * a query plan 
 */
enum op_dst_desc_t {
    MQ_CONSTRUCT_MAT,
    MQ_UPDATE_MAT
};

/*
 *
 */
enum op_type_t {
    MQ_MUL,     //multiplication of specified matrices
    MQ_TWOCB,   //first edge two constants, checks single entry in GRAPH
    MQ_TWOCM    //intermediate edge two constants, checks single entry in GRAPH,
                //and row in specified matrix
};

/*
 * Meta data sufficient to construct a matrix to take part in the
 * multiplications sequence of a query plan.
 */
struct mat_meta_t {
    mat_content_t   mat_content;
    mat_transform_t mat_transform;
    // in case of constants, specifies a '1' loc
    dedge_t         one_location;    
    // edge label (typically scalar mul with m)
    edge_label_t   scalar;                     
    // corresponding edge in q graph
    qedge_t         mat_e;           
    
    mat_meta_t()
    {
        mat_content = MQ_NONE;
        mat_transform = MQ_NOTRANSFORM;
        one_location = NULLP;
        scalar = -1;
        mat_e = NULLP;
    }
    
    void
    set_vals
    (
        mat_content_t   mat_content,
        mat_transform_t mat_transform,
        dedge_t         one_location,
        edge_label_t    scalar,
        qedge_t         mat_e
    )
    {
        this->mat_content = mat_content;
        this->mat_transform = mat_transform;
        this->one_location = one_location;
        this->scalar = scalar;
        this->mat_e = mat_e;
    }
    
    string
    stringify()
    {
        char buffer[1024];
        string ms;
        switch(mat_content) {
            case MQ_GRAPH:
            {
                ms = "G";
                ms += ( mat_transform == MQ_TRANSPOSE ? ".T()" : "" );
            }
                break;
            case MQ_ZERO: // ZERO  means there is a one in one place
            {
                string s = "";
                if(scalar > 1) s = string("*") + to_string(scalar);
                sprintf(buffer,
                        "{1@(%d,%d)}%s",
                        one_location.first, one_location.second,
                        s.c_str());
                ms =  buffer;
            }
                break;
            case MQ_IDENTITY:
            {
                sprintf(buffer,
                        "I*%d",
                        scalar);
                ms =  buffer;
            }
                break;
            case MQ_STORED:
            {
                string T, D, S;
                T = ""; D = ""; S = "";
                if(mat_transform == MQ_TRANSPOSE) T = ".T()";
                if(mat_transform == MQ_DIAGONALIZE) D = ".D()";
                if(mat_transform == MQ_TRANS_DIAG) {T = ".T()"; D = ".D()";}
                if(scalar > 1) {S = "*"; S += to_string(scalar);}
                sprintf(buffer,
                        "m_(%d,%d)%s%s%s",
                        mat_e.first, mat_e.second,
                        T.c_str(),
                        D.c_str(),
                        S.c_str());
                ms =  buffer;
            }
                break;
        }
        return ms;
    }
};

/*
 *
 */
struct query_op_t {
    op_type_t type;
    // m1, m2 for multiplication
    mat_meta_t m1;
    mat_meta_t m2;
    // for TWOCB and TWOCM
    dnode_t         c1, c2;
    edge_label_t   edge_label; 
    bool            G_transpose;   
    // for TWOCM
    dnode_t         mc;
    // always needed to specify result matrix
    qedge_t         dst_mat_e;
    op_dst_desc_t   dst_desc;
    //
    bool trans_res;
    
    query_op_t()
    {
        trans_res = false;
    }
    
    string
    stringify()
    {
        char buffer[1024];
        string os;
        switch(type) {
            case MQ_TWOCB:
            {
                string T = ( G_transpose ? ".T()" : "" );
                sprintf(buffer,
                    "m_(%d,%d) = ( G%s[%d,%d] == %d ? {1@(%d,%d)} : ZERO )",
                    dst_mat_e.first, dst_mat_e.second,
                    T.c_str(), c1, c2, edge_label,
                    c1, c2);
            }
                break;
            case MQ_TWOCM:
            {
                string GT = ( G_transpose ? ".T()" : "" );
                string MT = ( m1.mat_transform == MQ_TRANS_DIAG ? ".T()" : "" );
                sprintf(buffer,
                    "m_(%d,%d) = ( G%s[%d,%d] == %d AND reduce(m_(%d,%d)%s[%d, :]) ? {1@(%d,%d)} : ZERO )",
                    dst_mat_e.first, dst_mat_e.second,
                    GT.c_str(), c1, c2, edge_label,
                    m1.mat_e.first, m1.mat_e.second, MT.c_str(),
                    mc,
                    c1, c2);
            }
                break;
            case MQ_MUL:
            {
                string m1s = m1.stringify();
                string m2s = m2.stringify();
                string spaces = "";
                spaces.insert(0, min(30-m1s.size(), 100-m1s.size()), ' ');
                if(!trans_res) {
                    sprintf(buffer,
                        "m_(%d,%d) = %s %sx %s",
                        dst_mat_e.first, dst_mat_e.second,
                        m1s.c_str(),
                        spaces.c_str(),
                        m2s.c_str());
                } else {
                    sprintf(buffer,
                        "m_(%d,%d) = [%s %sx %s].T()",
                        dst_mat_e.first, dst_mat_e.second,
                        m1s.c_str(),
                        spaces.c_str(),
                        m2s.c_str());
                }
            }
                break;
        }
        os = buffer;
        return os;
    }
};





class QPlan
{
    QGraph              qgraph;
    vector<query_op_t>  op_vec;
    vector<qedge_t >    qwalk;
    vector<qedge_t >    eval_seq;
    bool                csc_friendly;
    int                 num_vars;
    
    
    /*
     *
     */
    void
    dfs_dwalk_to_eval_seq()
    {
        vector<qedge_t> res, t_vec;
        for(auto it = qwalk.begin(); it != qwalk.end(); ++it) {
            qedge_t e = *it;
            qedge_t e_r = make_pair(e.second, e.first);
            if(find(qwalk.begin(), it, e_r) == it) {
                t_vec.push_back(e);
            } else {
                t_vec.push_back(e_r);
            }
        }
        for(auto it = t_vec.begin(); it != t_vec.end(); ++it) {
            qedge_t e = *it;
            if(find(it+1, t_vec.end(), e) == t_vec.end()) {
                res.push_back(e);
            }
        }
        reverse(res.begin(), res.end());
        eval_seq = res;
    }
    
public:
    /*
     *
     */
    void
    print_qgraph()
    {
        qgraph.print();
    }
    
    /*
     *
     */
    void
    print_qwalk()
    {
       printf("DDFS walk:\n    ");
        for(auto p : qwalk) {
            printf("(%d,%d) ", p.first, p.second);
        }
        printf("\n");
    }
    
    /*
     *
     */
    void
    print_qplan()
    {
        printf("Query plan: ");
        if(csc_friendly)
            printf("(CSC friendly)\n");
        else
            printf("\n");
        for(auto op : op_vec) {
            string opstr = op.stringify();
            printf("%s\n", opstr.c_str());
        }
    }
    
    /*
     *
     */
    vector<query_op_t>
    get_op_vec()
    {
        return op_vec;
    }
    
    /*
     *
     */
    QGraph
    get_qgraph()
    {
        return qgraph;
    }
    
    /*
     *
     */
    int
    get_num_var()
    {
        return num_vars;
    }
    
    /*
     *
     */
    vector<qedge_t>
    get_eval_seq()
    {
        return eval_seq;
    }
    
    /*
     * 
     */
    bool
    is_csc_friend()
    {
        return csc_friendly;
    }
    
    /*
     *
     */
    vector<pair<qnode_t, qnode_t> >
    get_eq_var_pairs()
    {
        return qgraph.get_eq_var_pairs();
    }
    
    /*
     *
     */
    int
    get_num_actual_vars()
    {
        return qgraph.get_num_actual_vars();
    }
    
    /*
     *
     */
    string
    get_var_name
    (qnode_t v)
    {
        return qgraph.get_var_name(v);
    }
    
    /*
     *
     */
    vector<pair<string, qnode_t> >
    get_proj_vars() {
        return qgraph.get_proj_vars();
    }
    
    /*
     *
     */
    void
    make_sqplan
    (const Query& query)
    {
        op_vec.clear();
        csc_friendly = false;
        
        int start_node_mtree = 0;
        int start_node_ddfs  = 0;
        QGraph qg;
        qg.init(query);
        qg.make_tree(start_node_mtree);

        //start_node_ddfs controls the query plan
        qwalk = qg.extract_dfs_dwalk(start_node_ddfs);
        dfs_dwalk_to_eval_seq();
 
        // Each edge of qwalk produces an operation (mxm).
        // Some ops involve multiplication with the input graph, and some don't,
        // called corrections.
        // Handling the first edge in qwalk is special.
        set<qedge_t > done_edges_set;
        qedge_t e = NULLP;
        qedge_t prev_e = NULLP;
        bool e_dir_match = false;   // means edge direction in qwalk matches
                                    // that in edge_map (actual edge in qgraph)
        edge_label_t edge_scalar = -1;       // label of the edge in q graph

        // do the first edge
        e = qwalk[0];
        done_edges_set.insert(e);
        qg.get_edge_label(edge_scalar, e_dir_match, e);
        query_op_t op;
        op.type = MQ_MUL;
        op.dst_mat_e = e;
        op.dst_desc = MQ_CONSTRUCT_MAT;
        if(e_dir_match) {
            op.m1.set_vals(MQ_GRAPH, MQ_NOTRANSFORM, NULLP, -1, NULLP);
            op.m2.set_vals(MQ_IDENTITY, MQ_NOTRANSFORM, NULLP, edge_scalar, NULLP);
        } else {
            op.m1.set_vals(MQ_IDENTITY, MQ_NOTRANSFORM, NULLP, edge_scalar, NULLP);
            op.m2.set_vals(MQ_GRAPH, MQ_TRANSPOSE, NULLP, -1, NULLP);
        }
        op_vec.push_back(op);

        prev_e = e;
        for(int i = 1; i < qwalk.size(); ++i) {
            e = qwalk[i];
            qg.get_edge_label(edge_scalar, e_dir_match, e);
            qedge_t e_rev = make_pair(e.second, e.first);
            if(done_edges_set.find(e_rev) != done_edges_set.end()) { // correction operation
                query_op_t op;
                op.type = MQ_MUL;
                op.dst_desc = MQ_UPDATE_MAT;
                op.dst_mat_e = e_rev;
                op.m1.set_vals(MQ_STORED, MQ_NOTRANSFORM, NULLP, -1, e_rev);
                op.m2.set_vals(MQ_STORED, MQ_DIAGONALIZE, NULLP, 1, prev_e);
                op_vec.push_back(op);
                prev_e = e_rev;
                continue;
            }
            query_op_t op;
            op.type = MQ_MUL;
            op.dst_desc = MQ_CONSTRUCT_MAT;
            op.dst_mat_e = e;
            mat_transform_t transform;
            transform = ( prev_e.first == e.first ? MQ_DIAGONALIZE : MQ_TRANS_DIAG );
            op.m1.set_vals(MQ_STORED, transform, NULLP, edge_scalar, prev_e);
            transform = ( e_dir_match ? MQ_NOTRANSFORM : MQ_TRANSPOSE);
            op.m2.set_vals(MQ_GRAPH, transform, NULLP, -1, NULLP);
            op_vec.push_back(op);        
            prev_e = e;
            done_edges_set.insert(e);
        }

        // the previous part assumes all nodes in q graph are variables (?x)
        // handling constants is below
        op = op_vec[0];
        e = op.dst_mat_e;
        bool v1_is_c, v2_is_c;
        dnode_t c1, c2;
        qg.get_node_labels(v1_is_c, v2_is_c, c1, c2, e);
        qg.get_edge_label(edge_scalar, e_dir_match, e);
        if(v1_is_c && v2_is_c) {
            op.type = MQ_TWOCB;
            op.c1 = c1;
            op.c2 = c2;
            op.edge_label = edge_scalar;
            op.G_transpose = !e_dir_match;
        } else if(v1_is_c) {
            mat_transform_t transform;
            op.m1.set_vals(MQ_ZERO, MQ_NOTRANSFORM, make_pair(c1, c1), edge_scalar, NULLP);
            transform = ( e_dir_match ? MQ_NOTRANSFORM : MQ_TRANSPOSE );
            op.m2.set_vals(MQ_GRAPH, transform, NULLP, -1, NULLP);
        } else if(v2_is_c) {
            mat_transform_t transform;
            transform = ( e_dir_match ? MQ_NOTRANSFORM : MQ_TRANSPOSE );
            op.m1.set_vals(MQ_GRAPH, transform, NULLP, -1, NULLP);
            op.m2.set_vals(MQ_ZERO, MQ_NOTRANSFORM, make_pair(c2, c2), edge_scalar, NULLP);
        }
        op_vec[0] = op;
        for(int i = 1; i < op_vec.size(); ++i) {
            if(op_vec[i].dst_desc == MQ_UPDATE_MAT) {
                //fix op (back edge in qwalk), doesn't affect const handling  
                continue;
            }
            op = op_vec[i];
            e = op.dst_mat_e;
            qg.get_node_labels(v1_is_c, v2_is_c, c1, c2, e);
            qg.get_edge_label(edge_scalar, e_dir_match, e);
            if(v1_is_c && v2_is_c) {
                op.type = MQ_TWOCM;
                op.c1 = c1;
                op.c2 = c2;
                op.mc = c1;
                op.edge_label = edge_scalar;
                op.G_transpose = !e_dir_match;
                op_vec[i] = op;
            } else if(v2_is_c) {
                query_op_t new_op;
                new_op.type = MQ_MUL;
                new_op.dst_desc = MQ_UPDATE_MAT;
                new_op.dst_mat_e = op.dst_mat_e;
                new_op.m1.set_vals(MQ_STORED, MQ_NOTRANSFORM, NULLP, -1, op.dst_mat_e);
                new_op.m2.set_vals(MQ_ZERO, MQ_NOTRANSFORM, make_pair(c2, c2), 1, NULLP);
                op_vec.insert(op_vec.begin()+i+1, new_op);
            }
        }
        qgraph = qg;
        num_vars = qg.get_num_nodes();
        
        //======================================================================
        // EXPERIMENTAL - cycle fix
        auto eqp_vec = qg.get_eq_var_pairs();
        map<qnode_t, qnode_t> sn_map;
        for(auto p : eqp_vec) {sn_map[p.first] = p.second;}
        for(int i = 0; i < op_vec.size(); ++i) {
            op = op_vec[i];
            e = op.dst_mat_e;
            if(sn_map.find(e.second) != sn_map.end()) {
                // e is a back edge, part of a cycle
                // e.second is a leaf
                qnode_t sn = e.second;
                qnode_t rn = sn_map[sn];
                qedge_t col_sel_e;
                for(auto o : op_vec) {
                    if(o.dst_mat_e.first == rn) {
                        col_sel_e = o.dst_mat_e;
                        break;
                    }
                }
                query_op_t cycle_op;
                cycle_op.type = MQ_MUL;
                cycle_op.dst_desc = MQ_UPDATE_MAT;
                cycle_op.dst_mat_e = e;
                cycle_op.m1.set_vals
                    (MQ_STORED, MQ_NOTRANSFORM, NULLP, -1, cycle_op.dst_mat_e);
                cycle_op.m2.set_vals
                    (MQ_STORED, MQ_DIAGONALIZE, NULLP, 1, col_sel_e);
                op_vec.insert(op_vec.begin()+i+1, cycle_op);
                ++i;
                cycle_op.dst_mat_e = col_sel_e;
                cycle_op.m1.set_vals
                    (MQ_STORED, MQ_TRANS_DIAG, NULLP, 1, e);
                cycle_op.m2.set_vals
                    (MQ_STORED, MQ_NOTRANSFORM, NULLP, -1, cycle_op.dst_mat_e);
                op_vec.insert(op_vec.begin()+i+1, cycle_op);
                ++i;
            }
        }
        //======================================================================
        
        
    }
    
    /*
     *
     */
    void
    make_sqplan_csc_friendly
    (const Query& query)
    {
        make_sqplan(query);
        csc_friendly = true;
        for(auto& op : op_vec) {
            switch(op.type) {
                case MQ_TWOCB:
                {
                    op.dst_mat_e = make_pair(op.dst_mat_e.second, op.dst_mat_e.first);
                    auto tmp = op.c1;
                    op.c1 = op.c2;
                    op.c2 = tmp;
                    op.G_transpose = !op.G_transpose;
                }
                break;
                case MQ_TWOCM:
                {
                    op.dst_mat_e = make_pair(op.dst_mat_e.second, op.dst_mat_e.first);
                    auto tmp = op.c1;
                    op.c1 = op.c2;
                    op.c2 = tmp;
                    op.G_transpose = !op.G_transpose;
                    op.m1.mat_e = make_pair(op.m1.mat_e.second, op.m1.mat_e.first);
                    op.mc = op.c2;
                    switch(op.m1.mat_transform) {
                        case MQ_TRANSPOSE:
                            op.m1.mat_transform = MQ_NOTRANSFORM;
                            break;
                        case MQ_TRANS_DIAG:
                            op.m1.mat_transform = MQ_DIAGONALIZE;
                            break;
                        case MQ_NOTRANSFORM:
                            op.m1.mat_transform = MQ_TRANSPOSE;
                    }
                }
                break;
                case MQ_MUL:
                {
                    if(op.m1.mat_content == MQ_GRAPH && op.m2.mat_content == MQ_IDENTITY) {
                        op.trans_res = true;
                        op.dst_mat_e = make_pair(op.dst_mat_e.second, op.dst_mat_e.first);
                        break;
                    }
                    op.dst_mat_e = make_pair(op.dst_mat_e.second, op.dst_mat_e.first);
                    auto tmp = op.m1;
                    op.m1 = op.m2;
                    op.m2 = tmp;
                    op.m1.mat_e = make_pair(op.m1.mat_e.second, op.m1.mat_e.first);
                    op.m2.mat_e = make_pair(op.m2.mat_e.second, op.m2.mat_e.first);
                    switch(op.m1.mat_transform) {
                        case MQ_TRANSPOSE:
                            op.m1.mat_transform = MQ_NOTRANSFORM;
                            break;
                        case MQ_DIAGONALIZE:
                            op.m1.mat_transform = MQ_TRANS_DIAG;
                            break;
                        case MQ_TRANS_DIAG:
                            op.m1.mat_transform = MQ_DIAGONALIZE;
                            break;
                        case MQ_NOTRANSFORM:
                            if(op.m1.mat_e != op.dst_mat_e) // ==== EXP cycle fix
                                op.m1.mat_transform = MQ_TRANSPOSE;
                            break;
                    }
                    switch(op.m2.mat_transform) {
                        case MQ_TRANSPOSE:
                            if(op.m2.mat_content == MQ_GRAPH)
                                op.m2.mat_transform = MQ_NOTRANSFORM;
                            break;
                        case MQ_DIAGONALIZE:
                            op.m2.mat_transform = MQ_TRANS_DIAG;
                            break;
                        case MQ_TRANS_DIAG:
                            op.m2.mat_transform = MQ_DIAGONALIZE;
                            break;
                        case MQ_NOTRANSFORM:
                            if(op.m2.mat_content == MQ_GRAPH)
                                op.m2.mat_transform = MQ_TRANSPOSE;
                            break;
                    }
                    op.m1.one_location = make_pair(op.m1.one_location.second, op.m1.one_location.first);
                    op.m2.one_location = make_pair(op.m2.one_location.second, op.m2.one_location.first);
                }
                break;
            }
        }
        
    } 
};


#endif	/* MQ_QPLAN_H */

