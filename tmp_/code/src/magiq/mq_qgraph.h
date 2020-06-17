/* 
 * File:   mq_qgraph.h
 * Author: fuad
 *
 * Created on December 2, 2017, 10:29 AM
 */

#ifndef MQ_QGRAPH_H
#define	MQ_QGRAPH_H
#include <vector>
#include <map>
#include <string>
#include <tuple>
#include <stack>

#include "mq_types.h"

//Query parse includes
#include "query_parser.h"

//Global logger
#include "mq_simple_logger.h"

using namespace std;

/*
 * Query graph is a directed, labeled small graph.
 * Query graph could be cyclic (the undirected version of it)
 * IMP Different nodes cannot have the same label
 */
class QGraph
{
    // label to internal node index
    map<string, qnode_t>            node_label_map;
    
    // internal node index to label
    map<qnode_t, string>            node_label_map_r;
    
    // edge node pair to edge label
    map<qedge_t, edge_label_t>      edge_map;
    
    // directed adj list
    map<qnode_t, vector<qnode_t> >  d_adj_list;
    
    // undirected adj list
    map<qnode_t, vector<qnode_t> >  u_adj_list;
    
    // tree produced from query graph
    map<qnode_t, vector<qnode_t> >  t_adj_list;
    
    // to keep track of replicated vars
    // in cyclic query graphs
    vector<pair<qnode_t, qnode_t> > equal_vars_vec;                                      
    
    // projection variables
    vector<pair<string, qnode_t> >  proj_vars;

    /*
     *
     */
    void
    i_extract_dfs_dwalk_visitor
    (
        qnode_t             v,
        set<qnode_t>&       visited,
        vector<qedge_t >&   qwalk,
        bool                klbe = false // keep leaf back edge
    )
    {
        auto sorted_adj_list = t_adj_list[v];
        sort(sorted_adj_list.begin(), sorted_adj_list.end());
        for(auto nbr : sorted_adj_list) {
            if(visited.find(nbr) == visited.end()) {
                qwalk.push_back(make_pair(v, nbr));
                visited.insert(nbr);
                i_extract_dfs_dwalk_visitor(nbr, visited, qwalk, klbe);
                if(klbe){
                    qwalk.push_back(make_pair(nbr, v));
                } else {
                    if(t_adj_list[nbr].size() > 1) {
                        qwalk.push_back(make_pair(nbr, v));
                    }
                }

            }
        }
    }
    
public:
    /*
     * Returns the string label of a qnode (from the parsed query).
     */
    string
    get_var_name
    (qnode_t var_index)
    {
        return node_label_map_r[var_index];
    }
    
    /*
     * Returns the equal node indexes. Created node index first of each pair.
     */
    vector<pair<qnode_t, qnode_t> >
    get_equal_vars()
    {
        return equal_vars_vec;
    }
    
    /*
     * Returns the the number of variables of the parsed query.
     */
    int
    get_num_actual_vars()
    {
        return d_adj_list.size();
    }
    
    /*
     * Returns the number of nodes in the query tree.
     */
    int
    get_num_nodes()
    {
        return t_adj_list.size();
    }
    
    /*
     *
     */
    vector<pair<qnode_t, qnode_t> >
    get_eq_var_pairs()
    {
        return equal_vars_vec;
    }
    
    /*
     * Fills qgraph structures from a parsed query.
     */
    void
    init
    (const Query& query)
    {    
        // produce triples
        typedef tuple<string, string, string> triple_t;
        vector<triple_t> triples_vec;
        for (int i = 0; i < query.nodes.size(); i++) {
            string s, p, o;
            s = query.nodes.at(i).row.at(0);
            p = query.nodes.at(i).row.at(1);
            o = query.nodes.at(i).row.at(2);
            triples_vec.push_back(make_tuple(s, o, p));
        }
        for(int i = 0; i < triples_vec.size(); ++i) {
            string s, o;
            s = get<0>(triples_vec[i]);
            o = get<1>(triples_vec[i]);
            for(string ss : {s, o}) {
                 node_label_map[ss] = -1;
            }
        }
        int var_i = 0;
        for(auto& p : node_label_map) {
            p.second = var_i++;
            node_label_map_r[p.second] = p.first;
            d_adj_list[p.second] = vector<int>();
            u_adj_list[p.second] = vector<int>();
        }
        for(int i = 0; i < triples_vec.size(); ++i) {
            string s, o, p;
            s = get<0>(triples_vec[i]);
            o = get<1>(triples_vec[i]);
            p = get<2>(triples_vec[i]);
            int v1, v2;
            v1 = node_label_map[s];
            v2 = node_label_map[o];
            // edge map
            edge_map[make_pair(v1, v2)] = stoi(p);
            // adj lists
            d_adj_list[v1].push_back(v2);
            u_adj_list[v1].push_back(v2);
            u_adj_list[v2].push_back(v1);
        }
        
        // fill projection variables
        if(query.projections.size() > 0) {
            if(query.projections[0] == "*") {
                for(auto p : node_label_map) {
                    proj_vars.push_back(p);
                }
            } else {
                for(auto v : query.projections) {
                    proj_vars.push_back(make_pair(v, node_label_map[v]));
                }
            }
        }
    }
    
    /*
     *
     */
    void
    print()
    {
        printf("Node to qnode:\n    ");
        for(auto p : node_label_map) {
            printf("(%s->%d) ", p.first.c_str(), p.second);
        }
        printf("\n");
        if(equal_vars_vec.size() > 0) {
            printf("Equal pairs:\n    ");
            for(auto p : equal_vars_vec) {
                printf("(%d->%d) ", p.first, p.second);
            }
        } else {
            printf("No cycles");
        }
        printf("\n");
        printf("Edges:\n    ");
        int cnt = 1;
        for(auto p : edge_map) {
            if(cnt++%6 == 0) printf("\n    ");
            printf("{(%d,%d)->%d} ", p.first.first, p.first.second, p.second);
        }
        printf("\n");
        printf("Directed adjacency list:\n    ");
        for(auto v : d_adj_list) {
            printf("%d(", v.first);
            for(int i = 0; i < v.second.size(); ++i) {
                auto nbr = v.second[i];
                printf("%d", nbr);
                if(i != v.second.size()-1)
                    printf(",");
            }
            printf(")  ");
        }
        printf("\n");
        printf("Tree adjacency list:\n    ");
        for(auto v : t_adj_list) {
            printf("%d(", v.first);
            for(int i = 0; i < v.second.size(); ++i) {
                auto nbr = v.second[i];
                printf("%d", nbr);
                if(i != v.second.size()-1)
                    printf(",");
            }
            printf(")  ");
        }
        printf("\n");
    }
    
    /*
     * Checks if undirected version of graph has at least one cycle
     */
    bool
    is_tree()
    {
        // use the graph representation in u_adj_list to find cycles
        map<int, bool> visited;
        for(auto p : u_adj_list) {
            visited[p.first] = false;
        }
        int s = 0; //start bfs from node 0
        queue<int> Q;
        visited[s] = true;
        Q.push(s);
        while(!Q.empty()) {
            int v = Q.front(); Q.pop();
            int num_visited = 0;
            for(auto nbr : u_adj_list[v]) {
                if(!visited[nbr]) {
                    Q.push(nbr);
                    visited[nbr] = true;
                } else {
                    ++num_visited;
                }
            }
            if(num_visited > 1)
                return false;
        }
        return true;
    }
    
    /*
     * Updates internal structures to convert query graph to tree
     * , cycles are broken by adding a shadow node
     * fills: t_adj_list, equal_vars_vec
     * updates: node_label_map, node_label_map_r
     */
    void
    make_tree(int start_node = 0)
    {
        if(start_node >= u_adj_list.size()) start_node = 0;
        if(is_tree()) {
            t_adj_list = u_adj_list;
            return;
        }
        for(auto p : u_adj_list) {
            t_adj_list[p.first] = vector<int>();
        }

        set<int> visited;
        map<int, int> parent;
        parent[start_node] = -1;
        equal_vars_vec.clear();
        int next_node = u_adj_list.size();
        stack<int> S;
        S.push(start_node);
        while(!S.empty()) {
            int v = S.top(); S.pop();
            if(visited.find(v) == visited.end()) {
                visited.insert(v);
                vector<int> nbr_vec = u_adj_list[v];
                sort(nbr_vec.begin(), nbr_vec.end());
                reverse(nbr_vec.begin(), nbr_vec.end());
                for(auto nbr : nbr_vec) {
                    if(visited.find(nbr) == visited.end()) {
                        parent[nbr] = v;
                        S.push(nbr);
                    } else if(nbr != parent[v]) {
                        // add edge
                        t_adj_list[v].push_back(next_node);
                        t_adj_list[next_node].push_back(v);
                        equal_vars_vec.push_back(make_pair(next_node, nbr));
                        // fix corresponding edge in edge_map
                        pair<int, int> e(v, nbr);
                        pair<int, int> ne(v, next_node);
                        if(edge_map.find(e) == edge_map.end()) {
                            e = make_pair(nbr, v);
                            ne = make_pair(next_node, v);
                        }
                        auto it = edge_map.find(e);
                        auto val = it->second;
                        edge_map.erase(it);
                        edge_map.insert(make_pair(ne, val));
                        ++next_node;
                    }
                }
            }
        }
        for(auto p : parent) {
            int k = p.first; int v = p.second;
            if(v != -1) {
                t_adj_list[v].push_back(k);
                t_adj_list[k].push_back(v);
            }
        }
        // make necessary corrections
        // fix node_label_map_r (add mappings of new nodes)
        for(auto p : equal_vars_vec) {
            int new_var_i = p.first;
            int old_var_i = p.second;
            node_label_map_r[new_var_i] = node_label_map_r[old_var_i];
        }
    }
    
    /*
     * TODO describe this precisely
     * Extracts a dfs walk where an edge appears twice if it's not
     * an edge of a leaf. the 'back edges' are there to do corrections.
     */
    vector<qedge_t >
    extract_dfs_dwalk
    (qnode_t start_node = 0)
    {
        if(start_node >= t_adj_list.size()) start_node = 0;
        vector<qedge_t > qwalk;
        set<qnode_t> visited;
        visited.insert(start_node);
        i_extract_dfs_dwalk_visitor(start_node, visited, qwalk);
        return qwalk;
    }
    /*
     * Same as extract_dfs_dwalk but keeps leaf back edge
     */
    vector<qedge_t >
    extract_dfs_dwalk_klbe
    (qnode_t start_node = 0)
    {
        if(start_node >= t_adj_list.size()) start_node = 0;
        vector<qedge_t > qwalk;
        set<qnode_t> visited;
        visited.insert(start_node);
        bool klbe = true; // keep leaf back edge
        i_extract_dfs_dwalk_visitor(start_node, visited, qwalk, klbe);
        return qwalk;
    }
    
    /*
     * Returns integral value for edge label e from edge map.
     * If the actual edge direction doesn't match direction with e,
     * output.second is set to false
     */
    void
    get_edge_label
    (
        edge_label_t&  label,
        bool&           dir_match,
        qedge_t         e
    )
    {
        if(edge_map.find(e) == edge_map.end()) {
            dir_match = false;
            pair<int, int> ee = make_pair(e.second, e.first);
            label = edge_map[ee];
        } else {
            dir_match = true;
            label = edge_map[e];
        }
    }
    
    /*
     * If e.first corresponds to a constant returns the integral label (node id)
     * in c1, and true in v1_is_c.. same for 2
     */
    void
    get_node_labels
    (
        bool&   v1_is_c,
        bool&   v2_is_c,
        int&    c1,
        int&    c2,
        qedge_t e)
    {
        string s1, s2;
        s1 = node_label_map_r[e.first];
        s2 = node_label_map_r[e.second];
        if(s1[0] == '?') {
            v1_is_c = false;
            c1 = -1;
        } else {
            v1_is_c = true;
            c1 = stoi(s1);
        }
        if(s2[0] == '?') {
            v2_is_c = false;
            c2 = -1;
        } else {
            v2_is_c = true;
            c2 = stoi(s2);
        }
    }
    
    /*
     *
     */
    vector<pair<string, qnode_t> >
    get_proj_vars() {
        return proj_vars;
    }
    
    /*
     * Returns the predicates involved in q
     */
    vector<edge_label_t>
    get_predicates()
    {
        vector<edge_label_t> qpred;
        for(auto it : edge_map) {
            qpred.push_back(it.second);
        }
        sort(qpred.begin(), qpred.end());
        qpred.erase(unique(qpred.begin(), qpred.end()), qpred.end());
        return qpred;
    }
};


#endif	/* MQ_QGRAPH_H */

