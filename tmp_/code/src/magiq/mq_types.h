/* 
 * File:   mq_types.h
 * Author: fuad
 *
 * Created on December 2, 2017, 10:41 AM
 */

#ifndef MQ_TYPES_H
#define	MQ_TYPES_H

#include <utility>

using namespace std;

typedef int                     qnode_t;
typedef pair<qnode_t, qnode_t>  qedge_t;
typedef int                     edge_label_t;

typedef int                     dnode_t;
typedef pair<dnode_t, dnode_t>  dedge_t;

#define NULLP make_pair(-1, -1)



#endif	/* MQ_TYPES_H */

