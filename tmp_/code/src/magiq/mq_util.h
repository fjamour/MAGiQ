/* 
 * File:   mq_util.h
 * Author: fuad
 *
 * Created on December 2, 2017, 10:40 AM
 */

#ifndef MQ_UTIL_H
#define	MQ_UTIL_H

#include "query_parser.h"
#include <string>
#include <regex>
#include <algorithm>
using namespace std;
/*
 *
 */
void
read_query_file
(
    Query& query,
    string& qstr,
    string qpath

)
{
    string queryFile = qpath;
    ifstream queryIn(queryFile.c_str());
    string querystring, str;
    if (!queryIn) {
        throwException("Query file '" + queryFile
				+ "' was not found. Please try again!\n");
    }
    // read the query string
    querystring = "";
    int counter = 0;
    while (true) {
    	getline(queryIn,str);
        //clean str from \t
        str.erase(remove(str.begin(), str.end(), '\t'), str.end());
        if(str!=query_delim){
            querystring+=str;
            if (!queryIn.good())
		break;
            querystring+='\n';
        }
    }
    query = parse_sparql_query(querystring);
    
    char buffer[5000];
    sprintf(buffer, "query string: %s\n", querystring.c_str());
    qstr += buffer;
    sprintf(buffer, "query.type: %d\n", query.type);
    qstr += buffer;
    sprintf(buffer, "query projs ");
    qstr += buffer;
    for (unsigned i = 0; i < query.projections.size(); i++) {
        sprintf(buffer, "%s, ", query.projections[i].c_str());
        qstr += buffer;
    }
    sprintf(buffer, "\n");
    qstr += buffer;
    sprintf(buffer, "query #triples %d\n", (int)query.nodes.size());
    qstr += buffer;
    sprintf(buffer, "query #variables %d\n", (int)query.variables.size());
    qstr += buffer;
    for (unsigned i = 0; i < query.nodes.size(); i++) {
        sprintf(buffer, "query tripe # %d ", i);
        qstr += buffer;
        for (unsigned j = 0; j < query.nodes.at(i).row.size(); j++) {
            sprintf(buffer, "%s, ", query.nodes.at(i).row.at(j).c_str());
            qstr += buffer;
        }
        sprintf(buffer, "\n");
        qstr += buffer;
    }
}

/*
 *
 */
string
get_qname
(string qpath)
{
    string  out;
    regex   re(".*/(.*)\\..*");
    smatch  match;
    regex_search(qpath, match, re);
    out = match.str(1);
    return out;
}

/*
 * given a graph path of the form PATH/dataset_name/encoded.nt
 * , return dataset_name
 */
string
get_dataset_prefix
(string path)
{
    string  out;
    regex   re(".*/(.*)/.*\\..*");
    smatch  match;
    regex_search(path, match, re);
    out = match.str(1);
    return out;
}


#endif	/* MQ_UTIL_H */

