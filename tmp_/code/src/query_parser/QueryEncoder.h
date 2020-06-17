/* 
 * File:   QueryEncoder.h
 * Author: abdelai
 *
 * Created on November 13, 2017, 3:10 PM
 */

#ifndef QUERYENCODER_H
#define	QUERYENCODER_H

#include <cstdlib>
#include "utils.h"
#include <iostream>
#include <dirent.h>
#include <errno.h>
#include "SPARQLLexer.hpp"
#include "SPARQLParser.hpp"
#include "utils.h"
using namespace std;


int encoded_queries = 0;

void  encode_query(SPARQLParser & parser, ofstream &stream, boost::unordered_map<string, long long>& predicate_map, boost::unordered_map<string, long long>& subj_map, long long max_predicate, long long max_verts) {
	boost::unordered_map<string, long long>::iterator it;

	stream<<"SELECT ";
	for(unsigned i = 0 ; i < parser.projection.size() ;i++){
		stream<<"?"<<parser.getVariableName(parser.projection[i])<<" ";
	}
	stream<<"WHERE {"<<endl;

	for(unsigned i= 0; i < parser.patterns.patterns.size(); i++){
		if(parser.patterns.patterns[i].subject.type == SPARQLParser::Element::Variable){
			stream<<"\t?"<<parser.getVariableName(parser.patterns.patterns[i].subject.id)<<" ";
		}
		else if(parser.patterns.patterns[i].subject.type == SPARQLParser::Element::IRI){
			it = subj_map.find(parser.patterns.patterns[i].subject.value);
			if (it == subj_map.end()) {
				stream<<"\t<"<<max_verts<<"> ";
			}
			else
				stream<<"\t<"<<it->second<<"> ";
		}
		else{
			it = subj_map.find(parser.patterns.patterns[i].subject.value);
			if (it == subj_map.end()) {
				stream<<"\t\""<<max_verts<<"\" ";
			}
			else
				stream<<"\t\""<<it->second<<"\" ";
		}

		it = predicate_map.find(parser.patterns.patterns[i].predicate.value);
		if (it == predicate_map.end()) {
			stream<<"<"<<max_predicate<<"> ";
		}
		else
			stream<<"<"<<it->second<<"> ";

		if(parser.patterns.patterns[i].object.type == SPARQLParser::Element::Variable){
			stream<<"?"<<parser.getVariableName(parser.patterns.patterns[i].object.id)<<" ."<<endl;
		}
		else if(parser.patterns.patterns[i].object.type == SPARQLParser::Element::IRI){
			it = subj_map.find(parser.patterns.patterns[i].object.value);
			if (it == subj_map.end()) {
				stream<<"<"<<max_verts<<"> ."<<endl;
			}
			else
				stream<<"<"<<it->second<<"> ."<<endl;
		}
		else{
			it = subj_map.find(parser.patterns.patterns[i].object.value);
			if (it == subj_map.end()) {
				stream<<"\""<<max_verts<<"\" ."<<endl;
			}
			else
				stream<<"\""<<it->second<<"\" ."<<endl;
		}
	}
	stream<<"}\n#EOQ#"<<endl;
}

void load_encode_queries(string queryFile, boost::unordered_map<string, long long>& predicate_map, boost::unordered_map<string, long long>& subj_map, long long max_predicate, long long max_verts) {
	ifstream queryIn(queryFile.c_str());
	string outputFileName = queryFile+"_encoded";
	ofstream queryOut(outputFileName.c_str());
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
		if(str!=query_delim){
			querystring+=str;
			if (!queryIn.good())
				break;
			querystring+='\n';
		}
		else if ((str == query_delim || queryIn.eof()) && !querystring.empty()) {
			/*if(counter>2802102*/
			cout<<"query#"<<counter<<"\n"<<querystring<<endl;
			counter++;
			SPARQLLexer lexer(querystring);
			SPARQLParser parser(lexer);
			try {
				parser.parse();
			} catch (const SPARQLParser::ParserException& e) {
				cerr << "parse error: " << e.message << endl;
				return;
			}
			if(parser.variableCount != 0){
                            if(parser.patterns.patterns.size() > 0){
                                encode_query(parser, queryOut, predicate_map, subj_map, max_predicate, max_verts);
                                encoded_queries++;
                            }
			}
			querystring = "";
		}
	}
	queryIn.close();
	queryOut.close();
}


int main_q_encoder(int argc, const char** argv) {
	if (argc < 4) {
		throwException("Usage: QueryLoadEncoder <predicate_map_file> <subj_map_file> <query_files_directory>");
	}
	string line, directory = string(argv[3]);
	vector<string> files = vector<string>();
	boost::unordered_map<string, long long> verts_map, preds_map;
	long long id, max_predicate = numeric_limits<long long>::min(), max_verts = numeric_limits<long long>::min();
	char value[1000000];
	FILE * pFile;

	getdir(directory,files);
	cout<<"Will encode "<<(files.size())<<" query files."<<endl;

	//Reading predicates dictionary
	pFile = fopen(argv[1], "r");
	while(fscanf(pFile, "%lld %[^\r\n] ", &id, value)==2){
		if(id>max_predicate)
			max_predicate = id+1;
		preds_map[value] = id;
	}
	fclose(pFile);

	id = -1;
	value[0] = '\0';

	//Reading vertices dictionary
	pFile = fopen (argv[2],"r");
	while(fscanf(pFile, "%lld %[^\r\n] ", &id, value)==2){
		if(id>max_verts){
			max_verts = id+1;
		}
		verts_map[value] = id;
	}
	fclose(pFile);
	cout<<"Number of vertices: "<<max_verts<<", Number of Predicates: "<<max_predicate<<endl;
	for(unsigned q = 0 ; q < files.size(); q++){
		cout << "reading parsing and encoding queries from file: " << files[q] << endl;
		load_encode_queries(directory+files[q], preds_map, verts_map, max_predicate, max_verts);

		cout << "done" << endl;
		cout.flush();

	}
	
	
	return 0;

}

#endif	/* QUERYENCODER_H */

