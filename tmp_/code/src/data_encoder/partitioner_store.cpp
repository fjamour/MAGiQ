//============================================================================
// Name        : RDFDataEncoder
// Version     :
// Copyright   : KAUST-Infocloud
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "partitioner_store.h"
#include "utils.h"
#include "query.h"
#include <sys/types.h>





partitioner_store::partitioner_store() {
}

partitioner_store::~partitioner_store() {
	/*for (unsigned i = 0; i < rdf_data.size(); i++)
		delete rdf_data[i];*/
}

void partitioner_store::remove_double_edges(string input_file_name, string output_file_name)
{
    string subject,predicate,object,objectSubType;
    Type::ID objectType;
    ifstream fin(input_file_name.c_str());
    TurtleParser parser(fin);
    print_to_screen("Re-writing the encoded file with header and removing double edges: " + input_file_name);
    map<pair<string, string>, ll> edge_map;
    ll numLinesRead = 0;
    ll numRemovedEdges = 0;
    try {
            while (true) {
                    try {
                            if (!parser.parse(subject,predicate,object,objectType,objectSubType))
                                    break;
                    } catch (const TurtleParser::Exception& e) {
                            cerr << e.message << endl;
                            // recover...
                            while (fin.get()!='\n') ;
                            continue;
                    }
                    
                    ll pid = atoi(predicate.c_str());
                    pair<string, string> e = make_pair(subject, object);
                    if(edge_map.find(e) == edge_map.end()) {

                        edge_map[e] = pid;
                        
                       // cout<<"inserting edge first: s = "<<sid<<", o = "<<oid<<", preds = "<<pid<<endl;
                    }
                    else{
                        auto it = edge_map.find(e);
                        ll val = it->second;
//                        if(pid!=val)
                        {
                        
                            if(val > pid)
                            {
                                //cout<<"double edge exists: s = "<<subject<<", o = "<<object<<", preds = "<<val<<", "<<pid<<endl;
                                //edge_map.erase(it);
                                //cout<<"will use the pred: "<<pid<<endl;    
                                edge_map[e] = pid; 
                                numRemovedEdges++;
                            }
                        }

                    }
                    numLinesRead++;
                    if(numLinesRead%1000000==0)
                        cout<<"Finished Reading "<<numLinesRead<<endl;
            }
            
    }
    catch (const TurtleParser::Exception&) {
            return ;
    }
    fin.close();
    ofstream mapStream(output_file_name.c_str());
    
    mapStream<<edge_map.size()<<endl; //num_records
    mapStream<<so_map.size()<<endl;
    
    for(auto& p : edge_map) {
        pair<string, string> key = p.first;
        mapStream<<"<"<<key.first<<"> <"<<p.second<<"> <"<<key.second<<"> .\n";
    }   
    mapStream.close();
    std::remove(input_file_name.c_str()); // delete file
    cout<<"Number of removed edges: "<<numRemovedEdges<<endl;
    print_to_screen("Done with Re-writing the encoded file with header and no double edges");
}
void partitioner_store::load_encode_rdf_data(string input_dir, string output_file_name, bool fix_double_edges) {
	print_to_screen(part_string);
	Profiler profiler;
	profiler.startTimer("load_rdf_data");
	boost::unordered_map<string, ll>::iterator map_it;
	ll so_id = 0;
	ll predicate_id = 2;
	int64_t num_rec = 0;
	triple tmp_triple;
	ofstream ofs;
	string input_file;
	Type::ID objectType;
	ll sid, pid, oid;
	vector<triple> tmp_data;
	string subject,predicate,object,objectSubType;
	vector<string> files = vector<string>();

        //////
//        map<pair<int, int>, int> edge_map;
//        int numEdgesRemoved = 0;
        
        std::string out_file_name = string(output_file_name);
        std::string tmp_file_name = string(output_file_name+".tmp");
        
	getdir(input_dir,files);
	ofs.open(tmp_file_name.c_str(), ofstream::out);
	for (unsigned int file_id = 0; file_id < files.size();file_id++) {
		input_file = string(input_dir+"/"+files[file_id]);
		ifstream fin(input_file.c_str());
		TurtleParser parser(fin);
		print_to_screen("Reading triples from file: " + input_file);
		try {
			while (true) {
				try {
					if (!parser.parse(subject,predicate,object,objectType,objectSubType))
						break;
				} catch (const TurtleParser::Exception& e) {
					cerr << e.message << endl;
					// recover...
					while (fin.get()!='\n') ;
					continue;
				}

                                
				//lookup subject
				map_it = so_map.find(subject);
				if (map_it == so_map.end()) {
					so_map[subject] = so_id;
					sid = so_id;
					so_id++;
                                        //cout<<num_rec<<": increasing so_id subj: "<<so_id<<endl;
				}
				else{
					sid = map_it->second;
				}

				//lookup predicate
				map_it = predicate_map.find(predicate);
				if (map_it == predicate_map.end()) 
                                {
					predicate_map[predicate] = predicate_id;
					pid = predicate_id;
					predicate_id++;
				}
				else{
					pid = map_it->second;
				}

				//lookup object
				for(unsigned i = 0 ; i < object.length() ;i++){
					if(object[i] == '\n' || object[i] == '\r')
						object [i] = ' ';
				}
                              
                                //MAGIC: remove douple qoutes from object
                                //object = object.replace('"','');
                                if(objectType == 1)
                                {
                                    object.erase(remove(object.begin(), object.end(), '\"' ),  object.end());
                                    objectType = (Type::ID)0;//URI;
                                    //cout<<object<<endl;
                                }
                                
				map_it = so_map.find(object);
				if (map_it == so_map.end()) 
                                {
					/*if(objectType == 1){
                                            so_map["\""+object+"\""] = so_id; //ToDo: treat all as URIs
					}
					else
                                            so_map[object] = so_id;
                                        */
                                        //MAGIC: treat all objects as URIs
                                        so_map[object] = so_id;
                                    
					oid = so_id;
					so_id++;
                                       // cout<<num_rec<<": increasing so_id obj: "<<so_id<<endl;
				}
				else{
					oid = map_it->second;
				}
//
//
//                                ////////////////////////////
//                                pair<int, int> e = make_pair(sid, oid);
//                                if(edge_map.find(e) == edge_map.end()) {
//                                    //e = make_pair(sid, oid);
//                                    ///edge_map[e] = pid;
//                                    edge_map.insert(make_pair(e, pid)); 
//                                   // cout<<"inserting edge first: s = "<<sid<<", o = "<<oid<<", preds = "<<pid<<endl;
//                                }
//                                else{
//                                    auto it = edge_map.find(e);
//                                    auto val = it->second;
//                                        
//                                    if(pid!=val){
//                                        cout<<"double edge exists: s = "<<sid<<", o = "<<oid<<", preds = "<<val<<", "<<pid<<endl;
//                                        if(pid<val)
//                                        {
//                                            edge_map.erase(it);
//                                            edge_map.insert(make_pair(e, pid)); 
//                                        }
//                                        else{
//                                            numEdgesRemoved++;
//                                            continue;
//                                        }
//                                            
//                                    
//                                    }
//                                        
//                                }
                                ////////////////////////
                                
				tmp_triple = triple(sid, pid, oid, objectType);
				tmp_data.push_back(tmp_triple);
				num_rec++;
				if (num_rec % 1000000 == 0) {
					this->dump_encoded_data(ofs, tmp_data);
					tmp_data.clear();
					cout<<"Finished "<<num_rec<<endl;
//                                        cout<<"numEdgesRemoved: "<<numEdgesRemoved<<endl;
				}
			}
			this->dump_encoded_data(ofs, tmp_data);
			tmp_data.clear();
			fin.close();
		}catch (const TurtleParser::Exception&) {
			return ;
		}
	}
	ofs.close();


	total_data_size = num_rec;
	print_to_screen(part_string+"\nTotal number of triples: " + toString(this->total_data_size) + " records");
        print_to_screen("Number of nodes: "+toString(so_map.size()));
	profiler.pauseTimer("load_rdf_data");
	print_to_screen("Done with data encoding in " + toString(profiler.readPeriod("load_rdf_data")) + " sec");
	profiler.clearTimer("load_rdf_data");
	this->dump_dictionaries(output_file_name);
        
        cout<<so_id<<endl;
        if(fix_double_edges == true){
            remove_double_edges(tmp_file_name,out_file_name);//, num_rec, so_map.size());
        }
        else {
            //TODO HERE
            print_to_screen("Re-writing the encoded file with header ....");
            std::ifstream infile(tmp_file_name.c_str());
            ofs.open(out_file_name.c_str(), ofstream::out);
            ofs<<num_rec<<endl;
            ofs<<so_map.size()<<endl;
            std::string line;
            while (std::getline(infile, line))
            {
                ofs<<line<<endl;
            }
            ofs.close();
            std::remove(tmp_file_name.c_str()); // delete file
            print_to_screen("Done with Re-writing the encoded file with header");        }

}

void partitioner_store::dump_dictionaries(string file_name) {
	print_to_screen(part_string);
	Profiler profiler;
	profiler.startTimer("dump_dictionaries");
	print_to_screen("Dumping Dictionaries!");
	dump_map(so_map, file_name+"verts_map.dic", true);
	dump_map(predicate_map, file_name+"predicate_map.dic", true);
	profiler.pauseTimer("dump_dictionaries");
	print_to_screen("Done with dumping dictionaries in " + toString(profiler.readPeriod("dump_dictionaries")) + " sec");
	profiler.clearTimer("dump_dictionaries");

}

void partitioner_store::dump_encoded_data(ofstream &output_stream, vector<triple> & data){
	print_to_screen(part_string);
	Profiler profiler;
	profiler.startTimer("dump_encoded_data");

	for(unsigned i = 0 ; i < data.size(); i++){
		output_stream<<data[i].print()<<endl;
	}

	profiler.pauseTimer("dump_encoded_data");
	print_to_screen("Done with dump_encoded_data in "+toString(profiler.readPeriod("dump_encoded_data"))+" sec");
	profiler.clearTimer("dump_encoded_data");
}
