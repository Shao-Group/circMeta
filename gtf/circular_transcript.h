/*
(c) 2023 by Tasfia Zahin, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __GTF_CIRCULAR_TRANSCRIPT_H__
#define __GTF_CIRCULAR_TRANSCRIPT_H__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>

using namespace std;

class circular_transcript
{
public:
    circular_transcript();
    int write(ostream &fout) const;
    int print(int index);
    ~circular_transcript();
public:
    int sid;
    
	string id;
    string chrm;
	string source;
    string feature;

    int32_t start;
	int32_t end;
	double score;
    char strand;
    int frame;

	string gene_id;
	string transcript_id;
    int coverage;

    int exon_count;
    vector<int32_t> intron_chain;
    
    //for candidate transcripts whose bridging will be decided in meta assembly
    int frag_bridged_type; // 0 default, 1 if both circ and reg bridged, 2 if circ frag bridged but reg unbridged, 3 if reg bridged but circ frag unbridged
    bool h1_suppl;
    vector<vector<int32_t>> intermediate_chains;
};

#endif
