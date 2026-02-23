/*
(c) 2023 by Tasfia Zahin, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <map>
#include "circular_transcript.h"

circular_transcript::circular_transcript()
{
    sid = -1;
    
    id = "";
    chrm = "";
    source = "";
    feature = "";

    start = 0;
	end = 0;
    score = 0;
    strand = '.';
    frame = -1;

    gene_id = "";
    transcript_id = "";
    coverage = 0;

    exon_count = 0;
    intron_chain.clear();

    frag_bridged_type = 0;
    h1_suppl = true;
    intermediate_chains.clear();
}

int circular_transcript::write(ostream &fout) const
{
    fout.precision(4);
	fout<<fixed;
    
    fout<<chrm.c_str()<<"\t";
    fout<<source.c_str()<<"\t";
    fout<<feature.c_str()<<"\t";
    fout<<start + 1<<"\t";
    fout<<end<<"\t";
	fout<<score<<"\t";							
	fout<<strand<<"\t";							
	fout<<".\t";								            
	fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
	fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
    fout<<"cov \""<<coverage<<"\";"<<endl;
    
    vector<int32_t> exon_chain;
    exon_chain.push_back(start);
    exon_chain.insert(exon_chain.end(), intron_chain.begin(), intron_chain.end());
    exon_chain.push_back(end);

    // if( intron_chain.size() > 0 && intron_chain[0] == 152305224 && intron_chain[intron_chain.size()-1] == 152313562)
    // {
    //     printf("got full chain, exon count %d, read name:%s\n",exon_count,transcript_id);
    // }

    int cnt = 0;
    for(int i=0;i<exon_chain.size();i+=2)
    {
        fout<<chrm.c_str()<<"\t";
        fout<<source.c_str()<<"\t";
        fout<<"exon\t";

        fout<<exon_chain[i] + 1<<"\t";           
        fout<<exon_chain[i+1]<<"\t";

        fout<<score<<"\t";							
        fout<<strand<<"\t";							
        fout<<".\t";								           
        fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
        fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
        fout<<"exon \""<<++cnt<<"\"; "<<endl;
    }

    return 0;
}

int circular_transcript::print(int index)
{
    printf("circRNA_id:%s, chrm:%s, start:%d, end:%d, score:%d, strand:%c\n",id.c_str(), chrm.c_str(), start, end, score, strand);
    return 0;

}

circular_transcript::~circular_transcript()
{
}

