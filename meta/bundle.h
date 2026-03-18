/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __COMBINED_BUNDLE_H__
#define __COMBINED_BUNDLE_H__

#include <vector>
#include <stdint.h>
#include <iostream>
#include <string>

#include "parameters.h"
#include "bundle_base.h"
#include "sample_profile.h"
#include "circular_transcript.h"
#include "region.h"
#include "htslib/faidx.h"
#include "graph_builder.h"
#include <mutex>

using namespace std;

class bundle : public bundle_base
{
public:
	bundle(const parameters &cfg, const sample_profile &sp);
	bundle(const parameters &cfg, const sample_profile &sp, bundle_base &&bb);
	bundle(const bundle &cb) = default;
	bundle(bundle &&cb) = default;

public:
	string gid;
	const parameters &cfg;
	const sample_profile &sp;
	int num_combined;

	vector<circular_transcript> circ_trsts; //storing circular transcripts
	map<string, circular_transcript> unbridged_candidate_trsts; //storing partially unbridged circs

public:
	int set_gid(int instance, int subindex);
	int set_gid(int rid, int gid, int instance, int subindex);
	int copy_meta_information(const bundle &bb);
	int combine(const bundle &bb, bool combine_map);
	int print(int k);
	int bridge();
	
	int set_supplementaries();
	int build_supplementaries(); 		// setting hit supple
	int get_more_chimeric_reads(faidx_t *fai);	// for getting more chimeric reads for junction mapping and creating fake supplementaries
	int build_fake_circ_fragments();
	int set_chimeric_cigar_positions(); //setting h.first_pos/second_pos etc for getting back splice positions using cigars
	int build_circ_fragments();	//for building r2 and r1.supple frags
	int fix_alignment_boundaries();
	int bridge_circ(); //bridging regular and circ fragments for circrna
	int bridge_circ_optimized();
	int create_fake_supple(int fr_index, int32_t soft_len, int32_t pos1, int32_t pos2, int soft_clip_side);
	string get_fasta_seq(faidx_t *fai, int32_t pos1, int32_t pos2);
	double get_Jaccard(string s, string t);
	double get_Jaccard_correct(string s, string t);
};

#endif
