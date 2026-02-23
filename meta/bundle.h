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
	
	int build_supplementaries(); 		// setting hit supple
	int set_chimeric_cigar_positions(); //setting h.first_pos/second_pos etc for getting back splice positions using cigars
	int build_circ_fragments();	//for building r2 and r1.supple frags
	int bridge_circ(); //bridging regular and circ fragments for circrna
};

#endif
