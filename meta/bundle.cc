/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "bundle.h"
#include "config.h"
#include "essential.h"
#include "graph_builder.h"
#include "graph_cluster.h"
#include "bridge_solver.h"

#include <sstream>
#include <algorithm>

bundle::bundle(const parameters &c, const sample_profile &s)
	: cfg(c), sp(s)
{
	num_combined = 0;
	circ_trsts.clear();
}

bundle::bundle(const parameters &c, const sample_profile &s, bundle_base &&bb)
	: cfg(c), sp(s), bundle_base(bb)
{
	num_combined = 0;
	circ_trsts.clear();
}

int bundle::set_gid(int instance, int subindex)
{
	char name[10240];
	sprintf(name, "instance.%d.%d", instance, subindex);
	gid = name;
	return 0;
}

int bundle::set_gid(int rid, int g, int instance, int subindex)
{
	char name[10240];
	sprintf(name, "instance.%d.%d.%d.%d", rid, g, instance, subindex);
	gid = name;
	return 0;
}

int bundle::copy_meta_information(const bundle &bb)
{
	chrm = bb.chrm;
	strand = bb.strand;
	tid = bb.tid;
	lpos = bb.lpos;
	rpos = bb.rpos;
	return 0;
}

int bundle::set_chimeric_cigar_positions()
{
	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];
		if(h.suppl == NULL) continue;

		int32_t p;
		int32_t q;

		p = h.pos;
		int match_index = 0;;

		for(int j=0;j<h.n_cigar;j++)
		{
			if(h.cigar_vector[j].first == 'M')
			{
				match_index = j;
				break;
			}
		}
		for(int j=match_index-1;j>=0;j--)
		{
			p -= h.cigar_vector[j].second; //subtracting cigars before match index
		}
		
		q = h.suppl->pos;
		match_index = 0;;

		for(int j=0;j<h.suppl->n_cigar;j++)
		{
			if(h.suppl->cigar_vector[j].first == 'M')
			{
				match_index = j;
				break;
			}
		}
		for(int j=match_index-1;j>=0;j--)
		{
			q -= h.suppl->cigar_vector[j].second; //subtracting cigars before match index
		}
		

		int32_t x = p;
		int32_t diff_cigar1 = 10000;
		int32_t diff_cigar2 = 10000;
		int best_pos_flag = 0;

		for(int j=0;j<h.n_cigar-1;j++)
		{
			int32_t y = q;
			pair<char, int32_t> hp_cigar1 = h.cigar_vector[j];
			pair<char, int32_t> hp_cigar2 = h.cigar_vector[j+1];

			for(int k=0;k<h.suppl->n_cigar-1;k++)
			{

				pair<char, int32_t> hs_cigar1 = h.suppl->cigar_vector[k];
				pair<char, int32_t> hs_cigar2 = h.suppl->cigar_vector[k+1];

				//printf("check %d,%c\n",hp_cigar1.second,hp_cigar1.first);

				int32_t hp_cigar1_len = hp_cigar1.second;
				int32_t hp_cigar2_len = hp_cigar2.second;
				int32_t hs_cigar1_len = hs_cigar1.second;
				int32_t hs_cigar2_len = hs_cigar2.second;

				if(((hp_cigar1.first == 'S' || hp_cigar1.first == 'H') && hp_cigar2.first == 'M') && (hs_cigar1.first == 'M' && (hs_cigar2.first == 'S' || hs_cigar2.first == 'H')))
				{
					if(hp_cigar2.first == 'M' && j+3 < h.cigar_vector.size() && h.cigar_vector[j+3].first == 'M')
					{
						//&& h.cigar_vector[j+2].first == 'N' can be I or D as well
						for(int p=j+3;p<h.cigar_vector.size();p+=2) //traverse all Ms with 1 gap in the middle ex:SMNMNM/SIMIMIM/SMDMDMD
						{
							if(h.cigar_vector[p].first == 'M')
							{
								hp_cigar2_len += h.cigar_vector[p].second;
							}
						}			
					}

					if(hs_cigar1.first == 'M' && k-2 >= 0 && h.suppl->cigar_vector[k-2].first == 'M')
					{
						//&& h.suppl->cigar_vector[k-1].first == 'N' can be I or D as well
						for(int p=k-2;p>=0;p-=2)
						{
							if(h.suppl->cigar_vector[p].first == 'M')
							{
								hs_cigar1_len += h.suppl->cigar_vector[p].second;
							}
						}
					}
					
					diff_cigar1 = abs(hs_cigar1_len - hp_cigar1_len);
					diff_cigar2 = abs(hs_cigar2_len - hp_cigar2_len);
				
				}
				else if(((hs_cigar1.first == 'S' || hs_cigar1.first == 'H') && hs_cigar2.first == 'M') && (hp_cigar1.first == 'M' && (hp_cigar2.first == 'S' || hp_cigar2.first == 'H')))
				{
					if(hs_cigar2.first == 'M' && k+3 < h.suppl->cigar_vector.size() && h.suppl->cigar_vector[k+3].first == 'M')
					{
						//&& h.suppl->cigar_vector[k+2].first == 'N' can be I or D as well
						for(int p=k+3;p<h.suppl->cigar_vector.size();p+=2)
						{
							if(h.suppl->cigar_vector[p].first == 'M')
							{
								hs_cigar2_len += h.suppl->cigar_vector[p].second;
							}
						}			
					}

					if(hp_cigar1.first == 'M' && j-2 >=0 && h.cigar_vector[j-2].first == 'M')
					{
						//&& h.cigar_vector[j-1].first == 'N' can be I or D as well
						for(int p=j-2;p>=0;p-=2)
						{
							if(h.cigar_vector[p].first == 'M')
							{
								hp_cigar1_len += h.cigar_vector[p].second;
							}
						}
					}

					diff_cigar1 = abs(hs_cigar1_len - hp_cigar1_len);
					diff_cigar2 = abs(hs_cigar2_len - hp_cigar2_len);					
				}

				if(diff_cigar1 < 20 && diff_cigar2 < 20) //setting diff max 20 between complementing cigars
				{
					//printf("Found best positions\n");
					best_pos_flag = 1;

					//set first,second and third positions of prim and suppl
					h.first_pos = x;
					h.second_pos = x + hp_cigar1.second;
					h.third_pos = x + hp_cigar1.second + hp_cigar2.second;

					h.suppl->first_pos = y;
					h.suppl->second_pos = y + hs_cigar1.second;
					h.suppl->third_pos = y + hs_cigar1.second + hs_cigar2.second;

					h.left_cigar = hp_cigar1.first;
					h.left_cigar_len = hp_cigar1_len;
					h.right_cigar = hp_cigar2.first;
					h.right_cigar_len = hp_cigar2_len;

					h.suppl->left_cigar = hs_cigar1.first;
					h.suppl->left_cigar_len = hs_cigar1_len;
					h.suppl->right_cigar = hs_cigar2.first;
					h.suppl->right_cigar_len = hs_cigar2_len;

					break;

				}
				y += hs_cigar1.second;

			}
			if(best_pos_flag == 1) break;

			x += hp_cigar1.second;
		}
	}

	return 0;
}

int bundle::build_supplementaries()
{
	int max_index = hits.size() + 1;
	if(max_index > 1000000) max_index = 1000000;

	vector< vector<int> > vv;
	vv.resize(max_index);

	printf("hits size:%d\n",hits.size());

	// first build index
	for(int i = 0; i < hits.size(); i++)
	{
		const hit &h = hits[i];

		if(h.hid < 0) continue;
		if((h.flag & 0x800) == 0) continue; //skip non supplementary

		int k = (h.get_qhash() % max_index + (h.flag & 0x40) + (h.flag & 0x80)) % max_index;
		vv[k].push_back(i);
	}

	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];

		if(h.hid < 0) continue;
		if((h.flag & 0x800) >= 1) continue;  // skip supplemetary

		int k = (h.get_qhash() % max_index + (h.flag & 0x40) + (h.flag & 0x80)) % max_index;

		for(int j = 0; j < vv[k].size(); j++)
        {
            int u = vv[k][j];
			hit &z = hits[u];
			
            if(z.qname != h.qname) continue;
			if(z.mpos != h.mpos) continue; //primary and supple both have same mpos (mate pos)

            // TODO check 0x40 and 0x80 are the same for primary and supple
            if(((z.flag & 0x40) != (h.flag & 0x40)) || ((z.flag & 0x80) != (h.flag & 0x80))) continue;
			
        	h.suppl = &z;
			h.suppl_index = u;
        	break; //Taking the first supplementary read
        }
	}
	return 0;
}

int bundle::build_circ_fragments()
{
	for(int i = 0; i < frgs.size(); i++)
	{
		int idx1 = frgs[i][0];
		int idx2 = frgs[i][1];
		hit &h1 = hits[idx1];
		hit &h2 = hits[idx2];

		if(h1.suppl != NULL)
		{
			hit *h1_supple = h1.suppl;

			//use |prim.S/H + suppl.S/H - read-length| <= a threshold as a criteria for discarding cases
			int32_t len_HS = 0;

			if(h1.left_cigar == 'H' || h1.left_cigar == 'S')
			{
				len_HS += h1.left_cigar_len;
			}
			else if(h1.right_cigar == 'H' || h1.right_cigar == 'S')
			{
				len_HS += h1.right_cigar_len;
			}
			if(h1.suppl->left_cigar == 'H' || h1.suppl->left_cigar == 'S')
			{
				len_HS += h1.suppl->left_cigar_len;
			}
			else if(h1.suppl->right_cigar == 'H' || h1.suppl->right_cigar == 'S')
			{
				len_HS += h1.suppl->right_cigar_len;
			}

			// printf("len_HS = %d, read_length = %d\n",len_HS, cfg.read_length);

			if(abs(len_HS - cfg.read_length) > 5)
			{
				printf("read length criteria unsatisfied h1s.\n");
				continue;
			}

			// if(h1.pos < h2.pos && h2.pos < h1.suppl->pos)
			if(h1.second_pos <= h2.pos && h1.second_pos <= h1_supple->pos && h1_supple->second_pos >= h1.rpos && h1_supple->second_pos >= h2.rpos)
			{
				circ_frgs.push_back(AI3({idx2, h1.suppl_index, 0}));
			}
		}

		if(h2.suppl != NULL)
		{
			hit *h2_supple = h2.suppl;

			if(h2.qname == "E00512:127:HJNF3ALXX:3:1105:24211:69309")
			{
				printf("meta circRNA candidate read entered build circ frags:%s\n",h2.qname.c_str());
			}

			//use |prim.S/H + suppl.S/H - read-length| <= a threshold as a criteria for discarding cases
			int32_t len_HS = 0;

			if(h2.left_cigar == 'H' || h2.left_cigar == 'S')
			{
				len_HS += h2.left_cigar_len;
			}
			else if(h2.right_cigar == 'H' || h2.right_cigar == 'S')
			{
				len_HS += h2.right_cigar_len;
			}
			if(h2.suppl->left_cigar == 'H' || h2.suppl->left_cigar == 'S')
			{
				len_HS += h2.suppl->left_cigar_len;
			}
			else if(h2.suppl->right_cigar == 'H' || h2.suppl->right_cigar == 'S')
			{
				len_HS += h2.suppl->right_cigar_len;
			}

			if(abs(len_HS - cfg.read_length) > 5) //here 100 is the estimated read length, replace this with any related exisiting parameter
			{
				printf("read length criteria unsatisfied h2s.\n");
				continue;
			}

			if(h2.qname == "E00512:127:HJNF3ALXX:3:1105:24211:69309")
			{
				printf("h2.second pos:%d, h1.rpos:%d\n", h2.second_pos, h1.rpos);
			}

			// if(h2.suppl->pos < h1.pos && h1.pos < h2.pos)
			if(h2_supple->second_pos <= h1.pos && h2_supple->second_pos <= h2.pos && h2.second_pos >= h2_supple->rpos && h2.second_pos >= h1.rpos)
			{
				if(h2.qname == "E00512:127:HJNF3ALXX:3:1105:24211:69309")
				{
					printf("meta circRNA candidate read pushed in circ frags:%s\n",h2.qname.c_str());
				}
				circ_frgs.push_back(AI3({h2.suppl_index, idx1, 0}));
			}
		}
	}

	return 0;
}

int bundle::fix_alignment_boundaries()
{
	for(int k = 0; k < frgs.size(); k++)
	{
		int idx1 = frgs[k][0];
		int idx2 = frgs[k][1];
		hit &h1 = hits[idx1];
		hit &h2 = hits[idx2];

		if(h1.suppl != NULL)
		{
			hit *h1_supple = h1.suppl;

			//fixing boundaries for h1 and its suppl
			if(h1.pos > h2.pos && h1.pos <= h1_supple->pos && h1_supple->rpos >= h1.rpos && h1_supple->rpos >= h2.rpos && h1.pos - h2.pos <= cfg.alignment_boundary_error)
			{
				int32_t diff = h1.pos - h2.pos;

				int first_M_len = 0;
				for(int i=0;i<h2.cigar_vector.size();i++)
				{
					if(h2.cigar_vector[i].first == 'M')
					{
						first_M_len = h2.cigar_vector[i].second;
						break;
					}
				}

				if(first_M_len > cfg.alignment_boundary_error)
				{
					h2.pos = h2.pos + diff;
				}
			}
			else if(h1.pos <= h2.pos && h1.pos <= h1_supple->pos && h1_supple->rpos >= h1.rpos && h1_supple->rpos < h2.rpos && h2.rpos - h1_supple->rpos <= cfg.alignment_boundary_error)
			{
				int32_t diff = h2.rpos - h1_supple->rpos;

				int last_M_len = 0;
				for(int i=0;i<h2.cigar_vector.size();i++)
				{
					if(h2.cigar_vector[i].first == 'M')
					{
						last_M_len = h2.cigar_vector[i].second;
					}
				}

				if(last_M_len > cfg.alignment_boundary_error)
				{
					h2.rpos = h2.rpos - diff;
				}
			}
		}

		if(h2.suppl != NULL)
		{
			hit *h2_supple = h2.suppl;

			//fixing boundaries for h2 and its suppl
			if(h2_supple->pos <= h1.pos && h2_supple->pos <= h2.pos && h2.rpos >= h2_supple->rpos && h2.rpos < h1.rpos && h1.rpos - h2.rpos <= cfg.alignment_boundary_error)
			{
				int32_t diff = h1.rpos - h2.rpos;

				int last_M_len = 0;
				for(int i=0;i<h1.cigar_vector.size();i++)
				{
					if(h1.cigar_vector[i].first == 'M')
					{
						last_M_len = h1.cigar_vector[i].second;
					}
				}

				if(last_M_len > cfg.alignment_boundary_error)
				{
					h1.rpos = h1.rpos - diff;
				}
			}
			else if(h2_supple->pos > h1.pos && h2_supple->pos <= h2.pos && h2.rpos >= h2_supple->rpos && h2.rpos >= h1.rpos && h2_supple->pos - h1.pos <= cfg.alignment_boundary_error)
			{
				int32_t diff = h2_supple->pos - h1.pos;

				int first_M_len = 0;
				for(int i=0;i<h1.cigar_vector.size();i++)
				{
					if(h1.cigar_vector[i].first == 'M')
					{
						first_M_len = h1.cigar_vector[i].second;
						break;
					}
				}

				if(first_M_len > cfg.alignment_boundary_error)
				{
					h1.pos = h1.pos + diff;
				}
			}
		}
	} 
	return 0;
}

int bundle::bridge_circ_optimized()
{
	splice_graph gr;
	graph_builder gb(*this, cfg, sp);
	gb.build(gr);
	gr.build_vertex_index();

	graph_cluster gc(gr, *this, cfg.max_reads_partition_gap, false); //store hits true tasfia
	vector<pereads_cluster> vc;
	vector<pereads_cluster> vc_circ;

	gc.build_pereads_clusters(vc);
	gc.build_pereads_clusters_circ(vc_circ);

	printf("size of vc:%d\n",vc.size());
	printf("size of vc_circ:%d\n",vc_circ.size());

	vc.insert(vc.end(), vc_circ.begin(), vc_circ.end());

	printf("new size of vc:%d\n",vc.size());
	bridge_solver bs(gr, vc, cfg, sp.insertsize_low, sp.insertsize_high);
	// bridge_solver bs_circ(gr, vc_circ, cfg, sp.insertsize_low, sp.insertsize_high);

	assert(vc.size() == bs.opt.size());

	int cnt = 0;
	for(int k = 0; k < vc.size(); k++)
	{
		if(vc[k].is_circ == true) continue;
		if(bs.opt[k].type <= 0) continue; 
		cnt += update_bridges(vc[k].frlist, bs.opt[k].chain, bs.opt[k].strand);
	}

	printf("gid %s: total frags %lu, bridged frags = %d\n", gid.c_str(), frgs.size(), cnt);

	int non_empty_chain_vc_circ_count = 0;
	int bridge_cnt = 0;

	for(int k = 0; k < vc.size(); k++)
	{
		if(vc[k].is_circ == false) continue;
		if(hits[circ_frgs[vc[k].frlist[0]][0]].qname == "E00512:127:HJNF3ALXX:1:1105:24403:29630" && bs.opt[k].type <= 0)
		{
			printf("E00512:127:HJNF3ALXX:1:1105:24403:29630 path type negative\n");
		}
		if(bs.opt[k].type <= 0) continue;
		non_empty_chain_vc_circ_count += 1;
		bridge_cnt += update_bridges_circ(vc[k].frlist, bs.opt[k].chain, bs.opt[k].strand);
	}

	printf("gid %s: total circ frags %lu, bridged circ frags = %d\n", gid.c_str(), non_empty_chain_vc_circ_count, bridge_cnt);

	//join circ frags
	unordered_map<string, pair<int,int>> reg_index; 
	// key = qname; value = list of (cluster_index, frag_index_inside_cluster)
	
	//build qname index for frags in reg clusters that has a suppl
	for(int j = 0; j < vc.size(); j++)
	{
		if(vc[j].is_circ == true) continue;

		for(int k = 0; k < vc[j].frlist.size(); k++)
		{
			int reg_hit1_index = frgs[vc[j].frlist[k]][0];
			int reg_hit2_index = frgs[vc[j].frlist[k]][1];

			hit &h1 = hits[reg_hit1_index];
			hit &h2 = hits[reg_hit2_index];

			if(h1.suppl != NULL)
				reg_index[h1.qname] = {j,k};

			else if(h2.suppl != NULL)
				reg_index[h2.qname] = {j,k};
		}
	}

	for(int i = 0; i < vc.size(); i++)
	{
		pereads_cluster vc_circ = vc[i]; //one frag in each cluster
		if(vc_circ.is_circ == false) continue;

		int circ_hit1_index = circ_frgs[vc_circ.frlist[0]][0];
		int circ_hit2_index = circ_frgs[vc_circ.frlist[0]][1];

		hit &circ_h1 = hits[circ_hit1_index];
		hit &circ_h2 = hits[circ_hit2_index];

		vector<int32_t> vc_circ_bridge_chain = bs.opt[i].chain; //chain index and vc index corr to same peread and its bridged chain

		if((circ_h2.flag & 0x800) >= 1) //reg h1 has suppl part, as circ frag h2 has suppl not null
		{
			string circ_h2_suppl_qname = hits[circ_hit2_index].qname; 

			auto it1 = reg_index.find(circ_h2_suppl_qname);
			if(it1 != reg_index.end())
			{
				int j = it1->second.first;
				int k = it1->second.second;

				pereads_cluster vc_reg = vc[j]; //multiple frags in each cluster
				vector<int32_t> vc_reg_bridge_chain = bs.opt[j].chain;

				int reg_hit1_index = frgs[vc_reg.frlist[k]][0];
				int reg_hit2_index = frgs[vc_reg.frlist[k]][1];

				hit &h1 = hits[reg_hit1_index];
				hit &h2 = hits[reg_hit2_index];

				assert(h1.qname == h1.suppl->qname);
				assert(h1.suppl->hid == circ_h2.hid);

				printf("joined circRNA H1 has supple:\n");
				h1.print();
				h2.print();
				h1.suppl->print();

				vector<int32_t> x,y,z,final_intron_chain;
				merge_intron_chains(vc_reg.chain1, vc_reg_bridge_chain, x);
				merge_intron_chains(x, vc_reg.chain2, y);
				merge_intron_chains(y, vc_circ_bridge_chain, z);
				merge_intron_chains(z, vc_circ.chain2, final_intron_chain);

				printf("final intron chain H1 supple, read %s:\n",h1.qname.c_str());
				printv(final_intron_chain);

				//store circRNA
				circular_transcript circ;
				circ.sid = sp.sample_id;
				circ.chrm = chrm;
				circ.start = h1.pos;
				circ.end = h1.suppl->rpos;
				circ.id = chrm + ":" + tostring(circ.start) + "-" + tostring(circ.end);
				circ.source = "circMeta";
				circ.feature = "circRNA";
				circ.score = 1;
				circ.coverage = 1;
				circ.strand = strand;
				circ.gene_id = "gene_id";
				circ.transcript_id = h1.qname;
				circ.exon_count = final_intron_chain.size()/2 + 1;
				circ.intron_chain = final_intron_chain;

				circ.h1_suppl = true;
				circ.intermediate_chains.push_back(vc_reg.chain1);
				circ.intermediate_chains.push_back(vc_reg_bridge_chain);
				circ.intermediate_chains.push_back(vc_reg.chain2);
				circ.intermediate_chains.push_back(vc_circ_bridge_chain);
				circ.intermediate_chains.push_back(vc_circ.chain2);
				
				// check partially bridged
				int circ_frag_idx = vc_circ.frlist[0];
				int reg_frag_idx = vc_reg.frlist[k];
				
				if(circ_frgs[circ_frag_idx][2] > 0 && frgs[reg_frag_idx][2] > 0)
				{
					printf("circRNA H1 suppl formed, both bridged\n");
					circ.frag_bridged_type = 1;  //both bridged
					circ_trsts.push_back(circ);
				}
				else if(circ_frgs[circ_frag_idx][2] > 0 && frgs[reg_frag_idx][2] <= 0)
				{
					printf("circRNA H1 suppl formed, circ frag bridged but reg frag unbridged\n");
					circ.frag_bridged_type = 2; //circ frag unbridged but reg bridged
					unbridged_candidate_trsts[h1.qname] = circ;
					circ_frgs[circ_frag_idx][2] = 0;
				}
				else if(circ_frgs[circ_frag_idx][2] <= 0 && frgs[reg_frag_idx][2] > 0)
				{
					printf("circRNA H1 suppl formed, circ frag unbridged but reg frag bridged\n");
					circ.frag_bridged_type = 3; //reg bridged but circ unbridged
					unbridged_candidate_trsts[h1.qname] = circ;
					frgs[reg_frag_idx][2] = 0;
				}
				else
				{
					printf("circRNA H1 suppl formed, both unbridged\n");
				}
			}
		}

		else if((circ_h1.flag & 0x800) >= 1) //reg h2 has suppl part, as circ frag h1 has suppl not null
		{
			string circ_h1_suppl_qname = hits[circ_hit1_index].qname; 

			auto it2 = reg_index.find(circ_h1_suppl_qname);
			if(it2 != reg_index.end())
			{
				int j = it2->second.first;
				int k = it2->second.second;

				pereads_cluster vc_reg = vc[j]; //multiple frags in each cluster
				vector<int32_t> vc_reg_bridge_chain = bs.opt[j].chain;

				int reg_hit1_index = frgs[vc_reg.frlist[k]][0];
				int reg_hit2_index = frgs[vc_reg.frlist[k]][1];

				hit &h1 = hits[reg_hit1_index];
				hit &h2 = hits[reg_hit2_index];

				assert(h2.qname == h2.suppl->qname);
				assert(h2.suppl->hid == circ_h1.hid);

				printf("joined circRNA H2 has supple:\n");
				h2.suppl->print();
				h1.print();
				h2.print();

				vector<int32_t> x,y,z,final_intron_chain;
				merge_intron_chains(vc_circ.chain1, vc_circ_bridge_chain, x);
				merge_intron_chains(x, vc_reg.chain1, y);
				merge_intron_chains(y, vc_reg_bridge_chain, z);
				merge_intron_chains(z, vc_reg.chain2, final_intron_chain);

				printf("final intron chain H2 supple, read %s:\n",h2.qname.c_str());
				printv(final_intron_chain);

				if(h2.suppl->qname == "E00512:127:HJNF3ALXX:1:1105:24403:29630")
				{
					printf("printing circ bridge chain:");
					printv(vc_circ_bridge_chain);
					printf("printing reg bridge chain:");
					printv(vc_reg_bridge_chain);
				}

				//store circRNA
				circular_transcript circ;
				circ.sid = sp.sample_id;
				circ.chrm = chrm;
				circ.start = h2.suppl->pos;
				circ.end = h2.rpos;
				circ.id = chrm + ":" + tostring(circ.start) + "-" + tostring(circ.end);
				circ.source = "circMeta";
				circ.feature = "circRNA";
				circ.score = 1;
				circ.coverage = 1;
				circ.strand = strand;
				circ.gene_id = "gene_id";
				circ.transcript_id = h2.qname;
				circ.exon_count = final_intron_chain.size()/2 + 1;
				circ.intron_chain = final_intron_chain;

				circ.h1_suppl = false;
				circ.intermediate_chains.push_back(vc_circ.chain1);
				circ.intermediate_chains.push_back(vc_circ_bridge_chain);
				circ.intermediate_chains.push_back(vc_reg.chain1);
				circ.intermediate_chains.push_back(vc_reg_bridge_chain);
				circ.intermediate_chains.push_back(vc_reg.chain2);
				
				// check partially bridged
				int circ_frag_idx = vc_circ.frlist[0];
				int reg_frag_idx = vc_reg.frlist[k];

				if(circ_frgs[circ_frag_idx][2] > 0 && frgs[reg_frag_idx][2] > 0)
				{
					printf("circRNA H2 suppl formed, both bridged\n");
					circ.frag_bridged_type = 1;  //both bridged
					circ_trsts.push_back(circ);
				}
				else if(circ_frgs[circ_frag_idx][2] > 0 && frgs[reg_frag_idx][2] <= 0)
				{
					printf("circRNA H2 suppl formed, circ frag bridged but reg frag unbridged\n");
					circ.frag_bridged_type = 2; //circ frag unbridged but reg bridged
					unbridged_candidate_trsts[h1.qname] = circ;
					circ_frgs[circ_frag_idx][2] = 0;
				}
				else if(circ_frgs[circ_frag_idx][2] <= 0 && frgs[reg_frag_idx][2] > 0)
				{
					printf("circRNA H2 suppl formed, circ frag unbridged but reg frag bridged\n");
					circ.frag_bridged_type = 3; //reg bridged but circ frga unbridged
					unbridged_candidate_trsts[h1.qname] = circ;
					frgs[reg_frag_idx][2] = 0;
				}
				else
				{
					printf("circRNA H2 suppl formed, both unbridged\n");
				}
			}
		}
	}
	return 0;
}
int bundle::bridge_circ()
{
	splice_graph gr;
	graph_builder gb(*this, cfg, sp);
	gb.build(gr);
	gr.build_vertex_index();

	graph_cluster gc(gr, *this, cfg.max_reads_partition_gap, false); //store hits true tasfia
	vector<pereads_cluster> vc;
	vector<pereads_cluster> vc_circ;

	gc.build_pereads_clusters(vc);
	gc.build_pereads_clusters_circ(vc_circ);

	printf("size of vc:%d\n",vc.size());
	printf("size of vc_circ:%d\n",vc_circ.size());

	vc.insert(vc.end(), vc_circ.begin(), vc_circ.end());

	printf("new size of vc:%d\n",vc.size());
	bridge_solver bs(gr, vc, cfg, sp.insertsize_low, sp.insertsize_high);
	// bridge_solver bs_circ(gr, vc_circ, cfg, sp.insertsize_low, sp.insertsize_high);

	assert(vc.size() == bs.opt.size());

	int cnt = 0;
	for(int k = 0; k < vc.size(); k++)
	{
		if(vc[k].is_circ == true) continue;
		// if(bs.opt[k].type <= 0) continue; //negative for empty chain
		cnt += update_bridges(vc[k].frlist, bs.opt[k].chain, bs.opt[k].strand);
	}

	printf("gid %s: total frags %lu, bridged frags = %d\n", gid.c_str(), frgs.size(), cnt);

	int non_empty_chain_vc_circ_count = 0;
	int bridge_cnt = 0;

	for(int k = 0; k < vc.size(); k++)
	{
		if(vc[k].is_circ == false) continue;
		// if(bs.opt[k].type <= 0) continue; //negative for empty chain
		non_empty_chain_vc_circ_count += 1;
		bridge_cnt += update_bridges_circ(vc[k].frlist, bs.opt[k].chain, bs.opt[k].strand);
	}

	printf("gid %s: total circ frags %lu, bridged circ frags = %d\n", gid.c_str(), non_empty_chain_vc_circ_count, bridge_cnt);

	//join circ frags
	for(int i=0;i<vc.size();i++)
	{
		pereads_cluster vc_circ = vc[i]; //one frag in each cluster

		if(vc_circ.is_circ == false) continue;

		int circ_hit1_index = circ_frgs[vc_circ.frlist[0]][0];
		int circ_hit2_index = circ_frgs[vc_circ.frlist[0]][1];
		
		vector<int32_t> vc_circ_bridge_chain = bs.opt[i].chain; //chain index and vc index corr to same peread and its bridged chain

		bool partner_found = false;
		for(int j=0;j<vc.size();j++)
		{
			if(vc[j].is_circ == true) continue;
			// printf("entered regular loop\n");
			pereads_cluster vc_reg = vc[j]; //multiple frags in each cluster
			vector<int32_t> vc_reg_bridge_chain = bs.opt[j].chain;
			
			for(int k=0;k<vc_reg.frlist.size();k++)
			{
				// printf("entered vc_reg frlist\n");
				int reg_hit1_index = frgs[vc_reg.frlist[k]][0];
				int reg_hit2_index = frgs[vc_reg.frlist[k]][1];

				hit &h1 = hits[reg_hit1_index];
				hit &h2 = hits[reg_hit2_index];
				// printf("traverse vc_reg frlist\n");

				if(h1.suppl != NULL && h1.qname == hits[circ_hit2_index].qname)
				{
					// if(h1.suppl == NULL) continue;
					// if(h1.qname != vc_circ.hits2[0].qname) continue;
					printf("joined circRNA H1 has supple:\n");
					h1.print();
					h2.print();
					h1.suppl->print();
					// printf("joined circRNA bridge path H1 supple:\n");
					// printf("vc_reg_chain1:\n");
					// printv(vc_reg.chain1); // any merging function?
					// printf("vc_reg_bridge_chain:\n");
					// printv(vc_reg_bridge_chain);
					// printf("vc_reg_chain2:\n");
					// printv(vc_reg.chain2);
					// printf("vc_circ_bridge_chain:\n");
					// printv(vc_circ_bridge_chain);
					// printf("vc_circ_chain2:\n");
					// printv(vc_circ.chain2);

					vector<int32_t> x,y,z,final_intron_chain;
					merge_intron_chains(vc_reg.chain1, vc_reg_bridge_chain, x);
					merge_intron_chains(x, vc_reg.chain2, y);
					merge_intron_chains(y, vc_circ_bridge_chain, z);
					merge_intron_chains(z, vc_circ.chain2, final_intron_chain);

					printf("final intron chain H1 supple, read %s:\n",h1.qname.c_str());
					printv(final_intron_chain);

					//store circRNA
					circular_transcript circ;
					circ.sid = sp.sample_id;
					circ.chrm = chrm;
					circ.start = h1.pos;
					circ.end = h1.suppl->rpos;
					circ.id = chrm + ":" + tostring(circ.start) + "-" + tostring(circ.end);
					circ.source = "circMeta";
					circ.feature = "circRNA";
					circ.score = 1;
					circ.coverage = 1;
					circ.strand = strand;
					circ.gene_id = "gene_id";
					circ.transcript_id = h1.qname;
					circ.exon_count = final_intron_chain.size()/2 + 1;
					circ.intron_chain = final_intron_chain;

					circ.h1_suppl = true;
					circ.intermediate_chains.push_back(vc_reg.chain1);
					circ.intermediate_chains.push_back(vc_reg_bridge_chain);
					circ.intermediate_chains.push_back(vc_reg.chain2);
					circ.intermediate_chains.push_back(vc_circ_bridge_chain);
					circ.intermediate_chains.push_back(vc_circ.chain2);
					
					// check partially bridged
					int circ_frag_idx = vc_circ.frlist[0];
					int reg_frag_idx = vc_reg.frlist[k];
					
					if(circ_frgs[circ_frag_idx][2] > 0 && frgs[reg_frag_idx][2] > 0)
					{
						printf("circRNA H1 suppl formed, both bridged\n");
						circ.frag_bridged_type = 1;  //both bridged
						circ_trsts.push_back(circ);
					}
					else if(circ_frgs[circ_frag_idx][2] > 0 && frgs[reg_frag_idx][2] <= 0)
					{
						printf("circRNA H1 suppl formed, circ frag bridged but reg frag unbridged\n");
						circ.frag_bridged_type = 2; //circ frag unbridged but reg bridged
						unbridged_candidate_trsts[h1.qname] = circ;
						circ_frgs[circ_frag_idx][2] = 0;
					}
					else if(circ_frgs[circ_frag_idx][2] <= 0 && frgs[reg_frag_idx][2] > 0)
					{
						printf("circRNA H1 suppl formed, circ frag unbridged but reg frag bridged\n");
						circ.frag_bridged_type = 3; //reg bridged but circ unbridged
						unbridged_candidate_trsts[h1.qname] = circ;
						frgs[reg_frag_idx][2] = 0;
					}
					else
					{
						printf("circRNA H1 suppl formed, both unbridged\n");
					}
					partner_found = true;
					break;
				}
				else if(h2.suppl != NULL && h2.qname == hits[circ_hit1_index].qname)
				{
					// if(h2.suppl == NULL) continue;
					// if(h2.qname != vc_circ.hits1[0].qname) continue;
					printf("joined circRNA H2 has supple:\n");
					h2.suppl->print();
					h1.print();
					h2.print();
					// printf("joined circRNA bridge path H2 supple:\n");
					// printf("vc_circ_chain1:\n");
					// printv(vc_circ.chain1);
					// printf("vc_circ_bridge_chain:\n");
					// printv(vc_circ_bridge_chain);
					// printf("vc_reg_chain1:\n");
					// printv(vc_reg.chain1);
					// printf("vc_reg_bridge_chain:\n");
					// printv(vc_reg_bridge_chain);
					// printf("vc_reg_chain2:\n");
					// printv(vc_reg.chain2);

					vector<int32_t> x,y,z,final_intron_chain;
					merge_intron_chains(vc_circ.chain1, vc_circ_bridge_chain, x);
					merge_intron_chains(x, vc_reg.chain1, y);
					merge_intron_chains(y, vc_reg_bridge_chain, z);
					merge_intron_chains(z, vc_reg.chain2, final_intron_chain);

					printf("final intron chain H2 supple, read %s:\n",h2.qname.c_str());
					printv(final_intron_chain);

					//store circRNA
					circular_transcript circ;
					circ.sid = sp.sample_id;
					circ.chrm = chrm;
					circ.start = h2.suppl->pos;
					circ.end = h2.rpos;
					circ.id = chrm + ":" + tostring(circ.start) + "-" + tostring(circ.end);
					circ.source = "circMeta";
					circ.feature = "circRNA";
					circ.score = 1;
					circ.coverage = 1;
					circ.strand = strand;
					circ.gene_id = "gene_id";
					circ.transcript_id = h2.qname;
					circ.exon_count = final_intron_chain.size()/2 + 1;
					circ.intron_chain = final_intron_chain;

					circ.h1_suppl = false;
					circ.intermediate_chains.push_back(vc_circ.chain1);
					circ.intermediate_chains.push_back(vc_circ_bridge_chain);
					circ.intermediate_chains.push_back(vc_reg.chain1);
					circ.intermediate_chains.push_back(vc_reg_bridge_chain);
					circ.intermediate_chains.push_back(vc_reg.chain2);
					
					// check partially bridged
					int circ_frag_idx = vc_circ.frlist[0];
					int reg_frag_idx = vc_reg.frlist[k];

					if(circ_frgs[circ_frag_idx][2] > 0 && frgs[reg_frag_idx][2] > 0)
					{
						printf("circRNA H2 suppl formed, both bridged\n");
						circ.frag_bridged_type = 1;  //both bridged
						circ_trsts.push_back(circ);
					}
					else if(circ_frgs[circ_frag_idx][2] > 0 && frgs[reg_frag_idx][2] <= 0)
					{
						printf("circRNA H2 suppl formed, circ frag bridged but reg frag unbridged\n");
						circ.frag_bridged_type = 2; //circ frag unbridged but reg bridged
						unbridged_candidate_trsts[h1.qname] = circ;
						circ_frgs[circ_frag_idx][2] = 0;
					}
					else if(circ_frgs[circ_frag_idx][2] <= 0 && frgs[reg_frag_idx][2] > 0)
					{
						printf("circRNA H2 suppl formed, circ frag unbridged but reg frag bridged\n");
						circ.frag_bridged_type = 3; //reg bridged but circ frga unbridged
						unbridged_candidate_trsts[h1.qname] = circ;
						frgs[reg_frag_idx][2] = 0;
					}
					else
					{
						printf("circRNA H2 suppl formed, both unbridged\n");
					}
					partner_found = true;
					break;
				}
				if(partner_found == true) break;
			}
		}
	}
	return 0;
}

int bundle::bridge()
{
	/*
	int round = 0;
	while(round < 2)
	{
	*/
		splice_graph gr;
		graph_builder gb(*this, cfg, sp);
		gb.build(gr);
		gr.build_vertex_index();

		vector<pereads_cluster> vc;
		graph_cluster gc(gr, *this, cfg.max_reads_partition_gap, false);
		gc.build_pereads_clusters(vc);

		bridge_solver bs(gr, vc, cfg, sp.insertsize_low, sp.insertsize_high);

		int cnt = 0;
		assert(vc.size() == bs.opt.size());
		for(int k = 0; k < vc.size(); k++)
		{
			if(bs.opt[k].type <= 0) continue;
			cnt += update_bridges(vc[k].frlist, bs.opt[k].chain, bs.opt[k].strand);
		}

		printf("gid %s: total frags %lu, bridged frags = %d\n", gid.c_str(), frgs.size(), cnt);

	/*
		round++;
		if(cnt <= 0) break;
	}
	*/
	return 0;
}

int bundle::combine(const bundle &bb, bool combine_map)
{
	num_combined += bb.num_combined;
	assert(strand == bb.strand);
	assert(chrm == bb.chrm);
	assert(tid == bb.tid);
	if(lpos > bb.lpos) lpos = bb.lpos;
	if(rpos < bb.rpos) rpos = bb.rpos;
	hcst.add(bb.hcst);
	fcst.add(bb.fcst);
	//mmap.insert(mmap.end(), bb.mmap.begin(), bb.mmap.end());
	//imap.insert(imap.end(), bb.imap.begin(), bb.imap.end());
	if(combine_map) mmap += bb.mmap;
	if(combine_map) imap += bb.imap;
	//for(SIMI z = bb.mmap.begin(); z != bb.mmap.end(); z++) mmap += *z;
	//for(SIMI z = bb.imap.begin(); z != bb.imap.end(); z++) imap += *z;
	return 0;
}

int bundle::print(int index)
{
	printf("bundle %d: sample = %d, ", index, sp.sample_id);

	// statistic xs
	int n0 = 0, np = 0, nq = 0;
	for(int i = 0; i < hits.size(); i++)
	{
		if(hits[i].xs == '.') n0++;
		if(hits[i].xs == '+') np++;
		if(hits[i].xs == '-') nq++;
	}

	printf("tid = %d, range = %s:%d-%d, orient = %c, #hits = %lu, #frgs = %lu, +/-/. = %d / %d / %d\n",
			tid, chrm.c_str(), lpos, rpos, strand, hits.size(), frgs.size(), np, nq, n0);

	return 0;
}

/*
int bundle::count_unbridged_fragments()
{
	int c = 0;
	for(int j = 0; j < frgs.size(); j++)
	{
		if(frgs[j][2] <= 0) c++;
	}
	return c;
}
*/
