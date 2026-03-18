/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstring>
#include <cassert>
#include <cstdio>
#include <sstream>
#include <cmath>

#include "hit.h"
#include "util.h"
#include "constants.h"

hit::hit()
{
	hid = 0;
	rpos = 0;
	nh = 0;
	hi = 0;
	nm = 0;
	strand = '.';
	xs = '.';
	ts = '.';
	qname = "";

	suppl = NULL;
	suppl_index = -1;
	cigar_vector.clear();
	left_cigar = '.';					
	right_cigar = '.';					
	left_cigar_len = 0;
	right_cigar_len = 0;
	first_pos = 0;						
	second_pos = 0;
	third_pos = 0;

	seq = "";
	soft_left_clip_seqs.clear();
	soft_right_clip_seqs.clear();
	is_fake = false;
	fake_hit_index = -1;
	soft_clip_side = 0;
	has_fake_suppl = false;
	has_chimeric_suppl = false;
}

hit& hit::operator=(const hit &h)
{
	bam1_core_t::operator=(h);
	hid = h.hid;
	rpos = h.rpos;
	qname = h.qname;
	strand = h.strand;
	//spos = h.spos;
	xs = h.xs;
	ts = h.ts;
	hi = h.hi;
	nm = h.nm;
	nh = h.nh;

	suppl = h.suppl;
	suppl_index = h.suppl_index;
	cigar_vector = h.cigar_vector;
	left_cigar = h.left_cigar;
	right_cigar = h.right_cigar;
	left_cigar_len = h.left_cigar_len;
	right_cigar_len = h.right_cigar_len;
	first_pos = h.first_pos;
	second_pos = h.second_pos;
	third_pos = h.third_pos;
	seq  = h.seq;
	soft_left_clip_seqs = h.soft_left_clip_seqs;
	soft_right_clip_seqs = h.soft_right_clip_seqs;
	is_fake = h.is_fake;
	fake_hit_index = h.fake_hit_index;
	soft_clip_side = h.soft_clip_side;
	has_fake_suppl = h.has_fake_suppl;
	has_chimeric_suppl = h.has_chimeric_suppl;

	return *this;
}

hit::hit(const hit &h)
	:bam1_core_t(h)
{
	hid = h.hid;
	rpos = h.rpos;
	qname = h.qname;
	strand = h.strand;
	//spos = h.spos;
	xs = h.xs;
	ts = h.ts;
	hi = h.hi;
	nm = h.nm;
	nh = h.nh;

	suppl = h.suppl;
	suppl_index = h.suppl_index;
	cigar_vector = h.cigar_vector;
	left_cigar = h.left_cigar;
	right_cigar = h.right_cigar;
	left_cigar_len = h.left_cigar_len;
	right_cigar_len = h.right_cigar_len;
	first_pos = h.first_pos;
	second_pos = h.second_pos;
	third_pos = h.third_pos;
	seq  = h.seq;
	soft_left_clip_seqs = h.soft_left_clip_seqs;
	soft_right_clip_seqs = h.soft_right_clip_seqs;
	is_fake = h.is_fake;
	fake_hit_index = h.fake_hit_index;
	soft_clip_side = h.soft_clip_side;
	has_fake_suppl = h.has_fake_suppl;
	has_chimeric_suppl = h.has_chimeric_suppl;
}

hit::~hit()
{
}

hit::hit(bam1_t *b, int id)
	:bam1_core_t(b->core), hid(id)
{
	// fetch query name
	char buf[1024];
	char *qs = bam_get_qname(b);
	int l = strlen(qs);
	memcpy(buf, qs, l);
	buf[l] = '\0';
	qname = string(buf);

	// compute rpos
	rpos = pos + (int32_t)bam_cigar2rlen(n_cigar, bam_get_cigar(b));

	suppl = NULL;
	suppl_index = -1;
	left_cigar = '.';
	right_cigar = '.';
	left_cigar_len = 0;
	right_cigar_len = 0;
	first_pos = 0;
	second_pos = 0;
	third_pos = 0;
	seq = "";
	soft_left_clip_seqs.clear();
	soft_right_clip_seqs.clear();
	is_fake = false;
	fake_hit_index = -1;
	soft_clip_side = 0;
	has_fake_suppl = false;
	has_chimeric_suppl = false;

	int max_num_cigar = 10000; // TODO use from paemeters.h
	assert(n_cigar <= max_num_cigar);
	assert(n_cigar >= 1);
	uint32_t * cigar = bam_get_cigar(b);
	
	cigar_vector.clear();
	for(int k = 0; k < n_cigar; k++)
	{
		if(bam_cigar_op(cigar[k]) == BAM_CMATCH)
		{
			cigar_vector.push_back(pair<char, int32_t>('M',bam_cigar_oplen(cigar[k])));
		}
		else if(bam_cigar_op(cigar[k]) == BAM_CSOFT_CLIP)
		{
			cigar_vector.push_back(pair<char, int32_t>('S',bam_cigar_oplen(cigar[k])));
		}
		else if(bam_cigar_op(cigar[k]) == BAM_CHARD_CLIP)
		{
			cigar_vector.push_back(pair<char, int32_t>('H',bam_cigar_oplen(cigar[k])));
		}
		else if(bam_cigar_op(cigar[k]) == BAM_CREF_SKIP)
		{
			cigar_vector.push_back(pair<char, int32_t>('N',bam_cigar_oplen(cigar[k])));
		}
		else if(bam_cigar_op(cigar[k]) == BAM_CINS)
		{
			cigar_vector.push_back(pair<char, int32_t>('I',bam_cigar_oplen(cigar[k])));
		}
		else if(bam_cigar_op(cigar[k]) == BAM_CDEL)
		{
			cigar_vector.push_back(pair<char, int32_t>('D',bam_cigar_oplen(cigar[k])));
		}
		else
		{
			cigar_vector.push_back(pair<char, int32_t>('.',0));
		}
	}
	
	//putting a placeholder to avoid core dumped when cigar_vector[0] is accessed
	if(cigar_vector.size() == 0)
	{
		cigar_vector.push_back(pair<char, int32_t>('.',0));
	}

	if(cigar_vector[0].first == 'S' || cigar_vector[cigar_vector.size()-1].first == 'S')
	{
		set_seq(b);
	}
}

int hit::set_seq(bam1_t *b)
{
	uint32_t seq_len = b->core.l_qseq;
	l_qseq = seq_len;
	uint8_t *q = bam_get_seq(b); 
	vector<int> code;

	for(int i=0;i<seq_len;i++)
	{
		code.push_back(bam_seqi(q,i)); //gets nucleotide id
	}

	string hit_seq = convert_to_IUPAC(code);
	seq = hit_seq;

	set_soft_clip_seq_combo();

	return 0;
}

string hit::convert_to_IUPAC(vector<int> code)
{
	string seq = "";
	for(int i=0;i<code.size();i++)
	{
		if(code[i] == 1)
		{
			seq = seq + 'A';
		}
		else if(code[i] == 2)
		{
			seq = seq + 'C';
		}
		else if(code[i] == 4)
		{
			seq = seq + 'G';
		}
		else if(code[i] == 8)
		{
			seq = seq + 'T';
		}
		else if(code[i] == 15)
		{
			seq = seq + 'N';
		}
	}

	return seq;
}

int hit::set_soft_clip_seq_combo()
{
	if(cigar_vector[0].first == 'S')
	{
		int32_t len = cigar_vector[0].second;
		//index 0, extract start len bp
		string str0 = "";
		for(int i=1;i<len+1;i++)
		{
			str0 = str0 + seq[i];
		}
		soft_left_clip_seqs.push_back(str0);

		//index 1, extract start len bp rev comp
		soft_left_clip_seqs.push_back(get_reverse_complement(soft_left_clip_seqs[0]));
	}
	if(cigar_vector[cigar_vector.size()-1].first == 'S')
	{
		int32_t len = cigar_vector[cigar_vector.size()-1].second;

		//index 0, extract end len bp
		string str2 = "";
		for(int i=seq.size()-len;i<seq.size();i++)
		{
			str2 = str2 + seq[i];
		}
		soft_right_clip_seqs.push_back(str2);

		//index 1,extract end len bp rev comp
		soft_right_clip_seqs.push_back(get_reverse_complement(soft_right_clip_seqs[0]));
	}

	// printf("Printing four combos for read:%s\n",qname.c_str());
	// if(soft_left_clip_seqs.size() == 2)
	// {
	// 	printf("start: %s\n",soft_left_clip_seqs[0].c_str());
	// 	printf("start RC: %s\n",soft_left_clip_seqs[1].c_str());
	// }
	// if(soft_right_clip_seqs.size() == 2)
	// {
	// 	printf("end: %s\n",soft_right_clip_seqs[0].c_str());
	// 	printf("end RC: %s\n",soft_right_clip_seqs[1].c_str());
	// }

	return 0;
}

string hit::get_reverse_complement(string str)
{
	string out = "";
	for(int i=str.size()-1;i>=0;i--)
	{
		if(str[i] == 'A')
		{
			out = out + 'T';
		}
		else if(str[i] == 'T')
		{
			out = out + 'A';
		}
		else if(str[i] == 'C')
		{
			out = out + 'G';
		}
		else if(str[i] == 'G')
		{
			out = out + 'C';
		}
		else if(str[i] == 'N')
		{
			out = out + str[i];
		}
	}
	return out;
}

bool hit::contain_splices(bam1_t *b) const
{
	uint32_t *cigar = bam_get_cigar(b);
    for(int k = 0; k < n_cigar; k++)
	{
		if(bam_cigar_op(cigar[k]) == BAM_CREF_SKIP) return true;
	}
	return false;
}

vector<int32_t> hit::extract_splices(bam1_t *b) const
{
	vector<int32_t> spos;
	uint32_t *cigar = bam_get_cigar(b);
	int32_t p = pos;
	int32_t q = 0;
    for(int k = 0; k < n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
			p += bam_cigar_oplen(cigar[k]);

		if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
			q += bam_cigar_oplen(cigar[k]);

		if(k == 0 || k == n_cigar - 1) continue;
		if(bam_cigar_op(cigar[k]) != BAM_CREF_SKIP) continue;
		//if(bam_cigar_op(cigar[k-1]) != BAM_CMATCH) continue;
		//if(bam_cigar_op(cigar[k+1]) != BAM_CMATCH) continue;
		////if(bam_cigar_oplen(cigar[k-1]) < min_flank) continue;
		////if(bam_cigar_oplen(cigar[k+1]) < min_flank) continue;

		int32_t s = p - bam_cigar_oplen(cigar[k]);
		//spos.push_back(pack(s, p));
		spos.push_back(s);
		spos.push_back(p);
	}
	return spos;
}

int hit::set_tags(bam1_t *b)
{
	ts = '.';
	uint8_t *p0 = bam_aux_get(b, "ts");
	if(p0 && (*p0) == 'A') ts = bam_aux2A(p0);

	xs = '.';
	uint8_t *p1 = bam_aux_get(b, "XS");
	if(p1 && (*p1) == 'A') xs = bam_aux2A(p1);

	if(xs == '.' && ts != '.')
	{
		// convert ts to xs
		if((flag & 0x10) >= 1 && ts == '+') xs = '-';
		if((flag & 0x10) >= 1 && ts == '-') xs = '+';
		if((flag & 0x10) <= 0 && ts == '+') xs = '+';
		if((flag & 0x10) <= 0 && ts == '-') xs = '-';
	}

	hi = -1;
	uint8_t *p2 = bam_aux_get(b, "HI");
	if(p2) hi = bam_aux2i(p2);

	nh = -1;
	uint8_t *p3 = bam_aux_get(b, "NH");
	if(p3) nh = bam_aux2i(p3);

	nm = 0;
	uint8_t *p4 = bam_aux_get(b, "nM");
	if(p4) nm = bam_aux2i(p4);

	uint8_t *p5 = bam_aux_get(b, "NM");
	if(p5) nm = bam_aux2i(p5);

	return 0;
}

bool hit::get_concordance() const
{
	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) return true;		// F1R2
	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) return true;		// R1F2
	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) return true;		// F2R1
	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) return true;		// R2F1
	return false;
}

int hit::set_strand(int libtype)
{
	strand = '.';
	
	if(libtype == FR_FIRST && ((flag & 0x1) >= 1))
	{
		if((flag & 0x10) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '-';
		if((flag & 0x10) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '+';
		if((flag & 0x10) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '+';
		if((flag & 0x10) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '-';
	}

	if(libtype == FR_SECOND && ((flag & 0x1) >= 1))
	{
		if((flag & 0x10) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '+';
		if((flag & 0x10) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '-';
		if((flag & 0x10) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '-';
		if((flag & 0x10) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '+';
	}

	if(libtype == FR_FIRST && ((flag & 0x1) <= 0))
	{
		if((flag & 0x10) <= 0) strand = '-';
		if((flag & 0x10) >= 1) strand = '+';
	}

	if(libtype == FR_SECOND && ((flag & 0x1) <= 0))
	{
		if((flag & 0x10) <= 0) strand = '+';
		if((flag & 0x10) >= 1) strand = '-';
	}

	return 0;
}

bool hit::operator<(const hit &h) const
{
	if(qname < h.qname) return true;
	if(qname > h.qname) return false;
	if(hi != -1 && h.hi != -1 && hi < h.hi) return true;
	if(hi != -1 && h.hi != -1 && hi > h.hi) return false;
	return (pos < h.pos);
}

int hit::print() const
{
	// print basic information
	printf("Hit %s: tid = %d, hid = %d, [%d-%d), mpos = %d, flag = %d, quality = %d, strand = %c, xs = %c, ts = %c, isize = %d, hi = %d\n", 
			qname.c_str(), tid, hid, pos, rpos, mpos, flag, qual, strand, xs, ts, isize, hi);

	return 0;

	/*
	printf(" start position (%d - )\n", pos);
	for(int i = 0; i < spos.size() / 2; i++)
	{
		int32_t p1 = spos[i * 2 + 0];
		int32_t p2 = spos[i * 2 + 1];
		printf(" splice position (%d - %d)\n", p1, p2);
	}
	printf(" end position (%d - )\n", rpos);
	return 0;
	*/
}

size_t hit::get_qhash() const
{
	return string_hash(qname);
}

/*
int hit::get_aligned_intervals(vector<int64_t> &v) const
{
	v.clear();
	int32_t p1 = pos;
	for(int k = 0; k < spos.size() / 2; k++)
	{
		int32_t p2 = spos[k * 2 + 0];
		v.push_back(pack(p1, p2));
		p1 = spos[k * 2 + 1];
	}
	v.push_back(pack(p1, rpos));
	return 0;
}

size_t hit::get_phash() const
{
	vector<int32_t> v;
	v.push_back(pos);
	v.insert(v.end(), spos.begin(), spos.end());
	v.push_back(rpos);
	return vector_hash(v);
}

bool hit::equal(const hit &h) const
{
	if(strand != h.strand) return false;

	if(tid != h.tid) return false;
	if(mtid != h.mtid) return false;
	if(pos != h.pos) return false;
	if(mpos != h.mpos) return false;
	if(rpos != h.rpos) return false;
	if(spos != h.spos) return false;
	if(flag != h.flag) return false;
	if(isize != h.isize) return false;
	if(n_cigar != h.n_cigar) return false;
	if(l_qname != h.l_qname) return false;
	if(l_extranul != h.l_extranul) return false;
	return true;
}
*/
