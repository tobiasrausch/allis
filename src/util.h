/*
============================================================================
Allis: Allele-specific expression and open chromatin
============================================================================
Copyright (C) 2017 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#ifndef UTIL_H
#define UTIL_H

#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <htslib/sam.h>


namespace allis
{

  struct BiallelicVariant {
    int32_t pos;
    std::string ref;
    std::string alt;
    bool hap;

    BiallelicVariant(int32_t p) : pos(p), ref(""), alt(""), hap(0) {}
    BiallelicVariant(int32_t p, std::string r, std::string a, bool h) : pos(p), ref(r), alt(a), hap(h) {}
  };


  template<typename TRecord>
  struct SortVariants : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      return s1.pos < s2.pos;
    }
  };

  inline uint32_t
  lastAlignedPosition(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CDEL)) alen += bam_cigar_oplen(cigar[i]);
    return rec->core.pos + alen;
  }

  inline unsigned
  hash_string(const char *s) {
    unsigned h = *s;
    if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
    return h;
  }

  inline std::size_t
  hash_pair(bam1_t* rec) {
    std::size_t seed = 0;
    if (rec->core.flag & BAM_FREAD2) {
      boost::hash_combine(seed, rec->core.mtid);
      boost::hash_combine(seed, rec->core.mpos);
      boost::hash_combine(seed, rec->core.tid);
      boost::hash_combine(seed, rec->core.pos);
    } else {
      boost::hash_combine(seed, rec->core.tid);
      boost::hash_combine(seed, rec->core.pos);
      boost::hash_combine(seed, rec->core.mtid);
      boost::hash_combine(seed, rec->core.mpos);
    }
    return seed;
  }
  
  template<typename TVariants>
  inline bool
  _loadVariants(htsFile* ifile, hts_idx_t* bcfidx, bcf_hdr_t* hdr, std::string const& sample, std::string const& chrom, TVariants& pV) {
    typedef typename TVariants::value_type TVariant;
    
    int32_t sampleIndex = -1;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
      if (hdr->samples[i] == sample) sampleIndex = i;
    if (sampleIndex < 0) return false;
        
    // Genotypes
    int ngt = 0;
    int32_t* gt = NULL;

    // Collect het. bi-allelic variants for this chromosome
    int32_t chrid = bcf_hdr_name2id(hdr, chrom.c_str());
    if (chrid < 0) return false;
    hts_itr_t* itervcf = bcf_itr_querys(bcfidx, hdr, chrom.c_str());
    if (itervcf != NULL) {
      bcf1_t* rec = bcf_init1();
      while (bcf_itr_next(ifile, itervcf, rec) >= 0) {
	// Only bi-allelic variants
	if (rec->n_allele == 2) {
	  bcf_unpack(rec, BCF_UN_ALL);
	  bcf_get_genotypes(hdr, rec, &gt, &ngt);
	  if ((bcf_gt_allele(gt[sampleIndex*2]) != -1) && (bcf_gt_allele(gt[sampleIndex*2 + 1]) != -1) && (!bcf_gt_is_missing(gt[sampleIndex*2])) && (!bcf_gt_is_missing(gt[sampleIndex*2 + 1]))) {
	    int gt_type = bcf_gt_allele(gt[sampleIndex*2]) + bcf_gt_allele(gt[sampleIndex*2 + 1]);
	    if (gt_type == 1) {
	      std::vector<std::string> alleles;
	      for(std::size_t i = 0; i<rec->n_allele; ++i) alleles.push_back(std::string(rec->d.allele[i]));
	      //std::cerr << chrom << "\t" << (rec->pos + 1) << "\t" << std::string(alleles[0]) << "\t" << std::string(alleles[1]) << std::endl;
	      pV.push_back(TVariant(rec->pos, std::string(alleles[0]), std::string(alleles[1]), bcf_gt_allele(gt[sampleIndex*2])));
	    }
	  }
	}
      }
      bcf_destroy(rec);
      hts_itr_destroy(itervcf);
    }
    if (gt != NULL) free(gt);
    return true;
  }

  
  template<typename TVariants>
  inline bool
  _loadVariants(std::string const& sample, std::string const& chrom, std::string const& bcffile, TVariants& pV) {
    // Load BCF file
    htsFile* ifile = bcf_open(bcffile.c_str(), "r");
    hts_idx_t* bcfidx = bcf_index_load(bcffile.c_str());
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);

    bool retVal = _loadVariants(ifile, bcfidx, hdr, sample, chrom, pV);
    
    // Close BCF
    bcf_hdr_destroy(hdr);
    hts_idx_destroy(bcfidx);
    bcf_close(ifile);
    
    return retVal;
  }


}

#endif
