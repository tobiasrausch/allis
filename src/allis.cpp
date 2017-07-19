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

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#include <boost/math/distributions/binomial.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#ifdef OPENMP
#include <omp.h>
#endif

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include "util.h"
#include "infoparse.h"

using namespace allis;

struct Config {
  unsigned short minMapQual;
  unsigned short minBaseQual;
  std::string sample;
  boost::filesystem::path as;
  boost::filesystem::path genome;
  boost::filesystem::path h1bam;
  boost::filesystem::path h2bam;
  boost::filesystem::path bamfile;
  boost::filesystem::path vcffile;
};

template<typename TConfig>
inline int
phaseBamRun(TConfig& c) {

#ifdef PROFILE
  ProfilerStart("phasebam.prof");
#endif

  // Load bam files
  samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
  hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
  bam_hdr_t* hdr = sam_hdr_read(samfile);

  // Load bcf file
  htsFile* ibcffile = bcf_open(c.vcffile.string().c_str(), "r");
  hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
  bcf_hdr_t* bcfhdr = bcf_hdr_read(ibcffile);

  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Assign reads to haplotypes" << std::endl;
  boost::progress_display show_progress(hdr->n_targets);

  // Open output file
  samFile* h1bam = sam_open(c.h1bam.string().c_str(), "wb");
  if (sam_hdr_write(h1bam, hdr) != 0) {
    std::cerr << "Could not write ouptut file header!" << std::endl;
    return -1;
  }

  samFile* h2bam = sam_open(c.h2bam.string().c_str(), "wb");
  if (sam_hdr_write(h2bam, hdr) != 0) {
    std::cerr << "Could not write ouptut file header!" << std::endl;
    return -1;
  }

  // Allele support file
  boost::iostreams::filtering_ostream dataOut;
  dataOut.push(boost::iostreams::gzip_compressor());
  dataOut.push(boost::iostreams::file_sink(c.as.string().c_str(), std::ios_base::out | std::ios_base::binary));
  dataOut << "chr\tpos\tid\tref\talt\tdepth\trefsupport\taltsupport\tgt\tvaf\th1af\tpvalue\t1kgpaf" << std::endl;
  
  // Assign reads to SNPs
  uint32_t assignedReadsH1 = 0;
  uint32_t assignedReadsH2 = 0;
  uint32_t unassignedReads = 0;
  uint32_t ambiguousReads = 0;
  faidx_t* fai = fai_load(c.genome.string().c_str());
  for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
    std::string chrName(hdr->target_name[refIndex]);
    ++show_progress;

    // Load het. markers
    typedef std::vector<BiallelicVariant> TPhasedVariants;
    TPhasedVariants pv;
    if (!_loadVariants(ibcffile, bcfidx, bcfhdr, c.sample, chrName, pv)) continue;
    if (pv.empty()) continue;

    // Sort variants
    std::sort(pv.begin(), pv.end(), SortVariants<BiallelicVariant>());

    // Load reference
    int32_t seqlen = -1;
    char* seq = faidx_fetch_seq(fai, chrName.c_str(), 0, hdr->target_len[refIndex], &seqlen);    
    
    // Assign reads to haplotypes
    std::set<std::size_t> h1;
    std::set<std::size_t> h2;
    hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
    bam1_t* rec = bam_init1();
    while (sam_itr_next(samfile, iter, rec) >= 0) {
      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
      if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
      uint32_t hp1votes = 0;
      uint32_t hp2votes = 0;
      TPhasedVariants::const_iterator vIt = std::lower_bound(pv.begin(), pv.end(), BiallelicVariant(rec->core.pos), SortVariants<BiallelicVariant>());
      TPhasedVariants::const_iterator vItEnd = std::upper_bound(pv.begin(), pv.end(), BiallelicVariant(lastAlignedPosition(rec)), SortVariants<BiallelicVariant>());
      if (vIt != vItEnd) {
	// Get read sequence
	std::string sequence;
	sequence.resize(rec->core.l_qseq);
	uint8_t* seqptr = bam_get_seq(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	  
	// Parse CIGAR
	uint32_t* cigar = bam_get_cigar(rec);
	for(;vIt != vItEnd; ++vIt) {
	  int32_t gp = rec->core.pos; // Genomic position
	  int32_t sp = 0; // Sequence position
	  bool varFound = false;
	  for (std::size_t i = 0; ((i < rec->core.n_cigar) && (!varFound)); ++i) {
	    if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) gp += bam_cigar_oplen(cigar[i]);
	    else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	      //Nop
	    } else if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
	      if (gp + (int32_t) bam_cigar_oplen(cigar[i]) < vIt->pos) {
		gp += bam_cigar_oplen(cigar[i]);
		sp += bam_cigar_oplen(cigar[i]);
	      } else {
		for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
		  if (gp == vIt->pos) {
		    varFound = true;
		    // Check REF allele
		    if (vIt->ref == std::string(seq + gp, seq + gp + vIt->ref.size())) {
		      // Check ALT allele
		      if ((sp + vIt->alt.size() < sequence.size()) && (sp + vIt->ref.size() < sequence.size())) {
			if (vIt->ref.size() == vIt->alt.size()) {
			  // SNP
			  if ((sequence.substr(sp, vIt->alt.size()) == vIt->alt) && (sequence.substr(sp, vIt->ref.size()) != vIt->ref)) {
			    // ALT supporting read
			    if (vIt->hap) ++hp1votes;
			    else ++hp2votes;
			  } else if ((sequence.substr(sp, vIt->alt.size()) != vIt->alt) && (sequence.substr(sp, vIt->ref.size()) == vIt->ref)) {
			    // REF supporting read
			    if (vIt->hap) ++hp2votes;
			    else ++hp1votes;
			  }
			} else if (vIt->ref.size() < vIt->alt.size()) {
			  // Insertion
			  int32_t diff = vIt->alt.size() - vIt->ref.size();
			  std::string refProbe = vIt->ref + std::string(seq + gp + vIt->ref.size(), seq + gp + vIt->ref.size() + diff);
			  if ((sequence.substr(sp, vIt->alt.size()) == vIt->alt) && (sequence.substr(sp, vIt->alt.size()) != refProbe)) {
			    // ALT supporting read
			    if (vIt->hap) ++hp1votes;
			    else ++hp2votes;
			  } else if ((sequence.substr(sp, vIt->alt.size()) != vIt->alt) && (sequence.substr(sp, vIt->alt.size()) == refProbe)) {
			    // REF supporting read
			    if (vIt->hap) ++hp2votes;
			    else ++hp1votes;
			  }
			} else {
			  // Deletion
			  int32_t diff = vIt->ref.size() - vIt->alt.size();
			  std::string altProbe = vIt->alt + std::string(seq + gp + vIt->ref.size(), seq + gp + vIt->ref.size() + diff);
			  if ((sequence.substr(sp, vIt->ref.size()) == altProbe) && (sequence.substr(sp, vIt->ref.size()) != vIt->ref)) {
			    // ALT supporting read
			    if (vIt->hap) ++hp1votes;
			    else ++hp2votes;
			  } else if ((sequence.substr(sp, vIt->ref.size()) != altProbe) && (sequence.substr(sp, vIt->ref.size()) == vIt->ref)) {
			    // REF supporting read
			    if (vIt->hap) ++hp2votes;
			    else ++hp1votes;
			  }
			}
		      }
		    }
		  }
		}
	      }
	    } else {
	      std::cerr << "Unknown Cigar options" << std::endl;
	      return 1;
	    }
	  }
	}
	int32_t hp = 0;
	if (hp1votes > 2*hp2votes) hp = 1;
	else if (hp2votes > 2*hp1votes) hp = 2;
	if (hp) {
	  if (hp == 1) h1.insert(hash_pair(rec));
	  else h2.insert(hash_pair(rec));
	}
      }
    }
    bam_destroy1(rec);
    hts_itr_destroy(iter);

    // Fetch all pairs
    hts_itr_t* itr = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
    bam1_t* r = bam_init1();
    // Annotate REF and ALT support
    typedef std::vector<uint32_t> TAlleleSupport;
    TAlleleSupport ref(pv.size(), 0);
    TAlleleSupport alt(pv.size(), 0);
    while (sam_itr_next(samfile, itr, r) >= 0) {
      if (r->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
      if ((r->core.qual < c.minMapQual) || (r->core.tid<0)) continue;
      bool h1Found = false;
      bool h2Found = false;
      if (h1.find(hash_pair(r)) != h1.end()) h1Found = true;
      if (h2.find(hash_pair(r)) != h2.end()) h2Found = true;
      if ((h1Found) && (h2Found)) {
	// Inconsistent haplotype assignment for this pair
	//std::cout << "Read\t" << bam_get_qname(r) << "\t" <<  hdr->target_name[r->core.tid] << "\t" << r->core.pos << "\t" << hdr->target_name[r->core.mtid] << "\t" << r->core.mpos << std::endl;
	++ambiguousReads;
      } else if ((!h1Found) && (!h2Found)) {
	++unassignedReads;
      } else {
	// Valid haplotype assignment
	TPhasedVariants::const_iterator vIt = std::lower_bound(pv.begin(), pv.end(), BiallelicVariant(r->core.pos), SortVariants<BiallelicVariant>());
	TPhasedVariants::const_iterator vItEnd = std::upper_bound(pv.begin(), pv.end(), BiallelicVariant(lastAlignedPosition(r)), SortVariants<BiallelicVariant>());
	if (vIt != vItEnd) {
	  // Get read sequence
	  std::string sequence;
	  sequence.resize(r->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(r);
	  for (int32_t i = 0; i < r->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	  // Get base qualities
	  typedef std::vector<uint8_t> TQuality;
	  TQuality quality;
	  quality.resize(r->core.l_qseq);
	  uint8_t* qualptr = bam_get_qual(r);
	  for (int i = 0; i < r->core.l_qseq; ++i) quality[i] = qualptr[i];
	  
	  // Parse CIGAR
	  uint32_t* cigar = bam_get_cigar(r);
	  for(;vIt != vItEnd; ++vIt) {
	    int32_t gp = r->core.pos; // Genomic position
	    int32_t sp = 0; // Sequence position
	    bool varFound = false;
	    for (std::size_t i = 0; ((i < r->core.n_cigar) && (!varFound)); ++i) {
	      if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) gp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
		//Nop
	      } else if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
		if (gp + (int32_t) bam_cigar_oplen(cigar[i]) < vIt->pos) {
		  gp += bam_cigar_oplen(cigar[i]);
		  sp += bam_cigar_oplen(cigar[i]);
		} else {
		  for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
		    if (gp == vIt->pos) {
		      varFound = true;
		      if (quality[sp] >= c.minBaseQual) {
			// Check REF allele
			if (vIt->ref == std::string(seq + gp, seq + gp + vIt->ref.size())) {
			  // Check ALT allele
			  if ((sp + vIt->alt.size() < sequence.size()) && (sp + vIt->ref.size() < sequence.size())) {
			    if (vIt->ref.size() == vIt->alt.size()) {
			      // SNP
			      if ((sequence.substr(sp, vIt->alt.size()) == vIt->alt) && (sequence.substr(sp, vIt->ref.size()) != vIt->ref)) ++alt[vIt-pv.begin()];
			      else if ((sequence.substr(sp, vIt->alt.size()) != vIt->alt) && (sequence.substr(sp, vIt->ref.size()) == vIt->ref)) ++ref[vIt-pv.begin()];
			    }
			  } else if (vIt->ref.size() < vIt->alt.size()) {
			    // Insertion
			    int32_t diff = vIt->alt.size() - vIt->ref.size();
			    std::string refProbe = vIt->ref + std::string(seq + gp + vIt->ref.size(), seq + gp + vIt->ref.size() + diff);
			    if ((sequence.substr(sp, vIt->alt.size()) == vIt->alt) && (sequence.substr(sp, vIt->alt.size()) != refProbe)) ++alt[vIt-pv.begin()];
			    else if ((sequence.substr(sp, vIt->alt.size()) != vIt->alt) && (sequence.substr(sp, vIt->alt.size()) == refProbe)) ++ref[vIt-pv.begin()];
			  } else {
			    // Deletion
			    int32_t diff = vIt->ref.size() - vIt->alt.size();
			    std::string altProbe = vIt->alt + std::string(seq + gp + vIt->ref.size(), seq + gp + vIt->ref.size() + diff);
			    if ((sequence.substr(sp, vIt->ref.size()) == altProbe) && (sequence.substr(sp, vIt->ref.size()) != vIt->ref)) ++alt[vIt-pv.begin()];
			    else if ((sequence.substr(sp, vIt->ref.size()) != altProbe) && (sequence.substr(sp, vIt->ref.size()) == vIt->ref)) ++ref[vIt-pv.begin()];
			  }
			}
		      }
		    }
		  }
		}
	      }
	      else {
		std::cerr << "Unknown Cigar options" << std::endl;
		return 1;
	      }
	    }
	  }
	}
	// Output read to haplotype specific BAM
	if (h1Found) {
	  int32_t hp = 1;
	  ++assignedReadsH1;
	  bam_aux_append(r, "HP", 'i', 4, (uint8_t*)&hp);
	  if (!sam_write1(h1bam, hdr, r)) {
	    std::cerr << "Could not write to bam file!" << std::endl;
	    return -1;
	  }
	} else if (h2Found) {
	  int32_t hp = 2;
	  ++assignedReadsH2;
	  bam_aux_append(r, "HP", 'i', 4, (uint8_t*)&hp);
	  if (!sam_write1(h2bam, hdr, r)) {
	    std::cerr << "Could not write to bam file!" << std::endl;
	    return -1;
	  }
	}
      }
    }
    bam_destroy1(r);
    hts_itr_destroy(itr);
    if (seqlen) free(seq);

    // Output phased allele support
    hts_itr_t* itervcf = bcf_itr_querys(bcfidx, bcfhdr, chrName.c_str());
    if (itervcf != NULL) {
      bcf1_t* recvcf = bcf_init1();
      for (uint32_t i = 0; i<pv.size(); ++i) {
	// Fetch variant annotation from VCF
	int32_t itrRet = 0;
	do {
	  itrRet = bcf_itr_next(ibcffile, itervcf, recvcf);
	  if (itrRet >= 0) {
	    bcf_unpack(recvcf, BCF_UN_SHR);
	    std::vector<std::string> alleles;
	    for(std::size_t k = 0; k<recvcf->n_allele; ++k) alleles.push_back(std::string(recvcf->d.allele[k]));
	    if ((recvcf->pos == pv[i].pos) && (pv[i].ref == alleles[0]) && (pv[i].alt == alleles[1])) break;
	  } else {
	    std::cerr << "Error: Variant not found! " << chrName << ":" << (pv[i].pos + 1) << std::endl;
	    return 1;
	  }
	} while (itrRet >= 0);
	double af1kgp = fetch1kGPAF(bcfhdr, recvcf);
	uint32_t totalcov = ref[i] + alt[i];
	std::string hapstr;
	if (pv[i].hap) hapstr = "1|0";
	else hapstr = "0|1";
	if (totalcov > 0) {
	  double h1af = 0;
	  double vaf = (double) alt[i] / (double) totalcov;
	  if (pv[i].hap) h1af = (double) alt[i] / (double) totalcov;
	  else h1af = (double) ref[i] / (double) totalcov;
	  double pval = binomTest(alt[i], totalcov, 0.5);
	  dataOut << chrName << "\t" << (pv[i].pos + 1) << "\t" << recvcf->d.id << "\t" << pv[i].ref << "\t" << pv[i].alt << "\t" << totalcov << "\t" << ref[i] << "\t" << alt[i] << "\t" << hapstr << "\t" << vaf << "\t" << h1af << "\t" << pval << "\t" << af1kgp << std::endl;
	} else {
	  // No coverage
	  dataOut << chrName << "\t" << (pv[i].pos + 1) << "\t" << recvcf->d.id << "\t" << pv[i].ref << "\t" << pv[i].alt << "\t" << totalcov << "\t" << ref[i] << "\t" << alt[i] << "\t" << hapstr << "\tNA\tNA\tNA" << "\t" << af1kgp << std::endl;
	}
      }
      bcf_destroy(recvcf);
      hts_itr_destroy(itervcf);
    }
  }
  fai_destroy(fai);

  // Close bam
  bam_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  sam_close(samfile);

  // Close output BAMs
  sam_close(h1bam);
  sam_close(h2bam);
  bam_index_build(c.h1bam.string().c_str(), 0);
  bam_index_build(c.h2bam.string().c_str(), 0);

  // Close output allele file
  dataOut.pop();

  // Close BCF
  bcf_hdr_destroy(bcfhdr);
  hts_idx_destroy(bcfidx);
  bcf_close(ibcffile);
  
  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

  // Statistics
  uint64_t sumReads = assignedReadsH1 + assignedReadsH2 + unassignedReads + ambiguousReads;
  std::cout << "AssignedReadsH1=" << assignedReadsH1 << ", AssignedReadsH2=" << assignedReadsH2 << ", UnassignedReads=" << unassignedReads << ", AmbiguousReads=" << ambiguousReads << ", FractionReadsAssigned=" << (float) (assignedReadsH1 + assignedReadsH2) / (float) sumReads << ", FractionAmbiguousReads=" << (float) (ambiguousReads) / (float) (sumReads) << std::endl;


#ifdef PROFILE
  ProfilerStop();
#endif


  return 0;
}

int main(int argc, char **argv) {
  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("map-qual,m", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(10), "min. mapping quality")
    ("base-qual,b", boost::program_options::value<unsigned short>(&c.minBaseQual)->default_value(10), "min. base quality")
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "reference fasta file")
    ("hap1,p", boost::program_options::value<boost::filesystem::path>(&c.h1bam)->default_value("h1.bam"), "haplotype1 output file")
    ("hap2,q", boost::program_options::value<boost::filesystem::path>(&c.h2bam)->default_value("h2.bam"), "haplotype2 output file")
    ("sample,s", boost::program_options::value<std::string>(&c.sample)->default_value("NA12878"), "sample name")
    ("as,a", boost::program_options::value<boost::filesystem::path>(&c.as)->default_value("as.tsv.gz"), "allele-specific output file")
    ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input phased BCF file")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamfile), "input bam file")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome")) || (!vm.count("vcffile"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -g <ref.fa> -s NA12878 -v <snps.bcf> --as <as.tsv> --hap1 <h1.bam> --hap2 <h2.bam> <unphased.bam>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  }

  // Check input BAM file
  if (vm.count("input-file")) {
    if (!(boost::filesystem::exists(c.bamfile) && boost::filesystem::is_regular_file(c.bamfile) && boost::filesystem::file_size(c.bamfile))) {
      std::cerr << "Input BAM file is missing: " << c.bamfile.string() << std::endl;
      return 1;
    }
    samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
    if (samfile == NULL) {
      std::cerr << "Fail to open file " << c.bamfile.string() << std::endl;
      return 1;
    }
    hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
    if (idx == NULL) {
      std::cerr << "Fail to open index for " << c.bamfile.string() << std::endl;
      return 1;
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    if (hdr == NULL) {
      std::cerr << "Fail to open header for " << c.bamfile.string() << std::endl;
      return 1;
    }
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }
  
  // Check VCF/BCF file
  if (vm.count("vcffile")) {
    if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
      std::cerr << "Input SNP VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
      return 1;
    }
    htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
    if (ifile == NULL) {
      std::cerr << "Fail to open file " << c.vcffile.string() << std::endl;
      return 1;
    }
    hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
    if (bcfidx == NULL) {
      std::cerr << "Fail to open index file for " << c.vcffile.string() << std::endl;
      return 1;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);
    if (hdr == NULL) {
      std::cerr << "Fail to open header for " << c.vcffile.string() << std::endl;
      return 1;
    }
    bcf_hdr_destroy(hdr);
    hts_idx_destroy(bcfidx);
    bcf_close(ifile);
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  return phaseBamRun(c);
}
