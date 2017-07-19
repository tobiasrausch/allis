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

#ifndef INFOPARSE_H
#define INFOPARSE_H

#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <htslib/sam.h>


namespace allis
{

  inline double
  fetch1kGPAF(bcf_hdr_t* hdr, bcf1_t* rec) {
    int32_t naf = 0;
    float* af = NULL;
    double ret = 0;
    if (bcf_get_info_float(hdr, rec, "1kGP_AF", &af, &naf) > 0) ret = *af;
    if (af != NULL) free(af);
    return ret;
  }

  template<typename TColumnMap, typename TStrVec>
  inline void
  fetchCSQ(bcf_hdr_t* hdr, bcf1_t* rec, TColumnMap const& cmap, TStrVec& rtv) {
    int32_t ncsq = 0;
    char* csq = NULL;
    if (bcf_get_info_string(hdr, rec, "CSQ", &csq, &ncsq) > 0) {
      std::string vep = std::string(csq);
      TStrVec mgenes;
      boost::split(mgenes, vep, boost::is_any_of(std::string(",")));
      for(typename TStrVec::const_iterator mgIt = mgenes.begin(); mgIt != mgenes.end(); ++mgIt) {
	rtv.clear();
	boost::split(rtv, *mgIt, boost::is_any_of(std::string("|")));
	std::string canonical("NA");
	if ((!rtv.empty()) && (rtv[cmap.find("CANONICAL")->second].size())) canonical = rtv[cmap.find("CANONICAL")->second];
	if (canonical == "YES") break;
      }
    }
    if (csq != NULL) free(csq);
  }


  template<typename TMap>
  inline bool
  getCSQ(std::string const& header, TMap& cmap) {
    std::string delimiters("\n");
    typedef std::vector<std::string> TStrParts;
    TStrParts lines;
    boost::split(lines, header, boost::is_any_of(delimiters));
    TStrParts::const_iterator itH = lines.begin();
    TStrParts::const_iterator itHEnd = lines.end();
    bool foundCSQ = false;
    for(;itH!=itHEnd; ++itH) {
      if (itH->find("##INFO=<ID=CSQ,")==0) {
	foundCSQ = true;
	std::string delim(",");
	TStrParts keyval;
	boost::split(keyval, *itH, boost::is_any_of(delim));
	TStrParts::const_iterator itKV = keyval.begin();
	TStrParts::const_iterator itKVEnd = keyval.end();
	for(;itKV != itKVEnd; ++itKV) {
	  size_t sp = itKV->find("=");
	  if (sp != std::string::npos) {
	    std::string field = itKV->substr(0, sp);
	    if (field == "Description") {
	      std::string desc = itKV->substr(sp+1, itKV->size() - sp - 2);
	      size_t colon = desc.find(":");
	      if (colon != std::string::npos) {
		std::string format = desc.substr(colon+2);
		TStrParts columns;
		boost::split(columns, format, boost::is_any_of(std::string("|")));
		TStrParts::const_iterator itC = columns.begin();
		TStrParts::const_iterator itCEnd = columns.end();
		int32_t i = 0;
		for(;itC != itCEnd; ++itC, ++i) cmap.insert(std::make_pair(*itC, i));
	      }
	    }
	  }
	}
      }
    }
    return foundCSQ;
  }
  
}

#endif
