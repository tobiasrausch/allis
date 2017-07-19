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

}

#endif
