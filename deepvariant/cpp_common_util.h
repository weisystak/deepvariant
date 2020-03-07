#ifndef CPPCOMMONUTIL
#define CPPCOMMONUTIL

#include "deepvariant/protos/deepvariant.pb.h"
namespace nucleus {

nucleus::genomics::v1::Range make_range(std::string chrom, long long start ,long long end);
bool position_overlaps(std::string chrom, long long position, nucleus::genomics::v1::Range interval);
nucleus::genomics::v1::Range read_range(nucleus::genomics::v1::Read& read);
std::vector<nucleus::genomics::v1::Read> query(std::vector<nucleus::genomics::v1::Read>& reads, nucleus::genomics::v1::Range& region);
int find_max_overlapping(nucleus::genomics::v1::Range& query_range, std::vector<nucleus::genomics::v1::Range>& search_ranges);
bool ranges_overlap(nucleus::genomics::v1::Range& i1, nucleus::genomics::v1::Range& i2);
}
#endif
