#include "cpp_common_util.h"
#include <algorithm>

nucleus::genomics::v1::Range nucleus::make_range(std::string chrom, long long start ,long long end){

  nucleus::genomics::v1::Range range;
  range.set_reference_name(chrom); 
  range.set_start(start); 
  range.set_end(end);

  return range;
}

bool nucleus::position_overlaps(std::string chrom, long long position, nucleus::genomics::v1::Range interval){

  return (interval.reference_name()==chrom && position >= interval.start() && position < interval.end() );


}

nucleus::genomics::v1::Range nucleus::read_range(nucleus::genomics::v1::Read& read){

  long long start = read.alignment().position().position();
  long long sum=0;
  for ( int i=0; i<read.alignment().cigar_size(); i++){
    if (read.alignment().cigar(i).operation()==nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH ||
	    read.alignment().cigar(i).operation()==nucleus::genomics::v1::CigarUnit::SEQUENCE_MATCH ||
 	    read.alignment().cigar(i).operation()==nucleus::genomics::v1::CigarUnit::DELETE ||
	    read.alignment().cigar(i).operation()==nucleus::genomics::v1::CigarUnit::SKIP || 
		read.alignment().cigar(i).operation()==nucleus::genomics::v1::CigarUnit::SEQUENCE_MISMATCH )
      sum += read.alignment().cigar(i).operation_length();  
  }
  long long end = start+ sum;
  return nucleus::make_range(read.alignment().position().reference_name(), start, end);

  

}
bool nucleus::ranges_overlap(nucleus::genomics::v1::Range& i1, nucleus::genomics::v1::Range& i2){

  return (i1.reference_name() == i2.reference_name() && i1.end() > i2.start() &&
            i1.start() < i2.end());
}
std::vector<nucleus::genomics::v1::Read> nucleus::query(std::vector<nucleus::genomics::v1::Read>& reads, nucleus::genomics::v1::Range& region){


 std::vector<nucleus::genomics::v1::Read> a;
 for (auto i:reads){
   auto r=nucleus::read_range(i);
   if (nucleus::ranges_overlap(region, r))
     a.emplace_back(i);
 }
 return a;


}
long long overlap_len(nucleus::genomics::v1::Range& range1, nucleus::genomics::v1::Range& range2){
  long int a=0;
  if (range1.reference_name() != range2.reference_name())
    return a;
  return std::max(a, (std::min(range1.end(), range2.end()) - std::max(range1.start(), range2.start())));



}

int nucleus::find_max_overlapping(nucleus::genomics::v1::Range& query_range, std::vector<nucleus::genomics::v1::Range>& search_ranges){
  
  // if not overlapping, return -1
  std::vector<long int> overlaps;
  for (auto i:search_ranges)
    overlaps.emplace_back(overlap_len(query_range, i));
  int index = std::distance(overlaps.begin(), std::max_element(overlaps.begin(), overlaps.end()));
  if (overlaps[index]==0)
    return -1;
  else
    return index;

}


