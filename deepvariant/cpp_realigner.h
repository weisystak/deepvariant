#ifndef CPPREALIGNER
#define CPPREALIGNER


#include "deepvariant/protos/deepvariant.pb.h"
#include "cpp_common_util.h"
#include "third_party/nucleus/io/indexed_fasta_reader.h"
#include <unordered_map>
#include "deepvariant/allelecounter.h"
#include <vector>
#include <iostream>
#include <memory>
namespace nucleus{

std::vector<nucleus::genomics::v1::Range> select_windows(const learning::genomics::deepvariant::WindowSelectorOptions& config, nucleus::genomics::v1::Range region, std::unique_ptr<IndexedFastaReader> & fasta,
      std::vector<nucleus::genomics::v1::Read> & reads);


std::vector<nucleus::genomics::v1::Read> realign_reads(const learning::genomics::deepvariant::RealignerOptions& config, std::vector<nucleus::genomics::v1::Read> & reads, std::unique_ptr<IndexedFastaReader> & fasta, 
                   nucleus::genomics::v1::Range region);



}
#endif
