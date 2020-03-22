#include "cpp_realigner.h"
#include "deepvariant/realigner/window_selector.h"
#include "deepvariant/realigner/debruijn_graph.h"
#include "deepvariant/realigner/fast_pass_aligner.h"

namespace nucleus{
std::vector<nucleus::genomics::v1::Range> select_windows(const learning::genomics::deepvariant::WindowSelectorOptions& config, nucleus::genomics::v1::Range region, std::unique_ptr<IndexedFastaReader> & fasta,
	                std::vector<nucleus::genomics::v1::Read> & reads){
  learning::genomics::deepvariant::AlleleCounterOptions allelecounteroption;
  auto rqp = allelecounteroption.mutable_read_requirements();
  rqp -> set_min_mapping_quality(config.min_mapq());
  rqp -> set_min_base_quality(config.min_base_quality());
  //expand region
  long long new_start = (region.start()-config.region_expansion_in_bp())>0?(region.start()-config.region_expansion_in_bp()):0;
  long long new_end = region.end()+config.region_expansion_in_bp();
  auto ref_contig = fasta -> Contigs();

  std::unordered_map<std::string, nucleus::genomics::v1::ContigInfo> contig_dict;
  for (auto i:ref_contig)
    contig_dict.insert({i.name(), i});
  
  auto it = contig_dict.find(region.reference_name());
  new_end = (new_end- it -> second.n_bases())<0? new_end:it -> second.n_bases();
  //std::cout << new_start << std::endl << new_end << std::endl;
  auto expanded_region = make_range(region.reference_name(), new_start, new_end);
  auto allele_counter = std::unique_ptr<learning::genomics::deepvariant::AlleleCounter>(
        new learning::genomics::deepvariant::AlleleCounter(fasta.get(), expanded_region, allelecounteroption));

  for (auto i:reads)
  	allele_counter -> Add(i);
  auto scores_vec = AlleleCountLinearWindowSelectorCandidates(*allele_counter, config.window_selector_model().allele_count_linear_model());

  std::vector<long long> candidates;
  for (int i =0; i < scores_vec.size(); i++){
    if(scores_vec[i] > config.window_selector_model().allele_count_linear_model().decision_boundary()
    	&& position_overlaps(region.reference_name(), expanded_region.start()+i, region))
      candidates.push_back(expanded_region.start()+i);
  }
  std::sort(candidates.begin(), candidates.end());
  long long start_pos=-1;
  long long end_pos=-1;
  std::vector<nucleus::genomics::v1::Range> windows;
  //std::cout << candidates.size() << std::endl;
  for (auto i:candidates){
    if (start_pos<0){
      start_pos=i;
      end_pos=i;
    }
    else if(i > end_pos+2*config.min_windows_distance()) {
      windows.emplace_back(make_range(region.reference_name(), start_pos-config.min_windows_distance(), end_pos+config.min_windows_distance()));
      start_pos=i;
      end_pos=i;
    }
    else 
      end_pos=i; 
  }
  if (start_pos > -1)
  	windows.emplace_back(make_range(region.reference_name(), start_pos-config.min_windows_distance(), end_pos+config.min_windows_distance()));
  return windows;

}

std::vector<learning::genomics::deepvariant::CandidateHaplotypes> call_debruijn_graph(const learning::genomics::deepvariant::RealignerOptions& config, std::vector<nucleus::genomics::v1::Range>& windows, std::vector<nucleus::genomics::v1::Read> & reads, std::unique_ptr<IndexedFastaReader> & fasta){
  std::vector<learning::genomics::deepvariant::CandidateHaplotypes> windows_haplotypes;
  for (auto window:windows){
    if (window.end() - window.start() > config.ws_config().max_window_size())
	  continue;
	if (!fasta -> IsValidInterval(window))
	  continue;
	auto ref = fasta -> GetBases(window).ValueOrDie();
    auto window_reads = query(reads, window);
	auto graph = learning::genomics::deepvariant::DeBruijnGraph::Build(ref, window_reads, config.dbg_config());
	std::vector<string> tmp{ref};
	std::vector<string> candidate_haplotypes;
	if (!graph)
      candidate_haplotypes = tmp;
	else 
	  candidate_haplotypes = graph -> CandidateHaplotypes();
	if (candidate_haplotypes.size()>0 && (candidate_haplotypes!=tmp)){

	  windows_haplotypes.emplace_back(learning::genomics::deepvariant::CandidateHaplotypes());
      *windows_haplotypes.back().mutable_span() = window;
	  for (auto i:candidate_haplotypes)
	    windows_haplotypes.back().add_haplotypes(i);
	}
  }
  return windows_haplotypes;
}

std::vector<nucleus::genomics::v1::Read> call_fast_pass_aligner(std::vector<nucleus::genomics::v1::Read> & reads, std::unique_ptr<IndexedFastaReader> & fasta, nucleus::genomics::v1::Range & region,
                            nucleus::genomics::v1::Range& read_span, learning::genomics::deepvariant::AlignerOptions alnconfig,
                            learning::genomics::deepvariant::CandidateHaplotypes& CandidateHaplotypes){
  //UNDO
  int _REF_ALIGN_MARGIN = 20;
  long int zero =0;
  auto contig = region.reference_name();
  auto ref_start= std::max(
        zero,
        std::min(read_span.start(), region.start()) -
        _REF_ALIGN_MARGIN);
  auto ref_end = std::min(
        fasta -> Contig(contig).ValueOrDie()-> n_bases(),
        std::max(read_span.end(), region.end()) +
        _REF_ALIGN_MARGIN);
  auto ref_prefix = fasta -> GetBases(make_range(contig, ref_start, region.start())).ValueOrDie();
  auto ref = fasta -> GetBases(region).ValueOrDie();
  string ref_suffix;
  if (ref_end <= region.end())
    return reads;
  else
    ref_suffix = fasta -> GetBases(make_range(contig, region.end(), ref_end)).ValueOrDie();
  auto ref_seq = ref_prefix + ref + ref_suffix;
  auto fast_pass_realigner = learning::genomics::deepvariant::FastPassAligner();
  alnconfig.set_read_size(reads[0].aligned_sequence().size());
  fast_pass_realigner.set_options(alnconfig);
  fast_pass_realigner.set_reference(ref_seq);
  fast_pass_realigner.set_ref_start(contig, ref_start);
  std::vector<string> haplotypes;
  for (int i=0; i< CandidateHaplotypes.haplotypes_size(); i++){
    haplotypes.emplace_back(ref_prefix + CandidateHaplotypes.haplotypes(i) + ref_suffix);
  } 
  fast_pass_realigner.set_haplotypes(haplotypes);
  return *fast_pass_realigner.AlignReads(reads);

}

struct minmax{

  long int minspan;
  long int maxspan;
} ;

std::vector<nucleus::genomics::v1::Read> realign_reads(const learning::genomics::deepvariant::RealignerOptions& config, std::vector<nucleus::genomics::v1::Read> & reads, std::unique_ptr<IndexedFastaReader> & fasta, 
                   nucleus::genomics::v1::Range region){
  
  auto candidate_windows = select_windows(config.ws_config(), region, fasta, reads);
  //std::cout << candidate_windows[0].start() << " "<<candidate_windows[0].end() << std::endl;
  auto candidate_haplotypes = call_debruijn_graph(config, candidate_windows, reads, fasta);
  if (candidate_haplotypes.size() == 0)
    return reads;
  std::vector<nucleus::genomics::v1::Range> asregions;
  //for (auto i:candidate_haplotypes)
    //std::cout << i.span().start() << "  "<<i.span().end() << std::endl;
  for (auto i:candidate_haplotypes)
    asregions.emplace_back(i.span());
  std::vector<nucleus::genomics::v1::Read> unassigned_reads;
  std::vector<std::vector<nucleus::genomics::v1::Read>> assigned_reads(asregions.size());
  std::vector<minmax> asspan(asregions.size());
  for (auto& i : asspan){
    i.minspan = 2147483647;
    i.maxspan = -1;
  }
  for (auto read:reads){
    auto readrange = read_range(read);
    int window_i = find_max_overlapping(readrange, asregions);
    if (window_i == -1)
      unassigned_reads.emplace_back(read);
    else{
      assigned_reads[window_i].emplace_back(read);
      if (readrange.start()<asspan[window_i].minspan)
        asspan[window_i].minspan = readrange.start();
      if (readrange.end()>asspan[window_i].maxspan)
        asspan[window_i].maxspan = readrange.end();
    }
  }
  std::vector<nucleus::genomics::v1::Range> read_span;
  for (auto i : asspan){
    read_span.emplace_back(make_range(region.reference_name(), i.minspan, i.maxspan));
  }

  //std::cout << unassigned_reads.size() << std::endl;
  for (unsigned int i=0; i<assigned_reads.size(); i++){
    auto asreads =call_fast_pass_aligner(assigned_reads[i], fasta, asregions[i], read_span[i], config.aln_config(), candidate_haplotypes[i]);
    unassigned_reads.insert(unassigned_reads.end(), asreads.begin(), asreads.end());
  }

  return unassigned_reads;
 

  



}

}
