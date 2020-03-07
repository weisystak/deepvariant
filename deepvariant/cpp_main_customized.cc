#include <iostream>
#include <regex>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <tuple>
#include "third_party/libzmq/zmq.hpp"
#include "third_party/nucleus/io/indexed_fasta_reader.h"
#include "third_party/nucleus/io/sam_reader.h"
#include "third_party/nucleus/io/hts_verbose.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/variant_calling.h"
#include "deepvariant/pileup_image_native.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "cpp_main.h"
#include "cpp_common_util.h"   
#include "cpp_realigner.h"
#include <thread>
#include <boost/program_options.hpp>
namespace BPO = boost::program_options;

namespace nucleus {

void removeCharsFromString( std::string &str, std::string charsToRemove ) {
   for ( unsigned int i = 0; i < charsToRemove.size(); ++i ) {
      str.erase( std::remove(str.begin(), str.end(), charsToRemove[i]), str.end() );
   }
}



std::vector<nucleus::genomics::v1::Range> processing_regions_from_options(const learning::genomics::deepvariant::DeepVariantOptions& options) {

  StatusOr<std::unique_ptr<IndexedFastaReader>> fasta =
      IndexedFastaReader::FromFile(options.reference_filename(),
                                   options.reference_filename()+".fai");


  auto ref_contig = fasta.ValueOrDie() -> Contigs();
  nucleus::genomics::v1::SamReaderOptions samoptions;
  samoptions.set_random_seed(2928130004);
  samoptions.set_aux_field_handling(nucleus::genomics::v1::SamReaderOptions::SKIP_AUX_FIELDS);
  std::unique_ptr<SamReader> reader = std::move(SamReader::FromFile(options.reads_filename(), samoptions).ValueOrDie());
  auto header = reader->Header();
  std::cout << header.contigs_size() << std::endl;
  std::vector<nucleus::genomics::v1::ContigInfo> tmp;
  for (unsigned int i=0; i< ref_contig.size(); i++){
    bool a = true;
    for (auto j=0; j< options.exclude_contigs_size(); j++){
      if (ref_contig[i].name() == options.exclude_contigs(j)){
        a= false; 
        break;
      }
    }
    if (a)
      tmp.emplace_back(ref_contig[i]);

  }

  std::unordered_map<std::string, nucleus::genomics::v1::ContigInfo> contig_dict;
  for (auto i=0; i<header.contigs_size(); i++)
    contig_dict.insert({header.contigs(i).name(), header.contigs(i)});
  std::unordered_map<std::string, nucleus::genomics::v1::ContigInfo>::iterator it;
  std::vector<nucleus::genomics::v1::ContigInfo> common_contigs;

  for (unsigned int i=0; i< tmp.size(); i++){
    
    if ((it = contig_dict.find(tmp[i].name())) != contig_dict.end()){
      if (tmp[i].n_bases() == (it -> second).n_bases())
        common_contigs.emplace_back(tmp[i]);
    }

  }
  std::cout << "common contigs are: ";
  for (auto i : common_contigs)
    std::cout << i.name() << " ";
  std::cout << std::endl;

  std::vector<nucleus::genomics::v1::Range> rangeset;
  for (auto i : common_contigs){

    rangeset.emplace_back(make_range(i.name(),0,i.n_bases()));
  }
  
  // if --regions do intersection
  if (options.calling_regions_size() >0){
    std::vector<nucleus::genomics::v1::Range> regionsrangeset;
    std::vector<nucleus::genomics::v1::Range> finalrangeset;
    std::regex rgx("^(\\S+):([0-9,]+)-([0-9,]+)$");
    std::smatch matches;
    std::string s= options.calling_regions(0);
    //std::cout << s << std::endl;
    if(std::regex_search(s, matches, rgx)) {

      std::string start = matches[2];
      removeCharsFromString(start, ",");
      std::string end = matches[3];
      removeCharsFromString(end, ",");
      regionsrangeset.emplace_back(make_range(matches[1], stoll(start)-1, stoll(end)));
      std::unordered_map<std::string, nucleus::genomics::v1::Range> range_dict;
      for (auto i :regionsrangeset )
        range_dict.insert({i.reference_name(), i});
      std::unordered_map<std::string, nucleus::genomics::v1::Range>::iterator rit;
      for (auto i :rangeset ){
        if ((rit = range_dict.find(i.reference_name())) != range_dict.end())
          finalrangeset.emplace_back(make_range(i.reference_name(), std::max(i.start(),(rit->second).start()), std::min(i.end(),(rit->second).end())));
      }
      
      return finalrangeset;
      //std::cout << finalrangeset[0].start() << std::endl;

      
    } else {
        std::cout << "Could not parse --regions\n";
    }
  }

  return rangeset;
}





std::string extract_sample_name_from_sam_reader( std::string reads){


  nucleus::genomics::v1::SamReaderOptions options;
  options.set_random_seed(2928130004);
  options.set_aux_field_handling(nucleus::genomics::v1::SamReaderOptions::SKIP_AUX_FIELDS);
  std::unique_ptr<SamReader> reader = std::move(SamReader::FromFile(reads, options).ValueOrDie());
  auto header = reader->Header();
  auto id = header.read_groups(0).sample_id();
  return id;

}

std::string get_reference_bases(nucleus::genomics::v1::Variant& variant, 
                                const learning::genomics::deepvariant::PileupImageOptions& option, 
								std::unique_ptr<IndexedFastaReader> & fasta){
  int half_width = (option.width()-1)/2;
  int start = variant.start() - half_width;
  int end = start + option.width();
  auto region = make_range(variant.reference_name(), start, end);
  if(fasta -> IsValidInterval(region))
    return fasta -> GetBases(region).ValueOrDie();
  else
    return "";

}

}


void make_example(  learning::genomics::deepvariant::DeepVariantOptions options, std::vector<nucleus::genomics::v1::Range> finalrangeset, 
                     int hts_block_size, int task_id, int task_number){

  zmq::context_t context (1);
  zmq::socket_t socket (context, ZMQ_PUSH);
  std::string address="tcp://localhost:";
  address+=std::to_string(5558+task_id);
  socket.connect(address);



  //-----------------start processing---------------------------------------------------------
  nucleus::genomics::v1::SamReaderOptions samoptions;
  samoptions.set_downsample_fraction(options.downsample_fraction());
  samoptions.set_random_seed(609314161);
  samoptions.set_aux_field_handling(nucleus::genomics::v1::SamReaderOptions::SKIP_AUX_FIELDS);
  samoptions.set_hts_block_size(hts_block_size);
  auto srqp = samoptions.mutable_read_requirements();
  srqp -> set_min_base_quality(10);
  srqp -> set_min_mapping_quality(10);
  srqp -> set_min_base_quality_mode(nucleus::genomics::v1::ReadRequirements::ENFORCED_BY_CLIENT);
  
  std::unique_ptr<nucleus::SamReader> samreader = std::move(nucleus::SamReader::FromFile(options.reads_filename(), samoptions).ValueOrDie());
  std::srand(options.pic_options().random_seed());  
  
 // std::cout <<b[0].alignment().position().position() << std::endl;
  if (options.realigner_enabled())
    std::cout <<"realigner_enabled!" << std::endl;
  nucleus::StatusOr<std::unique_ptr<nucleus::IndexedFastaReader>> fasta =
      nucleus::IndexedFastaReader::FromFile(options.reference_filename(),
                                   options.reference_filename()+".fai");

  unsigned int max_reads=options.pic_options().height() - options.pic_options().reference_band_height();
  int examplecount=0;
  int rangecount=0;
  for (auto& finalrange:finalrangeset){
    auto finalreference = finalrange.reference_name();
    for (int finalindex=finalrange.start(); finalindex<finalrange.end(); finalindex+=options.allele_counter_options().partition_size()){
      if (rangecount%task_number==task_id){
        auto finalend = finalindex+1000<finalrange.end()?finalindex+1000:finalrange.end();
        auto region = nucleus::make_range(finalreference, finalindex, finalend);
        auto b= as_vector_reservoir_sample(samreader->Query(region), options.max_reads_per_partition());
        //std::cout <<b.size() << std::endl;

        auto in_memory_reads =nucleus::realign_reads(options.realigner_options(), b,  fasta.ValueOrDie(), region);
        //std::cout <<nnnnn.size() << std::endl;
        auto allele_counter = std::unique_ptr<learning::genomics::deepvariant::AlleleCounter>(
            new learning::genomics::deepvariant::AlleleCounter(fasta.ValueOrDie().get(), region, options.allele_counter_options()));
        for (unsigned int j=0; j< in_memory_reads.size(); j++){
          auto r=nucleus::read_range(in_memory_reads[j]);
          if (nucleus::ranges_overlap(region, r))
            allele_counter -> Add(in_memory_reads[j]);
        }
        learning::genomics::deepvariant::VariantCaller variant_caller(options.variant_caller_options());
        auto candidates = variant_caller.CallsFromAlleleCounter(*allele_counter);
        //std::cout <<candidates.size() << std::endl;
        learning::genomics::deepvariant::PileupImageEncoderNative encoder(options.pic_options());
        for (auto dv_call:candidates){
          auto variant=dv_call.variant();
          auto ref = nucleus::get_reference_bases(variant, options.pic_options(), fasta.ValueOrDie());
          //std::cout <<ref.size() << std::endl;
          if (ref.empty())
            continue;
          auto image_start_pos = variant.start() - (options.pic_options().width()-1)/2;
          auto query_start = variant.start() - options.pic_options().read_overlap_buffer_bp();
          auto query_end = variant.end() + options.pic_options().read_overlap_buffer_bp();
          auto queryregion = nucleus::make_range(variant.reference_name(), query_start, query_end);
          std::vector<std::vector<std::string>> alt_allele_combinations;
          for (unsigned int alti=0; alti<variant.alternate_bases_size()+1;alti++){
            std::string tmp="";
            if (alti>0)
              tmp=variant.alternate_bases(alti-1);
            for (unsigned int altj=alti+1; altj < variant.alternate_bases_size()+1; altj++){
              std::vector<std::string> tmpvector;
              if(!tmp.empty())
                tmpvector.push_back(tmp);
              tmpvector.push_back(variant.alternate_bases(altj-1));
              alt_allele_combinations.push_back(tmpvector);
            }
          }
          for (auto& alt:alt_allele_combinations){
            learning::genomics::deepvariant::CallVariantsOutput::AltAlleleIndices altindices;
            for (auto& altstr:alt){
              for (int alti=0; alti<variant.alternate_bases_size();alti++){
                if (altstr==variant.alternate_bases(alti))
                  altindices.add_indices(alti);
              } 
            }
            std::string Image;
            auto imgptr = &Image;
            auto refimgrow = encoder.EncodeReference(ref);
            std::string refrow;
            for (int ii=0; ii<221; ii++){
              refrow+=refimgrow->base[ii];
              refrow+=refimgrow->base_quality[ii];
              refrow+=refimgrow->mapping_quality[ii];
              refrow+=refimgrow->on_positive_strand[ii];
              refrow+=refimgrow->supports_alt[ii];
              refrow+=refimgrow->matches_ref[ii];
            }
            for (int tmpi=0; tmpi<options.pic_options().reference_band_height(); tmpi++)
              *imgptr+=refrow;


            //reads rows
            std::vector<learning::genomics::deepvariant::ImageRow> reads_rows;
            std::vector<std::tuple<int, int>> sortindex;

            unsigned int count=0;
            unsigned int random = 0;

            for (auto& read : in_memory_reads){
              auto r=nucleus::read_range(read);
              if (nucleus::ranges_overlap(queryregion, r)){
                if (count < max_reads){
                  auto ptr = encoder.EncodeRead(dv_call, ref, read, image_start_pos, alt);
                  if (ptr){
                    reads_rows.emplace_back(*ptr);
                    sortindex.emplace_back(std::make_tuple(read.alignment().position().position(), count)); 
                    count++;
                  }
                }
                else{
                  random = std::rand()%(count+1);
                  if (random < max_reads){
                    auto ptr = encoder.EncodeRead(dv_call, ref, read, image_start_pos, alt);
                    if(ptr){  
                      reads_rows[random]=*ptr;
                      sortindex[random]=std::make_tuple(read.alignment().position().position(), random);
                    }
                  }
                }
              }
            }
            std::sort(sortindex.begin(), sortindex.end());
            for (auto& si:sortindex){
              for (int jj=0; jj<221; jj++){
                *imgptr+=reads_rows[std::get<1>(si)].base[jj];
                *imgptr+=reads_rows[std::get<1>(si)].base_quality[jj];
                *imgptr+=reads_rows[std::get<1>(si)].mapping_quality[jj];
                *imgptr+=reads_rows[std::get<1>(si)].on_positive_strand[jj];
                *imgptr+=reads_rows[std::get<1>(si)].supports_alt[jj];
                *imgptr+=reads_rows[std::get<1>(si)].matches_ref[jj]; 
              }
            }
            auto n_missing_rows = options.pic_options().height() - reads_rows.size()-options.pic_options().reference_band_height();
            std::vector<unsigned char> emptyrow(n_missing_rows*6*221,0);
            imgptr->append(emptyrow.begin(), emptyrow.end());
            //std::cout <<n_missing_rows << std::endl;
            //std::cout <<(int)rows[0].base[0] << " ";
            //std::cout <<(int)rows[0].base_quality[0] << " ";
            //std::cout <<(int)rows[0].mapping_quality[0] << " ";
            //std::cout <<(int)rows[0].on_positive_strand[0] << " ";
            //std::cout <<(int)rows[0].supports_alt[0] <<" ";
            //std::cout <<(int)rows[0].matches_ref[0] << std::endl;
            std::string variantstring;
            variant.SerializeToString(&variantstring);
            std::string altindicesstring;
            altindices.SerializeToString(&altindicesstring);

            examplecount++;
            //zmq::message_t request;

            //  Wait for next request from client
            //socket.recv (&request);
            //std::cout << Image.alt_allele_indices()<< std::endl;

       
            //  Send reply back to client
            
            zmq::message_t reply1 (Image.size());
            memcpy (reply1.data (), Image.c_str(), Image.size());
            socket.send (reply1);
            zmq::message_t reply2 (variantstring.size());
            memcpy (reply2.data (), variantstring.c_str(), variantstring.size());
            socket.send (reply2);
            zmq::message_t reply3 (altindicesstring.size());
            memcpy (reply3.data (), altindicesstring.c_str(), altindicesstring.size());
            socket.send (reply3);
            
            //std::cout <<imgptr->size() << std::endl;

            std::cout <<examplecount <<"" << std::endl;
        
          }
        }
        rangecount++;
      }
      else{
        rangecount++;


      }

    }  
  }
  
  zmq::message_t reply1 (4);
  memcpy (reply1.data (), "Done", 4);
  socket.send (reply1);
  
}



int main( int argc, char** argv )
{


  bool realign_reads, write_run_info;
  int gvcf_gq_binsize, task, num_task, partition_size, max_reads_per_partition, hts_block_size, vsc_min_count_snps,
      vsc_min_count_indels, pileup_image_height, pileup_image_width, logging_every_n_candidates;
  float downsample_fraction, vsc_min_fraction_snps, vsc_min_fraction_indels, training_random_emit_ref_sites;
  std::string candidates, regions, exclude_regions, gvcf, confident_regions, truth_variants, multi_allelic_mode, 
              sample_name, hts_logging_level, labeler_algorithm, customized_classes_labeler_classes_list,
              customized_classes_labeler_info_field_name;

  bool ws_use_window_selector_model, emit_realigned_reads, use_fast_pass_aligner;
  std::string ws_window_selector_model, realigner_diagnostics;
  int ws_min_num_supporting_reads, ws_max_num_supporting_reads, ws_min_mapq, ws_min_base_quality, ws_min_windows_distance,
      ws_max_window_size, ws_region_expansion_in_bp, dbg_min_k, dbg_max_k, dbg_step_k, dbg_min_mapq, dbg_min_base_quality,
      dbg_min_edge_weight, dbg_max_num_paths, aln_match, aln_mismatch, aln_gap_open, aln_gap_extend, aln_k, max_num_mismatches,
      kmer_size;
  float aln_error_rate, realignment_similarity_threshold;
  int thread_number, total_thread, thread_offset;

  

  // setup program options description
  BPO::options_description bOptions( "MakeExamples Options" );
  bOptions.add_options() //remember add ; after last item
      ( "help", "Produce help message" )
      ( "ref", BPO::value<std::string>(),
        "Required. Genome reference to use. Must have an associated FAI index as well." 
        " Supports text or gzipped references. Should match the reference used to align "
        " the BAM file provided to --reads.")
      ( "reads", BPO::value<std::string>(), 
        "Required. Aligned, sorted, indexed BAM file containing the reads we want to call."
        " Should be aligned to a reference genome compatible with --ref." )
      ( "candidates", BPO::value<std::string>(&candidates)->default_value(""),
        "Candidate DeepVariantCalls in tfrecord format. For DEBUGGING.")
      ( "mode", BPO::value<std::string>(),
        "Mode to run. Must be one of calling or training")
      ( "regions", BPO::value<std::string>(&regions)->default_value(""),
        "Optional. Space-separated list of regions we want to process. Elements "
        "can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE files.")
      ( "exclude_regions", BPO::value<std::string>(&exclude_regions)->default_value(""),
        "Optional. Space-separated list of regions we want to exclude from "
        "processing. Elements can be region literals (e.g., chr20:10-20) or paths "
        "to BED/BEDPE files. Region exclusion happens after processing the "
        "--regions argument, so --region 20 --exclude_regions 20:100 does "
        "everything on chromosome 20 excluding base 100")
      ("gvcf", BPO::value<std::string>(&gvcf)->default_value(""),
       "Optional. Path where we should write gVCF records in TFRecord of Variant proto format.")
      ("gvcf_gq_binsize", BPO::value<int>(&gvcf_gq_binsize)->default_value(5),
       "Bin size in which to quantize gVCF genotype qualities. Larger bin size "
       "reduces the number of gVCF records at a loss of quality granularity.")
      ("confident_regions", BPO::value<std::string>(&confident_regions)->default_value(""),
       "Regions that we are confident are hom-ref or a variant in BED format. In "
       "BED or other equivalent format, sorted or unsorted. Contig names must "
       "match those of the reference genome.")
      ("truth_variants", BPO::value<std::string>(&truth_variants)->default_value(""),
       "Tabix-indexed VCF file containing the truth variant calls for this labels "
       "which we use to label our examples.")
      ("task", BPO::value<int>(&task)->default_value(0),
       "Task ID of this task")
      ("num_task", BPO::value<int>(&num_task)->default_value(1),
       "Number of tasks")
      ("partition_size", BPO::value<int>(&partition_size)->default_value(1000),
       "The maximum number of basepairs we will allow in a region before splitting"
       "it into multiple smaller subregions.")
      ("max_reads_per_partition", BPO::value<int>(&max_reads_per_partition)->default_value(1500),
       "The maximum number of reads per partition that we consider before "
       "following processing such as sampling and realigner.")
      ("multi_allelic_mode", BPO::value<std::string>(&multi_allelic_mode)->default_value(""), 
       "How to handle multi-allelic candidate variants. For DEBUGGING")
      ("realign_reads", BPO::value<bool>(&realign_reads)->default_value(true),
       "If True, locally realign reads before calling variants.")
      ("write_run_info", BPO::value<bool>(&write_run_info)->default_value(true),
       "If True, write out a MakeExamplesRunInfo proto besides our examples in "
       "text_format.")
      ("downsample_fraction", BPO::value<float>(&downsample_fraction)->default_value(0.0),
       "If not 0.0 must be a value between 0.0 and 1.0. "
       "Reads will be kept (randomly) with a probability of downsample_fraction "
       "from the input BAM. This argument makes it easy to create examples as "
       "though the input BAM had less coverage.")
      ("sample_name", BPO::value<std::string>(&sample_name)->default_value(""),
       "Variant/DeepVariantCall protos. If not specified, will be inferred from "
       "the header information from --reads.")
      ("hts_logging_level", BPO::value<std::string>(&hts_logging_level)->default_value("HTS_LOG_WARNING"),
        "Sets the htslib logging threshold.")
      ("hts_block_size", BPO::value<int>(&hts_block_size)->default_value(134217728),
       "Sets the htslib block size. Zero or negative uses default htslib setting; "
       "larger values (e.g. 1M) may be beneficial for using remote files. "
       "Currently only applies to SAM/BAM reading.")
      ("vsc_min_count_snps", BPO::value<int>(&vsc_min_count_snps)->default_value(2),
       "SNP alleles occurring at least this many times in our "
       "AlleleCount will be advanced as candidates.")
      ("vsc_min_count_indels", BPO::value<int>(&vsc_min_count_indels)->default_value(2),
       "Indel alleles occurring at least this many times in our "
       "AlleleCount will be advanced as candidates.")
      ("vsc_min_fraction_snps", BPO::value<float>(&vsc_min_fraction_snps)->default_value(0.12),
       "SNP alleles occurring at least this fraction of all "
       "counts in our AlleleCount will be advanced as candidates.")
      ("vsc_min_fraction_indels", BPO::value<float>(&vsc_min_fraction_indels)->default_value(0.12),
       "Indel alleles occurring at least this fraction of all "
       "counts in our AlleleCount will be advanced as candidates.")
      ("training_random_emit_ref_sites", BPO::value<float>(&training_random_emit_ref_sites)->default_value(0.0),
       "If > 0, emit extra random reference examples with this probability.")
      ("pileup_image_height", BPO::value<int>(&pileup_image_height)->default_value(0),
        "Height for the pileup image. If 0, uses the default height")
      ("pileup_image_width", BPO::value<int>(&pileup_image_width)->default_value(0),
        "Width for the pileup image. If 0, uses the default width")
      ("labeler_algorithm", BPO::value<std::string>(&labeler_algorithm)->default_value("haplotype_labeler"),
       "Algorithm to use to label examples in training mode. Must be one of the "
       "LabelerAlgorithm enum values in the DeepVariantOptions proto.")
      ("customized_classes_labeler_classes_list", BPO::value<std::string>(&customized_classes_labeler_classes_list)->default_value(""),
       "A comma-separated list of strings that defines customized class labels "
       "for variants. This is only set when labeler_algorithm is "
       "customized_classes_labeler.")
      ("customized_classes_labeler_info_field_name", BPO::value<std::string>(&customized_classes_labeler_info_field_name)->default_value(""),
       "The name from the INFO field of VCF where we should get the customized "
       "class labels from. This is only set when labeler_algorithm is "
       "customized_classes_labeler.")
      ("logging_every_n_candidates", BPO::value<int>(&logging_every_n_candidates)->default_value(100),
       "Print out the log every n candidates. The smaller the number, the more "
       "frequent the logging information emits.")
      //realigner options-------------------------------------------------------------------------
      ("ws_use_window_selector_model", BPO::value<bool>(&ws_use_window_selector_model)->default_value(true),
       "Activate the use of window selector models.")
      ("ws_window_selector_model", BPO::value<std::string>(&ws_window_selector_model)->default_value(""),
       "Path to a text format proto of the window selector model to use.")
      ("ws_min_num_supporting_reads", BPO::value<int>(&ws_min_num_supporting_reads)->default_value(-1),
       "Minimum number of supporting reads to call a reference position for local assembly.")
      ("ws_max_num_supporting_reads", BPO::value<int>(&ws_max_num_supporting_reads)->default_value(-1),
       "Maximum number of supporting reads to call a reference position for local assembly.")
      ("ws_min_mapq", BPO::value<int>(&ws_min_mapq)->default_value(20), 
       "Minimum read alignment quality to consider in calling a reference "
       "position for local assembly.")
      ("ws_min_base_quality", BPO::value<int>(&ws_min_base_quality)->default_value(20),
       "Minimum base quality to consider in calling a reference position for "
       "local assembly.")
      ("ws_min_windows_distance", BPO::value<int>(&ws_min_windows_distance)->default_value(80),
       "Minimum distance between candidate windows for local assembly.")
      ("ws_max_window_size", BPO::value<int>(&ws_max_window_size)->default_value(1000),
       "Maximum window size to consider for local assembly. Large noisy regions "
       "are skipped for realignment.")
      ("ws_region_expansion_in_bp", BPO::value<int>(&ws_region_expansion_in_bp)->default_value(20),
       "Number of bases to expand the region when calculating windows; larger "
       "values add overhead but allow larger nearby events to contribute evidence "
       "for assembling an region even if they are not contained by the region.")
      ("dbg_min_k", BPO::value<int>(&dbg_min_k)->default_value(10),
       "Initial k-mer size to build the graph.")
      ("dbg_max_k", BPO::value<int>(&dbg_max_k)->default_value(101),
       "Maximum k-mer size. Larger k-mer size is used to resolve graph cycles.")
      ("dbg_step_k", BPO::value<int>(&dbg_step_k)->default_value(1),
       "Increment size for k to try in resolving graph cycles.")
      ("dbg_min_mapq", BPO::value<int>(&dbg_min_mapq)->default_value(14),
       "Minimum read alignment quality to consider in building the graph.")
      ("dbg_min_base_quality", BPO::value<int>(&dbg_min_base_quality)->default_value(15),
       "Minimum base quality in a k-mer sequence to consider in building the graph.")
      ("dbg_min_edge_weight", BPO::value<int>(&dbg_min_edge_weight)->default_value(2),
       "Minimum number of supporting reads to keep an edge.")
      ("dbg_max_num_paths", BPO::value<int>(&dbg_max_num_paths)->default_value(256),
       "Maximum number of paths within a graph to consider for realignment. "
       "Set max_num_paths to 0 to have unlimited number of paths.")
      ("aln_match", BPO::value<int>(&aln_match)->default_value(4),
       "Match score (expected to be a non-negative score).")
      ("aln_mismatch", BPO::value<int>(&aln_mismatch)->default_value(6),
       "Mismatch score (expected to be a non-negative score).")
      ("aln_gap_open", BPO::value<int>(&aln_gap_open)->default_value(8),
       "Gap open score (expected to be a non-negative score). "
       "Score for a gap of length g is -(gap_open + (g - 1) * gap_extend).")
      ("aln_gap_extend", BPO::value<int>(&aln_gap_extend)->default_value(2),
       "Gap extend score (expected to be a non-negative score). "
       "Score for a gap of length g is -(gap_open + (g - 1) * gap_extend).")
      ("aln_k", BPO::value<int>(&aln_k)->default_value(23),
       "k-mer size used to index target sequence.")
      ("aln_error_rate", BPO::value<float>(&aln_error_rate)->default_value(0.01),
       "Estimated sequencing error rate.")
      ("realigner_diagnostics", BPO::value<std::string>(&realigner_diagnostics)->default_value(""),
       "Root directory where the realigner should place diagnostic output (such as"
       " a dump of the DeBruijn graph, and a log of metrics reflecting the graph "
       "and  realignment to the haplotypes).  If empty, no diagnostics are output.")
      ("emit_realigned_reads", BPO::value<bool>(&emit_realigned_reads)->default_value(false),
       "If True, we will emit realigned reads if our realigner_diagnostics are "
       "also enabled.")
      ("use_fast_pass_aligner", BPO::value<bool>(&use_fast_pass_aligner)->default_value(true),
       "If True, fast_pass_aligner (improved performance) implementation is used ")
      ("max_num_mismatches", BPO::value<int>(&max_num_mismatches)->default_value(2),
       "Num of maximum allowed mismatches for quick read to "
       "haplotype alignment.")
      ("realignment_similarity_threshold", BPO::value<float>(&realignment_similarity_threshold)->default_value(0.16934),
       "Similarity threshold used in realigner in Smith-Waterman alignment.")
      ("kmer_size", BPO::value<int>(&kmer_size)->default_value(32),
       "K-mer size for fast pass alinger reads index.")
      ("thread_number", BPO::value<int>(&thread_number)->default_value(1),
       "Launch thread_number threads")
      ("total_thread", BPO::value<int>(&total_thread)->default_value(1),
       "Launched total threads in multiple node")
      ("thread_offset", BPO::value<int>(&thread_offset)->default_value(0),
       "Thread ID offset");






  // parse program options
  BPO::variables_map mVMap;
  BPO::store( BPO::parse_command_line( argc, argv, bOptions ), mVMap );
  BPO::notify( mVMap );
  if (!mVMap.count( "reads" )){
    std::cout << "Requires reads" << std::endl;
    return 1;
  }
  if (!mVMap.count( "mode" )){
    std::cout << "Requires mode" << std::endl;
    return 1;
  }
  if (!mVMap.count( "ref" )){
    std::cout << "Requires ref" << std::endl;
    return 1;
  }

  // output help message if required
  if( mVMap.count( "help" ) )
  {
    std::cout << bOptions << std::endl;
    return 1;
  }

  htsLogLevel log = HTS_LOG_WARNING;
  nucleus::HtsSetLogLevel(log);
  std::cout << log << std::endl;
  learning::genomics::deepvariant::DeepVariantOptions options;
  auto rqp = options.mutable_read_requirements();
  rqp -> set_min_base_quality(10);
  rqp -> set_min_mapping_quality(10);
  rqp -> set_min_base_quality_mode(nucleus::genomics::v1::ReadRequirements::ENFORCED_BY_CLIENT);
  auto poq = options.mutable_pic_options();
  auto poq_rqp = poq -> mutable_read_requirements();
  poq_rqp -> set_min_base_quality(10);
  poq_rqp -> set_min_mapping_quality(10);
  poq_rqp -> set_min_base_quality_mode(nucleus::genomics::v1::ReadRequirements::ENFORCED_BY_CLIENT);
  poq ->set_reference_band_height(5);
  poq ->set_base_color_offset_a_and_g(40);
  poq ->set_base_color_offset_t_and_c(30);
  poq ->set_base_color_stride(70);
  poq ->set_allele_supporting_read_alpha(1.0);
  poq ->set_allele_unsupporting_read_alpha(0.6);
  poq ->set_reference_matching_read_alpha(0.2);
  poq ->set_reference_mismatching_read_alpha(1.0);
  poq ->set_indel_anchoring_base_char("*");
  poq ->set_reference_alpha(0.4);
  poq ->set_reference_base_quality(60);
  poq ->set_positive_strand_color(70);
  poq ->set_negative_strand_color(240);
  poq ->set_base_quality_cap(40);
  poq ->set_mapping_quality_cap(60);
  poq ->set_height(100);
  poq ->set_width(221);
  poq ->set_num_channels(6);
  poq ->set_read_overlap_buffer_bp(5);
  poq ->set_multi_allelic_mode(learning::genomics::deepvariant::PileupImageOptions::ADD_HET_ALT_IMAGES);
  poq ->set_random_seed(2101079370);
  auto acp = options.mutable_allele_counter_options();
  acp -> set_partition_size(partition_size);
  auto acp_rqp = acp -> mutable_read_requirements();
  acp_rqp -> set_min_base_quality(10);
  acp_rqp -> set_min_mapping_quality(10);
  acp_rqp -> set_min_base_quality_mode(nucleus::genomics::v1::ReadRequirements::ENFORCED_BY_CLIENT);
  std::string options_sample_name;
  if (!sample_name.empty())
    options_sample_name = sample_name;
  else if (!mVMap["reads"].as<std::string>().empty()){
    auto tmp= nucleus::extract_sample_name_from_sam_reader(mVMap["reads"].as<std::string>());
    options_sample_name = tmp;
  }
  else 
    options_sample_name="UNKNOWN";

  auto vcp = options.mutable_variant_caller_options();
  vcp->set_min_count_snps(vsc_min_count_snps);
  vcp->set_min_count_indels(vsc_min_count_indels);
  vcp->set_min_fraction_snps(vsc_min_fraction_snps);
  vcp->set_min_fraction_indels(vsc_min_fraction_indels);
  vcp->set_random_seed(1400605801);
  vcp->set_sample_name(options_sample_name);
  vcp->set_p_error(0.001);
  vcp->set_max_gq(50);
  vcp->set_gq_resolution(gvcf_gq_binsize);
  vcp->set_ploidy(2);
  for (uint i=0; i<EXCLUDED_HUMAN_CONTIGS.size(); i++)
    options.add_exclude_contigs(EXCLUDED_HUMAN_CONTIGS[i]);
  options.set_random_seed(609314161);
  options.set_n_cores(1);
  options.set_task_id(task);
  options.set_num_shards(0);
  options.set_min_shared_contigs_basepairs(0.9);
  if (mVMap["mode"].as<std::string>() == "calling")
    options.set_mode(learning::genomics::deepvariant::DeepVariantOptions::CALLING);
  else if (mVMap["mode"].as<std::string>() == "training")
    options.set_mode(learning::genomics::deepvariant::DeepVariantOptions::TRAINING);
  else {
    std::cout << "mode must be calling or training" << std::endl;
    return 1;
  }
  options.set_labeler_algorithm(learning::genomics::deepvariant::DeepVariantOptions::HAPLOTYPE_LABELER);
  options.set_reference_filename(mVMap["ref"].as<std::string>());
  options.set_reads_filename(mVMap["reads"].as<std::string>());
  if (!confident_regions.empty())
    options.set_confident_regions_filename(confident_regions);
  if (!truth_variants.empty())
    options.set_truth_variants_filename(truth_variants);
  if (downsample_fraction != 0.0)
    options.set_downsample_fraction(downsample_fraction);
  if (pileup_image_height)
    poq ->set_height(pileup_image_height);
  if (pileup_image_width)
    poq ->set_width(pileup_image_width);
  options.set_num_shards(num_task);
  if (write_run_info)
    options.set_run_info_filename("tmp1111111111.run_info.pbtxt");
  if (!regions.empty())
    options.add_calling_regions(regions);
  options.add_exclude_calling_regions(exclude_regions);
  options.set_realigner_enabled(realign_reads);
  options.set_max_reads_per_partition(max_reads_per_partition);
  auto rap = options.mutable_realigner_options();
  auto wsp = rap -> mutable_ws_config();
  auto wsmp = wsp -> mutable_window_selector_model();
  wsmp -> set_model_type(learning::genomics::deepvariant::WindowSelectorModel::ALLELE_COUNT_LINEAR);
  auto wsmpacp = wsmp -> mutable_allele_count_linear_model();
  wsmpacp -> set_bias(-0.683379);
  wsmpacp -> set_coeff_soft_clip(2.997000);
  wsmpacp -> set_coeff_substitution(-0.086644);
  wsmpacp -> set_coeff_insertion(2.493585);
  wsmpacp -> set_coeff_deletion(1.795914);
  wsmpacp -> set_coeff_reference(-0.059787);
  wsmpacp -> set_decision_boundary(3);
  wsp -> set_min_mapq(ws_min_mapq);
  wsp -> set_min_base_quality(ws_min_base_quality);
  wsp -> set_min_windows_distance(ws_min_windows_distance);
  wsp -> set_max_window_size(ws_max_window_size);
  wsp -> set_region_expansion_in_bp(ws_region_expansion_in_bp);
  auto dbgp = rap -> mutable_dbg_config();
  dbgp -> set_min_k(dbg_min_k);
  dbgp -> set_max_k(dbg_max_k);
  dbgp -> set_step_k(dbg_step_k);
  dbgp -> set_min_mapq(dbg_min_mapq);
  dbgp -> set_min_base_quality(dbg_min_base_quality);
  dbgp -> set_min_edge_weight(dbg_min_edge_weight);
  dbgp -> set_max_num_paths(dbg_max_num_paths);
  auto alnp = rap -> mutable_aln_config();
  alnp -> set_match(aln_match);
  alnp -> set_mismatch(aln_mismatch);
  alnp -> set_gap_open(aln_gap_open);
  alnp -> set_gap_extend(aln_gap_extend);
  alnp -> set_k(aln_k);
  alnp -> set_error_rate(aln_error_rate);
  alnp -> set_max_num_of_mismatches(max_num_mismatches);
  alnp -> set_realignment_similarity_threshold(realignment_similarity_threshold);
  alnp -> set_kmer_size(kmer_size);
  auto diagp = rap -> mutable_diagnostics();
  diagp -> set_enabled(bool(!realigner_diagnostics.empty()));
  diagp -> set_output_root(realigner_diagnostics);
  diagp -> set_emit_realigned_reads(emit_realigned_reads);

  auto finalrangeset = nucleus::processing_regions_from_options(options);
  std::cout << finalrangeset[0].start() << std::endl;
  std::thread t[thread_number];
  //Launch a group of threads
  for (int i = 0; i < thread_number; ++i){
    t[i] = std::thread(make_example, options, finalrangeset, hts_block_size, i+thread_offset, total_thread);
  }
  for (int i = 0; i < thread_number; ++i) 
    t[i].join();



  
  
}





  



