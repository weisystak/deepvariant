THREAD=20
TOTAL_THREAD=20
OFFSET=0
# LOCATION=/home/yang

BASE="/disk1/case-study"
MODEL_VERSION="0.7.2"
MODEL_NAME="DeepVariant-inception_v3-${MODEL_VERSION}+data-wgs_standard"
MODEL_HTTP_DIR="https://storage.googleapis.com/deepvariant/models/DeepVariant/${MODEL_VERSION}/${MODEL_NAME}"

INPUT_DIR="${BASE}/input"
MODELS_DIR="${INPUT_DIR}/models"
MODEL="${MODELS_DIR}/model.ckpt"
DATA_DIR="${INPUT_DIR}/data"
REF="${DATA_DIR}/hs37d5.fa.gz"
BAM="${DATA_DIR}/HG002_NIST_150bp_50x.bam"
TRUTH_VCF="${DATA_DIR}/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz"
TRUTH_BED="${DATA_DIR}/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"

OUTPUT_DIR="${BASE}/output"
EXAMPLES="${OUTPUT_DIR}/HG002.examples.tfrecord@${N_SHARDS}.gz"
GVCF_TFRECORDS="${OUTPUT_DIR}/HG002.gvcf.tfrecord@${N_SHARDS}.gz"
CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/HG002.cvo.tfrecord.gz"
OUTPUT_VCF="${OUTPUT_DIR}/HG002.output.vcf.gz"
OUTPUT_GVCF="${OUTPUT_DIR}/HG002.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

# export PYTHONPATH=${LOCATION}/runfiles:${LOCATION}/runfiles/protobuf_archive/python:${LOCATION}/runfiles/protobuf_archive:${LOCATION}/runfiles/absl_py:${LOCATION}/runfiles/six_archive:${LOCATION}/runfiles/com_google_deepvariant:${LOCATION}/runfiles/org_tensorflow_slim
# source gpu/bin/activate


time python ~/deepvariant/deepvariant/call_variants_2new_customed.py \
    --outfile "${CALL_VARIANTS_OUTPUT}" \
    --checkpoint "${MODEL}" \
    --process_number 20 &
#time python deepvariant/deepvariant/call_variants_multiip.py --outfile /home/yang/quickstart-output/call_variants_output.tfrecord1.gz --checkpoint /home/yang/DeepVariant-inception_v3-0.7.0+data-wgs_standard/model.ckpt --process_number 16 --gpu 0 --ip 10.10.10.52 --portoffset 0&
# time ~/deepvariant/bazel-bin/deepvariant/program \
#     --reads "${BAM}" \
#     --mode calling \
#     --ref "${REF}" \
#     --thread_number $THREAD \
#     --total_thread $TOTAL_THREAD \
#     --thread_offset $OFFSET
