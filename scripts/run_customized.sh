THREAD=20
TOTAL_THREAD=20
OFFSET=0
# LOCATION=/home/yang

BASE="${HOME}/case-study"
BIN_VERSION="0.7.2"
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

N_SHARDS="64"

OUTPUT_DIR="${BASE}/output"
EXAMPLES="${OUTPUT_DIR}/HG002.examples.tfrecord@${N_SHARDS}.gz"
GVCF_TFRECORDS="${OUTPUT_DIR}/HG002.gvcf.tfrecord@${N_SHARDS}.gz"
CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/HG002.cvo.tfrecord.gz"
OUTPUT_VCF="${OUTPUT_DIR}/HG002.output.vcf.gz"
OUTPUT_GVCF="${OUTPUT_DIR}/HG002.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

# # Create local directory structure
# mkdir -p "${OUTPUT_DIR}"
# mkdir -p "${DATA_DIR}"
# mkdir -p "${MODELS_DIR}"
# mkdir -p "${LOG_DIR}"



# # Download models, and test data
# Copy the model files to your local disk.
# aria2c -c -x10 -s10 -d "${MODELS_DIR}" "${MODEL_HTTP_DIR}"/model.ckpt.data-00000-of-00001
# aria2c -c -x10 -s10 -d "${MODELS_DIR}" "${MODEL_HTTP_DIR}"/model.ckpt.index
# aria2c -c -x10 -s10 -d "${MODELS_DIR}" "${MODEL_HTTP_DIR}"/model.ckpt.meta

# Copy the data
# aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed -d "${DATA_DIR}"
# aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz -d "${DATA_DIR}"
# aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi -d "${DATA_DIR}"
# aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_NIST_150bp_50x.bam -d "${DATA_DIR}"
# aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/HG002_NIST_150bp_50x.bam.bai -d "${DATA_DIR}"
# aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz -d "${DATA_DIR}"
# aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz.fai -d "${DATA_DIR}"
# aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gz.gzi -d "${DATA_DIR}"
# aria2c -c -x10 -s10 https://storage.googleapis.com/deepvariant/case-study-testdata/hs37d5.fa.gzi -d "${DATA_DIR}"
# export PYTHONPATH=${LOCATION}/runfiles:${LOCATION}/runfiles/protobuf_archive/python:${LOCATION}/runfiles/protobuf_archive:${LOCATION}/runfiles/absl_py:${LOCATION}/runfiles/six_archive:${LOCATION}/runfiles/com_google_deepvariant:${LOCATION}/runfiles/org_tensorflow_slim
# source gpu/bin/activate

echo "${CALL_VARIANTS_OUTPUT}"
echo "${MODEL}"
time  python ~/deepvariant/deepvariant/call_variants_2new_customed.py \
    --outfile "${CALL_VARIANTS_OUTPUT}" \
    --checkpoint "${MODEL}" \
    --process_number 20 \
    --gpu 0 \
    --ip 127.0.0.1  &
    #--regions "chr20:10,000,000-10,010,000"
#time python deepvariant/deepvariant/call_variants_multiip.py --outfile /home/yang/quickstart-output/call_variants_output.tfrecord1.gz --checkpoint /home/yang/DeepVariant-inception_v3-0.7.0+data-wgs_standard/model.ckpt --process_number 16 --gpu 0 --ip 10.10.10.52 --portoffset 0&
time ~/deepvariant/bazel-bin/deepvariant/program \
    --reads "${BAM}" \
    --mode calling \
    --ref "${REF}" \
    --thread_number $THREAD \
    --total_thread $TOTAL_THREAD \
    --thread_offset $OFFSET


## Run `postprocess_variants`, with gVCFs.
# echo "Start running postprocess_variants (with gVCFs)...Log will be in the terminal and also to ${LOG_DIR}/postprocess_variants.withGVCF.log."
# ( time sudo docker run \
#     -v "${BASE}":"${BASE}" \
#     gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
#     /opt/deepvariant/bin/postprocess_variants \
#     --ref "${REF}" \
#     --infile "${CALL_VARIANTS_OUTPUT}" \
#     --outfile "${OUTPUT_VCF}" \
#     --nonvariant_site_tfrecord_path "${GVCF_TFRECORDS}" \
#     --gvcf_outfile "${OUTPUT_GVCF}"
# ) 2>&1 | tee "${LOG_DIR}/postprocess_variants.withGVCF.log"
# echo "Done."
# echo


## Evaluation: run hap.py
# echo "Start evaluation with hap.py..."
# UNCOMPRESSED_REF="${OUTPUT_DIR}/hs37d5.fa"

# hap.py cannot read the compressed fa, so uncompress
# into a writable directory and index it.
# zcat <"${REF}" >"${UNCOMPRESSED_REF}"
# samtools faidx "${UNCOMPRESSED_REF}"
# 
# sudo docker pull pkrusche/hap.py
# ( sudo docker run -i \
# -v "${DATA_DIR}:${DATA_DIR}" \
# -v "${OUTPUT_DIR}:${OUTPUT_DIR}" \
# pkrusche/hap.py /opt/hap.py/bin/hap.py \
#   "${TRUTH_VCF}" \
#   "${OUTPUT_VCF}" \
#   -f "${TRUTH_BED}" \
#   -r "${UNCOMPRESSED_REF}" \
#   -o "${OUTPUT_DIR}/happy.output" \
#   --engine=vcfeval
# ) 2>&1 | tee "${LOG_DIR}/happy.log"
# echo "Done."
