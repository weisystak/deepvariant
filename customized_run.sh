THREAD=20
TOTAL_THREAD=20
OFFSET=0
LOCATION=/home/yang
export PYTHONPATH=${LOCATION}/runfiles:${LOCATION}/runfiles/protobuf_archive/python:${LOCATION}/runfiles/protobuf_archive:${LOCATION}/runfiles/absl_py:${LOCATION}/runfiles/six_archive:${LOCATION}/runfiles/com_google_deepvariant:${LOCATION}/runfiles/org_tensorflow_slim
source gpu/bin/activate
time python deepvariant/deepvariant/call_variants_2new_customed.py --outfile /home/yang/quickstart-output/call_variants_output.tfrecord0.gz --checkpoint /home/yang/DeepVariant-inception_v3-0.7.0+data-wgs_standard/model.ckpt --process_number 20 &
#time python deepvariant/deepvariant/call_variants_multiip.py --outfile /home/yang/quickstart-output/call_variants_output.tfrecord1.gz --checkpoint /home/yang/DeepVariant-inception_v3-0.7.0+data-wgs_standard/model.ckpt --process_number 16 --gpu 0 --ip 10.10.10.52 --portoffset 0&
time deepvariant/bazel-bin/deepvariant/program --reads /home/paslab/workspace/d97922039/deepvariant_testcases/whole_genome/case-study/input/data/HG002_NIST_150bp_50x.bam --mode calling --ref /home/paslab/workspace/d97922039/deepvariant_testcases/whole_genome/case-study/input/data/hs37d5.fa.gz --thread_number $THREAD --total_thread $TOTAL_THREAD --thread_offset $OFFSET
