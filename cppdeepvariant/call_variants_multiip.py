# Copyright 2017 Google Inc.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
"""Code for calling variants with a trained DeepVariant model."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import time



from absl import flags
from absl import logging
import numpy as np
import tensorflow as tf
import six
import zmq


from tensorflow.python.eager import context
from tensorflow.python.framework import ops
from tensorflow.python.framework import random_seed
from tensorflow.python.training import training
from tensorflow.python.estimator import model_fn as model_fn_lib


from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import errors
from third_party.nucleus.util import io_utils
from third_party.nucleus.util import proto_utils
from third_party.nucleus.util import variant_utils
from deepvariant import data_providers
from deepvariant import logging_level
from deepvariant import modeling
from deepvariant import tf_utils
from deepvariant.protos import deepvariant_pb2

_ALLOW_EXECUTION_HARDWARE = [
    'auto',  # Default, no validation.
    'cpu',  # Don't use accelerators, even if available.
    'accelerator',  # Must be hardware acceleration or an error will be raised.
]

# The number of digits past the decimal point that genotype likelihoods are
# rounded to, for numerical stability.
_GL_PRECISION = 10

# This number is estimated by the following logic:
# CPU run is roughly 0.2 sec per 100.
# 15000 examples will take about 30secs to print each line.
_LOG_EVERY_N = 15000

FLAGS = flags.FLAGS

flags.DEFINE_integer(
    'process_number', None,
    'launch parallel process_number processes')
flags.DEFINE_integer(
    'portoffset', 0,
    'portoffset')
flags.DEFINE_string(
    'gpu', None,
    'use which gpu')
flags.DEFINE_string(
    'ip', None,
    'listen on which ip')

flags.DEFINE_string(
    'outfile', None,
    'Required. Destination path where we will write output candidate variants '
    'with additional likelihood information in TFRecord format of '
    'CallVariantsOutput protos.')
flags.DEFINE_string(
    'checkpoint', None,
    'Required. Path to the TensorFlow model checkpoint to use to evaluate '
    'candidate variant calls.')
flags.DEFINE_integer(
    'batch_size', 512,
    'Number of candidate variant tensors to batch together during inference. '
    'Larger batches use more memory but are more computational efficient.')
flags.DEFINE_integer('max_batches', None,
                     'Max. batches to evaluate. Defaults to all.')
flags.DEFINE_integer('num_readers', 8,
                     'Number of parallel readers to create for examples.')
flags.DEFINE_string('model_name', 'inception_v3',
                    'The name of the model architecture of --checkpoint.')
flags.DEFINE_boolean('include_debug_info', False,
                     'If true, include extra debug info in the output.')
flags.DEFINE_boolean(
    'debugging_true_label_mode', False,
    'If true, read the true labels from examples and add to '
    'output. Note that the program will crash if the input '
    'examples do not have the label field. '
    'When true, this will also fill everything when '
    '--include_debug_info is set to true.')
flags.DEFINE_string(
    'execution_hardware', 'auto',
    'When in cpu mode, call_variants will not place any ops on the GPU, even '
    'if one is available. In accelerator mode call_variants validates that at '
    'least some hardware accelerator (GPU/TPU) was available for us. This '
    'option is primarily for QA purposes to allow users to validate their '
    'accelerator environment is correctly configured. In auto mode, the '
    'default, op placement is entirely left up to TensorFlow.  In tpu mode, '
    'use and require TPU.')

# Cloud TPU Cluster Resolvers
flags.DEFINE_string(
    'gcp_project', None,
    'Project name for the Cloud TPU-enabled project. If not specified, we '
    'will attempt to automatically detect the GCE project from metadata.')
flags.DEFINE_string(
    'tpu_zone', None,
    'GCE zone where the Cloud TPU is located in. If not specified, we '
    'will attempt to automatically detect the GCE project from metadata.')
flags.DEFINE_string(
    'tpu_name',
    None,
    'Name of the Cloud TPU for Cluster Resolvers. You must specify either '
    'this flag or --master. An empty value corresponds to no Cloud TPU. See '
    'https://www.tensorflow.org/api_docs/python/tf/contrib/cluster_resolver/TPUClusterResolver'  # pylint: disable=line-too-long
)

flags.DEFINE_string(
    'master', None,
    'GRPC URL of the master (e.g. grpc://ip.address.of.tpu:8470). You '
    'must specify either this flag or --tpu_name.')

flags.DEFINE_boolean('use_tpu', False, 'Use tpu if available.')

flags.DEFINE_string(
    'kmp_blocktime', '0',
    'Value to set the KMP_BLOCKTIME environment variable to for efficient MKL '
    'inference. See https://www.tensorflow.org/performance/performance_guide '
    'for more information. The default value is 0, which provides the best '
    'performance in our tests. Set this flag to "" to not set the variable.')




class ExecutionHardwareError(Exception):
  pass


def prepare_inputs(source_path, use_tpu=False, num_readers=None):
  """Return a tf.data input_fn from the source_path.

  Args:
    source_path: Path to a TFRecord file containing deepvariant tf.Example
      protos.
    use_tpu: boolean.  Use the tpu code path.
    num_readers: int > 0 or None. Number of parallel readers to use to read
      examples from source_path. If None, uses FLAGS.num_readers instead.

  Returns:
    A tf input_fn yielding batches of image, encoded_variant,
    encoded_alt_allele_indices.

    The image is a [batch_size, height, width, channel] tensor. The
    encoded_variants is a tf.string or tpu-encoded tensor containing a
    serialized Variant proto describing the variant call associated with
    image. The encoded_alt_allele_indices is a tf.string or tpu-encoded
    tensor containing a serialized CallVariantsOutput.AltAlleleIndices proto
    containing the alternate alleles indices used as "alt" when constructing
    the image.
  """
  if not num_readers:
    num_readers = FLAGS.num_readers

  return data_providers.get_input_fn_from_filespec(
      input_file_spec=source_path,
      mode=tf.estimator.ModeKeys.PREDICT,
      use_tpu=use_tpu,
      input_read_threads=num_readers,
      debugging_true_label_mode=FLAGS.debugging_true_label_mode,
  )


def round_gls(gls, precision=None):
  """Returns genotype likelihoods rounded to the desired precision level.

  Args:
    gls: A list of floats. The input genotype likelihoods at any precision.
    precision: Positive int. The number of places past the decimal point to
      round to. If None, no rounding is performed.

  Returns:
    A list of floats rounded to the desired precision.

  Raises:
    ValueError: The input gls do not sum to nearly 1.
  """
  if abs(sum(gls) - 1) > 1e-6:
    raise ValueError(
        'Invalid genotype likelihoods do not sum to one: sum({}) = {}'.format(
            gls, sum(gls)))
  if precision is None:
    return gls

  min_ix = 0
  min_gl = gls[0]
  for ix, gl in enumerate(gls):
    if gl < min_gl:
      min_gl = gl
      min_ix = ix

  rounded_gls = [round(gl, precision) for gl in gls]
  rounded_gls[min_ix] = max(
      0.0,
      round(1 - sum(rounded_gls[:min_ix] + rounded_gls[min_ix + 1:]),
            precision))
  return rounded_gls


def write_variant_call(writer, prediction, use_tpu):
  """Write the variant call based on prediction.

  Args:
    writer: A object with a write() function that will be called for each
      encoded_variant and genotype likelihoods.
    prediction: A [3] tensor of floats. These are the predicted
      genotype likelihoods (p00, p0x, pxx) for some alt allele x, in the same
      order as encoded_variants.
      use_tpu: bool.  Decode the tpu specific encoding of prediction.

  Returns:
    The return status from writer.
  """
  encoded_variant = prediction['variant']
  if use_tpu:
    encoded_variant = tf_utils.int_tensor_to_string(encoded_variant)

  encoded_alt_allele_indices = prediction['alt_allele_indices']
  if use_tpu:
    encoded_alt_allele_indices = tf_utils.int_tensor_to_string(
        encoded_alt_allele_indices)

  rounded_gls = round_gls(prediction['probabilities'], precision=_GL_PRECISION)

  # Write it out.
  true_labels = prediction['label'] if FLAGS.debugging_true_label_mode else None
  cvo = _create_cvo_proto(encoded_variant, rounded_gls,
                          encoded_alt_allele_indices, true_labels)
  return writer.write(cvo.SerializeToString())


def _create_cvo_proto(encoded_variant,
                      gls,
                      encoded_alt_allele_indices,
                      true_labels=None):
  """Returns a CallVariantsOutput proto from the relevant input information."""
  variant = variants_pb2.Variant.FromString(encoded_variant)
  alt_allele_indices = (
      deepvariant_pb2.CallVariantsOutput.AltAlleleIndices.FromString(
          encoded_alt_allele_indices))
  debug_info = None
  if FLAGS.include_debug_info or FLAGS.debugging_true_label_mode:
    debug_info = deepvariant_pb2.CallVariantsOutput.DebugInfo(
        has_insertion=variant_utils.has_insertion(variant),
        has_deletion=variant_utils.has_deletion(variant),
        is_snp=variant_utils.is_snp(variant),
        predicted_label=np.argmax(gls),
        true_label=true_labels,
    )
  call_variants_output = deepvariant_pb2.CallVariantsOutput(
      variant=variant,
      alt_allele_indices=alt_allele_indices,
      genotype_probabilities=gls,
      debug_info=debug_info)
  return call_variants_output


def call_variants(
                  checkpoint_path,
                  model,
                  output_file,
                  execution_hardware='auto',
                  batch_size=16,
                  max_batches=None,
                  use_tpu=False,
                  master=''):
  """Main driver of call_variants."""
  os.environ["CUDA_VISIBLE_DEVICES"] = FLAGS.gpu
  if FLAGS.kmp_blocktime:
    os.environ['KMP_BLOCKTIME'] = FLAGS.kmp_blocktime
    logging.info('Set KMP_BLOCKTIME to %s', os.environ['KMP_BLOCKTIME'])

  # Read a single TFExample to make sure we're not loading an older version.
 

  # Check accelerator status.
  if execution_hardware not in _ALLOW_EXECUTION_HARDWARE:
    raise ValueError(
        'Unexpected execution_hardware={} value. Allowed values are {}'.format(
            execution_hardware, ','.join(_ALLOW_EXECUTION_HARDWARE)))
  init_op = tf.group(tf.global_variables_initializer(),
                     tf.local_variables_initializer())
  device_count = {'GPU': 0, 'TPU': 0} if execution_hardware == 'cpu' else {}
  config = tf.ConfigProto(device_count=device_count)
  with tf.Session(config=config) as sess:
    sess.run(init_op)
    if execution_hardware == 'accelerator':
      if not any(dev.device_type != 'CPU' for dev in sess.list_devices()):
        raise ExecutionHardwareError(
            'execution_hardware is set to accelerator, but no accelerator '
            'was found')
    # redacted
    # sess.list_devices here doesn't return the correct answer. That can only
    # work later, after the device (on the other VM) has been initialized,
    # which is generally not yet.

  


  # Instantiate the prediction "stream", and select the EMA values from
  # the model.
  if checkpoint_path is None:
    # Unit tests use this branch.
    predict_hooks = []
  else:
    predict_hooks = [h(checkpoint_path) for h in model.session_predict_hooks()]
  

  def _create_global_step(graph):
    """Creates the global step tensor in graph.
    The global step tensor must be an integer type with name 'global_step' and
    be added to the collection `tf.GraphKeys.GLOBAL_STEP`.
    Args:
      graph: The graph in which to create the global step tensor.
    Returns:
      The global step `tf.Tensor`.
    """
    return training.create_global_step(graph)

  def _create_and_assert_global_step(graph):
    """Creates and asserts properties of the global step.
    Args:
      graph: The graph in which to create the global step tensor.
    Returns:
      The global step `tf.Tensor`.
    """
    step = _create_global_step(graph)
    assert step == training.get_global_step()
    assert step.dtype.is_integer
    return step
  with context.graph_mode():
      hooks = predict_hooks
      # Check that model has been trained.
      if not checkpoint_path:
        print ('required checkpoint_path')
        return
      with ops.Graph().as_default() as g:
        random_seed.set_random_seed(1234) #uncertain
        _create_and_assert_global_step(g)
      
        a = {}
        a['image'] = tf.placeholder(tf.uint8, shape=(None,100,221,6))
    
        a['variant'] = tf.placeholder(tf.string, shape = (None,))
        a['alt_allele_indices'] = tf.placeholder(tf.string, shape = (None,))
        estimator_spec = model.model_fn(
            a, None, model_fn_lib.ModeKeys.PREDICT, None)
        predictions = estimator_spec.predictions

        # Call to warm_start has to be after model_fn is called.
        zmqcontext = zmq.Context()
        consumer_receiver=[]
        address="tcp://"+FLAGS.ip+":"
        port =5558
        for i in range(FLAGS.process_number):
          consumer_receiver.append(zmqcontext.socket(zmq.PULL))
          consumer_receiver[i].bind(address+str(port+i+FLAGS.portoffset))
        address="tcp://10.10.9.51:"
        for i in range(8):
          consumer_receiver.append(zmqcontext.socket(zmq.PULL))
          consumer_receiver[-1].bind(address+str(port+i+16))

            

        all_hooks = hooks
        logging.info('Writing calls to %s', output_file)
        writer, _ = io_utils.make_proto_writer(output_file)
        with writer:
          with training.MonitoredSession(
              session_creator=training.ChiefSessionCreator(
                  checkpoint_filename_with_path=checkpoint_path),
              hooks=all_hooks) as mon_sess:
            count = 0
            index = 0
            while not mon_sess.should_stop():
              x=[]
              y=[]
              z=[]
              start = time.time()
              for i in range(512):
                work1 = consumer_receiver[index].recv()
                
                if work1 != "Done":
                  work1image=np.fromstring(work1, dtype=np.uint8).reshape((100, 6,221),order='C').swapaxes(1,2)
                  work2 = consumer_receiver[index].recv()
                  work3 = consumer_receiver[index].recv()
                  #workimage=np.ones((100,6,221), dtype=np.uint8)
                  x.append(work1image)
                  y.append(work2)
                  z.append(work3)
                  index+=1
                  index%=len(consumer_receiver)
                else:
                  count+=1

                  if count == FLAGS.process_number+8:
                    preds_evaluated = mon_sess.run(predictions, feed_dict = {a['image'] : x,
                                                                       a['variant'] : y,
                                                                       a['alt_allele_indices'] : z})
                    for i in range(len(x)):
                      p= {key: value[i]
                      for key, value in six.iteritems(preds_evaluated)}
                      write_variant_call(writer, p, use_tpu)
                    return 
                  del consumer_receiver[index]
                  index%=len(consumer_receiver)
              print ('recv {}'.format(time.time()-start))
              start=time.time()


              preds_evaluated = mon_sess.run(predictions, feed_dict = {a['image'] : x,
                                                                       a['variant'] : y,
                                                                       a['alt_allele_indices'] : z})
              print ('inference {}'.format(time.time()-start))

              for i in range(len(x)):
                p= {key: value[i]
                    for key, value in six.iteritems(preds_evaluated)}
                write_variant_call(writer, p, use_tpu)
              print ("Done a batch")
  logging.info('Done evaluating variants')
    
def main(argv=()):
  with errors.clean_commandline_error_exit():
    
    proto_utils.uses_fast_cpp_protos_or_die()

    logging_level.set_from_flag()

    if FLAGS.use_tpu:
      master = tf_utils.resolve_master(FLAGS.master, FLAGS.tpu_name,
                                       FLAGS.tpu_zone, FLAGS.gcp_project)
    else:
      master = ''

    model = modeling.get_model(FLAGS.model_name)
    call_variants(
        checkpoint_path=FLAGS.checkpoint,
        model=model,
        execution_hardware=FLAGS.execution_hardware,
        output_file=FLAGS.outfile,
        max_batches=FLAGS.max_batches,
        batch_size=FLAGS.batch_size,
        master=master,
        use_tpu=FLAGS.use_tpu,
    )


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'process_number',
      'outfile',
      'checkpoint',
  ])
  
  tf.app.run()


