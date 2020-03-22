import tensorflow as tf
zero_out_module = tf.load_op_library('./zero_out.so')
with tf.Session('') as sess:
  x = tf.Variable([[[1, 2], [3, 7]],[[1, 2], [3, 7]]], name='x', dtype=tf.int32)
  sess.run(tf.global_variables_initializer())
  a = zero_out_module.zero_out(x)
  a*=2
  for i in range(10):
    print sess.run(a)

