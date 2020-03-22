#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File: zmq_ops.py

import sys
import tensorflow as tf
import struct
import numpy as np
import os
import addressbook_pb2

address_book = addressbook_pb2.AddressBook()



_zmq_ops = tf.load_op_library('./zmq_pull_op.so')


class ZMQPullSocket(object):
    def __init__(self, end_point, num_port):

       

        self._zmq_handle = _zmq_ops.zmq_connection(
            end_point, num_port)


    def pull(self):
        return _zmq_ops.zmq_pull(
            self._zmq_handle)

sock = ZMQPullSocket('127.0.0.1',1)

sess = tf.Session()
count=1;
while(1):
  b=sock.pull()
  a = tf.reshape(b[0],[512,100,221,6])
  print a
  c=b[3]*2
  a={'ff': c, 'gg': b[0]}  
  #print sess.run(a)
  #print count
  count+=1
