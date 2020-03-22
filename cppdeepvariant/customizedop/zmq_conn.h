//File: zmq_conn.h
//Author: Yuxin Wu <ppwwyyxxc@gmail.com>

#pragma once

#include <string>
#include <iostream>
#include <vector>


#include <tensorflow/core/framework/tensor_shape.h>
#include <tensorflow/core/lib/gtl/inlined_vector.h>
#include <tensorflow/core/lib/strings/strcat.h>
#include <tensorflow/core/platform/mutex.h>
#include "zmq.hpp"


#ifdef __GNUC__
#ifndef __clang__
#if ((__GNUC__ <= 5) && (__GNUC_MINOR__ <= 3))
#error "GCC >= 5.3 is required!"
#endif
#endif  // clang
#endif  // gnuc

struct ZMQSocketDef {
  std::string end_point;
  int num_port;
  int portoffset;
};


class ZMQConnection : public tensorflow::ResourceBase {
  public:
  explicit ZMQConnection(const ZMQSocketDef& def):
    def_{def}, ctx_{1} {
      std::string address="tcp://"+def_.end_point+":";
      for (int i=0; i<def_.num_port; i++){
        vector_sock.emplace_back(ctx_, ZMQ_PULL);
        vector_sock[i].bind(address+std::to_string(5558+i+def_.portoffset));
      }     
  }
  std::string DebugString() { return "123"; }
  std::string DebugString() const { return "456"; }
  public:
    ZMQSocketDef def_;
    zmq::context_t ctx_;
    std::vector<zmq::socket_t> vector_sock;
};



