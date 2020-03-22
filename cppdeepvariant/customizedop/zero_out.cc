#include "tensorflow/core/framework/op.h"
#include "tensorflow/core/framework/shape_inference.h"
#include "tensorflow/core/framework/op_kernel.h"
#include <cuda_runtime.h>
#include <iostream>
#include "zmq.hpp"

using namespace tensorflow;
int v=0;
class ZeroOutOp : public OpKernel {
 public:
  explicit ZeroOutOp(OpKernelConstruction* context) : OpKernel(context) {}

  void Compute(OpKernelContext* context) override {
    // Grab the input tensor
    // Create an output tensor
    Tensor* output_tensor = NULL;
    OP_REQUIRES_OK(context, context->allocate_output(0, {8,},
                                                     &output_tensor));
    auto output_flat = output_tensor->flat<int32>();

    // Set all but the first element of the output tensor to 0.
    const int N = 8;
    static int v=8;
    int a[10]={1,2,3,4,5,6,7,v};
    cudaMemcpy(output_flat.data()+1, a,4*(N-1), cudaMemcpyHostToDevice);
    v+=1;
    //std::cout << a[0] << std::endl;
    // Preserve the first input value if possible.
  }
};

REGISTER_OP("ZeroOut")
    .Input("to_zero: int32")
    .Output("zeroed: int32")
    .SetShapeFn([](::tensorflow::shape_inference::InferenceContext* c) {
      c->set_output(0, c->input(0));
      return Status::OK();
    });
REGISTER_KERNEL_BUILDER(Name("ZeroOut").Device(DEVICE_GPU), ZeroOutOp);

