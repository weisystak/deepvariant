
#include <string>
#include <memory>

#include <tensorflow/core/framework/op.h>
#include <tensorflow/core/framework/op_kernel.h>
#include <tensorflow/core/framework/resource_op_kernel.h>
#include <tensorflow/core/framework/resource_mgr.h>
#include <tensorflow/core/framework/common_shape_fns.h>
#include <iostream>
#include "zmq_conn.h"
#include <stdlib.h>
using namespace std;
using namespace tensorflow;



// An op to create zmq connection as a resource.
// Use ResourceOpKernel to ensure singleton construction.
class ZMQConnectionHandleOp : public ResourceOpKernel<ZMQConnection> {
  public:
    explicit ZMQConnectionHandleOp(OpKernelConstruction* ctx)
        : ResourceOpKernel<ZMQConnection>(ctx) {}

  private:
    Status CreateResource(ZMQConnection** ret) override EXCLUSIVE_LOCKS_REQUIRED(mu_) {
      const NodeDef& ndef = def();
      ZMQSocketDef sockdef;
      TF_RETURN_IF_ERROR(GetNodeAttr(ndef, "num_port", &sockdef.num_port));
      TF_RETURN_IF_ERROR(GetNodeAttr(ndef, "end_point", &sockdef.end_point));
      TF_RETURN_IF_ERROR(GetNodeAttr(ndef, "portoffset", &sockdef.portoffset));

      *ret = new ZMQConnection(sockdef);
      cout << "Handle OK" << endl;
      return Status::OK();
    }
};


class ZMQPullOp: public AsyncOpKernel {
 public:
  explicit ZMQPullOp(OpKernelConstruction* context) : AsyncOpKernel(context) {
  }

  void ComputeAsync(OpKernelContext* ctx, DoneCallback done) override {

  
    ZMQConnection* conn = nullptr;
    OP_REQUIRES_OK_ASYNC(
        ctx, LookupResource(ctx, HandleFromInput(ctx, 0), &conn), done);

    Tensor* output_0 = nullptr;
    Tensor* output_1 = nullptr;
    Tensor* output_2 = nullptr;
    Tensor* output_3 = nullptr;
    Tensor* output_4 = nullptr;

    OP_REQUIRES_OK_ASYNC(ctx, ctx->allocate_output(0, {512,100,221,6}, &output_0), done);
    OP_REQUIRES_OK_ASYNC(ctx, ctx->allocate_output(1, {512}, &output_1), done);
    OP_REQUIRES_OK_ASYNC(ctx, ctx->allocate_output(2, {512}, &output_2), done);
    OP_REQUIRES_OK_ASYNC(ctx, ctx->allocate_output(3, {1}, &output_3), done); 
    OP_REQUIRES_OK_ASYNC(ctx, ctx->allocate_output(4, {1}, &output_4), done); 

    
    auto o0= output_0->flat<uint8>().data();
    auto o1= output_1->flat<string>();
    auto o2= output_2->flat<string>();
    auto o3= output_3->flat<int32>();
    auto o4= output_4->flat<int32>();
 
    o3(0)=512;
    o4(0)=0;
    static int index=0;
    static int count=0;
    
    int oindex=0;
    for (int i=0; i<512; i++){
      zmq::message_t reply;
       
      conn->vector_sock[index].recv (&reply);
      if (reply.size() !=4){
        memcpy(o0, reply.data(), reply.size() );
        o0 += reply.size();
        zmq::message_t reply1;
        conn->vector_sock[index].recv (&reply1);
        o1(oindex)=std::string(static_cast<char*>(reply1.data()), reply1.size());
        zmq::message_t reply2;
        conn->vector_sock[index].recv (&reply2);
        o2(oindex)=std::string(static_cast<char*>(reply2.data()), reply2.size());
        oindex+=1;
        index+=1;
        index%=conn->vector_sock.size();

      }
      else{
        count+=1;
        if (count== conn->def_.num_port){
          o4(0)=1;
          
          break;
        }
        else{
          cout << "erase" <<endl;
          conn->vector_sock.erase (conn->vector_sock.begin()+index);
          index%=conn->vector_sock.size();
        }

      }
    
    }

     //std::string rpl = std::string(static_cast<char*>(tlist.message.data()), tlist.message.size());
     //cout << rpl << endl;
     //az(0) = rpl;
      //OP_REQUIRES_OK_ASYNC(ctx, ctx->allocate_output(1, {2,3}, &output), done);
     //auto bz= output->flat<int32>();
     //for (int i=0; i<6 ; i++)
       //bz(i) =i;
     //memcpy(output->flat<string>().data(), tlist.message.data(), tlist.message.size() );
    o3(0)=oindex;
    done();
  }
 
};


REGISTER_KERNEL_BUILDER(Name("ZMQPull").Device(DEVICE_CPU), ZMQPullOp);
REGISTER_KERNEL_BUILDER(Name("ZMQConnection").Device(DEVICE_CPU), ZMQConnectionHandleOp);


REGISTER_OP("ZMQPull")
    .Input("handle: resource")
    .Output("output_0: uint8")
    .Output("output_1: string")
    .Output("output_2: string")
    .Output("output_3: int32")
    .Output("output_4: int32")
    .SetShapeFn(shape_inference::UnknownShape)

    .SetIsStateful();
   


REGISTER_OP("ZMQConnection")
    .Output("handle: resource")
    .Attr("end_point: string")
    .Attr("num_port: int")
    .Attr("portoffset: int")
    .Attr("container: string = ''")
    .Attr("shared_name: string = ''")
    .SetIsStateful()
    .SetShapeFn(shape_inference::ScalarShape);
    