module OnelabReplyContractChecks
using Test,Tessella,Sockets
const A=Tessella.API
function malformed(operation,mode)
 listener=Sockets.listen(Sockets.localhost,0);port=getsockname(listener)[2]
 send(peer,kind,payload="")=(write(peer,htol(Int32(kind)));write(peer,htol(Int32(ncodeunits(payload))));write(peer,payload);flush(peer))
 server=@async begin
  peer=accept(listener)
  try
   while isopen(peer)
    kind=Int(ltoh(read(peer,Int32)));count=Int(ltoh(read(peer,Int32)));payload=String(read(peer,count))
    kind==2 && break
    if kind==24
     wire=A._OnelabWireParameter(:number,"different");wire.values_f=[9.]
     send(peer,23,A._onelab_wire_encode(wire))
    elseif kind==25
     if mode==:end
      send(peer,2,"unexpected STOP")
     else
      if split(payload,'\0')[2]=="number"
       wire=A._OnelabWireParameter(:string,"wrong kind");wire.values_s=["x"]
       send(peer,25,A._onelab_wire_encode(wire))
      end
      send(peer,26)
     end
    end
   end
  catch error
   # A client that rejects the first bad list record can close before the
   # fixture finishes writing its trailing terminator.
   error isa EOFError || error isa Base.IOError || rethrow()
  finally
   close(peer)
  end
 end
 c=A.onelab.client_connect("127.0.0.1",port)
 try
  operation(c)
 finally
  A.onelab.client_disconnect(c);wait(server);close(listener)
 end
end
@testset "ONELAB replies identify the requested record and list completion" begin
 @test_throws ArgumentError malformed(c->A.onelab.client_get_number(c,"wanted"),:name)
 @test_throws ArgumentError malformed(c->A.onelab.client_define_number(c,"wanted";value=2.),:name)
 @test_throws ArgumentError malformed(c->A.onelab.client_get_names(c),:end)
 @test_throws ArgumentError malformed(c->A.onelab.client_get_names(c),:kind)
end


end
