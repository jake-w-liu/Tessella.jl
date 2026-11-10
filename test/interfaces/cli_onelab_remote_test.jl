module CliOnelabRemoteChecks
using Test,Tessella,Sockets,Random
const A=Tessella.API
const C=Tessella.CLI
function send(peer,kind,payload="")
 write(peer,htol(Int32(kind)));write(peer,htol(Int32(ncodeunits(payload))));write(peer,payload);flush(peer)
end
function recv(peer)
 kind=Int(ltoh(read(peer,Int32)));count=Int(ltoh(read(peer,Int32)))
 0<=count<=64*1024*1024 || error("fixture byte count")
 return kind,String(read(peer,count))
end
function cli_case(transport,service;queue_early=false,bad_open=false)
 mktempdir() do directory
  input=joinpath(directory,"initial.geo");output=joinpath(directory,"rpc.msh")
  write(input,"n=DefineNumber(3,Name \"Input/N\",Min 1,Max 10,Label \"Value\");"*
   "s=DefineString(\"default\",Name \"Input/S\");"*
   "SetNumber(\"Echo/N\",GetNumber(\"Input/N\"));Point(1)={0,0,0,1};\n"*
   "Macro later\nPoint(newp)={n,1,0,1}; SetNumber(\"Macro/result\",n);\nReturn\n"*
   "Printf(StrCat(\"ACTION=\",OnelabAction));\n")
  socket_path=Sys.iswindows() ? "\\\\.\\pipe\\tessella-"*randstring(12) : joinpath(directory,"ipc.sock")
  listener=transport==:tcp ? Sockets.listen(Sockets.localhost,0) : Sockets.listen(socket_path)
  address=transport==:tcp ? "127.0.0.1:"*string(Sockets.getsockname(listener)[2]) : socket_path
  name=service ? "GmshRemote" : "worker"
  parameters=Dict{Tuple{Symbol,String},Any}();messages=Tuple{Int,String}[]
  function number(name,values)
   p=A._OnelabWireParameter(:number,name);p.values_f=values;p.clients["observer"]=0;parameters[(:number,name)]=p
  end
  function text(name,values)
   p=A._OnelabWireParameter(:string,name);p.values_s=values;parameters[(:string,name)]=p
  end
  number("Input/N",[7.,8.]);text("Input/S",["λ","tail"]);text(name*"/Action",["check"])
  rpc="n+=1; Call later; Point(3)={1,0,0,1}; Line(1)={1,3}; Transfinite Curve{1}=3; Mesh 1; "*
   "Save \""*replace(output,'\\'=>'/')*"\"; Printf(\"RPC_DONE\"); Exit;"
  peer_ref=Ref{Any}(nothing)
  server=@async begin
   peer=nothing;command_sent=false;bad_sent=false
   try
    peer=Sockets.accept(listener);peer_ref[]=peer
    while true
     kind,payload=recv(peer);push!(messages,(kind,payload))
     if kind==2;break
     elseif kind==24
      q=A._onelab_wire_decode(payload,"fixture");p=get(parameters,(q.kind,q.name),nothing)
      p===nothing ? send(peer,29,q.name) : send(peer,23,A._onelab_wire_encode(p)*"\0")
     elseif kind==25
      if service && queue_early && !command_sent
       command_sent=true;send(peer,21,rpc)
      end
      q=A._onelab_wire_decode(payload,"fixture")
      for p in values(parameters);p.kind===q.kind && send(peer,25,A._onelab_wire_encode(p)*"\0");end
      send(peer,26)
     elseif kind in (23,32)
      p=A._onelab_wire_decode(payload,"fixture");parameters[(p.kind,p.name)]=p
      p.name=="Input/N" && p.label=="Value" && (p.values_f=[11.,12.])
      if service && !command_sent && p.name=="Gmsh/Model dimension"
       if bad_open && !bad_sent
        bad_sent=true;send(peer,33,joinpath(directory,"missing.geo"))
       elseif !bad_open
        command_sent=true;send(peer,21,rpc)
       end
      end
     elseif kind==12 && bad_open && bad_sent && !command_sent
      command_sent=true;send(peer,21,rpc)
     end
    end
   catch error
    error isa EOFError || error isa Base.IOError || rethrow()
   finally
    peer===nothing || close(peer)
   end
  end
  try
   args=service ? [input,"-setnumber","unused","2","-socket",address] :
    [input,"-onelab",name,address,"-parse_and_exit"]
   @test C.main(args)==input
  finally
   peer_ref[]===nothing || close(peer_ref[])
   close(listener);wait(server)
  end
  @test !isempty(messages) && first(messages)[1]==1 && last(messages)[1]==2
  @test parameters[(:number,"Echo/N")].values_f==[11.]
  @test parameters[(:string,"Input/S")].values_s==["λ","tail"]
  @test parameters[(:string,"Gmsh/Model name")].values_s==[input]
  @test parameters[(:string,"Gmsh/Model absolute path")].values_s==[dirname(abspath(input))*"/"]
  @test any(message->message==(10,"ACTION=check"),messages)
  @test isempty(A._ONELAB_ACTIVE_REMOTE) && isempty(A._ONELAB_ACTIVE_ENGINES)
  if service
   bad_open && @test any(message->message[1]==12 && occursin("missing.geo",message[2]),messages)
   @test parameters[(:number,"Macro/result")].values_f==[8.]
   @test any(message->message==(10,"RPC_DONE"),messages)
   mesh=Tessella.Elements.read_mixed_msh(output)
   @test Tessella.validate(mesh).ok
   @test any(row->mesh.coords[:,row]==[8.,1.,0.],axes(mesh.coords,2))
   @test sum(size(block.nodes,2) for block in mesh.blocks if block.msh==1)==2
  end
 end
end
@testset "CLI native remote batch and socket RPC" begin
 selected,ordinary=C._onelab_cli_flags(["-setstring","literal","-socket","-socket","address"])
 @test selected.name=="GmshRemote" && selected.address=="address"
 @test ordinary==["-setstring","literal","-socket"]
 selected,ordinary=C._onelab_cli_flags(["-setstring","literal","-2","-socket","address"])
 @test !selected.batch && ordinary==["-setstring","literal","-2"]
 @test_throws ArgumentError C._onelab_cli_flags(["-socket","x","-onelab","y"])
 @test_throws ArgumentError C._onelab_cli_flags(["-socket"])
 cli_case(:tcp,false)
 cli_case(:pipe,true)
 cli_case(:tcp,true;queue_early=true)
 cli_case(:tcp,true;bad_open=true)
end

function idle_case(reset)
 listener=Sockets.listen(Sockets.localhost,0);peer_ref=Ref{Any}(nothing)
 address="127.0.0.1:"*string(Sockets.getsockname(listener)[2])
 ready=Channel{Nothing}(1)
 server=@async try
  peer=Sockets.accept(listener);peer_ref[]=peer
  while true
   kind,payload=recv(peer)
   if kind==24
    p=A._onelab_wire_decode(payload,"fixture");send(peer,29,p.name)
   elseif kind==25
    send(peer,26)
   elseif kind==23
    p=A._onelab_wire_decode(payload,"fixture")
    p.name=="Gmsh/Number of physical groups" && put!(ready,nothing)
   elseif kind==2
    break
   end
  end
 catch error
  error isa EOFError || error isa Base.IOError || rethrow()
 end
 job=@async try C.main(["-socket",address]) catch e;e end
 @test timedwait(()->isready(ready),30.;pollint=.01)==:ok
 emergency=Timer(10.) do timer
  peer_ref[]===nothing || close(peer_ref[])
 end
 try
  if reset
   A.initialize();A.onelab.set_number("new session",[42.])
  else
   close(peer_ref[])
  end
  @test timedwait(()->istaskdone(job),5.;pollint=.01)==:ok
  result=fetch(job)
  @test reset ? result isa ArgumentError : result==""
  if reset
   @test A.CURRENT[]!==nothing
   @test A.onelab.get_number("new session")==[42.]
  end
  @test isempty(A._ONELAB_ACTIVE_REMOTE)
 finally
  close(emergency);peer_ref[]===nothing || close(peer_ref[])
  close(listener);wait(server);A.finalize()
 end
end
@testset "Remote CLI EOF and reset preserve the replacement API session" begin
 idle_case(false);idle_case(true)
end
@testset "GEO file namespace publishes only with its model" begin
 A.initialize();A.option("General.Terminal",0)
 try
  mktempdir() do directory
   path=joinpath(directory,"namespace.geo")
   write(path,"n=7; Point(1)={0,0,0,1};\nMacro later\nPoint(newp)={n,0,0,1};\nReturn\n")
   @test A.open(path)===nothing
   @test A.parser.get_number("n")==[7.]
   before=A.model.get_entities()
   write(path,"n=99; Point(99)={99,0,0,1}; Exit 7;\n")
   @test_throws ArgumentError A.open(path)
   @test A.parser.get_number("n")==[7.] && A.model.get_entities()==before
   write(path,"n+=1; Call later;\n")
   @test A.parser.parse(path)===nothing
   @test A.model.get_entities()==[(0,1),(0,2)]
   @test A.model.get_value(0,2,Float64[])==[8.,0.,0.]
  end
 finally
  A.finalize()
 end
end

@testset "Remote CLI delegates genuine solver processes to its controller" begin
 A.initialize();A.option("General.Terminal",0)
 quote_arg=Sys.iswindows() ? A._onelab_windows_quote : Base.shell_escape_posixly
 quote_geo_arg(value)=Sys.iswindows() ? A._onelab_windows_quote(String(value)) :
  "\""*replace(String(value),'\\'=>"\\\\",'"'=>"\\\"")*"\""
 try
  mktempdir() do directory
   result=joinpath(directory,"solver.txt")
   solver=joinpath(@__DIR__,"..","fixtures","onelab_native_client.jl")
   solver_command=join(quote_geo_arg.([String(Base.julia_cmd().exec[1]),
    "--startup-file=no","--compile=min","-O0",solver,"child",result])," ")
   source=joinpath(directory,"solver.geo")
   write(source,"OnelabRun(\"native worker\", '"*solver_command*"');\n"*
    "OnelabRun(\"native worker\");\n"*
    "SetNumber(\"Remote/result\",GetNumber(\"Output/child\"));\n")
   client=joinpath(@__DIR__,"..","fixtures","onelab_remote_cli_client.jl")
   command=join(quote_arg.(String[Base.julia_cmd().exec[1],"--startup-file=no",
    "--compiled-modules=no","--compile=min","-O0","--project="*dirname(Base.active_project()),
    client,source,"-parse_and_exit"])," ")
   @test A.onelab.run("remote mesher",command;timeout=180.,connect_timeout=120.)===nothing
   @test readlines(result)==["native worker","native worker"]
   @test A.onelab.get_number("Output/child")==[9.]
   @test A.onelab.get_number("Remote/result")==[9.]
   @test isempty(A._ONELAB_ACTIVE_ENGINES)
  end
 finally
  A.finalize()
 end
end

end
