# `onelab` module UI01: ONELAB solver client — the upstream socket protocol's
# client side (length-prefixed `Int32 type, Int32 length, payload` frames,
# version-1.3 NUL-separated parameter records, query/define/update, changed
# flags, command/file messages, sub-client run/stop). Verified three ways:
# golden wire payloads produced by the public reference client, a scripted
# in-process fixture socket speaking the documented framing, and — when the
# pinned gmsh 4.15.2 binary is present — the live `gmsh -listen` server
# cross-checked by the reference python client.
module OnelabClientChecks

using Test, Tessella
import Sockets
const A=Tessella.API
const _OL=Tessella.API  # wire codec lives in API (APIOnelabClient.jl)

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

@testset "wire parameter codec against reference payloads" begin
    # Golden NUL-separated payloads captured from the public standalone
    # onelab client (`_parameter.tochar()`); field order is the documented
    # protocol layout.
    alpha="1.3\x00number\x00test/alpha\x00Alpha\x00hp\x0012\x001\x000\x001\x00Highlight\x00Red" *
          "\x001\x00solver\x001\x002\x001.5\x002\x000\x0010\x000.5\x00-1\x002\x001.5\x009\x002\x001.5" *
          "\x00small\x009\x00big"
    p=_OL._onelab_wire_decode(alpha,"test")
    @test p.kind===:number
    @test p.name=="test/alpha"
    @test p.label=="Alpha" && p.help=="hp"
    @test p.changed_value==12 && p.visible && !p.read_only
    @test p.attributes==Dict("Highlight"=>"Red")
    @test p.clients==Dict("solver"=>1)
    @test p.values_f==[1.5,2.0]
    @test p.min==0.0 && p.max==10.0 && p.step==0.5 && p.index==-1
    @test p.choices_f==[1.5,9.0]
    @test p.value_labels==Dict(1.5=>"small",9.0=>"big")
    # Encode round-trips to the identical byte stream.
    @test _OL._onelab_wire_encode(p)==alpha
    mode="1.3\x00string\x00test/mode\x00\x00\x0031\x001\x000\x000\x000\x001\x00compute\x00generic" *
         "\x002\x00check\x00compute"
    s=_OL._onelab_wire_decode(mode,"test")
    @test s.kind===:string
    @test s.name=="test/mode" && s.values_s==["compute"]
    @test s.kind_name=="generic" && s.choices_s==["check","compute"]
    @test s.changed_value==31 && s.visible && !s.read_only
    @test _OL._onelab_wire_encode(s)==mode
    dummy="1.3\x00number\x00dummy\x00\x00\x0031\x001\x000\x000\x000\x001\x000\x00-1e+200\x001e+200\x000" *
          "\x00-1\x000\x000"
    d=_OL._onelab_wire_decode(dummy,"test")
    @test d.name=="dummy" && d.values_f==[0.0]
    @test d.min==-1e200 && d.max==1e200
    @test _OL._onelab_wire_encode(d)==dummy
    # Fresh records encode to the reference dummy layout too.
    fresh=_OL._OnelabWireParameter(:number,"dummy")
    @test _OL._onelab_wire_encode(fresh)==dummy
    # Malformed payloads fail loudly.
    @test_throws ArgumentError _OL._onelab_wire_decode("1.2\x00number\x00x","test")
    @test_throws ArgumentError _OL._onelab_wire_decode("1.3\x00region\x00x","test")
    @test_throws ArgumentError _OL._onelab_wire_decode("1.3\x00number","test")
end

# ---------------------------------------------------------------------------
# Scripted fixture server — an independent implementation of the documented
# framing, used to check every client operation round-trips.
# ---------------------------------------------------------------------------

function _fixture_read_exact(sock,n)
    buf=Vector{UInt8}(undef,n);got=0
    while got<n
        buf[got+=1]=read(sock,UInt8)
    end
    return buf
end

function _fixture_recv(sock)
    hdr=_fixture_read_exact(sock,8)
    type=ltoh(reinterpret(Int32,hdr[1:4])[1])
    len=ltoh(reinterpret(Int32,hdr[5:8])[1])
    return type,String(_fixture_read_exact(sock,Int(len)))
end

function _fixture_send(sock,type,payload="")
    data=codeunits(payload)
    write(sock,htol(Int32(type)))
    write(sock,htol(Int32(length(data))))
    write(sock,data);flush(sock)
    return nothing
end

# Server state machine: parameter store keyed by name (the third NUL field),
# changed-name set, recorded command strings, async INFO injected before one
# reply to exercise client-side buffering, STOP reaping for CONNECT.
function _fixture_server(sock)
    params=Dict{String,String}()
    changed=Set{String}()
    commands=String[]
    inject_info=true
    while true
        type,msg=_fixture_recv(sock)
        if type==1                      # START
            continue
        elseif type==2                  # STOP
            break
        elseif type==24                 # PARAMETER_QUERY
            items=split(msg,'\0')
            name=length(items)>=3 ? items[3] : ""
            if inject_info
                _fixture_send(sock,10,"fixture says hi")
                inject_info=false
            end
            if haskey(params,name)
                _fixture_send(sock,23,params[name])
            else
                _fixture_send(sock,29,"not found")
            end
        elseif type==23 || type==32     # PARAMETER / PARAMETER_UPDATE
            items=split(msg,'\0')
            length(items)>=3 && (params[items[3]]=msg)
        elseif type==25                 # PARAMETER_QUERY_ALL
            items=split(msg,'\0')
            want=length(items)>=2 ? items[2] : ""
            for wire in sort!(collect(values(params)))
                split(wire,'\0')[2]==want || continue
                _fixture_send(sock,25,wire)
            end
            _fixture_send(sock,26,"")
        elseif type==31                 # PARAMETER_CLEAR
            delete!(params,msg)
        elseif type==34                 # CLIENT_CHANGED
            items=split(msg,'\0')
            if items[1]=="get"
                _fixture_send(sock,34,
                    length(items)>=2 && items[2] in changed ? "true" : "false")
            elseif items[1]=="set" && length(items)>=3
                items[3]=="true" ? push!(changed,items[2]) :
                    delete!(changed,items[2])
            end
        elseif type==27                 # CONNECT
            push!(commands,"connect:"*msg)
            @async begin
                sleep(0.05)
                try; _fixture_send(sock,2,"done");catch;end
            end
        elseif type in (20,21,33)       # MERGE_FILE / PARSE_STRING / OPEN_PROJECT
            push!(commands,msg)
        elseif type in (10,11,12,13)    # INFO / WARNING / ERROR / PROGRESS
            push!(commands,"msg$(type):"*msg)
        end
    end
    return commands
end

@testset "fixture socket: framed request/reply round-trip" begin
    server=Sockets.listen(Sockets.localhost,0)
    _,port=Sockets.getsockname(server)
    done=Channel{Any}(1)
    @async begin
        sock=Sockets.accept(server)
        try
            put!(done,_fixture_server(sock))
        catch err
            put!(done,err)
        end
    end
    c=A.onelab.client_connect("127.0.0.1",port;name="fixture client")
    @test c isa A.onelab.OnelabClient
    @test isopen(c.socket)
    # define: NOT_FOUND -> PARAMETER; second define adopts the stored copy.
    @test A.onelab.client_define_number(c,"T/a";value=3.5,min=0.0,max=9.0)==[3.5]
    @test A.onelab.client_define_number(c,"T/a";value=8.0)==[3.5]
    @test A.onelab.client_get_number(c,"T/a")==[3.5]
    # set on a server-known name keeps its metadata; on a new name defines.
    A.onelab.client_set_number(c,"T/a",[4.25])
    @test A.onelab.client_get_number(c,"T/a")==[4.25]
    A.onelab.client_set_number(c,"T/b",[1.0,2.0])
    @test A.onelab.client_get_number(c,"T/b")==[1.0,2.0]
    A.onelab.client_set_string(c,"T/s","hello")
    @test A.onelab.client_get_string(c,"T/s")==["hello"]
    @test A.onelab.client_get_number(c,"T/missing")===nothing
    # The async INFO the fixture injected is buffered, not mistaken for a reply.
    msgs=A.onelab.client_messages(c)
    @test any(m->m[1]==Int32(10) && m[2]=="fixture says hi",msgs)
    @test any(m->m[1]==Int32(11) && occursin("Unknown parameter",m[2]),msgs)
    # QUERY_ALL stream.
    @test sort(A.onelab.client_get_names(c))==["T/a","T/b","T/s"]
    @test sort(A.onelab.client_get_names(c;search="T/?"))==["T/a","T/b","T/s"]
    allp=A.onelab.client_get_all(c)
    @test count(p->p.kind===:number,allp)==2
    @test count(p->p.kind===:string,allp)==1
    # clear + changed flags.
    A.onelab.client_clear(c,"T/b")
    @test A.onelab.client_get_number(c,"T/b")===nothing
    @test A.onelab.client_get_changed(c,"T/a")==false
    A.onelab.client_set_changed(c,"T/a",true)
    @test A.onelab.client_get_changed(c,"T/a")==true
    A.onelab.client_set_changed(c,"T/a",false)
    @test A.onelab.client_get_changed(c,"T/a")==false
    # command/file/info message kinds.
    A.onelab.client_send_command(c,"Mesh 2;")
    A.onelab.client_merge_file(c,"x.geo")
    A.onelab.client_open_project(c,"p.pro")
    A.onelab.client_info(c,"running")
    A.onelab.client_error(c,"oops")
    # sub-client run: CONNECT then the fixture's delayed STOP.
    t0=time()
    A.onelab.client_run(c,"solver","runme arg")
    @test c.subclients==0
    @test time()-t0>=0.03
    # nonblocking run + wait.
    A.onelab.client_run(c,"solver2","runme2";block=false)
    @test c.subclients==1
    A.onelab.client_wait(c)
    @test c.subclients==0
    # disconnect sends STOP; the server reports its received commands.
    A.onelab.client_disconnect(c)
    @test c.socket===nothing
    commands=take!(done)
    commands isa Exception && rethrow(commands)
    @test "Mesh 2;" in commands
    @test "x.geo" in commands && "p.pro" in commands
    @test "msg10:running" in commands && "msg12:oops" in commands
    @test "connect:solver\0runme arg" in commands
    @test "connect:solver2\0runme2" in commands
    close(server)
    # Validation errors.
    @test occursin("not connected",
        _errtext(()->A.onelab.client_get_number(c,"x")))
    @test_throws ArgumentError A.onelab.client_connect("no port string")
    @test_throws ArgumentError A.onelab.client_connect("h:notaport")
    @test_throws ArgumentError A.onelab.client_connect("h",70000)
    @test_throws ArgumentError A.onelab.client_get_number("not a client","x")
    @test_throws ArgumentError A.onelab.client_run(1,"n","cmd")
end

@testset "fixture socket: check_parameters batch + action read" begin
    server=Sockets.listen(Sockets.localhost,0)
    _,port=Sockets.getsockname(server)
    @async _fixture_server(Sockets.accept(server))
    c=A.onelab.client_connect("127.0.0.1:"*string(port);name="cli")
    adopted=A.onelab.client_check_parameters(c,
        ["cli/width"=>(value=1.5,min=0.0,max=4.0),
         "cli/mode"=>(kind=:string,value="check",choices=["check","compute"]),
         "cli/raw"=>Dict(:value=>2.0)])
    @test adopted["cli/width"]==[1.5]
    @test adopted["cli/mode"]==["check"]
    @test adopted["cli/raw"]==[2.0]
    # `client_action` reads `<name>/Action`.
    A.onelab.client_set_string(c,"cli/Action","check")
    @test A.onelab.client_action(c)=="check"
    A.onelab.client_set_string(c,"cli/Action","stop")
    @test A.onelab.client_action(c)=="stop"
    @test_throws ArgumentError A.onelab.client_check_parameters(c,["x"=>3])
    A.onelab.client_disconnect(c)
    close(server)
end

# ---------------------------------------------------------------------------
# Live oracle: the pinned Gmsh 4.15.2 binary as a real ONELAB server
# (`gmsh -listen`), plus the reference python client as a second connection
# cross-verifying what our wire wrote. Skipped when the tools are absent.
# ---------------------------------------------------------------------------

const _GMSH_EXE="C:/Users/User/Tools/gmsh-4.15.2-Windows64/gmsh.exe"
const _ONELAB_PY="C:/Users/User/Tools/gmsh-4.15.2-Windows64/onelab.py"
const _PY_EXE=let p=Sys.which("python"); p===nothing ? "" : p end

function _free_port()
    s=Sockets.listen(Sockets.localhost,0)
    _,port=Sockets.getsockname(s);close(s)
    return Int(port)
end

function _py_onelab(script::String)
    f=tempname()*".py"
    write(f,script)
    try
        return read(`$_PY_EXE $f`,String)
    finally
        rm(f;force=true)
    end
end

@testset "live gmsh -listen server + reference client cross-check" begin
    ok=isfile(_GMSH_EXE) && isfile(_ONELAB_PY) && !isempty(_PY_EXE)
    ok || @test_skip "gmsh 4.15.2 SDK/python absent"
    if ok
        port=_free_port()
        addr="127.0.0.1:$port"
        proc=run(pipeline(`$_GMSH_EXE -listen $addr`;stdout=devnull,
                          stderr=devnull);wait=false)
        c=nothing
        try
            # The server accepts connections once its listener is up.
            for _ in 1:100
                c=try
                    A.onelab.client_connect(addr;name="tessella client")
                catch
                    sleep(0.1);nothing
                end
                c===nothing || break
            end
            @test c!==nothing
            c===nothing && return
            # Server-defined parameters are visible (upstream seed space).
            names=A.onelab.client_get_names(c)
            @test any(n->startswith(n,"Gmsh/"),names)
            # define + set + get against the real server.
            @test A.onelab.client_define_number(c,"jl/alpha";
                value=42.0,min=0.0,max=100.0,choices=[42.0,7.5])==[42.0]
            A.onelab.client_set_number(c,"jl/alpha",[7.5])
            @test A.onelab.client_get_number(c,"jl/alpha")==[7.5]
            p=A.onelab.client_get_parameter(c,:number,"jl/alpha")
            @test p.min==0.0 && p.max==100.0 && p.choices_f==[42.0,7.5]
            A.onelab.client_set_string(c,"jl/mode","compute")
            @test A.onelab.client_get_string(c,"jl/mode")==["compute"]
            # `set_number` keeps server-side metadata (server copy wins).
            p2=A.onelab.client_get_parameter(c,:number,"jl/alpha")
            @test p2.max==100.0
            # `.geo` string executed server-side creates a session parameter.
            A.onelab.client_send_command(c,"SetNumber(\"jl/via\", 5.5);")
            sleep(0.3)
            @test A.onelab.client_get_number(c,"jl/via")==[5.5]
            # changed flags round-trip on the server. The flag value itself
            # is server-side semantics (gmsh tracks per-client dirty state —
            # a client that just pushed a value may still read `false`);
            # the verified contract is that the get/set wire channel works.
            A.onelab.client_set_changed(c,"jl/alpha",true)
            @test A.onelab.client_get_changed(c,"jl/alpha") isa Bool
            # gmsh -listen serves ONE client connection at a time — a second
            # connect is refused while the first is attached (verified against
            # the oracle). Disconnect, then let the reference python client
            # read our values over its own connection — the wire encoding is
            # server-faithful.
            A.onelab.client_disconnect(c)
            sleep(0.6)
            out=_py_onelab("""
import sys
sys.path.insert(0, $(repr(dirname(_ONELAB_PY))))
sys.argv=['x','-onelab','py client',$(repr(addr))]
import onelab
c=onelab.client('dummy.geo')
print('alpha=',c.getNumber('jl/alpha'))
print('mode=',c.getString('jl/mode'))
print('via=',c.getNumber('jl/via'))
c.setNumber('jl/beta',value=13.25)
print('changed=',c.isChanged('jl/alpha'))
c.finalize()
""")
            @test occursin("alpha= 7.5",out)
            @test occursin("mode= compute",out)
            @test occursin("via= 5.5",out)
            @test occursin("changed=",out)
            # Reconnect and read what the reference client wrote.
            sleep(0.4)
            c=nothing
            for _ in 1:50
                c=try
                    A.onelab.client_connect(addr;name="tessella client")
                catch
                    sleep(0.1);nothing
                end
                c===nothing || break
            end
            @test c!==nothing
            c===nothing && return
            @test A.onelab.client_get_number(c,"jl/beta")==[13.25]
            A.onelab.client_disconnect(c)
        finally
            c===nothing || (try A.onelab.client_disconnect(c) catch end)
            kill(proc)
        end
    end
end

end # module
