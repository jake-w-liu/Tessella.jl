# Independent stdlib-only process fixture: no Tessella or external solver load.
using Sockets

function frame(socket,kind,payload="")
    write(socket,htol(Int32(kind)));write(socket,htol(Int32(ncodeunits(payload))))
    write(socket,payload);flush(socket)
end
function reply(socket)
    kind=Int(ltoh(read(socket,Int32)));count=Int(ltoh(read(socket,Int32)))
    0<=count<=1024*1024 || error("invalid response length")
    return kind,String(read(socket,count))
end
function number(name,values;label="",attributes=String[])
    join(vcat(["1.3","number",name,label,"fixture help","31","1","0",
        string(length(attributes)÷2)],attributes,["0",string(length(values))],
        string.(values),["-10","100","0.25","-1","2","1","2","1","1","one"]),'\0')
end
function text(name,values)
    join(vcat(["1.3","string",name,"","","31","1","0","0","0",string(length(values))],
        values,["generic","0"]),'\0')
end
function query(socket,name;kind="number")
    payload=kind=="number" ? number(name,[0.]) : text(name,[""])
    frame(socket,24,payload);return reply(socket)
end

function main()
mode=ARGS[1];output=ARGS[2]
position=findfirst(==("-onelab"),ARGS)
position===nothing && error("missing client arguments")
name=ARGS[position+1];host,port=split(ARGS[position+2],':')
mode=="nostart" && (write(output,string(getpid()));exit(0))
mode=="fail" && exit(7)
if mode=="relay"
    write(output,ARGS[position+2])
    while !isfile(output*".done");sleep(.01);end
    exit(0)
end
socket=connect(host,parse(Int,port))
frame(socket,1,string(getpid()))
if mode=="basic"
    kind,payload=query(socket,"Input/seed")
    kind==23 && occursin("Input/seed",payload) || error("missing seed")
    frame(socket,23,number("Output/λ",[4.,5.];label="result λ",attributes=["Highlight","Blue"]))
    frame(socket,32,text("Output/text",["hello λ","second"]))
    frame(socket,34,"set\0"*name*"\0false");frame(socket,34,"get\0"*name)
    reply(socket)==(34,"false") || error("client changed flag")
    frame(socket,25,"1.3\0number")
    count=0
    while true
        kind,payload=reply(socket)
        kind==26 && break
        kind==25 || error("list reply");count+=1
    end
    frame(socket,10,"fixture ready λ")
    write(output,join([name,string(count)],"\n"))
elseif mode=="nested"
    command=ARGS[3]
    frame(socket,27,"nested child\0"*command)
    reply(socket)==(2,"nested child") || error("child completion")
    # Sequential reuse must not consume the simultaneous-client budget.
    frame(socket,27,"nested child\0"*command)
    reply(socket)==(2,"nested child") || error("second child completion")
    write(output,"nested complete")
elseif mode=="child"
    frame(socket,23,number("Output/child",[9.]))
    open(output,"a") do io;println(io,name);end
elseif mode=="parse"
    frame(socket,21,ARGS[3])
    # A query is a barrier after the command, not a timing assumption.
    query(socket,"Input/seed")
    write(output,"parsed")
elseif mode=="file"
    frame(socket,parse(Int,ARGS[3]),ARGS[4])
    query(socket,"Input/seed")
    write(output,"file complete")
elseif mode=="preprocess"
    frame(socket,28,name*"\0"*ARGS[3])
    write(output,string(reply(socket)))
elseif mode=="metadata"
    records=String[]
    for (kind,key) in (("number","Controls/N"),("string","Controls/S"))
        message,payload=query(socket,key;kind)
        message==23 || error("missing metadata parameter")
        push!(records,payload)
    end
    write(output,join(records,'\n'))
elseif mode=="noread"
    # Request the full string list and then never read: the server's response
    # write backpressures. Exercise the finite watchdog and cleanup paths.
    frame(socket,25,"1.3\0string")
    sleep(120)
elseif mode=="hang"
    write(output,string(getpid()));sleep(120)
elseif mode=="tree"
    descendant=output*".descendant"
    script="write("*repr(descendant)*",string(getpid()));sleep(120)"
    run(`$(Base.julia_cmd()) --startup-file=no --compile=min -O0 -e $script`;wait=false)
    write(output,string(getpid()));sleep(120)
elseif mode=="badlength"
    write(socket,htol(Int32(23)));write(socket,htol(typemax(Int32)));flush(socket)
    sleep(120)
elseif mode=="partial"
    write(socket,htol(Int32(23)));write(socket,htol(Int32(1024*1024)))
    write(socket,"short");flush(socket);close(socket);exit(0)
elseif mode=="nostop"
    write(output,string(getpid()));close(socket);exit(0)
elseif mode=="fast"
    write(output,name)
elseif mode=="signal"
    frame(socket,2,"Goodbye!");close(socket)
    ccall(:raise,Cint,(Cint,),Base.SIGTERM)
    exit(99)
else
    error("unknown fixture mode $mode")
end
frame(socket,2,"Goodbye!");close(socket)
end
main()
