module BdfLargeChecks
using Test,Tessella
const F=Tessella.IO

function largecard(name,values)
    lines=String[]
    for (chunk,start) in enumerate(1:4:length(values))
        tag=chunk==1 ? name*"*" : "*C$(chunk-1)"
        fields=values[start:min(start+3,length(values))]
        body=rpad(tag,8)*join(lpad(string(value),16) for value in fields)
        if start+4<=length(values)
            body=rpad(body,72)*rpad("*C$chunk",8)
        end
        push!(lines,body)
    end
    return join(lines,"\n")*"\n"
end

@testset "BDF standard sixteen-column fields retain geometry and identifiers" begin
    mktempdir() do dir
        path=joinpath(dir,"large.bdf")
        ids=[10_000_000_001,10_000_000_002,10_000_000_003]
        coords=[1.25 2.25 1.25;2.5 2.5 3.5;3.75 3.75 3.75]
        nodes=join(largecard("GRID",Any[id,"",coords[:,i]...])
                   for (i,id) in enumerate(ids))
        write(path,nodes*largecard("CTRIA3",[23,7,ids...])*"ENDDATA\n")
        mesh=F.read_bdf(path)
        @test mesh.node_tags==ids
        @test mesh.mesh.coords==coords
        @test length(mesh.mesh.blocks)==1 && mesh.mesh.blocks[1].msh==2
        @test mesh.mesh.blocks[1].nodes==reshape(Int32[1,2,3],3,1)
        @test mesh.mesh.blocks[1].tags==[7]
        @test mesh.element_tags==[[23]]

        # A large first line always has four data slots, even when the
        # physical line ends early. Y then defaults to zero; Z is on *C1.
        write(path,rpad("GRID*",8)*lpad("1",16)*lpad("",16)*lpad("1.25",16)*
                   "\n"*rpad("*C1",8)*lpad("3.75",16)*"\nENDDATA\n")
        @test F.read_bdf(path).mesh.coords==reshape([1.25,0.,3.75],3,1)

        # Large connectivity is checked across several continuation lines.
        grid=join(largecard("GRID",Any[i,"",Float64(i),0.,0.]) for i in 1:20)
        for (name,msh,arity) in (("CTETRA",11,10),("CPENTA",18,15),
                                 ("CHEXA",17,20),("CPYRA",19,13))
            write(path,grid*largecard(name,[17,9,collect(1:arity)...])*"ENDDATA\n")
            mesh=F.read_bdf(path)
            @test length(mesh.mesh.blocks)==1 && mesh.mesh.blocks[1].msh==msh
            @test mesh.mesh.blocks[1].nodes==reshape(Int32.(1:arity),arity,1)
            @test mesh.mesh.blocks[1].tags==[9]
        end
    end
end

@testset "BDF explicit higher-order nodes cannot silently become linear" begin
    mktempdir() do dir
        path=joinpath(dir,"solid.bdf")
        grid=join("GRID,$i,,$i.,0.,0.\n" for i in 1:20)
        for (name,linear,low,high) in (("CTETRA",4,4,10),("CPENTA",6,6,15),
                                      ("CHEXA",5,8,20),("CPYRA",7,5,13))
            for tail in (fill("",high-low),fill("0",high-low))
                values=[string.(1:low);tail]
                write(path,grid*"$name,1,7,"*join(values,",")*"\nENDDATA\n")
                @test F.read_bdf(path).mesh.blocks[1].msh==linear
            end
            for final in ("999","nonnumeric","-1")
                values=string.(1:high)
                values[end]=final
                write(path,grid*"$name,1,7,"*join(values,",")*"\nENDDATA\n")
                @test_throws ArgumentError F.read_bdf(path)
            end
            values=[string.(1:low);fill("",high-low)]
            values[end]=string(high)
            write(path,grid*"$name,1,7,"*join(values,",")*"\nENDDATA\n")
            @test_throws ArgumentError F.read_bdf(path)
        end
    end
end
end
