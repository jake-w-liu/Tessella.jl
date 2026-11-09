module SurfaceOrderChecks
using Test,Tessella
const F=Tessella.IO
const E=Tessella.Elements

function specimen(msh)
    coords=E.lagrange_nodes(msh)
    # Reverse file-side node IDs to check selection of local corners rather
    # than the first rows of the coordinate table.
    indices=reverse(Int32.(1:size(coords,2)))
    return E.MixedMesh(coords,[E.ElementBlock(msh,reshape(indices,:,1))])
end

@testset "OFF and PLY2 project all catalog surface orders to corners" begin
    mktempdir() do dir
        for family in (:tri,:qua)
            types=unique(E.msh_type(family,p;serendipity=incomplete)
                         for p in 1:10 for incomplete in (false,true))
            for msh in types
                mesh=specimen(msh)
                corners=family===:tri ? 3 : 4
                expected=Int.(mesh.blocks[1].nodes[1:corners,1])
                off=joinpath(dir,"order.off")
                @test F.write_off(off,mesh)==off
                lines=readlines(off)
                @test parse.(Int,split(lines[1])[2:4])==[size(mesh.coords,2),1,0]
                @test parse.(Int,split(lines[end]))==[corners;expected.-1]
                _,flat,records=F._read_off(off,"surface order regression")
                @test length(records)==1 && records[1].msh==(family===:tri ? 2 : 3)
                # OFF deliberately uses six-significant-digit %g rows.
                @test isapprox(reshape(flat,3,:)[:,records[1].nodes],
                               mesh.coords[:,expected];rtol=5e-6,atol=1e-12)

                ply=joinpath(dir,"order.ply2")
                @test F.write_ply2(ply,mesh)==ply
                lines=readlines(ply)
                @test parse(Int,lines[1])==size(mesh.coords,2)
                @test parse(Int,lines[2])==(family===:tri ? 1 : 0)
                if family===:tri
                    @test parse.(Int,split(lines[end]))==[3;expected.-1]
                    _,flat,records=F._read_ply2(ply,"surface order regression")
                    @test length(records)==1 && records[1].msh==2
                    @test reshape(flat,3,:)[:,records[1].nodes]≈mesh.coords[:,expected]
                else
                    @test length(lines)==2+size(mesh.coords,2)
                end
            end
        end
        # Type 31 is a fifth-order tetrahedron, not a triangle. Neither deck
        # may create a spurious surface cell from its 56 interpolation nodes.
        mesh=specimen(31)
        off=joinpath(dir,"volume.off")
        F.write_off(off,mesh)
        @test parse(Int,split(first(readlines(off)))[3])==0
        ply=joinpath(dir,"volume.ply2")
        F.write_ply2(ply,mesh)
        @test parse(Int,readlines(ply)[2])==0
        empty=E.MixedMesh(zeros(3,0),E.ElementBlock[])
        F.write_off(off,empty);F.write_ply2(ply,empty)
        @test read(off,String)=="OFF 0 0 0\n"
        @test read(ply,String)=="0\n0\n"
    end
end

@testset "surface decks validate mutable sources before replacing files" begin
    mktempdir() do dir
        for (ext,writer) in (("off",F.write_off),("ply2",F.write_ply2))
            path=joinpath(dir,"held."*ext)
            for corrupt in (m->(m.blocks[1].nodes[1,1]=99),
                            m->(m.blocks[1].nodes[1,1]=0),
                            m->(m.coords[1,1]=NaN),
                            m->pop!(m.blocks[1].tags))
                mesh=specimen(9)
                corrupt(mesh)
                write(path,"held output")
                before=readdir(dir)
                @test_throws ArgumentError writer(path,mesh)
                @test read(path,String)=="held output"
                @test readdir(dir)==before
            end
            # Order-zero triangles have one interpolation node and cannot
            # stand in for a three-corner polygon.
            write(path,"held output")
            @test_throws ArgumentError writer(path,specimen(E.msh_type(:tri,0)))
            @test read(path,String)=="held output"
        end
        off=joinpath(dir,"held.off")
        @test_throws ArgumentError F.write_off(off,specimen(E.msh_type(:qua,0)))
        @test read(off,String)=="held output"
    end
end

@testset "public OFF and PLY2 export accepts higher-order surfaces" begin
    A=Tessella.API
    mktempdir() do dir
        A.initialize()
        try
            A.option("General.Terminal",0)
            A.option("Mesh.SaveAll",1)
            A.model.add("surface-orders")
            A.model.add_discrete_entity(2,1)
            tri=E.lagrange_nodes(E.msh_type(:tri,6))
            quad=E.lagrange_nodes(E.msh_type(:qua,6))
            quad[1,:].+=3
            coords=hcat(tri,quad)
            nt=size(tri,2);nn=size(coords,2)
            A.mesh.add_nodes(2,1,collect(1:nn),vec(coords))
            A.mesh.add_elements_by_type(1,E.msh_type(:tri,6),[1],collect(1:nt))
            A.mesh.add_elements_by_type(1,E.msh_type(:qua,6),[2],collect(nt+1:nn))
            for ext in ("off","ply2")
                A.model.set_current("surface-orders")
                path=joinpath(dir,"orders."*ext)
                @test A.write(path)===nothing
                @test A.open(path)===nothing
                @test length(A.mesh.get_nodes()[1])==(ext=="off" ? 7 : 3)
                types,tags,_=A.mesh.get_elements()
                @test Set(types)==(ext=="off" ? Set([2,3]) : Set([2]))
                @test sum(length,tags)==(ext=="off" ? 2 : 1)
            end
        finally
            A.finalize()
        end
    end
end
end
