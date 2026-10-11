module QueuedHomologyTests
using Tessella,Test,TOML,SHA
const A=Tessella.API
const PRIMARY=joinpath(@__DIR__,"..","artifacts","queued_homology_native_v1.toml")
const NATIVE=TOML.parsefile(PRIMARY)
bits(values)=string.(reinterpret(UInt64,Float64.(values));base=16,pad=16)
cell_rows(types,tags,nodes)=[Dict("type"=>Int(kind),"tags"=>string.(labels),"nodes"=>string.(connected))
    for (kind,labels,connected) in zip(types,tags,nodes)]

function state()
    nodes,xyz,uv=A.mesh.get_nodes();details=Any[];entities=Any[]
    for node in nodes
        co,pa,d,t=A.mesh.get_node(node)
        push!(details,Dict("uid"=>string(node),"owner"=>[d,t],"xyz_bits"=>bits(co),"uv_bits"=>bits(pa)))
    end
    for (d,t) in A.model.get_entities()
        en,co,pa=A.mesh.get_nodes(d,t,false,true)
        push!(entities,Dict("owner"=>[d,t],"nodes"=>string.(en),"xyz_bits"=>bits(co),"uv_bits"=>bits(pa),
            "cells"=>cell_rows(A.mesh.get_elements(d,t)...)))
    end
    return Dict("nodes"=>string.(nodes),"xyz_bits"=>bits(xyz),"uv_bits"=>bits(uv),"node_details"=>details,
        "cells"=>cell_rows(A.mesh.get_elements()...),"entities"=>entities,
        "groups"=>[Dict("owner"=>[d,t],"name"=>A.model.get_physical_name(d,t),
            "entities"=>A.model.get_entities_for_physical_group(d,t)) for (d,t) in A.model.get_physical_groups()],
        "max_node"=>string(A.mesh.get_max_node_tag()),"max_element"=>string(A.mesh.get_max_element_tag()),
        "last_error"=>A.logger.get_last_error())
end

function check(expected)
    actual=state()
    for key in keys(expected)
        @test actual[key]==expected[key]
    end
end

@testset "API generation executes queued homology and repeated native requests" begin
    @test bytes2hex(sha256(replace(read(PRIMARY,String),"\r\n"=>"\n")))==
        "c385f369dd150ac58efddea2fa5f7137c4bb623b82c7af65ade75be595a8dcc7"
    @test NATIVE["gmsh_version"]=="4.15.2" && length(NATIVE["cases"])==4
    for case in NATIVE["cases"]
        @testset "dim$(case["dimension"])-OnlyEmpty$(case["only_empty"])" begin
            dimension=case["dimension"]
            A.initialize()
            try
                A.option("General.AbortOnError",0);A.option("Mesh.Renumber",0)
                A.option("Mesh.MeshOnlyEmpty",case["only_empty"])
                points=dimension==1 ? ((0.,0.),(1.,0.)) : ((0.,0.),(1.,0.),(1.,1.),(0.,1.))
                for (tag,(x,y)) in enumerate(points);A.model.add_point(x,y,0.;tag,meshSize=1.);end
                ends=dimension==1 ? ((1,2),) : ((1,2),(2,3),(3,4),(4,1))
                for (tag,(a,b)) in enumerate(ends);A.model.add_line(a,b;tag);A.mesh.set_transfinite_curve(tag,2);end
                if dimension==2
                    A.model.add_curve_loop([1,2,3,4];tag=1);A.model.add_plane_surface([1];tag=1)
                    A.mesh.set_transfinite_surface(1)
                end
                A.model.add_physical_group(dimension,[1];tag=11)
                A.mesh.generate(dimension);check(case["before"])
                A.mesh.add_homology_request("Homology",[11],Int[],collect(0:dimension))
                A.mesh.generate(dimension);check(case["after"])
                A.mesh.generate(dimension);check(case["repeated"])
            finally A.finalize() end
        end
    end
end
end
