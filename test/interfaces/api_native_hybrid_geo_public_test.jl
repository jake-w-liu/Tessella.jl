module NativeHybridGeoPublicChecks
using Test,TOML,SHA,Tessella
const A=Tessella.API
const M=Tessella.Model
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const native=TOML.parsefile(joinpath(ROOT,"test/artifacts/native_hybrid_geo_import_4152_v16.toml"))
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]
function check_current(case)
    @test A.model.get_entities()==Tuple{Int,Int}[Tuple(Int.(p["owner"])) for p in case["final"]]
    @test (A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())==
        (UInt64(case["final_node_max"]),UInt64(case["final_element_max"]))
    for product in case["final"]
        dim,tag=product["owner"]
        uid,xyz,uv=A.mesh.get_nodes(dim,tag,false,true)
        @test uid==UInt64.(product["node_tags"])
        @test isequal(xyz,Float64[x for p in product["xyz_hex"] for x in hexpoint(p)])
        @test isequal(uv,hexpoint(product["uv_hex"]))
        cells=get(product,"cell",[])
        kinds,labels,nodes=A.mesh.get_elements(dim,tag)
        @test kinds==Int32[cell["msh"] for cell in cells]
        @test labels==[UInt64.(cell["labels"]) for cell in cells]
        @test nodes==[UInt64[node for corners in cell["nodes"] for node in corners] for cell in cells]
    end
    uid,xyz,uv=A.mesh.get_nodes(-1,-1,false,true)
    @test uid==UInt64[node for p in case["final"] for node in p["node_tags"]]
    @test isequal(xyz,Float64[x for p in case["final"] for point in p["xyz_hex"] for x in hexpoint(point)])
    @test isempty(uv)
end
@testset "literal complete GEO Mesh2/Delete/Merge/Mesh3 native hybrid histories" begin
    @test bytes2hex(sha256(read(joinpath(ROOT,native["primary"]))))==native["primary_sha256"]
    sourcefile=joinpath(ROOT,"test/data/native_hybrid_import_events_4152_v1/native_hybrid_import_source_4152_v18.geo")
    @test bytes2hex(sha256(read(sourcefile)))=="477c7c9cc9deb9803376d0c5197c432b6a4bb333d48af12afaaca5f99f2fccec"
    @test length(native["case"])==9
    for case in native["case"]
        @testset "$(case["name"])" begin
            println("HYBRID_GEO_PUBLIC_CASE ",case["name"]);flush(stdout)
            A.initialize()
            try
                A.option("General.Verbosity",0.)
                mktempdir() do directory
                    meshfile=joinpath(directory,"input.msh")
                    cp(joinpath(ROOT,"test/data/native_hybrid_import_events_4152_v1",case["name"]*".msh"),meshfile)
                    geo=joinpath(directory,"workflow.geo")
                    body=read(sourcefile,String)*(case["warm"] ? "\nMesh 2; Delete Meshes;\n" : "")*
                        "\nMerge "*repr(replace(meshfile,'\\'=>'/'))*";\nMesh.MeshOnlyEmpty=1; Mesh 3;\n"
                    write(geo,body)
                    opened=A.open_geo!(geo;mesh_dim=0)
                    @test opened.mesh!==nothing
                    check_current(case)
                    before_nodes=A.mesh.get_nodes();before_cells=A.mesh.get_elements()
                    clocks=(A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
                    A._with_model() do m
                        authority=m.meshing.mesh_authority
                        @test authority!==nothing && authority.mesh!==nothing
                        @test Tessella.MeshTypes.validate(authority.mesh).ok
                        @test (3,Int32(1)) in authority.done
                        @test length(authority.surface_constructor_inputs)==6
                        for (owner,constructor) in authority.surface_constructor_inputs
                            @test owner[1]==2
                            current=M._model_mesh_product(authority,2,Int(owner[2]);caller="GEO public native constructor")
                            @test M._model_surface_constructor_binding(current,constructor,"GEO public native constructor")===nothing
                        end
                    end
                    @test isequal(A.mesh.get_nodes(),before_nodes)
                    @test A.mesh.get_elements()==before_cells
                    @test (A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())==clocks
                    println("HYBRID_GEO_PUBLIC_FINAL ",case["name"]," clocks=",clocks);flush(stdout)
                end
            finally
                A.finalize()
            end
        end
    end
end
end
