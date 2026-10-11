module GeoOnlyEmptyHigherTests
using Test,Tessella
using Tessella.Elements:MixedMesh
const API=Tessella.API
const GEO=Tessella.GeoExec
const IO=Tessella.Elements
if !isdefined(parentmodule(@__MODULE__),:QuadTriNoNewB4FreeStripCertificates)
    Base.include(parentmodule(@__MODULE__),joinpath(@__DIR__,"quadtri_nonew_b4_free_strip_certificates.jl"))
end
const F=getfield(parentmodule(@__MODULE__),:QuadTriNoNewB4FreeStripCertificates)

const SQUARE="""
Mesh.Renumber=0; Mesh.MeshOnlyEmpty=1;
Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
Transfinite Curve{1:4}=3; Transfinite Surface{1};
"""

# MSH 2.2 records public labels and elementary cell owners directly; it has no
# per-node owner/UV section. Read those literal records independently rather
# than treating a storage-oriented reader's absent entity_data as their loss.
function saved_payload(path)
    lines=readlines(path);node_marker=findfirst(==("\$Nodes"),lines)
    element_marker=findfirst(==("\$Elements"),lines)
    node_marker!==nothing && element_marker!==nothing || error("missing MSH payload")
    node_count=parse(Int,lines[node_marker+1]);cell_count=parse(Int,lines[element_marker+1])
    coordinates=Dict{UInt64,NTuple{3,Float64}}()
    for line in lines[node_marker+2:node_marker+1+node_count]
        values=split(line);tag=parse(UInt64,values[1])
        haskey(coordinates,tag) && error("duplicate saved node label")
        coordinates[tag]=Tuple(parse.(Float64,values[2:4]))
    end
    cells=NamedTuple[];seen=Set{UInt64}()
    for line in lines[element_marker+2:element_marker+1+cell_count]
        values=parse.(UInt64,split(line));label=values[1];msh=Int(values[2]);ntags=Int(values[3])
        label in seen && error("duplicate saved element label");push!(seen,label)
        ntags>=2 || error("saved native cell lacks elementary identity")
        connection=Tuple(values[4+ntags:end])
        all(node->haskey(coordinates,node),connection) || error("unresolved saved node label")
        push!(cells,(;label,msh,owner=Int(values[5]),connection))
    end
    return (;coordinates,cells)
end

counts(mesh::NamedTuple)=Dict(msh=>count(cell->cell.msh==msh,mesh.cells)
    for msh in unique(cell.msh for cell in mesh.cells))
nodes(mesh::NamedTuple)=mesh.coordinates
connections(mesh::NamedTuple,dimension,tag)=Dict(msh=>[cell.connection for cell in mesh.cells
    if cell.msh==msh && cell.owner==tag] for msh in unique(cell.msh for cell in mesh.cells)
    if Tessella.Elements.msh_dimension(msh)==dimension)

function counts(mesh)
    out=Dict{Int,Int}()
    for block in mesh.blocks
        out[Int(block.msh)]=get(out,Int(block.msh),0)+size(block.nodes,2)
    end
    return out
end
nodes(mesh)=Dict(tag=>Tuple(mesh.coords[:,row])
    for (row,tag) in enumerate(mesh.entity_data.external_node_tags))
function connections(mesh,dimension,tag)
    data=mesh.entity_data
    return Dict(Int(block.msh)=>[Tuple(data.external_node_tags[cell])
        for (column,cell) in enumerate(eachcol(block.nodes))
        if data.block_entities[bi][column]==tag]
        for (bi,block) in enumerate(mesh.blocks)
        if Tessella.Elements.msh_dimension(block.msh)==dimension)
end

@testset "GEO actual higher products preserve generation and mesh operations" begin
    # The detached model adapter runs without initializing an API session.
    current=API.CURRENT[];cache=API.LAST_MESH[];classification=API.LAST_MESH_CLASS[]
    mktempdir() do directory
        source=joinpath(directory,"actual.geo")
        write(source,SQUARE*"""
        Mesh 2; Save "first.msh";
        Mesh 2; Save "second.msh";
        RefineMesh; Save "refined.msh";
        SetOrder 2; Save "quadratic.msh";
        SetOrder 1; RenumberMeshNodes; RenumberMeshElements; Save "renumbered.msh";
        Delete Meshes; Mesh 2;
        """)
        execution=GEO.execute_geo(source)
        first=saved_payload(joinpath(directory,"first.msh"))
        second=saved_payload(joinpath(directory,"second.msh"))
        refined=saved_payload(joinpath(directory,"refined.msh"))
        quadratic=saved_payload(joinpath(directory,"quadratic.msh"))
        renumbered=saved_payload(joinpath(directory,"renumbered.msh"))
        @test counts(first)==Dict(1=>8,2=>8,15=>4)
        @test counts(second)==counts(first)
        @test nodes(first)==nodes(second)
        @test connections(first,2,1)==connections(second,2,1)
        # MSH2.2 serialization densely numbers the written elements. The
        # generation factory identities are checked on actual products below.
        @test Set(cell.label for cell in second.cells)==Set(UInt64.(1:20))
        @test counts(refined)==Dict(1=>16,2=>32,15=>4)
        @test counts(quadratic)==Dict(8=>16,9=>32,15=>4)
        @test counts(renumbered)==counts(refined)
        @test Set(keys(renumbered.coordinates))==Set(UInt64.(1:length(renumbered.coordinates)))
        @test Set(cell.label for cell in renumbered.cells)==Set(UInt64.(1:52))
        @test execution.mesh isa MixedMesh && counts(execution.mesh)==counts(first)
        @test length(GEO.geo_entity_mesh(execution,2,1).blocks)==1
        @test API.CURRENT[]===current && API.LAST_MESH[]===cache && API.LAST_MESH_CLASS[]===classification
    end
end

@testset "GEO generation retains Point identities and renews native cells" begin
    mktempdir() do directory
        first_path=joinpath(directory,"first.geo")
        second_path=joinpath(directory,"second.geo")
        write(first_path,SQUARE*"Mesh 2;\n")
        write(second_path,SQUARE*"Mesh 2; Mesh 2;\n")
        first=GEO.execute_geo(first_path).mesh
        second=GEO.execute_geo(second_path).mesh
        first_data=first.entity_data;second_data=second.entity_data
        old_max=maximum(Iterators.flatten(first_data.external_element_tags))
        old_points=Dict(first_data.block_entities[bi][column]=>
            (first_data.external_element_tags[bi][column],
             first_data.external_node_tags[block.nodes[1,column]])
            for (bi,block) in enumerate(first.blocks) if block.msh==15
            for column in axes(block.nodes,2))
        new_points=Dict(second_data.block_entities[bi][column]=>
            (second_data.external_element_tags[bi][column],
             second_data.external_node_tags[block.nodes[1,column]])
            for (bi,block) in enumerate(second.blocks) if block.msh==15
            for column in axes(block.nodes,2))
        @test new_points==old_points
        @test minimum(second_data.external_element_tags[bi][column]
            for (bi,block) in enumerate(second.blocks) if block.msh!=15
            for column in axes(block.nodes,2))>old_max
        @test nodes(first)==nodes(second)
        @test connections(first,2,1)==connections(second,2,1)
    end
end

@testset "GEO completed Curve authority survives higher and zero phases" begin
    mktempdir() do directory
        path=joinpath(directory,"curve_phase.geo")
        write(path,SQUARE*"Mesh 1; Save \"curves.msh\"; Mesh 2; Mesh 0; Save \"all.msh\";\n")
        execution=GEO.execute_geo(path)
        curves=saved_payload(joinpath(directory,"curves.msh"))
        actual=saved_payload(joinpath(directory,"all.msh"))
        @test counts(curves)==Dict(1=>8,15=>4)
        @test counts(actual)==Dict(1=>8,2=>8,15=>4)
        actual_nodes=nodes(actual)
        @test all(pair->actual_nodes[first(pair)]==last(pair),nodes(curves))
        @test connections(actual,1,1)==connections(curves,1,1)
        @test counts(GEO.geo_entity_mesh(execution,1,1))==Dict(1=>2)
        @test counts(GEO.geo_entity_mesh(execution,2,1))==Dict(2=>8)
    end
end

@testset "GEO higher extrusion consumes real Source and copied-cap products" begin
    fixture=F.fixture("geo-actual-cap";strip_length=5,layers=:one)
    mktempdir() do directory
        path=joinpath(directory,"actual_cap.geo")
        write(path,fixture.source*"\nMesh.Renumber=0; Mesh.MeshOnlyEmpty=1; Mesh 2; Mesh 3;\n")
        execution=GEO.execute_geo(path)
        sweep=Int.(execution.lists["sweep"])
        source=GEO.geo_entity_mesh(execution,2,fixture.surface)
        cap=GEO.geo_entity_mesh(execution,2,sweep[1])
        volume=GEO.geo_entity_mesh(execution,3,sweep[2])
        @test source isa MixedMesh && cap isa MixedMesh && volume isa MixedMesh
        @test counts(source)==Dict(3=>5) && counts(cap)==Dict(2=>10)
        certificate=F.certify(volume,fixture,source;native_top=false)
        @test certificate.total==1 && certificate.ncell>0
        @test Tessella.Model._model_mesh_status(execution.model,execution.model.meshing.mesh_authority)==3
    end
end

@testset "GEO completed recombination publishes real entity products and factory labels" begin
    # Executed Gmsh4.15.2 primary84ade001...: Tri9/10 reserve temporary11;
    # angle2 creates Quad12. Every preceding angle rejection reserves one.
    square=replace(SQUARE,"Mesh.MeshOnlyEmpty=1"=>"Mesh.MeshOnlyEmpty=0",
        "Transfinite Curve{1:4}=3"=>"Transfinite Curve{1:4}=2")*
        "Mesh.RecombinationAlgorithm=0; Mesh 2;\n"
    for rejected in 0:2
        mktempdir() do directory
            path=joinpath(directory,"recombine.geo")
            write(path,square*repeat("Recombine Surface{1}=0.5;RecombineMesh;\n",rejected)*
                "Recombine Surface{1}=2;RecombineMesh;\n")
            execution=GEO.execute_geo(path)
            whole=execution.mesh;surface=GEO.geo_entity_mesh(execution,2,1)
            @test counts(surface)==Dict(3=>1)
            @test counts(whole)==Dict(1=>4,3=>1,15=>4)
            @test only(only(surface.entity_data.external_element_tags))==12+rejected
            @test surface.entity_data.external_node_tags[only(surface.blocks).nodes[:,1]]==UInt64[2,3,4,1]
            for tag in 1:4
                point=GEO.geo_entity_mesh(execution,0,tag)
                curve=GEO.geo_entity_mesh(execution,1,tag)
                @test only(only(point.entity_data.external_element_tags))==tag
                @test only(point.entity_data.external_node_tags)==tag
                @test only(only(curve.entity_data.external_element_tags))==4+tag
            end
            authority=execution.model.meshing.mesh_authority
            @test authority.mesh===whole && authority.carrier===whole
            @test allunique(authority.private_rank)
        end
    end
end
end
