module PublicCompactPrismFunctionalChecks
using Test, Tessella, LinearAlgebra, TOML

const E = Tessella.Elements
const PRISM_SOURCE = raw"""
Point(1)={0,0,0,0.3}; Point(2)={1,0,0,0.3}; Point(3)={0,1,0,0.3};
Point(4)={0,0,1,0.3}; Point(5)={1,0,1,0.3}; Point(6)={0,1,1,0.3};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,1};
Line(4)={4,5}; Line(5)={5,6}; Line(6)={6,4};
Line(7)={1,4}; Line(8)={2,5}; Line(9)={3,6};
Curve Loop(1)={1,2,3}; Curve Loop(2)={4,5,6};
Curve Loop(3)={1,8,-4,-7}; Curve Loop(4)={2,9,-5,-8};
Curve Loop(5)={3,7,-6,-9};
Plane Surface(1)={1}; Plane Surface(2)={2};
Plane Surface(3)={3}; Plane Surface(4)={4}; Plane Surface(5)={5};
Surface Loop(1)={1,2,3,4,5}; Volume(1)={1};
Mesh.TransfiniteTri=1;
Transfinite Curve{:}=4;
Transfinite Surface{:};
Recombine Surface{:};
Transfinite Volume{1};
"""

function used_rows(mesh)
    used=falses(size(mesh.coords,2))
    for block in mesh.blocks
        used[vec(block.nodes)].=true
    end
    return used
end

function cell_counts(mesh)
    result=Dict{Int,Int}()
    for block in mesh.blocks
        result[block.msh]=get(result,block.msh,0)+size(block.nodes,2)
    end
    return filter(pair->last(pair)>0,result)
end

function polygon_area(points)
    abs(sum(points[1,i]*points[2,mod1(i+1,size(points,2))]-
            points[2,i]*points[1,mod1(i+1,size(points,2))]
            for i in axes(points,2)))/2
end

# Independent containment for this convex planar fixture. This checks
# geometric conformity; production identity certification remains UID based.
function contains_triangle(polygon,triangle)
    tol=2e-10
    for point in eachcol(triangle)
        signed=Float64[]
        for i in axes(polygon,2)
            a=polygon[:,i];b=polygon[:,mod1(i+1,size(polygon,2))]
            push!(signed,(b[1]-a[1])*(point[2]-a[2])-
                         (b[2]-a[2])*(point[1]-a[1]))
        end
        (all(value->value>=-tol,signed) || all(value->value<=tol,signed)) || return false
    end
    return true
end

function prism_exterior(mesh)
    incidence=Dict{Tuple,Int}()
    # MSH Prism6: two triangular ends and three cyclic rectangular sides.
    faces=((1,2,3),(4,5,6),(1,2,5,4),(2,3,6,5),(3,1,4,6))
    for block in mesh.blocks
        block.msh==6 || continue
        for cell in eachcol(block.nodes),face in faces
            key=Tuple(sort(Int[cell[index] for index in face]))
            incidence[key]=get(incidence,key,0)+1
        end
    end
    return incidence,[mesh.coords[:,collect(face)] for (face,count) in incidence if count==1]
end

function independent_prism_volumes(mesh)
    result=Float64[]
    # Three disjoint tetrahedra in the ordinary Prism6 reference cell.
    partition=((1,2,3,4),(2,3,4,5),(3,4,5,6))
    for block in mesh.blocks
        block.msh==6 || continue
        for cell in eachcol(block.nodes)
            volume=0.0
            for local_tet in partition
                a,b,c,d=(mesh.coords[:,cell[index]] for index in local_tet)
                volume+=abs(det(hcat(b-a,c-a,d-a)))/6
            end
            push!(result,volume)
        end
    end
    return result
end

println("PUBLIC_COMPACT_BEGIN");flush(stdout)
@testset "public compact prism shifted cap coverage and orphan retention" begin
    execution=mktemp() do path,stream
        write(stream,PRISM_SOURCE);close(stream)
        Tessella.execute_geo(path;mesh_dim=3)
    end
    mesh=execution.mesh
    @test mesh isa E.MixedMesh
    @test Tessella.Model.model_entities(execution.model,3)==[(3,1)]
    @test cell_counts(mesh)==Dict(1=>27,2=>6,3=>33,6=>27)
    @test Tessella.validate(mesh).ok
    @test size(mesh.coords,2)==46 && count(used_rows(mesh))==44

    volume=Tessella.geo_entity_mesh(execution,3,1;geometry=true)
    volumes=independent_prism_volumes(volume)
    @test length(volumes)==27 && all(>(0),volumes)
    @test sum(volumes)≈0.5 atol=1e-10
    incidence,exterior=prism_exterior(volume)
    @test all(count->count in (1,2),values(incidence))

    for (surface,z) in ((1,0.0),(2,1.0))
        cap=Tessella.geo_entity_mesh(execution,2,surface;geometry=true)
        @test cell_counts(cap)==Dict(2=>3,3=>3)
        polygons=[cap.coords[:,collect(cell)] for block in cap.blocks for cell in eachcol(block.nodes)]
        triangles=[face for face in exterior if size(face,2)==3 &&
            all(value->isapprox(value,z;atol=1e-10),face[3,:])]
        @test length(triangles)==9
        @test sum(polygon_area,polygons)≈0.5 atol=1e-10
        @test all(triangle->count(polygon->contains_triangle(polygon,triangle),polygons)==1,triangles)
        covered=[sum(polygon_area(triangle) for triangle in triangles if
                     contains_triangle(polygon,triangle);init=0.0) for polygon in polygons]
        @test all(isapprox.(covered,polygon_area.(polygons);atol=1e-10))
    end

    # The public product retains actual Region-owned unused records. A
    # serialization round trip must preserve their IDs, owners and XYZ words.
    actual=execution.model.meshing.mesh_authority.mesh
    data=actual.entity_data
    @test cell_counts(actual)==Dict(15=>6,1=>27,2=>6,3=>33,6=>27)
    orphan_rows=findall(.!used_rows(actual))
    @test length(orphan_rows)==2 && all(row->data.node_entities[row]==(3,Int32(1)),orphan_rows)
    orphan_tags=data.external_node_tags[orphan_rows]
    orphan_xyz=copy(actual.coords[:,orphan_rows])
    @test allunique(data.external_node_tags) && all(isfinite,orphan_xyz)
    mktempdir() do directory
        path=joinpath(directory,"compact-prism.msh")
        E.write_mixed_msh(path,actual;version=4.1)
        recovered=E.read_mixed_msh(path)
        recovered_data=recovered.entity_data
        recovered_rows=[findfirst(==(uid),recovered_data.external_node_tags) for uid in orphan_tags]
        @test all(row->row!==nothing,recovered_rows)
        rows=Int[row for row in recovered_rows]
        @test reinterpret(UInt64,vec(recovered.coords[:,rows]))==reinterpret(UInt64,vec(orphan_xyz))
        @test recovered_data.node_entities[rows]==data.node_entities[orphan_rows]
        @test all(row->!used_rows(recovered)[row],rows)
    end

    # This existing opaque native capture is unrecombined. It supports only
    # allocation/orphan counts here, not matching recombined cell identities.
    artifact=joinpath(dirname(pathof(Tessella)),"..","test","artifacts",
                      "native_compact_transfinite_orphan_births_4152_v1.toml")
    native=only(filter(case->case["name"]=="compact_n3",TOML.parsefile(artifact)["case"]))
    @test size(actual.coords,2)==length(native["node_tags"])==46
    @test length(orphan_rows)==length(native["orphan_tags"])==2
    println("PUBLIC_COMPACT_COUNTS ",cell_counts(mesh)," nodes=",size(mesh.coords,2),
            " orphan_tags=",orphan_tags," volume=",sum(volumes));flush(stdout)
end
println("PUBLIC_COMPACT_OK");flush(stdout)

end
