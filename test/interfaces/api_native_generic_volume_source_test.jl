module GenericVolumeNativeSourceControls
using Tessella,Test,TOML,SHA
const A=Tessella.API
const _TOPOLOGY_API=A
const CASES=Dict{String,Any}[]
# Exact SDK4.15.2 public states before setSize and two generate(3) calls.
# The extra Point record remains queryable while only its first UID is consumed.
const ARTIFACT=joinpath(@__DIR__,"..","artifacts","api_generic_volume_gmsh4152")
const PRIMARY=joinpath(ARTIFACT,"primary.json")
const STATES=joinpath(ARTIFACT,"states.toml")
@assert bytes2hex(sha256(read(PRIMARY)))=="18ee8486acccec3c470f490c8aee7d97c504c2d60e511a90ccccd9e3c1b1e635"
@assert bytes2hex(sha256(read(STATES)))=="c793f9ac27b0a22f5388373e34a70b51bba50611ccab20954e9733438687b3ab"
@assert bytes2hex(sha256(read(joinpath(ARTIFACT,"capture.py"))))=="d954dd4f0e370f9633663d99a8de570ab746d81c32581c0ffacba8e2acdabbc0"
const GENERIC_VOLUME_PRIMARY=TOML.parsefile(STATES)["cases"]
function _add_api_topology_tetrahedron()
    for (tag,x,y,z) in ((10,0.0,0.0,0.0),(2,1.0,0.0,0.0),
                        (7,0.0,1.0,0.0),(5,0.0,0.0,1.0))
        _TOPOLOGY_API.model.add_point(x,y,z;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((8,10,2),(3,2,7),(11,7,10),(6,10,5),(12,2,5),(4,7,5))
        _TOPOLOGY_API.model.add_line(first_point,last_point;tag=tag)
    end
    for (tag,curves) in ((21,[8,3,11]),(22,[8,12,-6]),
                         (23,[3,4,-12]),(24,[11,6,-4]))
        _TOPOLOGY_API.model.add_curve_loop(curves;tag=tag)
        _TOPOLOGY_API.model.add_plane_surface([tag];tag=tag)
    end
    _TOPOLOGY_API.model.add_surface_loop([21,-22,23,-24];tag=30)
    _TOPOLOGY_API.model.add_volume([30];tag=40)
    return nothing
end


_topology_bits(values)=[string(reinterpret(UInt64,Float64(value));base=16,pad=16) for value in values]
function topology_snapshot()
    entities=Dict{String,Any}[]
    for (dim,tag) in A.model.get_entities()
        ids,xyz,uv=A.mesh.get_nodes(dim,tag)
        types,labels,cells=A.mesh.get_elements(dim,tag)
        push!(entities,Dict("dimension"=>dim,"tag"=>tag,"type"=>A.model.get_type(dim,tag),
            "bounds_bits"=>_topology_bits(A.model.get_bounding_box(dim,tag)),
            "node_tags"=>collect(ids),"xyz_bits"=>_topology_bits(xyz),"param_bits"=>_topology_bits(uv),
            "types"=>collect(types),"element_tags"=>[collect(v) for v in labels],"connectivity"=>[collect(v) for v in cells]))
    end
    ids,xyz,uv=A.mesh.get_nodes();types,labels,cells=A.mesh.get_elements()
    scalar=Dict{String,Any}[]
    for tag in ids
        p,u,dim,entity=A.mesh.get_node(tag)
        push!(scalar,Dict("tag"=>tag,"xyz_bits"=>_topology_bits(p),"params_bits"=>_topology_bits(u),"owner"=>[dim,entity]))
    end
    Dict("entities"=>entities,"global_node_tags"=>collect(ids),"global_xyz_bits"=>_topology_bits(xyz),
        "global_param_bits"=>_topology_bits(uv),"global_element_types"=>collect(types),
        "global_element_tags"=>[collect(v) for v in labels],"global_connectivity"=>[collect(v) for v in cells],
        "scalar_nodes"=>scalar,"max_node_tag"=>A.mesh.get_max_node_tag(),"max_element_tag"=>A.mesh.get_max_element_tag(),
        "volume_boundary"=>[collect(pair) for pair in A.model.get_boundary([(3,40)],false,true,false)])
end
function generic_volume_certificate(seeded)
    tags,flat,_=A.mesh.get_nodes()
    xyz=Dict(tag=>Rational{BigInt}.(flat[3i-2:3i]) for (i,tag) in enumerate(tags))
    labels,cells=A.mesh.get_elements_by_type(4,40)
    faces=Dict{NTuple{3,UInt64},Int}()
    edges=Set{Tuple{UInt64,UInt64}}()
    volume=zero(Rational{BigInt})
    for cell in eachcol(reshape(cells,4,:))
        p,q,r,s=[xyz[tag] for tag in cell];u=q-p;v=r-p;w=s-p
        determinant=u[1]*(v[2]*w[3]-v[3]*w[2])-u[2]*(v[1]*w[3]-v[3]*w[1])+u[3]*(v[1]*w[2]-v[2]*w[1])
        @test determinant>0
        volume+=determinant/6
        for i in 1:4,j in i+1:4;push!(edges,minmax(cell[i],cell[j]));end
        for pattern in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
            key=Tuple(sort!(collect(cell[collect(pattern)])))
            faces[key]=get(faces,key,0)+1
        end
    end
    @test length(labels)==4
    @test volume==(seeded ? 7//48 : 1//6)
    _,triangles=A.mesh.get_elements_by_type(2)
    @test Set(face for (face,count) in faces if count==1)==Set(Tuple(sort!(collect(cell))) for cell in eachcol(reshape(triangles,3,:)))
    @test all(value->value in (1,2),values(faces))
    _,lines=A.mesh.get_elements_by_type(1)
    @test all(line->minmax(line...) in edges,eachcol(reshape(lines,2,:)))
    if seeded
        point_tags,point_xyz,_=A.mesh.get_nodes(0,10)
        @test point_xyz[1:3]==[.125,0.,0.]
        @test first(point_tags) in cells
        if length(point_tags)==2
            @test point_xyz[4:6]==[.2,0.,0.]
            @test !(last(point_tags) in cells)
            @test A.mesh.get_elements(0,10)[2]==[UInt64[401,402]]
        end
    end
    return Dict("volume"=>string(volume),"tet_count"=>length(labels))
end
function generic_volume_source_control(reference)
    A.initialize()
    try
        A.option("Mesh.MeshOnlyEmpty",reference["onlyempty"])
        _add_api_topology_tetrahedron()
        seeded=reference["seeded"]
        if seeded
            extra=reference["extra_Point_node"]
            tags=extra ? [101,102] : [101]
            xyz=extra ? [.125,0.,0.,.2,0.,0.] : [.125,0.,0.]
            labels=extra ? [401,402] : [401]
            A.mesh.add_nodes(0,10,tags,xyz)
            A.mesh.add_elements_by_type(10,15,labels,tags)
            A.mesh.generate(1)
        end
        A.mesh.generate(3)
        before=topology_snapshot()
        for field in keys(reference["before_set_size"])
            @test before[field]==reference["before_set_size"][field]
        end
        generic_volume_certificate(seeded)
        A.mesh.set_size([(0,10),(0,2),(0,7),(0,5)],.35)
        @test topology_snapshot()==before
        entry=Dict{String,Any}("name"=>"generic_Volume_Point_source","seeded"=>seeded,
            "extra"=>reference["extra_Point_node"],"onlyempty"=>reference["onlyempty"],"before"=>before,
            "regenerated"=>Dict{String,Any}[])
        push!(CASES,entry)
        for iteration in 1:2
            A.mesh.generate(3)
            actual=topology_snapshot()
            expected=reference["regenerated"][iteration]["state"]
            for field in keys(expected)
                @test actual[field]==expected[field]
            end
            science=generic_volume_certificate(seeded)
            push!(entry["regenerated"],Dict("actual"=>actual,"science"=>science))
        end
        println("GENERIC_VOLUME_NATIVE seeded=",seeded," extra=",reference["extra_Point_node"],
            " OE=",reference["onlyempty"]," max=",A.mesh.get_max_node_tag(),"/",A.mesh.get_max_element_tag())
    finally
        A.finalize()
    end
end

@testset "Six generic Volume setter and actual Point consumers" begin
    @test length(GENERIC_VOLUME_PRIMARY)==6
    for reference in GENERIC_VOLUME_PRIMARY
        @testset "Volume seed$(reference["seeded"]) extra$(reference["extra_Point_node"]) OE$(reference["onlyempty"])" begin
            generic_volume_source_control(reference)
        end
    end
end
empty!(CASES)
A.finalize()
end
