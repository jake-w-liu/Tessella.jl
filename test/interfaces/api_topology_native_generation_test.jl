using Test, Tessella, TOML, SHA

module TopologyNativeGenerationCertificate
using Test, Tessella, TOML, SHA
const A=Tessella.API
const RESULT=Ref{Any}(nothing)
const ARTIFACT=joinpath(@__DIR__,"..","artifacts","api_topology_gmsh4152")
const PRIMARY=joinpath(ARTIFACT,"primary.json")
const SNAPSHOT=joinpath(ARTIFACT,"native_h1_snapshot.toml")
@assert bytes2hex(sha256(read(PRIMARY)))=="9237b93365b86c030fbb7c97707bb1422fc3b2b9e473e64c97cc540190163723"
@assert bytes2hex(sha256(read(SNAPSHOT)))=="4b187457365e816a34405786e1730aec9f2b89b8f220083480daeaf17579a828"
const EXPECTED=TOML.parsefile(SNAPSHOT)

function _add_api_topology_tetrahedron()
    for (tag,x,y,z) in ((10,0.0,0.0,0.0),(2,1.0,0.0,0.0),
                        (7,0.0,1.0,0.0),(5,0.0,0.0,1.0))
        A.model.add_point(x,y,z;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((8,10,2),(3,2,7),(11,7,10),(6,10,5),(12,2,5),(4,7,5))
        A.model.add_line(first_point,last_point;tag=tag)
    end
    for (tag,curves) in ((21,[8,3,11]),(22,[8,12,-6]),
                         (23,[3,4,-12]),(24,[11,6,-4]))
        A.model.add_curve_loop(curves;tag=tag)
        A.model.add_plane_surface([tag];tag=tag)
    end
    A.model.add_surface_loop([21,-22,23,-24];tag=30)
    A.model.add_volume([30];tag=40)
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


@testset "Gmsh4152 native h1 tetrahedron records and exact volume" begin
    A.finalize()
    A.initialize()
    try
        _add_api_topology_tetrahedron()
        generated=A.mesh.generate(3)
        actual=topology_snapshot()
        # Three fresh SDK4.15.2 h1 captures agree on every field, including
        # public node order, raw connectivity, owner parameters and factories.
        for field in keys(EXPECTED)
            @test actual[field]==EXPECTED[field]
        end
        @test Tessella.validate(generated).ok
        nodes,flat,_=A.mesh.get_nodes();xyz=Dict(tag=>Rational{BigInt}.(flat[3i-2:3i]) for (i,tag) in enumerate(nodes))
        labels,cells=A.mesh.get_elements_by_type(4,40)
        volume=zero(Rational{BigInt});faces=Dict{NTuple{3,UInt64},Int}()
        patterns=((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        for cell in eachcol(reshape(cells,4,:))
            p,q,r,s=[xyz[tag] for tag in cell];u=q-p;v=r-p;w=s-p
            determinant=u[1]*(v[2]*w[3]-v[3]*w[2])-u[2]*(v[1]*w[3]-v[3]*w[1])+u[3]*(v[1]*w[2]-v[2]*w[1])
            @test determinant>0;volume+=determinant/6
            for pattern in patterns
                key=Tuple(sort!(collect(cell[collect(pattern)])));faces[key]=get(faces,key,0)+1
            end
        end
        @test volume==1//6
        _,triangles=A.mesh.get_elements_by_type(2)
        @test Set(face for (face,count) in faces if count==1)==Set(Tuple(sort!(collect(cell))) for cell in eachcol(reshape(triangles,3,:)))
        @test all(value->value in (1,2),values(faces))

        crc=Tessella.mesh_crc(generated)
        @test crc.n_nodes==7
        @test crc.n_tets==4
        @test crc.n_boundary_faces==10
        @test Tessella.mesh_crc(A.mesh.get())==crc
        @test crc.sha=="357e66dce7c727fc08cba1f1d0e9e6970a18ffd0e9f78aee3326b227c250a094"
        # The returned facade uses a separate storage order. Certify its
        # represented coordinates independently of the public entity catalog.
        carrier_volume=zero(Rational{BigInt})
        for tet in eachcol(generated.tets)
            p,q,r,s=[Rational{BigInt}.(generated.coords[:,row]) for row in tet]
            u=q-p;v=r-p;w=s-p
            determinant=u[1]*(v[2]*w[3]-v[3]*w[2])-u[2]*(v[1]*w[3]-v[3]*w[1])+u[3]*(v[1]*w[2]-v[2]*w[1])
            @test determinant>0
            carrier_volume+=determinant/6
        end
        @test carrier_volume==1//6
        io=IOBuffer();TOML.print(io,actual;sorted=true)
        RESULT[]=(mesh_crc=crc,public_sha256=bytes2hex(sha256(take!(io))))
        println("TOPOLOGY_NATIVE_CERT mesh_crc=",crc.sha," public_crc=",RESULT[].public_sha256)
    finally
        A.finalize()
    end
end
end
