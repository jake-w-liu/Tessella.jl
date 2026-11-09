module OCCImportBREPNaturalChecks
using Test,Tessella
using LinearAlgebra: cross,det,dot,norm
const A=Tessella.API
const M=Tessella.Model
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function session(f)
    A.initialize()
    try
        return f(A.model.occ)
    finally
        A.finalize()
    end
end

function surface_snapshot(tag=1)
    tags,xyz,parameters=A.mesh.get_nodes(2,tag)
    coordinates=reshape(xyz,3,:)
    rows=Dict(uid=>row for (row,uid) in enumerate(tags))
    families,_,groups=A.mesh.get_elements(2,tag)
    connectivity=reduce(vcat,groups;init=UInt64[])
    triangles=reshape(Int32[rows[uid] for uid in connectivity],3,:)
    incidence=Dict{Tuple{UInt64,UInt64},Int}()
    for i in 1:3:length(connectivity),pair in ((1,2),(2,3),(3,1))
        edge=minmax(connectivity[i+pair[1]-1],connectivity[i+pair[2]-1])
        incidence[edge]=get(incidence,edge,0)+1
    end
    surface=Tessella.Mesh(coordinates;tris=triangles)
    maximum_edge=maximum(norm(coordinates[:,rows[a]].-coordinates[:,rows[b]])
        for (a,b) in keys(incidence))
    return (;tags,xyz,parameters,coordinates,rows,families,connectivity,
        incidence,surface,maximum_edge)
end

@testset "BREP finite natural faces retain their actual Surface topology" begin
    # Public CAD writer fixtures, with the Face's wire list omitted. The
    # reference reader imports exactly one Surface and no lower entities.
    cases=(("sphere",4π,2,true),("torus",π^2,0,true),
        ("bspline",1.0158281995745169,1,false),("plane",1.0,1,false))
    for (kind,area,euler,closed) in cases,
            algorithm in (kind=="plane" ? (6,1,7) : (6,1))
        session() do occ
            out=occ.import_shapes(joinpath(FIXTURES,"brep_native_natural_$kind.brep"))
            occ.synchronize()
            @test out==[(2,1)]
            @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(0,0,1,0)
            @test occ.get_mass(2,1)≈area rtol=1e-10
            @test all(isfinite,A.model.get_bounding_box(2,1))
            A.option("Mesh.MeshSizeMin",1.)
            A.option("Mesh.MeshSizeMax",1.)
            A.option("Mesh.Algorithm",algorithm)
            A.mesh.generate(2)
            snapshot=surface_snapshot()
            @test snapshot.families==Int32[2]
            @test !isempty(snapshot.tags)
            @test length(snapshot.parameters)==2length(snapshot.tags)
            @test A.model.get_value(2,1,snapshot.parameters)≈snapshot.xyz atol=2e-12
            data=A.CURRENT[].meshing.mesh_authority.mesh.entity_data
            @test all(==((2,Int32(1))),data.node_entities)
            @test Tessella.validate(snapshot.surface).ok
            @test M._plc_self_intersection3(snapshot.surface)===nothing
            @test length(snapshot.tags)-length(snapshot.incidence)+
                size(snapshot.surface.tris,2)==euler
            if closed
                @test all(==(2),values(snapshot.incidence))
            else
                @test any(==(1),values(snapshot.incidence)) &&
                    all(count->count in (1,2),values(snapshot.incidence))
                degree=Dict{UInt64,Int}()
                for (edge,count) in snapshot.incidence
                    count==1 || continue
                    for uid in edge;degree[uid]=get(degree,uid,0)+1;end
                end
                @test all(==(2),values(degree))
                @test_throws ArgumentError M._validate_surface_loop(
                    A.CURRENT[],[1],"natural open Surface")
            end
            algorithm in (1,7) && @test snapshot.maximum_edge<=1+1e-12
        end
    end
end

@testset "Natural chart shifts retain the physical constant size field" begin
    session() do occ
        occ.import_shapes(joinpath(FIXTURES,"brep_native_natural_torus_shifted.brep"))
        occ.synchronize()
        lo,hi=A.model.get_parametrization_bounds(2,1)
        @test collect(lo)≈[π/3,0.]
        @test collect(hi)≈[π/3+2π,2π]
        @test occ.get_mass(2,1)≈π^2 rtol=1e-12
        A.option("Mesh.MeshSizeMin",1.)
        A.option("Mesh.MeshSizeMax",1.)
        A.option("Mesh.Algorithm",1)
        A.mesh.generate(2)
        snapshot=surface_snapshot()
        @test snapshot.maximum_edge<=1+1e-12
        @test snapshot.maximum_edge>.7
        @test length(snapshot.tags)-length(snapshot.incidence)+
            size(snapshot.surface.tris,2)==0
        @test all(==(2),values(snapshot.incidence))
        @test M._plc_self_intersection3(snapshot.surface)===nothing
        @test A.model.get_value(2,1,snapshot.parameters)≈snapshot.xyz atol=2e-12
    end
end

@testset "Closed natural BREP solids consume their actual published Face mesh" begin
    for (kind,cad_volume,algorithm,element_size) in (("sphere",4π/3,6,2.),("torus",π^2/8,1,1.))
        session() do occ
            mktempdir() do directory
                source=read(joinpath(FIXTURES,"brep_native_natural_$kind.brep"),String)
                # The generated table ends Face→Shell→Solid→Compound.
                # Select the retained Solid root instead of the isolated Face.
                source=replace(source,r"(?m)^\+4 0[ \t]*\r?$"=>"+2 0")
                path=joinpath(directory,"natural_solid.brep")
                write(path,source)
                @test occ.import_shapes(path)==[(3,1)]
                occ.synchronize()
                @test isempty(A.model.get_entities(0)) && isempty(A.model.get_entities(1))
                @test occ.get_mass(3,1)≈cad_volume rtol=1e-10
                @test all(x->abs(x)<1e-10,occ.get_center_of_mass(3,1))
                I=occ.get_matrix_of_inertia(3,1)
                if kind=="sphere"
                    @test I≈[8π/15,0.,0., 0.,8π/15,0., 0.,0.,8π/15] atol=1e-8
                else
                    @test I≈[0.6650417,0.,0., 0.,0.6650417,0.,
                        0.,0.,1.29153026] atol=1e-6
                end
                @test all(isfinite,A.model.get_bounding_box(3,1))
                A.option("Mesh.MeshSizeMin",element_size)
                A.option("Mesh.MeshSizeMax",element_size)
                A.option("Mesh.Algorithm",algorithm)
                A.mesh.generate(3)
                tags,xyz,_=A.mesh.get_nodes()
                coordinates=reshape(xyz,3,:)
                rows=Dict(uid=>row for (row,uid) in enumerate(tags))
                families,_,groups=A.mesh.get_elements(3)
                @test families==Int32[4]
                incidence=Dict{NTuple{3,UInt64},Int}()
                volumes=Float64[]
                for cells in groups,i in 1:4:length(cells)
                    tet=cells[i:i+3];p=coordinates[:,rows[tet[1]]]
                    push!(volumes,det(hcat((coordinates[:,rows[tet[j]]].-p for j in 2:4)...))/6)
                    for face in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
                        key=Tuple(sort(tet[collect(face)]))
                        incidence[key]=get(incidence,key,0)+1
                    end
                end
                @test !isempty(volumes) && all(>(0),volumes)
                @test 0<sum(volumes)<cad_volume
                _,_,groups=A.mesh.get_elements(2)
                facets=Set(Tuple(sort(cells[i:i+2]))
                    for cells in groups for i in 1:3:length(cells))
                @test facets==Set(face for (face,count) in incidence if count==1)
                surface_volume=sum(dot(coordinates[:,rows[cells[i]]],
                    cross(coordinates[:,rows[cells[i+1]]],coordinates[:,rows[cells[i+2]]]))/6
                    for cells in groups for i in 1:3:length(cells))
                @test sum(volumes)≈surface_volume rtol=2e-12
            end
        end
    end
end

@testset "Natural BREP copy and transform keep finite domains and held state" begin
    session() do occ
        occ.import_shapes(joinpath(FIXTURES,"brep_native_natural_sphere.brep"))
        occ.synchronize()
        uv=[.4,.2]
        original=A.model.get_value(2,1,uv)
        @test occ.copy([(2,1)])==[(2,2)]
        occ.translate([(2,2)],3.,2.,1.)
        occ.dilate([(2,2)],3.,2.,1.,2.,2.,2.)
        occ.synchronize()
        @test A.model.get_value(2,1,uv)≈original
        @test A.model.get_value(2,2,uv)≈2original.+[3.,2.,1.]
        @test occ.get_mass(2,2)≈16π rtol=1e-12
        @test collect(A.model.get_bounding_box(2,2))≈[1.,0.,-1.,5.,4.,3.]
        A.option("Mesh.MeshSizeMin",2.)
        A.option("Mesh.MeshSizeMax",2.)
        A.option("Mesh.Algorithm",6)
        A.mesh.generate(2)
        held=deepcopy(A.mesh.get_nodes())
        geometry=deepcopy(A.CURRENT[].surface_geometry)
        @test_throws ArgumentError occ.dilate([(2,2)],0.,0.,0.,2.,1.,1.)
        @test A.CURRENT[].surface_geometry==geometry
        @test A.mesh.get_nodes()==held
        @test A.model.get_value(2,1,uv)≈original
    end
end

@testset "Natural tensor reflection and bounded plane stretching retain source parameters" begin
    session() do occ
        occ.import_shapes(joinpath(FIXTURES,"brep_native_natural_bspline.brep"))
        occ.synchronize()
        uv=[.3,.7]
        original=A.model.get_value(2,1,uv)
        occ.mirror([(2,1)],1.,0.,0.,0.)
        occ.synchronize()
        @test A.model.get_value(2,1,uv)≈[-original[1],original[2],original[3]]
        @test occ.get_mass(2,1)≈1.0158281995745169 rtol=1e-10
        @test collect(A.model.get_bounding_box(2,1))≈
            [-1.0-1e-7,-1e-7,-1e-7,1e-7,1.0+1e-7,.3+1e-7]
        A.option("Mesh.MeshSizeMin",1.)
        A.option("Mesh.MeshSizeMax",1.)
        A.option("Mesh.Algorithm",1)
        A.mesh.generate(2)
        snapshot=surface_snapshot()
        @test A.model.get_value(2,1,snapshot.parameters)≈snapshot.xyz atol=2e-12
        @test snapshot.maximum_edge<=1+1e-12
    end
    session() do occ
        occ.import_shapes(joinpath(FIXTURES,"brep_native_natural_plane.brep"))
        occ.synchronize()
        occ.dilate([(2,1)],0.,0.,0.,2.,3.,4.)
        occ.synchronize()
        @test A.model.get_value(2,1,[.2,.3])≈[.4,.9,0.]
        @test occ.get_mass(2,1)≈6.
        @test collect(A.model.get_bounding_box(2,1))≈[0.,0.,0.,2.,3.,0.]
        A.option("Mesh.MeshSizeMin",1.)
        A.option("Mesh.MeshSizeMax",1.)
        A.option("Mesh.Algorithm",1)
        A.mesh.generate(2)
        snapshot=surface_snapshot()
        @test snapshot.maximum_edge<=1+1e-12
        @test A.model.get_value(2,1,snapshot.parameters)≈snapshot.xyz atol=2e-12
    end
end
end
