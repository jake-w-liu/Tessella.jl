# LAY01 bounded slice: two adjacent QuadTriNoNewVerts sweeps that AutoCoherence
# welded onto ONE literal lateral Surface. The shared Surface's corner/column
# incidence binds once in the operation scope and transports into both Region
# plans; the joint certificate checks opposite child orientation, matching
# wall diagonals, positive typed volumes and the closed exterior boundary.
# Everything outside the literal two-Region equal-translation pair keeps its
# planner diagnostic. The certified pair also classifies through the API
# projection: the bound Surface part's `model_to_mixed` metadata is valid for
# both owning Regions, so `execute_geo`/`mesh.generate(3)` emit it once.
using Test
using Tessella
using Tessella.Elements: MixedMesh, ElementBlock, _VOLUME_CELL_FACES, msh_dimension
using Tessella.Model: mesh_model_volume, mesh_model_surface,
    _extrude_cell_signed_volume, _extrude_quadtri_regions, _extrude_volume_pass,
    model_to_mixed
using Tessella.MeshTypes: Mesh, nnodes, validate

const _SL_CALLER="extrude_nonew_shared_lateral_test"

# Two transfinite unit squares welded onto the literal edge Curve(2): the
# 1×1 grids stay on the generic single-quad catalog, and the two identical
# translations make each lateral coincide so AutoCoherence merges Curve(2)'s
# products into one literal Surface.
function _sl_squares()
    return """
    Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
    Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
    Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
    Point(5)={2,0,0,1};Point(6)={2,1,0,1};
    Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
    Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
    Transfinite Curve{:}=2;
    Transfinite Surface{1};Transfinite Surface{2};
    Recombine Surface{1,2};
    """
end

function _sl_model(source;mesh_dim=0)
    mktempdir() do directory
        path=joinpath(directory,"nonew_shared.geo")
        write(path,source)
        return execute_geo(path;mesh_dim)
    end
end

_sl_blocks(mesh::Mesh)=[ElementBlock(4,mesh.tets)]
_sl_blocks(mesh::MixedMesh)=mesh.blocks

# Position-keyed facet sets: each volume cell face reduced to its sorted
# vertex-position tuple, counting incidences across the union.
function _sl_face_counts(meshes...)
    counts=Dict{Tuple,Int}()
    for mesh in meshes
        for block in _sl_blocks(mesh)
            block isa ElementBlock && haskey(_VOLUME_CELL_FACES,Int(block.msh)) ||
                continue
            for cell in axes(block.nodes,2),face in _VOLUME_CELL_FACES[Int(block.msh)]
                vertices=ntuple(k->Tuple(mesh.coords[:,block.nodes[face[k],cell]]),
                                length(face))
                key=Tuple(sort!(collect(vertices)))
                counts[key]=get(counts,key,0)+1
            end
        end
    end
    return counts
end

function _sl_surface_faces(part)
    faces=Set{Tuple}()
    blocks=part isa Mesh ? [ElementBlock(2,part.tris)] : part.blocks
    for block in blocks
        block isa ElementBlock && block.msh in (2,3) || continue
        for cell in axes(block.nodes,2)
            push!(faces,Tuple(sort!([Tuple(part.coords[:,block.nodes[k,cell]])
                                     for k in axes(block.nodes,1)])))
        end
    end
    return faces
end

function _sl_part_equal(a,b)
    a isa Mesh && b isa Mesh &&
        return a.coords==b.coords && a.segs==b.segs && a.tris==b.tris &&
            a.tets==b.tets && a.seg_tag==b.seg_tag && a.tri_tag==b.tri_tag &&
            a.tet_tag==b.tet_tag
    a isa MixedMesh && b isa MixedMesh || return false
    a.coords==b.coords && length(a.blocks)==length(b.blocks) || return false
    for (x,y) in zip(a.blocks,b.blocks)
        x isa ElementBlock && y isa ElementBlock &&
            x.msh==y.msh && x.nodes==y.nodes && x.tags==y.tags || return false
    end
    return true
end

# Outward-oriented face cycles, mirroring the incidence certificate's
# per-type tables — `_VOLUME_CELL_FACES` records vertex incidence only and
# does not orient all faces outwards.
const _SL_OUTWARD_FACES=Dict(
    4=>((1,3,2),(1,2,4),(2,3,4),(3,1,4)),
    5=>((1,4,3,2),(5,6,7,8),(1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8)),
    6=>((1,3,2),(4,5,6),(1,2,5,4),(2,3,6,5),(3,1,4,6)),
    7=>((1,4,3,2),(1,2,5),(2,3,5),(3,4,5),(4,1,5)))

# The outward cycle a volume cell emits on one shared wall facet, so
# opposite child orientation can be checked without a second registry.
function _sl_wall_cycles(mesh,x)
    cycles=Dict{Tuple,NTuple}()
    for block in _sl_blocks(mesh)
        block isa ElementBlock && haskey(_SL_OUTWARD_FACES,Int(block.msh)) ||
            continue
        for cell in axes(block.nodes,2),face in _SL_OUTWARD_FACES[Int(block.msh)]
            vertices=ntuple(k->Tuple(mesh.coords[:,block.nodes[face[k],cell]]),
                            length(face))
            all(v->v[1]==x,vertices) || continue
            cycles[Tuple(sort!(collect(vertices)))]=vertices
        end
    end
    return cycles
end

function _sl_reversed(a::NTuple{N,Tuple},b::NTuple{N,Tuple}) where N
    reversed=collect(reverse(a))
    target=collect(b)
    for shift in 0:N-1
        circshift(reversed,shift)==target && return true
    end
    return false
end

function _sl_positive(mesh)
    for block in _sl_blocks(mesh)
        block isa ElementBlock && msh_dimension(block.msh)==3 || continue
        for cell in axes(block.nodes,2)
            @test _extrude_cell_signed_volume(mesh.coords,block.nodes,cell,
                                              block.msh)>0
        end
    end
end

function _sl_shared_tag(execution)
    left=Int.(execution.lists["left"])
    right=Int.(execution.lists["right"])
    shared=intersect(left[3:end],right[3:end])
    @test length(shared)==1
    return only(shared)
end

# All dim-2 cells of `mesh` lying entirely on the `axis == value` plane —
# used on the merged product, whose blocks bucket cells by entity.
function _sl_dim2_faces_on(mesh,axis::Int,value::Float64)
    faces=Set{Tuple}()
    blocks=mesh isa Mesh ? [ElementBlock(2,mesh.tris)] : mesh.blocks
    for block in blocks
        block isa ElementBlock && block.msh in (2,3) || continue
        for cell in axes(block.nodes,2)
            vertices=ntuple(k->Tuple(mesh.coords[:,block.nodes[k,cell]]),
                            size(block.nodes,1))
            all(v->v[axis]==value,vertices) &&
                push!(faces,Tuple(sort!(collect(vertices))))
        end
    end
    return faces
end

# A plain untagged view of a bound part's dim-2 cells — the input contract a
# standalone `model_to_mixed` projection expects.
function _sl_plain_surface(part)
    cells=NTuple{3,Int32}[]
    for block in (part isa Mesh ? [ElementBlock(2,part.tris)] : part.blocks)
        block isa ElementBlock && block.msh==2 || continue
        for cell in axes(block.nodes,2)
            push!(cells,ntuple(k->block.nodes[k,cell],3))
        end
    end
    nodes=Matrix{Int32}(undef,3,length(cells))
    for (i,cell) in pairs(cells),k in 1:3
        nodes[k,i]=cell[k]
    end
    return Mesh(part.coords;tris=nodes)
end

@testset "NoNew literal shared lateral" begin
    @testset "bound pair emits one certified wall" for laterals in (false,true)
        flag=laterals ? " RecombLaterals" : ""
        execution=_sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts$flag;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts$flag;};
        """)
        model=execution.model
        s=_sl_shared_tag(execution)
        # One literal Surface, owned by both Regions.
        @test sort(_extrude_quadtri_regions(model,s,_SL_CALLER))==[1,2]
        pass=_extrude_volume_pass(model,_SL_CALLER)
        scope=pass.scope
        @test scope!==nothing
        left=pass.parts[1];right=pass.parts[2]
        bound=scope.surfaces[s]
        @test validate(left).ok && validate(right).ok && validate(bound).ok
        # One bound product: both Region plans emitted the same part and the
        # scope stores it once; the standalone pair lookup returns it too.
        @test _sl_part_equal(scope.volumes[1].surfaces[s],bound)
        @test _sl_part_equal(scope.volumes[2].surfaces[s],bound)
        standalone=mesh_model_surface(model,s)
        @test _sl_part_equal(standalone,bound)
        # The shared wall incidences each Region emits equal the bound
        # Surface's cells — the transported corner/column contract.
        wall=_sl_surface_faces(bound)
        @test wall==Set(keys(_sl_wall_cycles(left,1.0)))==
            Set(keys(_sl_wall_cycles(right,1.0)))
        # Opposite child incidences: every shared facet arrives once per
        # Region with reversed outward cycles.
        cycles_left=_sl_wall_cycles(left,1.0)
        cycles_right=_sl_wall_cycles(right,1.0)
        @test keys(cycles_left)==keys(cycles_right)
        for key in keys(cycles_left)
            @test _sl_reversed(cycles_left[key],cycles_right[key])
        end
        # Positive typed volume cells in both Regions.
        _sl_positive(left);_sl_positive(right)
        # Each standalone Region plan also emits a certified per-Region
        # product on the same model.
        @test validate(mesh_model_volume(model,1)).ok
        @test validate(mesh_model_volume(model,2)).ok
        # Joint boundary: the shared wall welds (incidence two); the rest of
        # both Regions' surfaces is exactly the external boundary.
        counts=_sl_face_counts(left,right)
        @test all(c in (1,2) for c in values(counts))
        external=Set(k for (k,c) in counts if c==1)
        interior=Set(k for (k,c) in counts if c==2)
        interior_wall=Set(k for k in interior if all(v->v[1]==1.0,k))
        @test wall ⊆ interior && interior_wall==wall
        exterior=Set{Tuple}()
        for volume in (1,2),shell in model.volumes[volume],
            signed in model.surface_loops[abs(shell)]
            tag=abs(signed)
            tag==s && continue
            union!(exterior,_sl_surface_faces(scope.surfaces[tag]))
        end
        @test exterior==external
    end

    @testset "shared laterals reject outside the literal pair" begin
        # Different layer records weld the same literal Surface but break the
        # identical-operation contract.
        execution=_sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{3};Recombine;QuadTriNoNewVerts;};
        """)
        s=_sl_shared_tag(execution)
        @test_throws r"shared-region planner" mesh_model_volume(execution.model,1)
        @test_throws r"shared-region planner" mesh_model_volume(execution.model,2)
        @test_throws r"shared-region planner" mesh_model_surface(execution.model,s)
        # A RecombLaterals asymmetry is still a different operation record.
        execution=_sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts RecombLaterals;};
        """)
        s=_sl_shared_tag(execution)
        @test_throws r"shared-region planner" mesh_model_volume(execution.model,2)
        # A non-NoNew neighbor keeps the original loud diagnostic.
        execution=_sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriAddVerts;};
        """)
        s=_sl_shared_tag(execution)
        @test_throws r"shared-region planner" mesh_model_volume(execution.model,1)
        @test_throws r"shared-region planner" mesh_model_surface(execution.model,s)
        # Boundary overrides attached to the shared Surface keep the
        # constrained-boundary diagnostic before any shared admission.
        execution=_sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """)
        s=_sl_shared_tag(execution)
        execution.model.meshing.transfinite_surfaces[s]=
            (arrangement=:left,corners=Int[])
        @test_throws r"constrained-boundary planner" mesh_model_volume(execution.model,1)
    end

    @testset "certified pair classifies through the API projection" begin
        # execute_geo(mesh_dim=3) routes through the shared classify path the
        # API generate uses — the bound lateral Surface must classify once
        # and merge into one welded wall.
        execution=_sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        s=_sl_shared_tag(execution)
        merged=execution.mesh
        @test validate(merged).ok
        # Each merged vertex has one coordinate — the wall nodes are shared
        # bitwise between the Surface part and both Region boundaries.
        @test allunique(Tuple(merged.coords[:,node]) for node in 1:nnodes(merged))
        part=geo_entity_mesh(execution,2,s)
        wall=_sl_surface_faces(part)
        @test !isempty(wall)
        # The shared wall's dim-2 cells are emitted exactly once, on x=1.
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
        # Both Regions' lateral boundaries reference the identical wall
        # incidences: each shared facet has joint incidence two.
        volume_faces=_sl_face_counts(geo_entity_mesh(execution,3,1),
                                     geo_entity_mesh(execution,3,2))
        @test all(count->count in (1,2),values(volume_faces))
        interior=Set(key for (key,count) in volume_faces if count==2)
        @test wall ⊆ interior
        @test Set(key for key in interior if all(v->v[1]==1.0,key))==wall
        # The certified wall part classifies as Surface[s] plus its boundary
        # closure — the same entity metadata the isolated planner produces.
        fresh=_sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """)
        model=fresh.model
        fs=_sl_shared_tag(fresh)
        scope=_extrude_volume_pass(model,_SL_CALLER).scope
        bound=scope.surfaces[fs]
        closure=Set{Tuple{Int,Int}}(((2,fs),))
        for loop in model.surfaces[fs],signed in model.loops[loop]
            curve=abs(signed)
            push!(closure,(1,curve))
            for point in model.curves[curve]
                push!(closure,(0,point))
            end
        end
        for owner_scope in (nothing,scope)
            projected=model_to_mixed(model,bound,fs;_extrude_scope=owner_scope)
            @test validate(projected).ok
            @test projected.coords==bound.coords
            @test _sl_surface_faces(projected)==_sl_surface_faces(bound)
            data=projected.entity_data
            @test data!==nothing
            @test all(owner->owner in closure,data.node_entities)
            @test issubset(Set(keys(data.entities)),closure)
        end
        # A scope that never planned this operation stays a loud refusal.
        empty_scope=Tessella.Model._ExtrudeNoNewScope(model,
            Dict{Int,Tessella.Model._ExtrudeNoNewCompletePlan}(),
            Dict{Int,Union{Mesh,MixedMesh}}())
        @test_throws r"absent from this mesh scope" model_to_mixed(
            model,bound,fs;_extrude_scope=empty_scope)
        # Uncertified multi-owner incidence still rejects at API level — the
        # planner diagnostic surfaces both from the projection and from a
        # full execute_geo(mesh_dim=3) request.
        bad=_sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{3};Recombine;QuadTriNoNewVerts;};
        """)
        bs=_sl_shared_tag(bad)
        @test_throws r"shared-region planner" model_to_mixed(
            bad.model,_sl_plain_surface(bound),bs)
        @test_throws r"shared-region planner" _sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{3};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
    end
end
