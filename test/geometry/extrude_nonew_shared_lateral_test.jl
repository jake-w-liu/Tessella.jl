# LAY01 bounded slice: adjacent QuadTriNoNewVerts sweeps that AutoCoherence
# welded onto literal lateral Surfaces. Each shared Surface's corner/column
# incidence binds once in the operation scope and transports into both Region
# plans; the joint certificate checks opposite child orientation, matching
# wall diagonals, positive typed volumes and the closed exterior boundary.
# A Region may bound several such walls — chains and corners of paired sweeps
# plan the whole connected component, carrying one transported diagonal
# contract per earlier partner. Everything outside the literal two-Region
# equal-translation pair keeps its planner diagnostic. The certified pairs
# also classify through the API projection: a bound Surface part's
# `model_to_mixed` metadata is valid for both owning Regions, so
# `execute_geo`/`mesh.generate(3)` emit it once.
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

# The same pair with a triangle-meshed right Source — the certified
# unrecombined-partner wall. The oracle (Gmsh 4.15.2) welds the shared
# lateral onto the partner's prism-subdivision triangles.
function _sl_squares_tri()
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
    Recombine Surface{1};
    """
end

# Three transfinite unit squares in a row — the middle sweep welds its left
# lateral onto the left pair's wall and its right lateral onto the right
# pair's wall, so one Region bounds two literal shared Surfaces.
function _sl_chain_squares()
    return """
    Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
    Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
    Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
    Point(5)={2,0,0,1};Point(6)={2,1,0,1};
    Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
    Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
    Point(7)={3,0,0,1};Point(8)={3,1,0,1};
    Line(8)={5,7};Line(9)={7,8};Line(10)={8,6};
    Curve Loop(3)={-6,8,9,10};Plane Surface(3)={3};
    Transfinite Curve{:}=2;
    Transfinite Surface{:};
    Recombine Surface{:};
    """
end

# Three transfinite unit squares in an L — the corner sweep bounds shared
# walls with both neighbors on adjacent source edges.
function _sl_lshape_squares()
    return """
    Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
    Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
    Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
    Point(5)={2,0,0,1};Point(6)={2,1,0,1};
    Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
    Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
    Point(7)={1,2,0,1};Point(8)={0,2,0,1};
    Line(8)={3,7};Line(9)={7,8};Line(10)={8,4};
    Curve Loop(3)={-3,8,9,10};Plane Surface(3)={3};
    Transfinite Curve{:}=2;
    Transfinite Surface{:};
    Recombine Surface{:};
    """
end

# Four transfinite unit squares in a 2×2 arrangement — the last sweep carries
# two transported diagonal contracts on the same plan.
function _sl_grid22_squares()
    return """
    Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
    Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
    Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
    Point(5)={2,0,0,1};Point(6)={2,1,0,1};
    Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
    Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
    Point(7)={1,2,0,1};Point(8)={0,2,0,1};
    Line(8)={3,7};Line(9)={7,8};Line(10)={8,4};
    Curve Loop(3)={-3,8,9,10};Plane Surface(3)={3};
    Point(9)={2,2,0,1};
    Line(11)={7,9};Line(12)={9,6};
    Curve Loop(4)={-7,-12,-11,-8};Plane Surface(4)={4};
    Transfinite Curve{:}=2;
    Transfinite Surface{:};
    Recombine Surface{:};
    """
end

# The literal Surface two named extrusion return-lists welded together.
function _sl_pair_tag(execution,left::AbstractString,right::AbstractString)
    shared=intersect(Int.(execution.lists[left])[3:end],
                     Int.(execution.lists[right])[3:end])
    @test length(shared)==1
    return only(shared)
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

# Element-type multiset over all blocks of a completed mesh.
_sl_msh_counts(mesh::Mesh)=Dict{Int,Int}(
    k=>v for (k,v) in (
        1=>size(mesh.segs,2),2=>size(mesh.tris,2),
        4=>size(mesh.tets,2)) if v>0)
function _sl_msh_counts(mesh::MixedMesh)
    counts=Dict{Int,Int}()
    for block in _sl_blocks(mesh)
        block isa ElementBlock && msh_dimension(block.msh)>=1 ||
            continue
        counts[Int(block.msh)]=
            get(counts,Int(block.msh),0)+size(block.nodes,2)
    end
    return counts
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

    @testset "three-Region chain welds both shared laterals" for laterals in (false,true)
        flag=laterals ? " RecombLaterals" : ""
        execution=_sl_model(_sl_chain_squares()*"""
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts$flag;};
        b[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts$flag;};
        c[]=Extrude{0,0,1}{Surface{3};Layers{2};Recombine;QuadTriNoNewVerts$flag;};
        """;mesh_dim=3)
        model=execution.model
        s12=_sl_pair_tag(execution,"a","b")
        s23=_sl_pair_tag(execution,"b","c")
        @test s12!=s23
        # Each wall keeps its literal two-Region pair; the middle Region
        # bounds both.
        @test sort(_extrude_quadtri_regions(model,s12,_SL_CALLER))==[1,2]
        @test sort(_extrude_quadtri_regions(model,s23,_SL_CALLER))==[2,3]
        merged=execution.mesh
        @test validate(merged).ok
        @test allunique(Tuple(merged.coords[:,node]) for node in 1:nnodes(merged))
        pass=_extrude_volume_pass(model,_SL_CALLER)
        scope=pass.scope
        @test scope!==nothing
        for (s,x,lo,hi) in ((s12,1.0,1,2),(s23,2.0,2,3))
            bound=scope.surfaces[s]
            @test validate(bound).ok
            wall=_sl_surface_faces(bound)
            # One bound product per wall: both Region plans emitted the same
            # part and the scope stores it once.
            @test _sl_part_equal(scope.volumes[lo].surfaces[s],bound)
            @test _sl_part_equal(scope.volumes[hi].surfaces[s],bound)
            # The standalone surface path replays the connected component and
            # returns the identical emission.
            @test _sl_part_equal(mesh_model_surface(model,s),bound)
            # The merged mesh buckets each shared wall under its own Surface
            # exactly once, on its plane.
            @test _sl_surface_faces(geo_entity_mesh(execution,2,s))==wall
            @test _sl_dim2_faces_on(merged,1,x)==wall
            # Opposite child incidences between the wall's two owners.
            cycles_lo=_sl_wall_cycles(pass.parts[lo],x)
            cycles_hi=_sl_wall_cycles(pass.parts[hi],x)
            @test Set(keys(cycles_lo))==Set(keys(cycles_hi))==wall
            for key in keys(cycles_lo)
                @test _sl_reversed(cycles_lo[key],cycles_hi[key])
            end
        end
        for tag in (1,2,3)
            _sl_positive(pass.parts[tag])
        end
        # Joint boundary across all three Regions: each wall welds exactly on
        # its plane (incidence two); the remainder is exactly the exterior
        # boundary.
        counts=_sl_face_counts(pass.parts[1],pass.parts[2],pass.parts[3])
        @test all(c in (1,2) for c in values(counts))
        interior=Set(k for (k,c) in counts if c==2)
        @test Set(k for k in interior if all(v->v[1]==1.0,k))==
            _sl_surface_faces(scope.surfaces[s12])
        @test Set(k for k in interior if all(v->v[1]==2.0,k))==
            _sl_surface_faces(scope.surfaces[s23])
        external=Set(k for (k,c) in counts if c==1)
        exterior=Set{Tuple}()
        for volume in (1,2,3),shell in model.volumes[volume],
            signed in model.surface_loops[abs(shell)]
            tag=abs(signed)
            tag in (s12,s23) && continue
            union!(exterior,_sl_surface_faces(scope.surfaces[tag]))
        end
        @test exterior==external
    end

    @testset "corner Region claims two shared laterals" begin
        # The L-arrangement's corner Region plans first and binds both walls;
        # each neighbor transports its own diagonal contract.
        execution=_sl_model(_sl_lshape_squares()*"""
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        b[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        c[]=Extrude{0,0,1}{Surface{3};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        model=execution.model
        @test validate(execution.mesh).ok
        sx=_sl_pair_tag(execution,"a","b")
        sy=_sl_pair_tag(execution,"a","c")
        scope=_extrude_volume_pass(model,_SL_CALLER).scope
        @test scope!==nothing
        for s in (sx,sy)
            @test _sl_part_equal(mesh_model_surface(model,s),scope.surfaces[s])
        end
        # Two contracts on one plan: the 2×2 corner Region must satisfy the
        # diagonal transported through each earlier partner.
        execution=_sl_model(_sl_grid22_squares()*"""
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        b[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        c[]=Extrude{0,0,1}{Surface{3};Layers{2};Recombine;QuadTriNoNewVerts;};
        d[]=Extrude{0,0,1}{Surface{4};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        model=execution.model
        merged=execution.mesh
        @test validate(merged).ok
        @test allunique(Tuple(merged.coords[:,node]) for node in 1:nnodes(merged))
        walls=(_sl_pair_tag(execution,"a","b"),_sl_pair_tag(execution,"a","c"),
               _sl_pair_tag(execution,"b","d"),_sl_pair_tag(execution,"c","d"))
        @test allunique(collect(walls))
        scope=_extrude_volume_pass(model,_SL_CALLER).scope
        @test scope!==nothing
        owners=Dict(walls[1]=>(1,2),walls[2]=>(1,3),
                    walls[3]=>(2,4),walls[4]=>(3,4))
        for (s,(lo,hi)) in owners
            bound=scope.surfaces[s]
            @test validate(bound).ok
            @test _sl_part_equal(scope.volumes[lo].surfaces[s],bound)
            @test _sl_part_equal(scope.volumes[hi].surfaces[s],bound)
            @test _sl_part_equal(mesh_model_surface(model,s),bound)
        end
        # All four walls weld exactly once across the grid's volume parts —
        # filtered by wall plane since interior incidence also counts each
        # Region's own cell interfaces.
        counts=_sl_face_counts((scope.volumes[tag].volume for tag in 1:4)...)
        @test all(c in (1,2) for c in values(counts))
        interior=Set(k for (k,c) in counts if c==2)
        @test Set(k for k in interior if all(v->v[1]==1.0,k))==
            union(_sl_surface_faces(scope.surfaces[walls[1]]),
                  _sl_surface_faces(scope.surfaces[walls[4]]))
        @test Set(k for k in interior if all(v->v[2]==1.0,k))==
            union(_sl_surface_faces(scope.surfaces[walls[2]]),
                  _sl_surface_faces(scope.surfaces[walls[3]]))
        for tag in 1:4
            _sl_positive(scope.volumes[tag].volume)
        end
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
        # An unrecombined partner whose own Source carries quadrangles keeps
        # a different loud diagnostic — upstream (Gmsh 4.15.2 oracle) rejects
        # that sweep with `Cannot extrude quadrangles without Recombine`
        # before the wall is ever welded. That is a distinguishable
        # producer diagnostic, not the shared-planner diagnostic and not a
        # nonconformity: a triangle-sourced partner welds the certified
        # split wall (the "unrecombined partner" testset below). The NoNew
        # owner's own part and the bound wall part still certify
        # standalone.
        execution=_sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};};
        """)
        s=_sl_shared_tag(execution)
        @test validate(mesh_model_volume(execution.model,1)).ok
        @test validate(mesh_model_surface(execution.model,s)).ok
        @test_throws r"Cannot extrude quadrangles without Recombine" mesh_model_volume(execution.model,2)
        @test_throws r"Cannot extrude quadrangles" _sl_model(
            _sl_squares()*"""
            left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
            right[]=Extrude{0,0,1}{Surface{2};Layers{2};};
            """;mesh_dim=3)
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

    @testset "three Region owners weld one shared lateral" begin
        # Gmsh 4.15.2 admits one literal Surface bounded by three identical
        # QuadTriNoNewVerts sweeps: two owners sit on the same wall side and
        # all three emit the same wall cells (2 tris on the interior layer,
        # the whole pyramid-base quad on the cap layer).
        execution=_sl_model(_sl_squares()*"""
        Point(7)={0.5,0.25,0,1};Point(8)={0.5,0.75,0,1};
        Line(8)={3,8};Line(9)={8,7};Line(10)={7,2};
        Curve Loop(3)={2,8,9,10};Plane Surface(3)={3};
        Transfinite Curve{8,9,10}=2;
        Transfinite Surface{3};Recombine Surface{3};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        b[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        c[]=Extrude{0,0,1}{Surface{3};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        model=execution.model
        merged=execution.mesh
        @test validate(merged).ok
        @test allunique(Tuple(merged.coords[:,node]) for node in 1:nnodes(merged))
        shared=reduce(intersect,(Int.(execution.lists[k])[3:end] for k in ("a","b","c")))
        @test length(shared)==1
        s=only(shared)
        @test sort(_extrude_quadtri_regions(model,s,_SL_CALLER))==[1,2,3]
        pass=_extrude_volume_pass(model,_SL_CALLER)
        scope=pass.scope
        bound=scope.surfaces[s]
        @test validate(bound).ok
        wall=_sl_surface_faces(bound)
        # Every owner emitted the bound part and sees it as its own wall
        # cells; same-side owners orient alike, the far side opposes.
        for owner in 1:3
            @test _sl_part_equal(scope.volumes[owner].surfaces[s],bound)
            @test Set(keys(_sl_wall_cycles(pass.parts[owner],1.0)))==wall
            _sl_positive(pass.parts[owner])
        end
        @test _sl_part_equal(mesh_model_surface(model,s),bound)
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
    end

    @testset "ordinary and AddVerts partners keep the whole-quad wall" begin
        # Gmsh 4.15.2 welds the literal lateral and keeps whole quadrangles
        # on it whenever a partner cannot emit the NoNewVerts split — an
        # ordinary `Layers+Recombine` sweep and a `QuadTriAddVerts` sweep
        # alike — while the NoNewVerts owner fans its layers so its own
        # boundary arrives on the same quads (pyramid bases), in either
        # extrusion order.
        cases=(
            # NoNewVerts first, ordinary Layers+Recombine second.
            ("""
            a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
            b[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;};
            """,("a","b"),1),
            # Ordinary partner extruded first — same welded whole-quad wall.
            ("""
            b[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;};
            a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
            """,("a","b"),2),
            # A QuadTriAddVerts partner also keeps its whole-quad wall.
            ("""
            a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
            b[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriAddVerts;};
            """,("a","b"),1),
        )
        for (body,(nonew_tag,other_tag),nonew) in cases
            execution=_sl_model(_sl_squares()*body;mesh_dim=3)
            model=execution.model
            merged=execution.mesh
            @test validate(merged).ok
            s=_sl_pair_tag(execution,nonew_tag,other_tag)
            @test sort(_extrude_quadtri_regions(model,s,_SL_CALLER))==[1,2]
            pass=_extrude_volume_pass(model,_SL_CALLER)
            scope=pass.scope
            bound=scope.surfaces[s]
            @test validate(bound).ok
            wall=_sl_surface_faces(bound)
            @test !isempty(wall)
            # The welded wall stays in whole quads — upstream's outcome.
            @test all(length(key)==4 for key in wall)
            @test _sl_dim2_faces_on(merged,1,1.0)==wall
            # The NoNewVerts owner's boundary arrives on the identical
            # whole-quad facets.
            @test Set(keys(_sl_wall_cycles(pass.parts[nonew],1.0)))==wall
            _sl_positive(pass.parts[nonew])
        end
    end

    @testset "unrecombined partner welds the certified split wall" begin
        # Gmsh 4.15.2 oracle (`_nonew_mixed` capture): the literal lateral
        # is emitted once as the partner's prism-subdivision triangles and
        # every owner arrives on the identical facets — the NoNewVerts
        # owner fans no centroids here, the welded wall is split. The
        # bound wall's diagonals seed the global `SubdivideExtrudedMesh`
        # edge set so the partner's phase-3 tetrahedra land on them.
        execution=_sl_model(_sl_squares_tri()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};};
        """;mesh_dim=3)
        model=execution.model
        merged=execution.mesh
        @test validate(merged).ok
        s=_sl_shared_tag(execution)
        @test sort(_extrude_quadtri_regions(model,s,_SL_CALLER))==[1,2]
        pass=_extrude_volume_pass(model,_SL_CALLER)
        scope=pass.scope
        @test scope!==nothing
        bound=scope.surfaces[s]
        @test validate(bound).ok && validate(pass.parts[1]).ok &&
            validate(pass.parts[2]).ok
        # Exactly the oracle's emitted wall: each layer quadrangle split on
        # the (bottom-first, top-opposite) diagonal.
        wall=_sl_surface_faces(bound)
        @test wall==Set((
            ((1.0,0.0,0.0),(1.0,0.0,0.5),(1.0,1.0,0.5)),
            ((1.0,0.0,0.0),(1.0,1.0,0.0),(1.0,1.0,0.5)),
            ((1.0,0.0,0.5),(1.0,0.0,1.0),(1.0,1.0,1.0)),
            ((1.0,0.0,0.5),(1.0,1.0,0.5),(1.0,1.0,1.0))))
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
        # Both Regions' wall incidences equal the bound Surface's cells —
        # the seeded subdivision lands the partner's prism faces on them.
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        # Opposite child incidences on every shared facet.
        left=_sl_wall_cycles(pass.parts[1],1.0)
        right=_sl_wall_cycles(pass.parts[2],1.0)
        @test keys(left)==keys(right)
        for key in keys(left)
            @test _sl_reversed(left[key],right[key])
        end
        # The shared wall welds (incidence two); the rest of both Regions'
        # surfaces is exactly the external boundary.
        counts=_sl_face_counts(pass.parts[1],pass.parts[2])
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
            # The NoNewVerts Region's surfaces live in the scope; the
            # partner's own laterals remesh against the same seeded edges.
            part=get(scope.surfaces,tag,nothing)
            part===nothing && (part=mesh_model_surface(
                model,tag;_extrude_edges=pass.edges))
            union!(exterior,_sl_surface_faces(part))
        end
        @test exterior==external
        # The standalone surface lookup replays the same bound wall.
        @test _sl_surface_faces(mesh_model_surface(model,s))==wall
    end

    @testset "stacked copied cap is the next sweep's Source" begin
        # Gmsh 4.15.2 welds the copied top cap literally as the stacked
        # sweep's Source — equal and unequal layer records alike. Tessella
        # emits the cap once through the top owner's plan and consumes it
        # bitwise; a different transform (bent stack) produces a degenerate
        # upstream volume and keeps the loud diagnostic.
        for layers in (2,3)
            execution=_sl_model("""
            Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
            Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
            Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
            Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
            Transfinite Curve{:}=2;
            Transfinite Surface{1};
            Recombine Surface{1};
            a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
            b[]=Extrude{0,0,1}{Surface{a[0]};Layers{$layers};Recombine;QuadTriNoNewVerts;};
            """;mesh_dim=3)
            model=execution.model
            merged=execution.mesh
            @test validate(merged).ok
            cap=Int(execution.lists["a"][1])
            @test sort(_extrude_quadtri_regions(model,cap,_SL_CALLER))==[1,2]
            pass=_extrude_volume_pass(model,_SL_CALLER)
            scope=pass.scope
            bound=scope.surfaces[cap]
            @test validate(bound).ok
            cap_faces=_sl_surface_faces(bound)
            @test !isempty(cap_faces)
            # The welded cap is exactly the facet set both volumes bound:
            # volume 1's top and volume 2's bottom coincide on z=1.
            @test _sl_dim2_faces_on(merged,3,1.0)==cap_faces
            _sl_positive(pass.parts[1]);_sl_positive(pass.parts[2])
        end
        # A bent stack (the second sweep turns +x) produces a degenerate
        # upstream volume — the loud diagnostic stays.
        @test_throws ArgumentError _sl_model("""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Transfinite Curve{:}=2;
        Transfinite Surface{1};
        Recombine Surface{1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        b[]=Extrude{1,0,0}{Surface{a[0]};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
    end

    @testset "non-generic NoNewVerts owner transports the split wall" begin
        # A two-quadrangle strip Source lands the left Region on the
        # `_ExtrudeNoNewQuadStripCatalog` family while the right Region is
        # the generic single-quad catalog. The Gmsh 4.15.2 oracle welds the
        # literal lateral between them (24 nodes; both Regions remesh the
        # wall) — its joint solver keeps layer two a whole quad while
        # splitting layer one. Tessella's strip catalog always splits its
        # boundary quadrangles, so the certified weld here is the strip's
        # own all-split wall transported as a contract to the generic
        # partner — a conforming decomposition that differs from the
        # oracle's per-layer mixed choice only inside the welded wall.
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={2,0,0,1};
        Point(3)={2,1,0,1}; Point(4)={0,1,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={3,0,0,1};Point(6)={3,1,0,1};
        Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
        Transfinite Curve{1,3}=3;
        Transfinite Curve{2,4,5,6,7}=2;
        Transfinite Surface{1};Transfinite Surface{2};
        Recombine Surface{:};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        model=execution.model
        merged=execution.mesh
        @test validate(merged).ok
        @test allunique(Tuple(merged.coords[:,node]) for node in 1:nnodes(merged))
        @test nnodes(merged)==24
        s=_sl_shared_tag(execution)
        @test sort(_extrude_quadtri_regions(model,s,_SL_CALLER))==[1,2]
        pass=_extrude_volume_pass(model,_SL_CALLER)
        scope=pass.scope
        @test scope!==nothing
        bound=scope.surfaces[s]
        @test validate(bound).ok && validate(pass.parts[1]).ok &&
            validate(pass.parts[2]).ok
        # The certified welded wall: the strip binder's per-layer diagonals
        # transported through the generic partner's forced chain contract —
        # identical bound part, identical incidences, opposite orientation.
        wall=_sl_surface_faces(bound)
        @test all(length(key)==3 for key in wall)
        @test Set(keys(_sl_wall_cycles(pass.parts[1],2.0)))==
            Set(keys(_sl_wall_cycles(pass.parts[2],2.0)))==wall
        left=_sl_wall_cycles(pass.parts[1],2.0)
        right=_sl_wall_cycles(pass.parts[2],2.0)
        for key in keys(left)
            @test _sl_reversed(left[key],right[key])
        end
        @test _sl_part_equal(scope.volumes[1].surfaces[s],bound)
        @test _sl_part_equal(scope.volumes[2].surfaces[s],bound)
        @test _sl_part_equal(mesh_model_surface(model,s),bound)
        @test _sl_dim2_faces_on(merged,1,2.0)==wall
        _sl_positive(pass.parts[1]);_sl_positive(pass.parts[2])
        counts=_sl_face_counts(pass.parts[1],pass.parts[2])
        @test all(c in (1,2) for c in values(counts))
        interior=Set(k for (k,c) in counts if c==2)
        @test Set(k for k in interior if all(v->v[1]==2.0,k))==wall
        # The generic-first order binds the other way: the generic plan's
        # all-split wall coincides bitwise with the strip's emission, so the
        # same certified weld assembles.
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={3,0,0,1};Point(6)={3,1,0,1};
        Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
        Transfinite Curve{5,7}=3;
        Transfinite Curve{1,2,3,4,6}=2;
        Transfinite Surface{1};Transfinite Surface{2};
        Recombine Surface{:};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==24
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        bound=pass.scope.surfaces[s]
        wall=_sl_surface_faces(bound)
        @test Set(keys(_sl_wall_cycles(pass.parts[1],1.0)))==
            Set(keys(_sl_wall_cycles(pass.parts[2],1.0)))==wall
        _sl_positive(pass.parts[1]);_sl_positive(pass.parts[2])
        # The same strip Source swept against an unrecombined partner
        # already welded the split wall through the seeded subdivision —
        # pin the exact emitted facets the oracle keeps (Gmsh 4.15.2
        # `r1_nongeneric_norecomb` capture: four triangles on the wall).
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={2,0,0,1};
        Point(3)={2,1,0,1}; Point(4)={0,1,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={3,0,0,1};Point(6)={3,1,0,1};
        Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
        Transfinite Curve{1,3}=3;
        Transfinite Curve{2,4,5,6,7}=2;
        Transfinite Surface{1};Transfinite Surface{2};
        Recombine Surface{1};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==24
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        bound=pass.scope.surfaces[s]
        wall=_sl_surface_faces(bound)
        @test wall==Set((
            ((2.0,0.0,0.0),(2.0,0.0,0.5),(2.0,1.0,0.5)),
            ((2.0,0.0,0.0),(2.0,1.0,0.0),(2.0,1.0,0.5)),
            ((2.0,0.0,0.5),(2.0,0.0,1.0),(2.0,1.0,1.0)),
            ((2.0,0.0,0.5),(2.0,1.0,0.5),(2.0,1.0,1.0))))
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],2.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        # A strip Source whose wall must stay whole — the partner is an
        # ordinary `Layers+Recombine` sweep — cannot fan inside the retained
        # catalog. The oracle welds that wall in whole quads (its strip
        # boundary cells adapt); the bounded catalog keeps the loud
        # retained-policy diagnostic instead.
        @test_throws r"policy differs from its retained catalog" _sl_model("""
        Point(1)={0,0,0,1}; Point(2)={2,0,0,1};
        Point(3)={2,1,0,1}; Point(4)={0,1,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={3,0,0,1};Point(6)={3,1,0,1};
        Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
        Transfinite Curve{1,3}=3;
        Transfinite Curve{2,4,5,6,7}=2;
        Transfinite Surface{1};Transfinite Surface{2};
        Recombine Surface{:};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;};
        """;mesh_dim=3)
        # Two strip Sources welded on one wall: upstream's QuadToTri solve
        # fails outright on the partner's lateral (37 errors, no volume 2
        # cells), so the loud planner diagnostic stays here as well — the
        # partner's three-quad Source matches no bounded catalog and the
        # boundary-category diagnosis now reaches it before the wall
        # contract is consulted.
        @test_throws r"planner" _sl_model("""
        Point(1)={0,0,0,1}; Point(2)={2,0,0,1};
        Point(3)={2,1,0,1}; Point(4)={0,1,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={4,0,0,1};Point(6)={4,1,0,1};
        Line(5)={3,5};Line(6)={5,6};Line(7)={6,4};
        Curve Loop(2)={5,6,7,-3};Plane Surface(2)={2};
        Transfinite Curve{1,3,5,7}=3;
        Transfinite Curve{2,4,6}=2;
        Transfinite Surface{:};
        Recombine Surface{:};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
    end

    @testset "reversed-orientation owners weld the shared lateral" begin
        # A reversed source loop makes its whole Region inside-out — every
        # lateral shell sign flips — while upstream still welds the shared
        # wall. The geometric side test (generatrix × sweep against each
        # Source's interior) admits exactly that reversal; the Gmsh 4.15.2
        # oracle emits 18 nodes and a welded shared wall for both loop
        # directions. Upstream's wall decomposition itself is
        # solver-dependent across runs (all-split, split+whole-quad and
        # all-quad walls all observed on the same fixture); the
        # deterministic emission pinned here is the all-split member of
        # that observed family.
        oracle_wall=Set((
            ((1.0,0.0,0.0),(1.0,0.0,0.5),(1.0,1.0,0.5)),
            ((1.0,0.0,0.0),(1.0,1.0,0.0),(1.0,1.0,0.5)),
            ((1.0,0.0,0.5),(1.0,0.0,1.0),(1.0,1.0,1.0)),
            ((1.0,0.0,0.5),(1.0,1.0,0.5),(1.0,1.0,1.0))))
        for (loop1,loop2) in (("{-4,-3,-2,-1}","{5,6,7,-2}"),
                             ("{1,2,3,4}","{2,-7,-6,-5}"))
            execution=_sl_model("""
            Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
            Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
            Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
            Curve Loop(1)=$loop1;Plane Surface(1)={1};
            Point(5)={2,0,0,1};Point(6)={2,1,0,1};
            Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
            Curve Loop(2)=$loop2;Plane Surface(2)={2};
            Transfinite Curve{:}=2;
            Transfinite Surface{1};Transfinite Surface{2};
            Recombine Surface{:};
            left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
            right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
            """;mesh_dim=3)
            model=execution.model
            merged=execution.mesh
            @test validate(merged).ok
            @test nnodes(merged)==18
            s=_sl_shared_tag(execution)
            @test sort(_extrude_quadtri_regions(model,s,_SL_CALLER))==[1,2]
            pass=_extrude_volume_pass(model,_SL_CALLER)
            bound=pass.scope.surfaces[s]
            @test validate(bound).ok
            wall=_sl_surface_faces(bound)
            @test wall==oracle_wall
            @test _sl_dim2_faces_on(merged,1,1.0)==wall
            left=_sl_wall_cycles(pass.parts[1],1.0)
            right=_sl_wall_cycles(pass.parts[2],1.0)
            @test Set(keys(left))==Set(keys(right))==wall
            for key in keys(left)
                @test _sl_reversed(left[key],right[key])
            end
            _sl_positive(pass.parts[1]);_sl_positive(pass.parts[2])
            counts=_sl_face_counts(pass.parts[1],pass.parts[2])
            @test Set(k for (k,c) in counts if c==2 &&
                all(v->v[1]==1.0,k))==wall
        end
        # The inside-out partner may also be an ordinary unrecombined
        # sweep — the seeded subdivision welds it the same way upstream
        # (Gmsh 4.15.2 `a_rev_loop_mixed` capture: 18 nodes).
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={2,0,0,1};Point(6)={2,1,0,1};
        Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={2,-7,-6,-5};Plane Surface(2)={2};
        Transfinite Curve{:}=2;
        Transfinite Surface{1};Transfinite Surface{2};
        Recombine Surface{1};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==18
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        wall=_sl_surface_faces(pass.scope.surfaces[s])
        @test wall==oracle_wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        # Mesh-level reversals keep the canonical plan and flip only the
        # published windings. `Reverse Curve{2}` on the generatrix and
        # `Reverse Surface{2}` on the Source both weld identically
        # upstream (Gmsh 4.15.2: 18 nodes, same observed wall family).
        for suffix in ("Reverse Curve{2};\n","Reverse Surface{2};\n")
            execution=_sl_model(_sl_squares()*"""
            left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
            right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
            """*suffix;mesh_dim=3)
            merged=execution.mesh
            @test validate(merged).ok
            @test nnodes(merged)==18
            s=_sl_shared_tag(execution)
            pass=_extrude_volume_pass(execution.model,_SL_CALLER)
            bound=pass.scope.surfaces[s]
            wall=_sl_surface_faces(bound)
            @test wall==oracle_wall
            @test _sl_dim2_faces_on(merged,1,1.0)==wall
            for tag in (1,2)
                @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
                _sl_positive(pass.parts[tag])
            end
        end
        # `Reverse Surface` on the shared wall itself flips the published
        # part's windings — the certificates keep canonical orientation —
        # while the welded facet set stays the oracle's (Gmsh 4.15.2
        # `a_rev_wall` capture: 18 nodes, same wall family).
        canonical=_sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        reversed=_sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        Reverse Surface{left[3]};
        """;mesh_dim=3)
        @test validate(reversed.mesh).ok
        @test nnodes(reversed.mesh)==18
        s=_sl_shared_tag(reversed)
        pass=_extrude_volume_pass(reversed.model,_SL_CALLER)
        bound=pass.scope.surfaces[s]
        wall=_sl_surface_faces(bound)
        @test wall==oracle_wall
        @test _sl_dim2_faces_on(reversed.mesh,1,1.0)==wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        # The stored scope part and the published surface mesh carry the
        # same flipped emission — classification sees one part.
        bound_c=_extrude_volume_pass(canonical.model,_SL_CALLER).scope.surfaces[
            _sl_shared_tag(canonical)]
        @test !_sl_part_equal(bound,bound_c)
        @test _sl_part_equal(mesh_model_surface(reversed.model,s),bound)
        @test _sl_surface_faces(geo_entity_mesh(reversed,2,s))==wall
        # A generatrix defined in the opposite direction (Line(2)={3,2}
        # with both loops re-ordered) welds the same wall — the solver
        # pins the mirror diagonal, which is the diagonal orientation the
        # oracle's own split-layer emission shows on this fixture.
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        Line(1)={1,2};Line(2)={3,2};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,-2,3,4};Plane Surface(1)={1};
        Point(5)={2,0,0,1};Point(6)={2,1,0,1};
        Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={5,6,7,2};Plane Surface(2)={2};
        Transfinite Curve{:}=2;
        Transfinite Surface{1};Transfinite Surface{2};
        Recombine Surface{:};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==18
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        wall=_sl_surface_faces(pass.scope.surfaces[s])
        @test wall==Set((
            ((1.0,0.0,0.0),(1.0,0.0,0.5),(1.0,1.0,0.0)),
            ((1.0,0.0,0.5),(1.0,0.0,1.0),(1.0,1.0,0.5)),
            ((1.0,0.0,0.5),(1.0,1.0,0.0),(1.0,1.0,0.5)),
            ((1.0,0.0,1.0),(1.0,1.0,0.5),(1.0,1.0,1.0))))
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
    end

    @testset "multi-cell and strip sources on the shared lateral" begin
        # Two-quad strip Sources on both sweeps: the 1-segment wall at
        # x=2 welds four triangles bitwise (Gmsh 4.15.2 `b_strip_short`
        # capture: 30 nodes, same four facets); the all-split emission is
        # the strip catalog's own boundary policy transported through the
        # bitwise shared-part dedup rather than a diagonal contract.
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={2,0,0,1};
        Point(3)={2,1,0,1}; Point(4)={0,1,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={4,0,0,1};Point(6)={4,1,0,1};
        Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
        Transfinite Curve{1,3,5,7}=3; Transfinite Curve{2,4,6}=2;
        Transfinite Surface{1};Transfinite Surface{2};
        Recombine Surface{:};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==30
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        bound=pass.scope.surfaces[s]
        wall=_sl_surface_faces(bound)
        @test wall==Set((
            ((2.0,0.0,0.0),(2.0,0.0,0.5),(2.0,1.0,0.5)),
            ((2.0,0.0,0.0),(2.0,1.0,0.0),(2.0,1.0,0.5)),
            ((2.0,0.0,0.5),(2.0,0.0,1.0),(2.0,1.0,1.0)),
            ((2.0,0.0,0.5),(2.0,1.0,0.5),(2.0,1.0,1.0))))
        @test _sl_dim2_faces_on(merged,1,2.0)==wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],2.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        # The 2-segment generatrix welds eight triangles — the oracle's
        # `b_strip_long` capture (27 nodes) including its transfinite
        # evaluation noise on the interior generatrix node.
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,2,0,1}; Point(4)={0,2,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={2,0,0,1};Point(6)={2,2,0,1};
        Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
        Transfinite Curve{1,3,5,7}=2; Transfinite Curve{2,4,6}=3;
        Transfinite Surface{1};Transfinite Surface{2};
        Recombine Surface{:};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==27
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        wall=_sl_surface_faces(pass.scope.surfaces[s])
        mid=0.9999999999973842
        @test wall==Set((
            ((1.0,0.0,0.0),(1.0,0.0,0.5),(1.0,mid,0.5)),
            ((1.0,0.0,0.0),(1.0,mid,0.0),(1.0,mid,0.5)),
            ((1.0,0.0,0.5),(1.0,0.0,1.0),(1.0,mid,1.0)),
            ((1.0,0.0,0.5),(1.0,mid,0.5),(1.0,mid,1.0)),
            ((1.0,mid,0.0),(1.0,mid,0.5),(1.0,2.0,0.5)),
            ((1.0,mid,0.0),(1.0,2.0,0.0),(1.0,2.0,0.5)),
            ((1.0,mid,0.5),(1.0,mid,1.0),(1.0,2.0,1.0)),
            ((1.0,mid,0.5),(1.0,2.0,0.5),(1.0,2.0,1.0))))
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        # Unstructured multi-cell Sources sweep each source triangle into
        # its own prism column through the general per-column planner:
        # boundary columns keep the shared wall whole while every owner
        # does (`b_multi_nonew_nonew` oracle capture: 63 nodes, a 4-quad
        # wall), and the per-segment bound part transports the identical
        # facets to the second owner bitwise. Node and cell multisets pin
        # the Gmsh 4.15.2 emission exactly.
        execution=_sl_model("""
        Point(1)={0,0,0,0.6}; Point(2)={1,0,0,0.6};
        Point(3)={1,1,0,0.6}; Point(4)={0,1,0,0.6};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={2,0,0,0.6};Point(6)={2,1,0,0.6};
        Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==63
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        bound=pass.scope.surfaces[s]
        wall=_sl_surface_faces(bound)
        mid=0.49999999999869194
        @test wall==Set((
            ((1.0,0.0,0.0),(1.0,0.0,0.5),(1.0,mid,0.0),(1.0,mid,0.5)),
            ((1.0,0.0,0.5),(1.0,0.0,1.0),(1.0,mid,0.5),(1.0,mid,1.0)),
            ((1.0,mid,0.0),(1.0,mid,0.5),(1.0,1.0,0.0),(1.0,1.0,0.5)),
            ((1.0,mid,0.5),(1.0,mid,1.0),(1.0,1.0,0.5),(1.0,1.0,1.0))))
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        left=_sl_wall_cycles(pass.parts[1],1.0)
        right=_sl_wall_cycles(pass.parts[2],1.0)
        for key in keys(left)
            @test _sl_reversed(left[key],right[key])
        end
        # The coarse unstructured pair welds the same way: each native
        # Source carries two triangles (the oracle fans its own sources
        # through an extra centroid — 24 nodes — while the native pair
        # stays at 18), the certified wall stays two whole quadrangles.
        execution=_sl_model("""
        Point(1)={0,0,0,3.0}; Point(2)={1,0,0,3.0};
        Point(3)={1,1,0,3.0}; Point(4)={0,1,0,3.0};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={2,0,0,3.0};Point(6)={2,1,0,3.0};
        Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==18
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        wall=_sl_surface_faces(pass.scope.surfaces[s])
        @test wall==Set((
            ((1.0,0.0,0.0),(1.0,0.0,0.5),(1.0,1.0,0.0),(1.0,1.0,0.5)),
            ((1.0,0.0,0.5),(1.0,0.0,1.0),(1.0,1.0,0.5),(1.0,1.0,1.0))))
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
    end

    @testset "unstructured multi-cell prism columns" begin
        # A multi-cell NoNew sweep against an ordinary unrecombined
        # partner (`b_multi_nonew_ord` oracle: 63 nodes): the partner's
        # prism-subdivision vote splits the shared wall, the NoNew owner
        # free-emits the same diagonal on every segment quadrangle, and
        # the seeded `SubdivideExtrudedMesh` edge set lands the partner's
        # tets on the identical facets — the oracle's 8-triangle wall
        # bitwise, including its mid-generatrix node.
        execution=_sl_model("""
        Point(1)={0,0,0,0.6}; Point(2)={1,0,0,0.6};
        Point(3)={1,1,0,0.6}; Point(4)={0,1,0,0.6};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={2,0,0,0.6};Point(6)={2,1,0,0.6};
        Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};};
        """;mesh_dim=3)
        model=execution.model
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==63
        @test _sl_msh_counts(merged)==Dict(1=>40,2=>112,4=>140,6=>4,7=>8)
        s=_sl_shared_tag(execution)
        @test sort(_extrude_quadtri_regions(model,s,_SL_CALLER))==[1,2]
        pass=_extrude_volume_pass(model,_SL_CALLER)
        bound=pass.scope.surfaces[s]
        wall=_sl_surface_faces(bound)
        mid=0.49999999999869194
        @test wall==Set((
            ((1.0,0.0,0.0),(1.0,0.0,0.5),(1.0,mid,0.5)),
            ((1.0,0.0,0.0),(1.0,mid,0.0),(1.0,mid,0.5)),
            ((1.0,0.0,0.5),(1.0,0.0,1.0),(1.0,mid,1.0)),
            ((1.0,0.0,0.5),(1.0,mid,0.5),(1.0,mid,1.0)),
            ((1.0,mid,0.0),(1.0,mid,0.5),(1.0,1.0,0.5)),
            ((1.0,mid,0.0),(1.0,1.0,0.0),(1.0,1.0,0.5)),
            ((1.0,mid,0.5),(1.0,mid,1.0),(1.0,1.0,1.0)),
            ((1.0,mid,0.5),(1.0,1.0,0.5),(1.0,1.0,1.0))))
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        left=_sl_wall_cycles(pass.parts[1],1.0)
        right=_sl_wall_cycles(pass.parts[2],1.0)
        for key in keys(left)
            @test _sl_reversed(left[key],right[key])
        end
        # The same pair on the coarse two-triangle Sources: the wall
        # carries four triangles through the partner's vote. The twist
        # fixpoint selects the mirror diagonal on this source — a valid
        # certified weld the partner reproduces bitwise, oriented the
        # opposite way from the oracle's own centroid-fanned sources.
        execution=_sl_model("""
        Point(1)={0,0,0,3.0}; Point(2)={1,0,0,3.0};
        Point(3)={1,1,0,3.0}; Point(4)={0,1,0,3.0};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={2,0,0,3.0};Point(6)={2,1,0,3.0};
        Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==21
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        wall=_sl_surface_faces(pass.scope.surfaces[s])
        @test wall==Set((
            ((1.0,0.0,0.0),(1.0,0.0,0.5),(1.0,1.0,0.0)),
            ((1.0,0.0,0.5),(1.0,0.0,1.0),(1.0,1.0,0.5)),
            ((1.0,0.0,0.5),(1.0,1.0,0.0),(1.0,1.0,0.5)),
            ((1.0,0.0,1.0),(1.0,1.0,0.5),(1.0,1.0,1.0))))
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        # The isolated sweep (`b_multi_solo` oracle: 36 nodes — the same
        # node/cell multiset Tessella emits): fourteen prism columns, all
        # boundary laterals split on the chain-preferred diagonal.
        execution=_sl_model("""
        Point(1)={0,0,0,0.6}; Point(2)={1,0,0,0.6};
        Point(3)={1,1,0,0.6}; Point(4)={0,1,0,0.6};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==36
        @test _sl_msh_counts(merged)==Dict(1=>24,2=>60,4=>56,6=>4,7=>8)
        mid=0.49999999999869194
        @test _sl_dim2_faces_on(merged,1,1.0)==Set((
            ((1.0,0.0,0.0),(1.0,0.0,0.5),(1.0,mid,0.5)),
            ((1.0,0.0,0.0),(1.0,mid,0.0),(1.0,mid,0.5)),
            ((1.0,0.0,0.5),(1.0,0.0,1.0),(1.0,mid,1.0)),
            ((1.0,0.0,0.5),(1.0,mid,0.5),(1.0,mid,1.0)),
            ((1.0,mid,0.0),(1.0,mid,0.5),(1.0,1.0,0.5)),
            ((1.0,mid,0.0),(1.0,1.0,0.0),(1.0,1.0,0.5)),
            ((1.0,mid,0.5),(1.0,mid,1.0),(1.0,1.0,1.0)),
            ((1.0,mid,0.5),(1.0,1.0,0.5),(1.0,1.0,1.0))))
        _sl_positive(merged)
        # A four-triangle dart Source (two boundary Steiner nodes split
        # its long edges — `b_dart_solo` oracle: 18 nodes, the identical
        # multiset): the reflex-corner triangulation sweeps four prism
        # columns through the same per-column planner.
        execution=_sl_model("""
        Point(1)={0,0,0,1.2}; Point(2)={2,0,0,1.2};
        Point(3)={1,0.3,0,1.2}; Point(4)={0,1,0,1.2};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==18
        @test _sl_msh_counts(merged)==Dict(1=>20,2=>32,4=>24)
        _sl_positive(merged)
        # The non-transfinite warped Source (BSpline generatrix bulging
        # to z=0.3): the former boundary-category residual — the planner
        # needs no planar assumption, and the 38-node oracle multiset
        # matches Tessella's emission exactly (`b_warped_nonew`).
        execution=_sl_model("""
        Point(1)={0,0,0,0.6}; Point(2)={1,0,0,0.6};
        Point(3)={1,1,0,0.6}; Point(4)={0,1,0,0.6};
        Point(5)={0.5,0.0,0.3,0.6};
        BSpline(1)={1,5,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Ruled Surface(1)={1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==38
        @test _sl_msh_counts(merged)==Dict(1=>24,2=>60,4=>56,6=>4,7=>8)
        _sl_positive(merged)
        # The minimal unstructured source — a single square carrying two
        # triangles: one column's shared interior edge plus its three
        # split boundary laterals is the twisted all-split combination
        # the fixpoint resolves. The oracle's default Frontal-Delaunay
        # surface mesher inserts a centroid into this fixture's Source
        # and fans it to four triangles (15 solo nodes, 24 nodes on the
        # pair), but that insertion is the 2D mesher's member choice, not
        # the sweep's: `Mesh.Algorithm=8` reproduces this native member
        # bitwise (the same 2-4 diagonal, 12 nodes / 12 tets), a
        # transfinite square never fans (12 nodes upstream too), and the
        # insertion is source-geometry dependent — a parallelogram,
        # trapezoid, skew and kite Source all stay two triangles there
        # and match this emission bitwise. The sweep's contract is the
        # producer's actual Source mesh, so the planner must not insert
        # sweep-time centroids.
        execution=_sl_model("""
        Point(1)={0,0,0,3.0}; Point(2)={1,0,0,3.0};
        Point(3)={1,1,0,3.0}; Point(4)={0,1,0,3.0};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==12
        @test _sl_msh_counts(merged)==Dict(1=>16,2=>20,4=>12)
        _sl_positive(merged)
        # `QuadTriNoNewVerts RecombLaterals` on the same unstructured
        # source retains every lateral quadrangle whole — the oracle's
        # `b_multi_recomblat` capture (36 nodes, 28 prisms, no splits).
        # The per-column planner forces every interior edge whole on this
        # record, so all fourteen columns land on the retained-prism
        # template instead of a twisted state the fixpoint could not
        # resolve.
        execution=_sl_model("""
        Point(1)={0,0,0,0.6}; Point(2)={1,0,0,0.6};
        Point(3)={1,1,0,0.6}; Point(4)={0,1,0,0.6};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts RecombLaterals;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==36
        @test _sl_msh_counts(merged)==Dict(1=>24,2=>28,3=>16,6=>28)
        _sl_positive(merged)
    end

    @testset "irregular boundary-curve counts" begin
        # LAY01 residual: Sources whose boundaries carry five, six, or
        # irregular curve counts reach the same per-column multi-cell
        # planner — every bounded catalog answers its own side count and
        # the actual Source mesh decides the rest. The oracle sweeps all
        # of these natively (`d_pentagon`, `d_hexagon`, `d_lshape`,
        # `d_lens` captures).
        #
        # A five-curve pentagon: the oracle's Frontal-Delaunay Source fan
        # (five triangles around one centroid) is the native mesh
        # bitwise, so the sweep's node/cell multiset matches the oracle
        # exactly — 18 nodes, thirty tetrahedra.
        execution=_sl_model("""
        Point(1)={0,0,0,3.0};Point(2)={1,0,0,3.0};Point(3)={1,1,0,3.0};
        Point(4)={0.5,1.25,0,3.0};Point(5)={0,1,0,3.0};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,5};Line(5)={5,1};
        Curve Loop(1)={1,2,3,4,5};Plane Surface(1)={1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==18
        @test _sl_msh_counts(merged)==Dict(1=>20,2=>30,4=>30)
        _sl_positive(merged)
        # A six-curve hexagon: both meshers triangulate it with two
        # interior nodes, so the multisets match bitwise — 24 nodes, 40
        # tetrahedra plus 4 pyramids. The pyramids are the interior-pair
        # rule at work: the centroid-to-centroid edge joins two interior
        # cells (no boundary edge on either), so its wall stays whole
        # and each incident column emits a pyramid on the shared base.
        execution=_sl_model("""
        Point(1)={0,0,0,3.0};Point(2)={1,0,0,3.0};Point(3)={2,0.5,0,3.0};
        Point(4)={2,1.5,0,3.0};Point(5)={1,2,0,3.0};Point(6)={0,1.5,0,3.0};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,5};Line(5)={5,6};Line(6)={6,1};
        Curve Loop(1)={1,2,3,4,5,6};Plane Surface(1)={1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==24
        @test _sl_msh_counts(merged)==Dict(1=>24,2=>40,4=>40,7=>4)
        _sl_positive(merged)
        # The L-shaped six-curve Source: the oracle's two-centroid fan
        # (eight triangles, interior nodes inside both arms, the corner
        # chord 1-4 kept whole between two interior cells) reaches 24
        # nodes; the native four-triangle all-corner triangulation has
        # no interior cells, splits every interior edge, and emits 18
        # nodes / 24 tetrahedra — a certified member inside the same
        # source-mesher family, recorded like the two-triangle square.
        execution=_sl_model("""
        Point(1)={0,0,0,3.0};Point(2)={2,0,0,3.0};Point(3)={2,1,0,3.0};
        Point(4)={1,1,0,3.0};Point(5)={1,2,0,3.0};Point(6)={0,2,0,3.0};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,5};Line(5)={5,6};Line(6)={6,1};
        Curve Loop(1)={1,2,3,4,5,6};Plane Surface(1)={1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==18
        @test _sl_msh_counts(merged)==Dict(1=>24,2=>32,4=>24)
        _sl_positive(merged)
        # A two-curve lens (one arc, one closing chord): the oracle
        # sweeps it to 14 nodes / 12 tetrahedra and the native emission
        # matches that multiset bitwise. A single closed generatrix
        # still fails upstream — its swept cap cannot be meshed — and
        # keeps the boundary-category diagnostic here.
        execution=_sl_model("""
        Point(1)={0,0,0,3.0};Point(2)={1,0,0,3.0};Point(3)={0.5,-0.5,0,3.0};
        Circle(1)={1,3,2};Line(2)={2,1};
        Curve Loop(1)={1,2};Plane Surface(1)={1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==14
        @test _sl_msh_counts(merged)==Dict(1=>12,2=>20,4=>12)
        _sl_positive(merged)
        @test_throws r"single-curve source boundaries" _sl_model("""
        Point(1)={1,0,0,3.0};Point(2)={0.7,0.7,0,3.0};Point(3)={0,1,0,3.0};
        Point(4)={-0.7,0.7,0,3.0};Point(5)={-1,0,0,3.0};
        Point(6)={-0.7,-0.7,0,3.0};Point(7)={0,-1,0,3.0};Point(8)={0.7,-0.7,0,3.0};
        BSpline(1)={1,2,3,4,5,6,7,8,1};
        Curve Loop(1)={1};Plane Surface(1)={1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
    end

    @testset "irregular shared-wall pairs" begin
        # A pentagon welded to a heptagon over two generatrices — the
        # same-owner multi-wall slice the joint certificate now groups
        # per owner set. The oracle (`d_pent_pair`) emits 36 nodes and
        # the identical element multiset: sixty tetrahedra, twelve
        # pyramids, and two whole-quadrangle walls (four quads) — both
        # NoNewVerts owners abstain, so the shared walls stay whole.
        execution=_sl_model("""
        Point(1)={0,0,0,3.0}; Point(2)={1,0,0,3.0};
        Point(3)={1.4,0.6,0,3.0}; Point(4)={0.7,1.2,0,3.0}; Point(5)={-0.2,0.7,0,3.0};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,5};Line(5)={5,1};
        Curve Loop(1)={1,2,3,4,5};Plane Surface(1)={1};
        Point(6)={2.2,0.2,0,3.0}; Point(7)={2.6,0.9,0,3.0};
        Point(8)={1.8,1.4,0,3.0}; Point(9)={0.8,1.5,0,3.0};
        Line(6)={2,6};Line(7)={6,7};Line(8)={7,8};Line(9)={8,9};Line(10)={9,4};
        Curve Loop(2)={6,7,8,9,10,-3,-2};Plane Surface(2)={2};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        model=execution.model
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==36
        @test _sl_msh_counts(merged)==Dict(1=>38,2=>60,3=>4,4=>60,7=>12)
        shared=intersect(Int.(execution.lists["left"])[3:end],
                         Int.(execution.lists["right"])[3:end])
        @test length(shared)==2
        pass=_extrude_volume_pass(model,_SL_CALLER)
        for s in shared
            @test sort(_extrude_quadtri_regions(model,s,_SL_CALLER))==[1,2]
            wall=_sl_surface_faces(pass.scope.surfaces[s])
            @test length(wall)==2
            @test all(length(key)==4 for key in wall)
        end
        for tag in (1,2)
            _sl_positive(pass.parts[tag])
        end
        # The same pair with an ordinary unrecombined partner welds the
        # two walls split (`d_pent_ord` oracle: 36 nodes, the identical
        # element multiset — all tetrahedra, the partner's prism columns
        # subdivided three apiece). The partner's vote forces the split,
        # the bound per-segment parts transport the contract, and the
        # seeded edge set reproduces the diagonal facets on the partner
        # side bitwise.
        execution=_sl_model("""
        Point(1)={0,0,0,3.0}; Point(2)={1,0,0,3.0};
        Point(3)={1.4,0.6,0,3.0}; Point(4)={0.7,1.2,0,3.0}; Point(5)={-0.2,0.7,0,3.0};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,5};Line(5)={5,1};
        Curve Loop(1)={1,2,3,4,5};Plane Surface(1)={1};
        Point(6)={2.2,0.2,0,3.0}; Point(7)={2.6,0.9,0,3.0};
        Point(8)={1.8,1.4,0,3.0}; Point(9)={0.8,1.5,0,3.0};
        Line(6)={2,6};Line(7)={6,7};Line(8)={7,8};Line(9)={8,9};Line(10)={9,4};
        Curve Loop(2)={6,7,8,9,10,-3,-2};Plane Surface(2)={2};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};};
        """;mesh_dim=3)
        model=execution.model
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==36
        @test _sl_msh_counts(merged)==Dict(1=>38,2=>68,4=>84)
        shared=intersect(Int.(execution.lists["left"])[3:end],
                         Int.(execution.lists["right"])[3:end])
        @test length(shared)==2
        pass=_extrude_volume_pass(model,_SL_CALLER)
        for s in shared
            @test sort(_extrude_quadtri_regions(model,s,_SL_CALLER))==[1,2]
            wall=_sl_surface_faces(pass.scope.surfaces[s])
            @test length(wall)==4
            @test all(length(key)==3 for key in wall)
        end
        for tag in (1,2)
            _sl_positive(pass.parts[tag])
        end
    end

    @testset "constrained and compound shared-wall boundaries" begin
        # `Recombine Surface{...}` applied to the generated wall — the same
        # outcome an end-of-file `Recombine Surface{:}` reaches — makes the
        # oracle keep the welded wall in whole quadrangles while every free
        # lateral stays split (Gmsh 4.15.2 `r2_recomb_*` captures: both
        # Regions fan, 2 quads on the wall). The certified owners fan their
        # layers and the bound part is the oracle's whole-quad pair.
        for suffix in (
                "Recombine Surface{left[3]};\n",
                "Recombine Surface{:};\n")
            execution=_sl_model(_sl_squares()*"""
            left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
            right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
            """*suffix;mesh_dim=3)
            model=execution.model
            merged=execution.mesh
            @test validate(merged).ok
            s=_sl_shared_tag(execution)
            pass=_extrude_volume_pass(model,_SL_CALLER)
            bound=pass.scope.surfaces[s]
            wall=_sl_surface_faces(bound)
            @test all(length(key)==4 for key in wall)
            @test length(wall)==2
            for tag in (1,2)
                @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
                _sl_positive(pass.parts[tag])
            end
        end
        # A periodic relation declared on the shared wall or its
        # generatrix resolves nothing geometrically — upstream reports the
        # failed connection and welds the wall anyway; the identical
        # certified weld assembles here.
        execution=_sl_model(_sl_squares()*"""
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        Periodic Curve{4}={-6} Translate{2,0,0};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        bound=pass.scope.surfaces[s]
        wall=_sl_surface_faces(bound)
        @test Set(keys(_sl_wall_cycles(pass.parts[1],1.0)))==
            Set(keys(_sl_wall_cycles(pass.parts[2],1.0)))==wall
        # Compound surfaces over the sweep Sources break upstream's
        # extruded-vertex search outright (23 errors, no volume 2 cells),
        # so the dedicated compound diagnostic stays.
        @test_throws r"compound boundary entities" _sl_model(
            _sl_squares()*"""
            Compound Surface{1,2};
            left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
            right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
            """;mesh_dim=3)
    end

    @testset "warped ruled and BSpline sources" begin
        # A BSpline generatrix bowing past its own patch corner folds the
        # ruled surface — `model_to_mixed`'s on-surface audit inverted an
        # exact corner node through the loose `XYZtoUV` variant, which
        # accepted a stalled Newton iterate and rejected the source as
        # 0.124 off-patch. The strict inverse (the `getParametrization`
        # contract: on-surface and XYZ convergence together) recovers the
        # corner and the whole extrude/shared-lateral certificate
        # assembles. Node and cell multisets pin the Gmsh 4.15.2 captures
        # for these fixtures.
        _sl_msh_counts(mesh::Mesh)=Dict{Int,Int}(
            k=>v for (k,v) in (
                1=>size(mesh.segs,2),2=>size(mesh.tris,2),
                4=>size(mesh.tets,2)) if v>0)
        _sl_msh_counts(mesh::MixedMesh)=begin
            counts=Dict{Int,Int}()
            for block in _sl_blocks(mesh)
                block isa ElementBlock && msh_dimension(block.msh)>=1 ||
                    continue
                counts[Int(block.msh)]=
                    get(counts,Int(block.msh),0)+size(block.nodes,2)
            end
            counts
        end
        # Ruled Sources on both sweeps (`c_ruled_both`): the folded
        # BSpline edge sits on the right Region's Source — the exact
        # corner (2,0,0) stalled the loose inverse at 0.124 residual.
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        Point(5)={0.5,0.0,0.3,1};
        BSpline(1)={1,5,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Ruled Surface(1)={1};
        Point(6)={2,0,0,1};Point(7)={2,1,0,1};
        Point(8)={2.5,0.0,0.25,1};
        BSpline(5)={2,8,6};Line(6)={6,7};Line(7)={7,3};
        Curve Loop(2)={5,6,7,-2};Ruled Surface(2)={2};
        Transfinite Curve{:}=2;
        Transfinite Surface{1};Transfinite Surface{2};
        Recombine Surface{:};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        model=execution.model
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==22
        @test _sl_msh_counts(merged)==Dict(1=>26,2=>32,3=>2,4=>20,7=>2)
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(model,_SL_CALLER)
        bound=pass.scope.surfaces[s]
        wall=_sl_surface_faces(bound)
        # Oracle wall facets at x=1 (the all-split member of the observed
        # solver-dependent family — the same deterministic emission the
        # planar-source certificates carry).
        @test wall==Set((
            ((1.0,0.0,0.0),(1.0,0.0,0.5),(1.0,1.0,0.5)),
            ((1.0,0.0,0.0),(1.0,1.0,0.0),(1.0,1.0,0.5)),
            ((1.0,0.0,0.5),(1.0,0.0,1.0),(1.0,1.0,1.0)),
            ((1.0,0.0,0.5),(1.0,1.0,0.5),(1.0,1.0,1.0))))
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        # The ruled Source part classifies through `model_to_mixed` —
        # the exact call site that rejected on the folded corner.
        source_part=pass.scope.surfaces[2]
        plain=MixedMesh(source_part.coords,
            ElementBlock[ElementBlock(b.msh,b.nodes)
                         for b in source_part.blocks
                         if b isa ElementBlock && msh_dimension(b.msh)==2])
        projected=model_to_mixed(model,plain,2)
        @test validate(projected).ok
        @test projected.coords==source_part.coords
        # And a node lifted off the patch still fails the audit loudly —
        # the strict inverse does not weaken the certificate.
        offcoords=copy(source_part.coords);offcoords[3,2]+=0.01
        @test_throws r"off Surface\[2\]" model_to_mixed(model,
            MixedMesh(offcoords,plain.blocks),2)
        # Ruled Source on the right of a planar partner (`c_flip_ruled`):
        # same folded patch, same oracle multiset — 20 nodes this time
        # (one fewer bulge column on the left's flat boundary).
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Point(5)={2,0,0,1};Point(6)={2,1,0,1};
        Point(7)={2.5,0.0,0.3,1};
        BSpline(5)={2,7,5};Line(6)={5,6};Line(7)={6,3};
        Curve Loop(2)={5,6,7,-2};Ruled Surface(2)={2};
        Transfinite Curve{:}=2;
        Transfinite Surface{1};Transfinite Surface{2};
        Recombine Surface{:};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==20
        @test _sl_msh_counts(merged)==Dict(1=>26,2=>32,3=>2,4=>20,7=>2)
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        wall=_sl_surface_faces(pass.scope.surfaces[s])
        @test wall==Set((
            ((1.0,0.0,0.0),(1.0,0.0,0.5),(1.0,1.0,0.5)),
            ((1.0,0.0,0.0),(1.0,1.0,0.0),(1.0,1.0,0.5)),
            ((1.0,0.0,0.5),(1.0,0.0,1.0),(1.0,1.0,1.0)),
            ((1.0,0.0,0.5),(1.0,1.0,0.5),(1.0,1.0,1.0))))
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        # Ruled NoNew Source against an ordinary unrecombined partner
        # (`c_ruled_ord`): the certified split wall welds onto the
        # partner's prism-subdivision triangles (20 nodes oracle).
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        Point(5)={0.5,0.0,0.3,1};
        BSpline(1)={1,5,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Ruled Surface(1)={1};
        Point(6)={2,0,0,1};Point(7)={2,1,0,1};
        Line(5)={2,6};Line(6)={6,7};Line(7)={7,3};
        Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
        Transfinite Curve{:}=2;
        Transfinite Surface{1};Transfinite Surface{2};
        Recombine Surface{1};
        left[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{2};};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==20
        @test _sl_msh_counts(merged)==Dict(1=>26,2=>34,3=>1,4=>22,7=>1)
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        wall=_sl_surface_faces(pass.scope.surfaces[s])
        @test wall==Set((
            ((1.0,0.0,0.0),(1.0,0.0,0.5),(1.0,1.0,0.5)),
            ((1.0,0.0,0.0),(1.0,1.0,0.0),(1.0,1.0,0.5)),
            ((1.0,0.0,0.5),(1.0,0.0,1.0),(1.0,1.0,1.0)),
            ((1.0,0.0,0.5),(1.0,1.0,0.5),(1.0,1.0,1.0))))
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        # Three-layer sweep on the same ruled/planar pair (`c_ruled3`):
        # per-layer wall facets pin all three intervals, 26 nodes oracle.
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        Point(5)={0.5,0.0,0.3,1};
        BSpline(1)={1,5,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Ruled Surface(1)={1};
        Point(6)={2,0,0,1};Point(7)={2,1,0,1};
        Line(5)={2,6};Line(6)={6,7};Line(7)={7,3};
        Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
        Transfinite Curve{:}=2;
        Transfinite Surface{1};Transfinite Surface{2};
        Recombine Surface{:};
        left[]=Extrude{0,0,1}{Surface{1};Layers{3};Recombine;QuadTriNoNewVerts;};
        right[]=Extrude{0,0,1}{Surface{2};Layers{3};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==26
        @test _sl_msh_counts(merged)==Dict(1=>32,2=>46,3=>2,4=>32,7=>2)
        s=_sl_shared_tag(execution)
        pass=_extrude_volume_pass(execution.model,_SL_CALLER)
        third=1/3
        wall=_sl_surface_faces(pass.scope.surfaces[s])
        @test wall==Set((
            ((1.0,0.0,0.0),(1.0,0.0,third),(1.0,1.0,third)),
            ((1.0,0.0,0.0),(1.0,1.0,0.0),(1.0,1.0,third)),
            ((1.0,0.0,third),(1.0,0.0,2*third),(1.0,1.0,2*third)),
            ((1.0,0.0,third),(1.0,1.0,third),(1.0,1.0,2*third)),
            ((1.0,0.0,2*third),(1.0,0.0,1.0),(1.0,1.0,1.0)),
            ((1.0,0.0,2*third),(1.0,1.0,2*third),(1.0,1.0,1.0))))
        @test _sl_dim2_faces_on(merged,1,1.0)==wall
        for tag in (1,2)
            @test Set(keys(_sl_wall_cycles(pass.parts[tag],1.0)))==wall
            _sl_positive(pass.parts[tag])
        end
        # The solo ruled NoNew sweep (`c_solo_ruled`): 14 nodes oracle,
        # the generic single-quad catalog emission.
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        Point(5)={0.5,0.0,0.3,1};
        BSpline(1)={1,5,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Ruled Surface(1)={1};
        Transfinite Curve{:}=2;
        Transfinite Surface{1};
        Recombine Surface{1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;QuadTriNoNewVerts;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==14
        @test _sl_msh_counts(merged)==Dict(1=>16,2=>18,3=>1,4=>10,7=>1)
        _sl_positive(merged)
        # The ordinary recombined sweep of the same warped Source — 14
        # nodes oracle (16 lines/10 quads/2 hexes).
        execution=_sl_model("""
        Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
        Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
        Point(5)={0.5,0.0,0.3,1};
        BSpline(1)={1,5,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Ruled Surface(1)={1};
        Transfinite Curve{:}=2;
        Transfinite Surface{1};
        Recombine Surface{1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};Recombine;};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==14
        @test _sl_msh_counts(merged)==Dict(1=>16,3=>10,5=>2)
        _sl_positive(merged)
        # The non-transfinite warped Source: the ruled mesher discretizes
        # the (u,v) domain and the sweep welds the generatrices — 38
        # nodes oracle, the full unstructured fill.
        execution=_sl_model("""
        Point(1)={0,0,0,0.6}; Point(2)={1,0,0,0.6};
        Point(3)={1,1,0,0.6}; Point(4)={0,1,0,0.6};
        Point(5)={0.5,0.0,0.3,0.6};
        BSpline(1)={1,5,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Ruled Surface(1)={1};
        a[]=Extrude{0,0,1}{Surface{1};Layers{2};};
        """;mesh_dim=3)
        merged=execution.mesh
        @test validate(merged).ok
        @test nnodes(merged)==38
        @test _sl_msh_counts(merged)==Dict(1=>24,2=>60,4=>84)
        _sl_positive(merged)
    end
end
