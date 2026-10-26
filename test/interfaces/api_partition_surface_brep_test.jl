module PartitionSurfaceBRepFunctionalChecks
# Functional checks for top-dimension-2 partition topology — the `t21.geo`
# shape: two adjacent Surfaces sharing a declared Curve, partitioned under
# `Mesh.PartitionCreateTopology=1`. Upstream `PartitionMesh` creates one
# partition Surface child per (surface, partition) pair, one partition
# Curve per connected run of cross-partition border edges and partition
# Points where the dim-1 partition-entity network changes; a declared
# Curve already separating partitions carries the interface role through
# its union partition vector instead of duplicating an entity.
using Test, Tessella
const A=Tessella.API
const E=Tessella.Elements

# Unit square pair sharing the x=1 edge; node/Point tag i sits at vertex i.
const CORNERS=[0. 1. 0. 1. 2. 2.;0. 0. 1. 1. 0. 1.;0. 0. 0. 0. 0. 0.]
# Curve tag, stored Line2 cell direction, declared signed Point bounds.
const CURVES=[(11,(1,2),[1,-2]),(12,(2,4),[2,-4]),(13,(4,3),[4,-3]),
              (14,(3,1),[3,-1]),(15,(2,5),[2,-5]),(16,(5,6),[5,-6]),
              (17,(6,4),[6,-4])]
# Surface tag, Tri3 columns into TRIANGLES, declared signed Curve boundary.
const TRIANGLES=UInt64[1 1 2 2;2 4 5 6;4 3 6 4]
const SURFACES=[(21,(1,2),[11,12,13,14]),(22,(3,4),[15,16,17,-12])]

function fixture!(;declared=true)
    A.finalize();A.initialize();A.model.add("surface BRep")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    if declared
        for tag in 1:6
            A.model.add_discrete_entity(0,tag)
        end
        for (tag,_,bounds) in CURVES
            A.model.add_discrete_entity(1,tag,bounds)
        end
    end
    for (tag,_,bounds) in SURFACES
        A.model.add_discrete_entity(2,tag,declared ? bounds : Int[])
    end
    A.model.add_physical_group(2,[21,22];tag=300,name="sheet")
    for tag in 1:6
        # Undeclared fixtures own border vertices straight on the Surfaces;
        # declared ones own them on the Point entities like the producer.
        dim,entity=declared ? (0,tag) : (2,tag<=4 ? 21 : 22)
        A.mesh.add_nodes(dim,entity,[tag],CORNERS[:,tag])
    end
    if declared
        uid=UInt64(201)
        for (tag,(a,b),_) in CURVES
            A.mesh.add_elements_by_type(tag,1,[uid],[a,b]);uid+=1
        end
    end
    uid=UInt64(101)
    for (tag,columns,_) in SURFACES
        A.mesh.add_elements_by_type(tag,2,collect(UInt64,uid:uid+length(columns)-1),
            vec(TRIANGLES[:,collect(columns)]));uid+=length(columns)
    end
    declared && for tag in 1:6
        A.mesh.add_elements_by_type(tag,15,[400+tag],[tag])
    end
    A.mesh.generate(2)
    A.option("Mesh.PartitionCreateTopology",1)
    return nothing
end
# Faces T1,T2 of Surface 21 split into partitions 1/2; U1,U2 of Surface 22
# likewise. Every declared Curve/Point cell follows one of its face
# carriers. Cross-partition border edges: the Surface-21 internal edge
# (1,4), the Surface-22 internal edge (2,6) and the shared declared Curve
# 12 on edge (2,4) — the last carried by its declared Curve child.
const TAGS=vcat(UInt64[101,102,103,104],UInt64(201):UInt64(207),UInt64(401):UInt64(406))
const PARTS=Int32[1,2,1,2, 1,1,2,2,1,1,2, 1,1,2,2,1,1]
partition!()=A.mesh.partition(2,TAGS,PARTS)
partition_undeclared!()=A.mesh.partition(2,UInt64[101,102,103,104],Int32[1,2,1,2])
children(dim,parent)=sort!([tag for (_,tag) in A.model.get_entities(dim)
    if A.model.get_parent(dim,tag)==(dim,parent)])
child(dim,parent)=only(children(dim,parent))
function child_at(dim,parent,part)
    matches=[tag for tag in children(dim,parent)
        if Int32.(A.model.get_partitions(dim,tag))==Int32[part]]
    return only(matches)
end
signed_bounds(dim,tag)=Tuple{Int,Int}[pair for pair in
    A.model.get_boundary([(dim,tag)],false,true,false)]
function family(dim,owner,kind)
    types,ids,nodes=A.mesh.get_elements(dim,owner);index=findfirst(==(kind),types)
    index===nothing && return UInt64[],zeros(UInt64,E.msh_num_nodes(kind),0)
    return ids[index],reshape(nodes[index],E.msh_num_nodes(kind),length(ids[index]))
end
partitions_of(dim,tag)=sort!(Int32.(A.model.get_partitions(dim,tag)))
snapshot()=(A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements(),
    A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag(),
    A.model.get_number_of_partitions())
function interface_children()
    # Partition Curves whose parent is a Surface child — the emitted
    # interface entities (declared Curve children keep their Curve parent;
    # unparented originals report (-1,-1)).
    return sort!([tag for (_,tag) in A.model.get_entities(1)
        if A.model.get_parent(1,tag)[1]==2])
end
point_children()=sort!([tag for (_,tag) in A.model.get_entities(0)
    if A.model.get_parent(0,tag)[1]>=0])

@testset "shared-curve surface pair partition topology" begin
    fixture!();partition!()
    @test A.model.get_number_of_partitions()==2
    # One dim-2 child per (surface, partition); parent/membership retained.
    s11,s12,s21_,s22_=child_at(2,21,1),child_at(2,21,2),child_at(2,22,1),child_at(2,22,2)
    @test partitions_of(2,s11)==Int32[1];@test partitions_of(2,s12)==Int32[2]
    @test partitions_of(2,s21_)==Int32[1];@test partitions_of(2,s22_)==Int32[2]
    @test A.model.get_parent(2,s11)==(2,21);@test A.model.get_parent(2,s12)==(2,21)
    @test A.model.get_parent(2,s21_)==(2,22);@test A.model.get_parent(2,s22_)==(2,22)
    @test Set(A.model.get_entities_for_physical_group(2,300))==
        Set([21,22,s11,s12,s21_,s22_])
    # Declared Curve children: the shared Curve 12 spans both partitions —
    # its child is the coincident interface record — the others stay single.
    c1,c2,c3,c4,c5,c6,c7=(child(1,c) for c in (11,12,13,14,15,16,17))
    @test partitions_of(1,c1)==Int32[1];@test partitions_of(1,c3)==Int32[2]
    @test partitions_of(1,c4)==Int32[2];@test partitions_of(1,c5)==Int32[1]
    @test partitions_of(1,c6)==Int32[1];@test partitions_of(1,c7)==Int32[2]
    @test partitions_of(1,c2)==Int32[1,2]
    @test A.model.get_parent(1,c2)==(1,12)
    # Declared Point bounds transport to the Point children.
    p1,p2,p3,p4,p5,p6=(child(0,p) for p in 1:6)
    @test signed_bounds(1,c2)==[(0,p2),(0,-p4)]
    @test signed_bounds(1,c1)==[(0,p1),(0,-p2)]
    # Two emitted interface Curves — the internal edges (1,4) of Surface 21
    # and (2,6) of Surface 22 — each parented on the source Surface of its
    # lower-partition carrier and carrying both partition indices.
    ia,ib=interface_children()
    @test length([ia,ib])==2
    pa,pb=A.model.get_parent(1,ia),A.model.get_parent(1,ib)
    @test sort!([pa[2],pb[2]])==[21,22] && pa[1]==2 && pb[1]==2
    @test partitions_of(1,ia)==Int32[1,2] && partitions_of(1,ib)==Int32[1,2]
    # Junctions reuse the declared Point children — the only dim-0
    # partition entities are the six curve-endpoint children.
    @test point_children()==[p1,p2,p3,p4,p5,p6]
    @test all(tag->A.model.get_parent(0,tag)[1]==0,point_children())
    @test partitions_of(0,p1)==Int32[1,2];@test partitions_of(0,p4)==Int32[1,2]
    @test partitions_of(0,p3)==Int32[2];@test partitions_of(0,p5)==Int32[1]
    # Interface entity bounds are its junction endpoints signed by stored
    # traversal: ia stores (4,1) (the lower-partition carrier T1's edge),
    # so its start binds +p4 and its end −p1; likewise ib stores (6,2).
    @test sort!(signed_bounds(1,ia))==sort!([(0,-p1),(0,p4)])
    @test sort!(signed_bounds(1,ib))==sort!([(0,-p2),(0,p6)])
    # Surface children bind their lower entities with actual-traversal
    # signs: the two sides of the shared Curve carry opposite signs.
    @test sort!(signed_bounds(2,s11))==sort!([(1,c1),(1,c2),(1,ia)])
    @test sort!(signed_bounds(2,s12))==sort!([(1,c3),(1,c4),(1,-ia)])
    @test sort!(signed_bounds(2,s21_))==sort!([(1,c5),(1,c6),(1,ib)])
    @test sort!(signed_bounds(2,s22_))==sort!([(1,-c2),(1,c7),(1,-ib)])
    # Cell ownership: original cells sit on their partition children.
    @test family(2,s11,2)[1]==UInt64[101]
    @test family(2,s12,2)[1]==UInt64[102]
    @test family(2,s21_,2)[1]==UInt64[103]
    @test family(2,s22_,2)[1]==UInt64[104]
    @test family(1,c2,1)[1]==UInt64[202]
    # The interface cells carry the stored (lower-partition) traversal.
    _,nodes=family(1,ia,1);@test Tuple(vec(nodes))==(4,1)
    _,nodesb=family(1,ib,1);@test Tuple(vec(nodesb))==(6,2)
    # Two emitted interface cells are the only additions to the 17 originals.
    @test sum(size(block.nodes,2) for block in A.mesh.get().blocks)==17+2
    @test Tessella.validate(A.mesh.get()).ok
    A.mesh.unpartition()
    @test A.model.get_entities()==sort!([(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),
        (1,11),(1,12),(1,13),(1,14),(1,15),(1,16),(1,17),(2,21),(2,22)])
    @test A.model.get_number_of_partitions()==0
end

@testset "surface partition topology disabled stays bare" begin
    # The producer still runs when ghosts are requested; with
    # PartitionCreateTopology=0 it emits ownership buckets only — no
    # interface entities and no boundary rows, matching the ghost-only
    # mode the volume-tier tests document for Gmsh 4.15.2.
    fixture!();A.option("Mesh.PartitionCreateTopology",0)
    A.option("Mesh.PartitionCreateGhostCells",1);partition!()
    @test A.model.get_number_of_partitions()==2
    @test isempty(interface_children())
    s11=child_at(2,21,1)
    @test isempty(signed_bounds(2,s11))
    @test isempty(signed_bounds(1,child(1,12)))
    @test partitions_of(2,s11)==Int32[1]
    @test partitions_of(1,child(1,12))==Int32[1,2]
    @test !isempty(A.mesh.get().partition_data.ghost_elements)
    @test Tessella.validate(A.mesh.get()).ok
    A.mesh.unpartition()
end

@testset "surface partition ghost cells" begin
    fixture!();A.option("Mesh.PartitionCreateGhostCells",1);partition!()
    owned=A.mesh.get()
    # Each of the four faces borders the other partition — all ghost.
    @test length(owned.partition_data.ghost_elements)==4
    @test all(row->length(row.ghost_partitions)==1,owned.partition_data.ghost_elements)
    # Ghost entities materialize at the top dimension, one per partition.
    @test length(owned.partition_data.ghost_entities)==2
    @test sort!([part for (_,part) in owned.partition_data.ghost_entities])==Int32[1,2]
    A.mesh.unpartition()
end

@testset "undeclared shared border surfaces partition" begin
    fixture!(declared=false);partition_undeclared!()
    @test A.model.get_number_of_partitions()==2
    # All three cross-partition edges emit interface Curves — the shared
    # border edge (2,4) included, since no declared Curve carries it.
    interfaces=interface_children()
    @test length(interfaces)==3
    for tag in interfaces
        @test partitions_of(1,tag)==Int32[1,2]
        @test A.model.get_parent(1,tag)[1]==2
        _,nodes=family(1,tag,1);@test size(nodes,1)==2
    end
    # Every interface endpoint shares a vertex with another dim-1 entity or
    # terminates a chain: all four chain ends earn a partition Point on the
    # existing node.
    points=point_children()
    @test length(points)==4
    for tag in points
        @test partitions_of(0,tag)==Int32[1,2]
        _,nodes=family(0,tag,15);@test size(nodes,2)==1
    end
    for tag in interfaces
        @test length(signed_bounds(1,tag))==2
    end
    @test sum(size(block.nodes,2) for block in A.mesh.get().blocks)==4+3+4
    @test Tessella.validate(A.mesh.get()).ok
    A.mesh.unpartition()
end

@testset "surface partition atomic rejections" begin
    # A Curve cell assigned a partition no carrier face has is rejected.
    fixture!();before=snapshot()
    parts=copy(PARTS);parts[5]=Int32(2)
    @test_throws ArgumentError A.mesh.partition(2,TAGS,parts)
    @test snapshot()==before

    # A nonmanifold edge (three faces on one edge) is rejected.
    A.finalize();A.initialize();A.model.add("bad")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    A.model.add_discrete_entity(2,31)
    A.mesh.add_nodes(2,31,collect(UInt64,1:5),
        vec([0. 1. 0. 0. 0.;0. 0. 1. 0. 0.;0. 0. 0. 1. -1.]))
    A.mesh.add_elements_by_type(31,2,[501,502,503],vec(UInt64[1 1 1;2 2 2;3 4 5]))
    A.mesh.generate(2);A.option("Mesh.PartitionCreateTopology",1)
    before=snapshot()
    @test_throws ArgumentError A.mesh.partition(2,UInt64[501,502,503],Int32[1,2,1])
    @test snapshot()==before

    # Faces sharing an edge in the same traversal direction are rejected.
    A.finalize();A.initialize();A.model.add("bad2")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    A.model.add_discrete_entity(2,31);A.model.add_discrete_entity(2,32)
    A.mesh.add_nodes(2,31,collect(UInt64,1:4),
        vec([0. 1. 0. 1.;0. 0. 1. 1.;0. 0. 0. 0.]))
    A.mesh.add_elements_by_type(31,2,[501],UInt64[1,2,4])
    A.mesh.add_elements_by_type(32,2,[502],UInt64[1,2,3])
    A.mesh.generate(2);A.option("Mesh.PartitionCreateTopology",1)
    before=snapshot()
    @test_throws ArgumentError A.mesh.partition(2,UInt64[501,502],Int32[1,2])
    @test snapshot()==before

    # A mesh with no Surface family still reports the original blocker.
    A.finalize();A.initialize();A.model.add("lines only")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    A.model.add_discrete_entity(0,1);A.model.add_discrete_entity(0,2)
    A.model.add_discrete_entity(1,11,[1,-2])
    A.mesh.add_nodes(0,1,[1],[0.,0.,0.]);A.mesh.add_nodes(0,2,[2],[1.,0.,0.])
    A.mesh.add_elements_by_type(11,1,[701],[1,2])
    A.mesh.add_elements_by_type(1,15,[801],[1]);A.mesh.add_elements_by_type(2,15,[802],[2])
    A.mesh.generate(2);A.option("Mesh.PartitionCreateTopology",1)
    before=snapshot()
    @test_throws ArgumentError A.mesh.partition(2,UInt64[701,801,802],Int32[1,1,2])
    @test snapshot()==before
end
end
