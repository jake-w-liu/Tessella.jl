module PartitionTailFunctionalChecks
# Partition-parity tail checks against Gmsh 4.15.2 public behaviour:
#
# * `Mesh.PartitionCreateTopology=0` keeps partition children as bare
#   ownership buckets — upstream `$PartitionedEntities` rows carry no
#   boundary entries for them — while `Mesh.PartitionCreateGhostCells=1`
#   still emits the ghost layer and physical memberships transport to
#   children regardless of the topology gate.
# * A single Region holding a Prism6 and a Hex8 across a conformal Quad4
#   facet partitions with a quad interface and mixed Tri3/Quad4 exterior
#   surfaces (mixed-element coverage beyond the Tet4 fixtures).
# * The `.geo` `PartitionMesh` frontend reaches the same producer: an
#   imported `$Periodic` link stays bound to its original slave/master
#   entities, and the ghost-only topology gate suppresses child boundaries
#   on that path as well.
using Test, Tessella
const A=Tessella.API
const E=Tessella.Elements

# Two Tet4 cells sharing one facet, with a fully declared
# Surface/Curve/Point BRep — exercises every boundary-bearing child kind.
const T4_CORNERS=[0. 1. 0. 0. 0.;0. 0. 1. 0. 0.;0. 0. 0. 1. -1.]
const T4_TETS=UInt64[11 11;12 13;13 12;14 15]
const T4_EXTERIOR=UInt64[11 12 13 11 13 12;12 13 11 13 12 11;14 14 14 15 15 15]
const T4_CURVES=[(41,[(11,12)],[1,-2]),(42,[(12,13)],[2,-3]),(43,[(13,11)],[3,-1]),
    (44,[(11,14)],[1,-4]),(45,[(12,14)],[2,-4]),(46,[(13,14)],[3,-4]),
    (47,[(11,15)],[1,-5]),(48,[(12,15)],[2,-5]),(49,[(13,15)],[3,-5])]
const T4_SURFACES=[(31,1,[41,45,-44]),(32,2,[42,46,-45]),(33,3,[43,44,-46]),
    (34,4,[-43,49,-47]),(35,5,[-42,48,-49]),(36,6,[-41,47,-48])]

function tet_brep_fixture!(topology,ghosts)
    A.finalize();A.initialize();A.model.add("partition tail tets")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    for tag in 1:5
        A.model.add_discrete_entity(0,tag)
    end
    for (tag,_,bounds) in T4_CURVES
        A.model.add_discrete_entity(1,tag,bounds)
    end
    for (tag,_,bounds) in T4_SURFACES
        A.model.add_discrete_entity(2,tag,bounds)
    end
    A.model.add_discrete_entity(3,11,[31,32,33,34,35,36])
    A.model.add_physical_group(1,[41,44];tag=210,name="weld")
    A.model.add_physical_group(0,[1];tag=110,name="ports")
    for tag in 1:5
        A.mesh.add_nodes(0,tag,[10+tag],T4_CORNERS[:,tag])
    end
    for (tag,column,_) in T4_SURFACES
        A.mesh.add_elements_by_type(tag,2,[630+column],T4_EXTERIOR[:,column])
    end
    uid=701
    for (tag,cells,_) in T4_CURVES
        for (a,b) in cells
            A.mesh.add_elements_by_type(tag,1,[uid],[a,b]);uid+=1
        end
    end
    for tag in 1:5
        A.mesh.add_elements_by_type(tag,15,[800+tag],[10+tag])
    end
    A.mesh.add_elements_by_type(11,4,[501,502],vec(T4_TETS))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",topology)
    A.option("Mesh.PartitionCreateGhostCells",ghosts)
end

partition_two_tets!()=A.mesh.partition(2,
    [631:636...,701:709...,801:805...,501,502],
    vcat(Int32[1,1,1,2,2,2],fill(Int32(1),9),fill(Int32(1),5),Int32[1,2]))

# One Region holding a Prism6 (nodes a b c / d e f) and a Hex8 on the
# mirrored side of its a-b-e-d quad facet. Node uids 11..20 run
# a,b,c,d,e,f,g,h,i,j; the Hex8 order [a,d,i,h, b,e,j,g] keeps the shared
# facet winding opposite to the Prism's canonical (1,2,5,4) face.
const MIX_COORDS=[0. 1. 0. 0. 0. 0. 1. 0. 0. 1.;
                  0. 0. 1. 0. 0. 1. -1. -1. -1. -1.;
                  0. 0. 0. 1. 1. 1. 0. 0. 1. 1.]
const MIX_PRISM=UInt64[11,12,13,14,15,16]
const MIX_HEX=UInt64[11,14,19,18,12,15,20,17]
# Exterior facets in outward canonical order; the a-b-e-d quad between the
# two cells stays undeclared so it surfaces only as the partition interface.
const MIX_SURFACES=[(31,2,[11,13,12]),(32,2,[14,15,16]),
    (33,3,[11,14,16,13]),(34,3,[12,13,16,15]),
    (35,3,[11,18,19,14]),(36,3,[11,12,17,18]),(37,3,[14,19,20,15]),
    (38,3,[19,18,17,20]),(39,3,[12,15,20,17])]
# First surface owning each corner node keeps the realistic classification.
const MIX_NODE_SURFACES=[31,31,31,32,32,32,36,36,38,38]

function mixed_fixture!(topology,ghosts;physical=true)
    A.finalize();A.initialize();A.model.add("partition tail mixed")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    for (tag,_,_) in MIX_SURFACES
        A.model.add_discrete_entity(2,tag)
    end
    A.model.add_discrete_entity(3,11,[surface for (surface,_,_) in MIX_SURFACES])
    physical && A.model.add_physical_group(3,[11];tag=310,name="bulk")
    physical && A.model.add_physical_group(2,[31,36];tag=220,name="caps")
    for (slot,tag) in enumerate(MIX_NODE_SURFACES)
        A.mesh.add_nodes(2,tag,[10+slot],MIX_COORDS[:,slot])
    end
    uid=631
    for (tag,msh,nodes) in MIX_SURFACES
        A.mesh.add_elements_by_type(tag,msh,[uid],nodes);uid+=1
    end
    A.mesh.add_elements_by_type(11,6,[601],MIX_PRISM)
    A.mesh.add_elements_by_type(11,5,[602],MIX_HEX)
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",topology)
    A.option("Mesh.PartitionCreateGhostCells",ghosts)
end

# Surface cells follow their single Tet4/Hex8 facet carrier; the interface
# facet's two sides take the owning cell partitions 1 (prism) and 2 (hex).
partition_mixed!()=A.mesh.partition(2,[631:639...,601,602],
    vcat(fill(Int32(1),4),fill(Int32(2),5),Int32[1,2]))

support(nodes)=Tuple(sort(collect(nodes)))

@testset "ghost-only lower-BRep partition keeps children boundary-free" begin
    try
        tet_brep_fixture!(0,1);partition_two_tets!()
        @test A.model.get_number_of_partitions()==2
        mesh=A._cached_mesh_locked("test")
        @test mesh.partition_data!==nothing && mesh.partition_data.num_partitions==2
        @test !isempty(mesh.partition_data.ghost_elements)
        entities=mesh.partition_data.entities
        # Children materialize for every populated (parent, partition) pair —
        # declared Region, Surface, Curve and Point parents alike — but stay
        # bare ownership buckets: Gmsh 4.15.2 writes no boundary rows for
        # partition entities when PartitionCreateTopology=0.
        dims=Set(d for (d,_) in keys(entities))
        @test dims==Set([0,1,2,3])
        @test all(isempty(entity.boundaries) for entity in values(entities))
        for (dim,tag) in keys(entities)
            @test isempty(A.model.get_boundary([(dim,tag)],false,true,false))
            @test isempty(A.model.get_boundary([(dim,tag)],false,false,false))
        end
        # The declared parents keep their own signed topology.
        @test length(A.model.get_boundary([(3,11)],false,true,false))==6
        @test length(A.model.get_boundary([(1,41)],false,true,false))==2
        # Ghost visibility and physical transport are unaffected by the
        # topology gate.
        weld=sort!([t for ((d,t),e) in entities if d==1 && 210 in e.physical_tags])
        expected=sort!([only(t for ((d,t),e) in entities
            if d==1 && e.parent==(1,parent)) for parent in (41,44)])
        @test weld==expected
        port=only(e for ((d,t),e) in entities if d==0 && e.parent==(0,1))
        @test 110 in port.physical_tags
        @test Tessella.validate(mesh).ok
    finally
        A.finalize()
    end
end

@testset "topology gate still emits declared child topology" begin
    try
        tet_brep_fixture!(1,1);partition_two_tets!()
        mesh=A._cached_mesh_locked("test")
        entities=mesh.partition_data.entities
        @test mesh.partition_data.num_partitions==2
        @test any(!isempty(entity.boundaries) for entity in values(entities))
        @test !isempty(mesh.partition_data.ghost_elements)
        # Region children keep their declared surfaces plus the interface,
        # and each produced curve child keeps its signed point bounds.
        region_parts=Dict(only(e.partitions)=>tag
            for ((d,tag),e) in entities if d==3)
        @test length(region_parts)==2
        for (part,tag) in region_parts
            boundary=A.model.get_boundary([(3,tag)],false,true,false)
            @test length(boundary)==4 && all(b->b[1]==2,boundary)
        end
        curve_child=only(t for ((d,t),e) in entities if d==1 && e.parent==(1,44))
        parent_boundary=A.model.get_boundary([(1,44)],false,true,false)
        child_boundary=A.model.get_boundary([(1,curve_child)],false,true,false)
        # The child keeps the parent's endpoint convention against the
        # produced Point children.
        @test length(child_boundary)==2 && all(b->b[1]==0,child_boundary)
        @test sign.(last.(child_boundary))==sign.(last.(parent_boundary))
        point_children=Dict(e.parent[2]=>t for ((d,t),e) in entities if d==0)
        @test Set(abs(last(b)) for b in child_boundary)==
            Set(point_children[abs(last(b))] for b in parent_boundary)
        @test Tessella.validate(mesh).ok
    finally
        A.finalize()
    end
end

@testset "mixed Prism6+Hex8 partition emits a Quad4 interface" begin
    try
        mixed_fixture!(1,0);partition_mixed!()
        @test A.model.get_number_of_partitions()==2
        mesh=A._cached_mesh_locked("test")
        entities=mesh.partition_data.entities
        # One interface entity: the shared quad facet between the parts.
        interfaces=[(tag,e) for ((d,tag),e) in entities
            if d==2 && e.parent[1]==3]
        @test length(interfaces)==1
        (interface_tag,interface)=only(interfaces)
        @test sort!(collect(interface.partitions))==Int32[1,2]
        types,tags,nodes=A.mesh.get_elements(2,interface_tag)
        @test types==Int32[3]
        @test support(vec(only(nodes)))==support([11,12,14,15])
        # Both Region children bind the interface with opposite signs; which
        # side carries the positive sign follows the first carrier's stored
        # winding, not the partition number.
        region_children=Dict(only(e.partitions)=>tag
            for ((d,tag),e) in entities if d==3)
        @test length(region_children)==2
        b1=A.model.get_boundary([(3,region_children[1])],false,true,false)
        b2=A.model.get_boundary([(3,region_children[2])],false,true,false)
        @test ((2,interface_tag) in b1 && (2,-interface_tag) in b2) ||
              ((2,-interface_tag) in b1 && (2,interface_tag) in b2)
        # Each child keeps its own exterior surfaces — mixed Tri3/Quad4 on
        # the prism side, Quad4 only on the hex side.
        prism_boundary=Set(abs(last(b)) for b in
            A.model.get_boundary([(3,region_children[1])],false,true,false))
        hex_boundary=Set(abs(last(b)) for b in
            A.model.get_boundary([(3,region_children[2])],false,true,false))
        prism_surfaces=Set(e.parent[2] for ((d,t),e) in entities
            if d==2 && e.parent[1]==2 && t in prism_boundary)
        hex_surfaces=Set(e.parent[2] for ((d,t),e) in entities
            if d==2 && e.parent[1]==2 && t in hex_boundary)
        @test prism_surfaces==Set([31,32,33,34])
        @test hex_surfaces==Set([35,36,37,38,39])
        # Physical memberships transport to the children of their parents.
        @test all(310 in entities[(3,tag)].physical_tags
                  for (_,tag) in region_children)
        @test 220 in entities[(2,only(t for ((d,t),e) in entities
            if d==2 && e.parent==(2,31)))].physical_tags
        @test Tessella.validate(mesh).ok
    finally
        A.finalize()
    end
end

@testset ".geo PartitionMesh transports imported periodic links" begin
    mktempdir() do directory
        try
            # Two disconnected Tet4 regions, fully bounded: the deterministic
            # `PartitionMesh` split lands each in its own partition, so the
            # appended $Periodic link crosses the partition seam.
            A.finalize();A.initialize();A.model.add("periodic pair")
            A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
            for tag in 41:44
                A.model.add_discrete_entity(2,tag)
            end
            for tag in 51:54
                A.model.add_discrete_entity(2,tag)
            end
            A.model.add_discrete_entity(3,11,[41,42,43,44])
            A.model.add_discrete_entity(3,22,[51,52,53,54])
            # Periodic-linked nodes classify onto the linked surfaces.
            A.mesh.add_nodes(2,41,[11],[0.,0.,0.])
            A.mesh.add_nodes(2,44,[12,13,14],[1.,0.,0., 0.,1.,0., 0.,0.,1.])
            A.mesh.add_nodes(2,51,[21],[3.,0.,0.])
            A.mesh.add_nodes(2,54,[22,23,24],[4.,3.,3., 0.,1.,0., 0.,0.,1.])
            for (tag,uid,nodes) in ((41,631,[11,13,12]),(42,632,[11,12,14]),
                    (43,633,[11,14,13]),(44,634,[14,12,13]),
                    (51,635,[21,23,22]),(52,636,[21,22,24]),
                    (53,637,[21,24,23]),(54,638,[24,22,23]))
                A.mesh.add_elements_by_type(tag,2,[uid],nodes)
            end
            A.mesh.add_elements_by_type(11,4,[501],[11,12,13,14])
            A.mesh.add_elements_by_type(22,4,[502],[21,22,23,24])
            A.mesh.generate(3)
            source=joinpath(directory,"pair.msh");A.write(source)
            # Slave surface 44 -> master surface 54, translation +3 in x.
            section=["\$Periodic","1","2 44 54","16",
                "1.0 0.0 0.0 3.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0",
                "3","14 24","12 22","13 23","\$EndPeriodic"]
            write(source,join(readlines(source),"\n")*"\n"*join(section,"\n")*"\n")
            out=joinpath(directory,"partitioned.msh")
            geo=joinpath(directory,"case.geo")
            write(geo,"Merge \"$(replace(source,'\\'=>'/'))\";\n" *
                "PartitionMesh(2);\nSave \"$(replace(out,'\\'=>'/'))\";\n")
            A.finalize()
            ex=Tessella.GeoExec.execute_geo(geo)
            mesh=ex.mesh
            @test mesh.partition_data!==nothing &&
                mesh.partition_data.num_partitions==2
            # The periodic-node correspondence survives across partitions
            # bound to the original entities — the Gmsh convention.
            @test length(mesh.periodic_links)==1
            link=only(mesh.periodic_links)
            @test (link.dim,link.slave_entity,link.master_entity)==(2,44,54)
            @test length(link.slave_nodes)==3 &&
                length(link.master_nodes)==3 && link.affine!==nothing
            @test Tessella.validate(mesh).ok
            text=read(out,String)
            @test occursin("\$Periodic",text) && occursin("\$PartitionedEntities",text)
        finally
            A.finalize()
        end
    end
end

@testset ".geo PartitionMesh topology gate suppresses child boundaries" begin
    mktempdir() do directory
        for (topology,ghosts) in ((0,1),(1,0))
            src=joinpath(directory,"case-$(topology)-$(ghosts).geo")
            write(src,"SetFactory(\"OpenCASCADE\");\nBox(1)={0,0,0,1,1,1};\n" *
                "Physical Volume(7) = {1};\nPhysical Surface(8) = {1,2};\n" *
                "Mesh 3;\n" *
                "Mesh.PartitionCreateTopology=$topology;\n" *
                "Mesh.PartitionCreateGhostCells=$ghosts;\n" *
                "PartitionMesh(2);\n")
            ex=Tessella.GeoExec.execute_geo(src)
            mesh=ex.mesh
            entities=mesh.partition_data.entities
            @test mesh.partition_data.num_partitions==2
            region_children=[e for ((d,_),e) in entities if d==3]
            @test length(region_children)==2
            # Physical groups transport to children under either gate.
            @test all(7 in e.physical_tags for e in region_children)
            @test any(8 in e.physical_tags for ((d,t),e) in entities
                if d==2 && e.parent==(2,1))
            if topology==0
                @test all(isempty(e.boundaries) for e in values(entities))
                @test !isempty(mesh.partition_data.ghost_elements)
            else
                @test any(e->!isempty(e.boundaries),
                    (e for ((d,_),e) in entities if d==3))
                @test isempty(mesh.partition_data.ghost_elements)
            end
            @test Tessella.validate(mesh).ok
        end
    end
end
end
