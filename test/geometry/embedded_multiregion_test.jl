# EMB01 residual: multi-Region embedded composition in the native hybrid
# batch planner. Two unit cubes share the x=1 interface Face (Surface 4):
# interior entities land in the strictly-containing Region regardless of how
# many Regions declare them; records on the shared interface retriangulate
# that Face's literal panels into certified conforming sub-facets serving
# both neighbors; uncertifiable incidence (non-node boundary crossing,
# sheet overlap on the interface, outside/exterior records) rejects loudly
# with a specific message. Verified against the Gmsh 4.15.2 oracle: node
# counts and interior-entity cell counts match exactly; on-interface cases
# keep node parity while the native conformal split serves both neighbors
# (upstream splits only the absorbing side into a non-conforming interface).
# Run: julia --project test/geometry/embedded_multiregion_test.jl
module EmbeddedMultiregionChecks
using Test, LinearAlgebra, Tessella
using Tessella.Elements: MixedMesh

# Two unit cubes sharing the x=1 face (Surface 4 is the shared interface).
const BASE = """
Point(1)={0,0,0,0.4}; Point(2)={1,0,0,0.4}; Point(3)={1,1,0,0.4}; Point(4)={0,1,0,0.4};
Point(5)={0,0,1,0.4}; Point(6)={1,0,1,0.4}; Point(7)={1,1,1,0.4}; Point(8)={0,1,1,0.4};
Point(9)={2,0,0,0.4}; Point(10)={2,1,0,0.4}; Point(11)={2,0,1,0.4}; Point(12)={2,1,1,0.4};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Line(5)={5,6}; Line(6)={6,7}; Line(7)={7,8}; Line(8)={8,5};
Line(9)={1,5}; Line(10)={2,6}; Line(11)={3,7}; Line(12)={4,8};
Line(13)={9,10}; Line(14)={10,12}; Line(15)={12,11}; Line(16)={11,9};
Line(21)={2,9}; Line(22)={10,3}; Line(23)={6,11}; Line(24)={12,7};
Line Loop(1)={1,2,3,4}; Line Loop(2)={5,6,7,8};
Line Loop(3)={1,10,-5,-9}; Line Loop(4)={2,11,-6,-10};
Line Loop(5)={3,12,-7,-11}; Line Loop(6)={4,9,-8,-12};
Line Loop(7)={21,13,22,-2}; Line Loop(8)={23,-15,24,-6};
Line Loop(9)={13,14,15,16}; Line Loop(10)={21,-16,-23,-10};
Line Loop(11)={-22,14,24,-11};
Plane Surface(1)={1}; Plane Surface(2)={2}; Plane Surface(3)={3};
Plane Surface(4)={4}; Plane Surface(5)={5}; Plane Surface(6)={6};
Plane Surface(7)={7}; Plane Surface(8)={8}; Plane Surface(9)={9};
Plane Surface(10)={10}; Plane Surface(11)={11};
Surface Loop(1)={1,2,3,4,5,6};
Volume(1)={1};
Surface Loop(2)={7,8,9,10,11,-4};
Volume(2)={2};
Recombine Surface{1,2,3,5,6,7,8,9,10,11};
"""

function run_volume(name,extra)
    path=joinpath(mktempdir(),name*".geo");write(path,BASE*extra)
    result=Tessella.execute_geo(path;mesh_dim=3)
    return result.model.meshing.mesh_authority.mesh
end

# (entity => tet columns) owner map plus the total mesh volume and
# positivity; Pyramid cells do not enter the owner counts but their two fan
# tets contribute to the volume exactly like the kernel fan facets.
function tet_report(mesh)
    data=mesh.entity_data
    owners=Dict{Int,Vector{Int32}}();volume=0.0;positive=true
    for (bi,block) in enumerate(mesh.blocks)
        block.msh in (4,7) || continue
        for ci in axes(block.nodes,2)
            cell=block.nodes[:,ci]
            if block.msh==4
                push!(get!(owners,Int(data.block_entities[bi][ci]),Int32[]),ci)
                p,q,r,s=(mesh.coords[:,i] for i in cell)
                det=LinearAlgebra.det(hcat(q-p,r-p,s-p))
                positive &= det>0; volume+=det/6
            else
                b1,b2,b3,b4,a=(mesh.coords[:,cell[i]] for i in 1:5)
                for (x,y,z) in ((b1,b2,b3),(b1,b3,b4))
                    det=LinearAlgebra.det(hcat(y-x,z-x,a-x))
                    positive &= det>0; volume+=det/6
                end
            end
        end
    end
    return owners,volume,positive
end

# Tets of one Region containing a global node row.
function region_tets_using(mesh,region,row)
    data=mesh.entity_data;out=0
    for (bi,block) in enumerate(mesh.blocks)
        block.msh==4 || continue
        for ci in axes(block.nodes,2)
            Int(data.block_entities[bi][ci])==region && row in block.nodes[:,ci] && (out+=1)
        end
    end
    return out
end

entity_nodes(mesh,dim,tag)=
    Int[row for (row,owner) in enumerate(mesh.entity_data.node_entities)
        if owner==(dim,Int32(tag))]

@testset "multi-Region embedded composition" begin
    @testset "interior Point lands in the containing Region however declared" begin
        for (name,extra) in (
            ("mr_pt_v1","Point(50)={0.5,0.5,0.5,0.3};\nPoint{50} In Volume{1};"),
            ("mr_pt_both","Point(50)={0.5,0.5,0.5,0.3};\nPoint{50} In Volume{1};\nPoint{50} In Volume{2};"),
            ("mr_pt_v2","Point(50)={0.5,0.5,0.5,0.3};\nPoint{50} In Volume{2};"))
            mesh=run_volume(name,extra)
            owners,volume,positive=tet_report(mesh)
            # Oracle m_pt_v1_only/c_pt_both: 449 nodes, 730 + 717 Tets.
            @test size(mesh.coords,2)==449
            @test length(get(owners,1,Int32[]))==730
            @test length(get(owners,2,Int32[]))==717
            @test volume≈2.0 atol=1e-9
            @test positive
            rows=entity_nodes(mesh,0,50)
            @test length(rows)==1 && Tuple(mesh.coords[:,only(rows)])==(0.5,0.5,0.5)
            @test region_tets_using(mesh,1,only(rows))>0
            @test region_tets_using(mesh,2,only(rows))==0
        end
    end
    @testset "Point on the shared interface splits its Face conformally" begin
        for (name,extra) in (
            ("mr_pt_iface","Point(50)={1.0,0.5,0.5,0.3};\nPoint{50} In Volume{1};"),
            ("mr_pt_iface2","Point(50)={1.0,0.5,0.5,0.3};\nPoint{50} In Volume{1};\nPoint{50} In Volume{2};"))
            mesh=run_volume(name,extra)
            owners,volume,positive=tet_report(mesh)
            # Oracle b_pt_shared/d_pt_shared_both: 450 nodes. Native keeps a
            # conforming interface so the node serves both neighbor products.
            @test size(mesh.coords,2)==450
            @test volume≈2.0 atol=1e-9
            @test positive
            rows=entity_nodes(mesh,0,50)
            @test length(rows)==1 && Tuple(mesh.coords[:,only(rows)])==(1.0,0.5,0.5)
            @test region_tets_using(mesh,1,only(rows))>0
            @test region_tets_using(mesh,2,only(rows))>0
        end
    end
    @testset "interior Curve deduplicates across Region declarations" begin
        for (name,extra) in (
            ("mr_cv_v1","Point(50)={0.3,0.5,0.5,0.3}; Point(51)={0.7,0.5,0.5,0.3};\nLine(50)={50,51};\nCurve{50} In Volume{1};"),
            ("mr_cv_both","Point(50)={0.3,0.5,0.5,0.3}; Point(51)={0.7,0.5,0.5,0.3};\nLine(50)={50,51};\nCurve{50} In Volume{1};\nCurve{50} In Volume{2};"))
            mesh=run_volume(name,extra)
            owners,volume,positive=tet_report(mesh)
            # Oracle h_curve_v1_only/h_curve_v1_both: 451 nodes, 733 + 717 Tets.
            @test size(mesh.coords,2)==451
            @test length(get(owners,1,Int32[]))==733
            @test length(get(owners,2,Int32[]))==717
            @test volume≈2.0 atol=1e-9
            @test positive
            rows=sort!(collect(Set(entity_nodes(mesh,0,50) ∪ entity_nodes(mesh,0,51) ∪
                                   entity_nodes(mesh,1,50)));by=row->mesh.coords[1,row])
            @test length(rows)>=2
            @test Tuple(mesh.coords[:,first(rows)])==(0.3,0.5,0.5)
            @test Tuple(mesh.coords[:,last(rows)])==(0.7,0.5,0.5)
            edges=Set{NTuple{2,Int32}}()
            for block in mesh.blocks
                block.msh==4 || continue
                for ci in axes(block.nodes,2)
                    tet=block.nodes[:,ci]
                    for (i,j) in ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
                        push!(edges,minmax(tet[i],tet[j]))
                    end
                end
            end
            @test all(index->minmax(rows[index],rows[index+1]) in edges,1:length(rows)-1)
        end
    end
    @testset "Curve ending on the shared interface is certified" begin
        for (name,extra) in (
            ("mr_cv_touch","Point(50)={0.5,0.5,0.5,0.3}; Point(51)={1.0,0.5,0.5,0.3};\nLine(50)={50,51};\nCurve{50} In Volume{1};"),
            ("mr_cv_touch_both","Point(50)={0.5,0.5,0.5,0.3}; Point(51)={1.0,0.5,0.5,0.3};\nLine(50)={50,51};\nCurve{50} In Volume{1};\nCurve{50} In Volume{2};"))
            mesh=run_volume(name,extra)
            owners,volume,positive=tet_report(mesh)
            # Oracle o_curve_touch_v1: 451 nodes (upstream absorbs into one
            # side; the native conformal split serves both neighbors).
            @test size(mesh.coords,2)==451
            @test volume≈2.0 atol=1e-9
            @test positive
            ends=entity_nodes(mesh,0,51)
            @test length(ends)==1 && Tuple(mesh.coords[:,only(ends)])==(1.0,0.5,0.5)
            @test region_tets_using(mesh,1,only(ends))>0
            @test region_tets_using(mesh,2,only(ends))>0
        end
    end
    @testset "Curve lying on the shared interface serves both Regions" begin
        mesh=run_volume("mr_cv_oniface",
            "Point(50)={1.0,0.3,0.5,0.3}; Point(51)={1.0,0.7,0.5,0.3};\nLine(50)={50,51};\nCurve{50} In Volume{1};")
        owners,volume,positive=tet_report(mesh)
        @test size(mesh.coords,2)==452
        @test volume≈2.0 atol=1e-9
        @test positive
        rows=sort!(collect(Set(entity_nodes(mesh,0,50) ∪ entity_nodes(mesh,0,51) ∪
                               entity_nodes(mesh,1,50)));by=row->mesh.coords[2,row])
        @test length(rows)>=2
        @test all(row->mesh.coords[1,row]==1.0,rows)
        @test all(row->region_tets_using(mesh,1,row)>0 &&
                      region_tets_using(mesh,2,row)>0,rows)
    end
    @testset "interior sheet deduplicates and keeps interior facets" begin
        extra="Point(50)={0.3,0.3,0.5,0.3}; Point(51)={0.7,0.3,0.5,0.3}; Point(52)={0.5,0.7,0.5,0.3};\n"*
              "Line(50)={50,51}; Line(51)={51,52}; Line(52)={52,50};\n"*
              "Line Loop(50)={50,51,52};\nPlane Surface(50)={50};\nSurface{50} In Volume{1};\nSurface{50} In Volume{2};"
        mesh=run_volume("mr_sheet",extra)
        owners,volume,positive=tet_report(mesh)
        @test volume≈2.0 atol=1e-9
        @test positive
        faces=Dict{NTuple{3,Int32},Int}()
        for (bi,block) in enumerate(mesh.blocks)
            block.msh==4 || continue
            for ci in axes(block.nodes,2)
                tet=block.nodes[:,ci]
                for skip in 1:4
                    key=Tuple(sort!(Int32[tet[i] for i in 1:4 if i!=skip]))
                    faces[key]=get(faces,key,0)+1
                end
            end
        end
        data=mesh.entity_data
        sheet=[]
        for (bi,block) in enumerate(mesh.blocks)
            block.msh==2 || continue
            for ci in axes(block.nodes,2)
                Int(data.block_entities[bi][ci])==50 &&
                    push!(sheet,Tuple(sort!(collect(block.nodes[:,ci]))))
            end
        end
        @test !isempty(sheet)
        @test all(face->get(faces,face,0)==2,sheet)
    end
    @testset "In Surface declarations on the shared interface" begin
        # Point{50} In Surface{4}: the shared Face product carries the node,
        # so both neighbor boundaries recover it — oracle
        # k_pt_in_shared_surf: 445 nodes, 715 + 712 Tets.
        mesh=run_volume("mr_s_pt","Point(50)={1.0,0.5,0.5,0.3};\nPoint{50} In Surface{4};")
        owners,volume,positive=tet_report(mesh)
        @test size(mesh.coords,2)==445
        @test length(get(owners,1,Int32[]))==715
        @test length(get(owners,2,Int32[]))==712
        @test volume≈2.0 atol=1e-9
            @test positive
        rows=entity_nodes(mesh,0,50)
        @test length(rows)==1 && Tuple(mesh.coords[:,only(rows)])==(1.0,0.5,0.5)
        @test region_tets_using(mesh,1,only(rows))>0
        @test region_tets_using(mesh,2,only(rows))>0
        # Curve{50} In Surface{4}: endpoints and internal nodes bind the
        # shared Face and serve both Regions.
        mesh=run_volume("mr_s_cv",
            "Point(50)={1.0,0.3,0.5,0.3}; Point(51)={1.0,0.7,0.5,0.3};\nLine(50)={50,51};\nCurve{50} In Surface{4};")
        owners,volume,positive=tet_report(mesh)
        @test size(mesh.coords,2)==452
        @test volume≈2.0 atol=1e-9
        @test positive
        rows=sort!(collect(Set(entity_nodes(mesh,0,50) ∪ entity_nodes(mesh,0,51) ∪
                               entity_nodes(mesh,1,50)));by=row->mesh.coords[2,row])
        @test length(rows)>=2
        @test all(row->mesh.coords[1,row]==1.0,rows)
        @test all(row->region_tets_using(mesh,1,row)>0 &&
                      region_tets_using(mesh,2,row)>0,rows)
        # In Surface + In Volume on the same point deduplicates.
        mesh=run_volume("mr_s_pt_bothv",
            "Point(50)={1.0,0.5,0.5,0.3};\nPoint{50} In Surface{4};\nPoint{50} In Volume{1};")
        owners,volume,positive=tet_report(mesh)
        @test size(mesh.coords,2)==445
        @test length(get(owners,1,Int32[]))==715
        @test length(get(owners,2,Int32[]))==712
    end
    @testset "uncertifiable incidence rejects with a specific message" begin
        crossing="Point(50)={0.7,0.5,0.5,0.3}; Point(51)={1.3,0.5,0.5,0.3};\n"*
                 "Line(50)={50,51};\nCurve{50} In Volume{1};"
        sheet_overlap="Point(50)={1.0,0.3,0.3,0.3}; Point(51)={1.0,0.7,0.3,0.3}; Point(52)={1.0,0.5,0.7,0.3};\n"*
                      "Line(50)={50,51}; Line(51)={51,52}; Line(52)={52,50};\n"*
                      "Line Loop(50)={50,51,52};\nPlane Surface(50)={50};\nSurface{50} In Volume{1};"
        sheet_overlap_both=sheet_overlap*"\nSurface{50} In Volume{2};"
        outside="Point(50)={3.0,0.5,0.5,0.3};\nPoint{50} In Volume{1};"
        exterior="Point(50)={0.0,0.5,0.5,0.3};\nPoint{50} In Volume{1};"
        for (name,extra,fragment) in (
            ("mr_bad_cross",crossing,"non-shared-boundary"),
            ("mr_bad_surf",sheet_overlap,"overlapping facets"),
            ("mr_bad_surf_both",sheet_overlap_both,"overlapping facets"),
            ("mr_bad_out",outside,"outside or on the hybrid batch boundary"),
            ("mr_bad_ext",exterior,"hybrid batch boundary"))
            path=joinpath(mktempdir(),name*".geo");write(path,BASE*extra)
            error=try
                Tessella.execute_geo(path;mesh_dim=3);nothing
            catch e
                sprint(showerror,e)
            end
            @test error!==nothing && occursin(fragment,error)
        end
    end
end
end # module
