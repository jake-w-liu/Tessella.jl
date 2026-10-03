using Test
using Tessella
using Tessella.MeshTypes: nnodes
using Tessella.Elements: ElementBlock, MixedMesh

function _mesh_identity_execute(source)
    mktemp() do path,io
        write(io,source);close(io)
        execute_geo(path)
    end
end

# Small, valid carrier cells isolate the ownership propagation from mesh
# generation. The embedded surface also belongs to a volume, so recovered
# curve nodes must retain one identity in segments, triangles and tetrahedra.
function _mesh_identity_common_carrier(;extra_point=false)
    model=GeoModel()
    add_box!(model,-1.,-2.,-1.,3.,4.,2.;tag=1)
    for (tag,x,y) in ((101,0.,0.),(102,1.,0.),(103,0.,-1.),
                      (104,1.,-1.),(105,1.,1.),(106,0.,1.),(110,0.5,0.))
        add_point!(model,x,y,0.;tag=tag)
    end
    add_line!(model,101,102;tag=101)
    for (tag,a,b) in ((102,103,104),(103,104,105),(104,105,106),(105,106,103))
        add_line!(model,a,b;tag=tag)
    end
    add_curve_loop!(model,[102,103,104,105];tag=101)
    add_plane_surface!(model,[101];tag=101)
    model.embeds[(2,101)]=[(0,110),(1,101)]
    model.embeds[(3,1)]=[(2,101)]
    model.curve_params[101]=[0.,0.25,0.5,1.]
    if extra_point
        Tessella.Model.add_discrete_nodes!(model,0,110,[1001,1002],
                                         [0.5,0,0,0.25,0,0])
    end
    curve=Mesh(Float64[0 0.25 0.5 1;0 0 0 0;0 0 0 0];
               segs=Int32[1 2 3;2 3 4])
    carrier=Mesh(Float64[0 0.25 0.5 1 0 0.25 0.5 1;
                        0 0 0 0 1 1 1 1;
                        0 0 0 0 0 0 0 0];
                 tris=Int32[1 1 2 2 3 3;2 6 3 7 4 8;6 5 7 6 8 7])
    volume=Mesh(hcat(carrier.coords,[0.5,0.5,0.75]);
                tets=vcat(carrier.tris,fill(Int32(9),1,6)))
    parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[
        (1,101,curve),(2,101,carrier),(3,1,volume)]
    append!(parts,Tessella.Model._model_point_mesh_parts(model,"test"))
    return model,parts
end

@testset "Recovered common-carrier point incidence" begin
    for extra_point in (false,true)
        model,parts=_mesh_identity_common_carrier(;extra_point)
        owners=Tessella.GeoExec._geo_mesh_part_node_entities(model,parts,"test")
        @test owners[1][1]==(0,101)
        @test owners[1][end]==(0,102)
        @test owners[1][3]==owners[2][3]==owners[3][3]==(0,110)
        @test owners[1][2]==owners[2][2]==owners[3][2]==(1,101)
        mesh,node_owners=Tessella.GeoExec._geo_merge_entity_meshes(parts;node_entities=owners)
        midpoint=findall(n->Tuple(mesh.coords[:,n])==(0.5,0.,0.),1:nnodes(mesh))
        quarter=findall(n->Tuple(mesh.coords[:,n])==(0.25,0.,0.),1:nnodes(mesh))
        @test length(midpoint)==1
        @test length(quarter)==(extra_point ? 2 : 1)
        node=only(midpoint)
        @test node_owners[node]==(0,110)
        @test node in mesh.segs && node in mesh.tris && node in mesh.tets
        @test validate(mesh).ok
    end
    # An orphan Point and a retained raw Curve have independent identities,
    # even at exactly matching coordinates in the same geometry.
    model,parts=_mesh_identity_common_carrier()
    model.embeds[(2,101)]=[(1,101)]
    owners=Tessella.GeoExec._geo_mesh_part_node_entities(model,parts,"test")
    @test owners[1][3]==owners[2][3]==owners[3][3]==(1,101)
    mesh,_=Tessella.GeoExec._geo_merge_entity_meshes(parts;node_entities=owners)
    @test count(n->Tuple(mesh.coords[:,n])==(0.5,0.,0.),1:nnodes(mesh))==2
    model,parts=_mesh_identity_common_carrier()
    empty!(model.curve_params)
    Tessella.Model.add_discrete_nodes!(model,1,101,[1001,1002,1003,1004],
        [0.,0,0,0.25,0,0,0.5,0,0,1.,0,0])
    owners=Tessella.GeoExec._geo_mesh_part_node_entities(model,parts,"test")
    @test owners[1][3]==(1,101)
    @test owners[2][3]==owners[3][3]==(0,110)
    # Unused carrier vertices cannot establish recovered curve incidence.
    model,parts=_mesh_identity_common_carrier()
    carrier=parts[2][3]
    parts[2]=(2,101,Mesh(carrier.coords;tris=carrier.tris[:,1:2]))
    owners=Tessella.GeoExec._geo_mesh_part_node_entities(model,parts,"test")
    @test owners[1][3]==(1,101)
    # Coincident competing embedded Points require an explicit diagnostic.
    model,parts=_mesh_identity_common_carrier()
    add_point!(model,0.5,0,0;tag=111)
    push!(model.embeds[(2,101)],(0,111))
    @test_throws r"preserving their separate constraint identities is not implemented" Tessella.GeoExec._geo_mesh_part_node_entities(model,parts,"test")
    model,parts=_mesh_identity_common_carrier()
    empty!(model.embeds)
    @test Tessella.GeoExec._geo_mesh_constraint_point_owners(model,parts,"test")===nothing
end

function _mesh_identity_discrete_intervals(count=2)
    model=GeoModel()
    for interval in 1:count
        first_point=2interval-1;last_point=2interval
        add_point!(model,0.,0,0;tag=first_point)
        add_point!(model,1.,0,0;tag=last_point)
        Tessella.Model.add_discrete_entity!(model,1,interval,
            [(0,first_point),(0,last_point)])
        node_offset=10interval
        Tessella.Model.add_discrete_nodes!(model,1,interval,
            node_offset.+[1,2,3],[0.5,0,0,1,0,0,0,0,0])
        Tessella.Model.add_discrete_elements!(model,1,interval,[1],
            [100interval.+[1,2]],[node_offset.+[3,1,1,2]])
    end
    return model
end

function _mesh_identity_homology_parts(model;curve_parts=nothing)
    parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[]
    append!(parts,curve_parts===nothing ?
        Tessella.Model._model_curve_mesh_parts(model,"test") : curve_parts)
    append!(parts,Tessella.Model._model_point_mesh_parts(model,"test"))
    owners=Tessella.GeoExec._geo_mesh_part_node_entities(model,parts,"test")
    mesh,node_owners=Tessella.GeoExec._geo_merge_entity_meshes(parts;node_entities=owners)
    cells=Tessella.GeoExec._geo_homology_cells(model,mesh,parts,"test";
        node_owner=node_owners,part_node_entities=owners)
    return mesh,node_owners,cells
end

@testset "Discrete coincident intervals retain independent homology" begin
    model=_mesh_identity_discrete_intervals()
    mesh,owners,cells=_mesh_identity_homology_parts(model)
    @test nnodes(mesh)==10
    for interval in 1:2
        @test length(cells[(1,interval)])==2
        expected=Set(((1,interval),))
        for (_,nodes) in cells[(1,interval)]
            @test Set(owners[Int(node)] for node in nodes) ⊆ expected
        end
        Tessella.Model.add_physical_group!(model,1,[interval];tag=interval+6)
        Tessella.Model.add_homology_request!(model;kind="Homology",
            domain_tags=[interval+6],dims=[0])
    end
    Tessella.Model.compute_homology!(model,cells)
    @test Set(values(model.physical_names))==Set(["H_0{7}1","H_0{8}1"])
end

@testset "Homology uses the emitted discrete curve mesh after remeshing" begin
    model=_mesh_identity_discrete_intervals(1)
    # The current emitted mesh has three subdivisions; the stored record
    # still contains its former two-element chain at the old midpoint.
    coords=Float64[0 1/3 2/3 1;0 0 0 0;0 0 0 0]
    current=Mesh(coords;segs=Int32[1 2 3;2 3 4])
    model.curve_params[1]=[0.,1/3,2/3,1.]
    mesh,owners,cells=_mesh_identity_homology_parts(model;
        curve_parts=[(1,1,current)])
    @test nnodes(mesh)==4
    @test length(cells[(1,1)])==3
    @test all(cell->all(node->1<=node<=nnodes(mesh),cell[2]),cells[(1,1)])
    Tessella.Model.add_physical_group!(model,1,[1];tag=7)
    Tessella.Model.add_homology_request!(model;kind="Homology",domain_tags=[7],dims=[0,1])
    Tessella.Model.compute_homology!(model,cells)
    @test Set(values(model.physical_names))==Set(["H_0{7}1"])
end

@testset "Homology resolves curve connectivity through shared point records" begin
    model=GeoModel()
    for (point,x) in ((1,0.),(2,1.))
        Tessella.Model.add_discrete_entity!(model,0,point)
        Tessella.Model.add_discrete_nodes!(model,0,point,[10+point],[x,0,0])
    end
    Tessella.Model.add_discrete_entity!(model,1,1,[(0,1),(0,2)])
    # Imported record connectivity can reference nodes owned by lower entities.
    # The public constructor currently requires local node classification.
    record=model.discrete[(1,1)]
    push!(record.element_types,Int32(1));push!(record.element_tags,Int32(101))
    push!(record.element_nodes,Int32[11,12])
    mesh,owners,cells=_mesh_identity_homology_parts(model;curve_parts=[])
    @test nnodes(mesh)==2
    @test length(cells[(1,1)])==1
    @test Set(owners[Int(node)] for node in cells[(1,1)][1][2])==Set([(0,1),(0,2)])
    Tessella.Model.add_physical_group!(model,1,[1];tag=7)
    Tessella.Model.add_homology_request!(model;kind="Homology",domain_tags=[7],dims=[0])
    Tessella.Model.compute_homology!(model,cells)
    @test Set(values(model.physical_names))==Set(["H_0{7}1"])
end

@testset "Discrete point records retain every distinct homology node" begin
    model=GeoModel()
    Tessella.Model.add_discrete_entity!(model,0,1)
    Tessella.Model.add_discrete_nodes!(model,0,1,[11,12],[0.,0,0,1,0,0])
    Tessella.Model.add_discrete_elements!(model,0,1,[15],[[101,102]],[[11,12]])
    mesh,owners,cells=_mesh_identity_homology_parts(model;curve_parts=[])
    @test nnodes(mesh)==2
    @test Set(values(owners))==Set([(0,1)])
    @test Set(only(nodes) for (_,nodes) in cells[(0,1)])==Set(Int32[1,2])
    Tessella.Model.add_physical_group!(model,0,[1];tag=7)
    Tessella.Model.add_homology_request!(model;kind="Homology",domain_tags=[7],dims=[0])
    Tessella.Model.compute_homology!(model,cells)
    @test Set(values(model.physical_names))==Set(["H_0{7}1","H_0{7}2"])
end

@testset "Point homology follows actual point cells rather than all stored vertices" begin
    for native in (true,false), explicit in (:none,:first,:both)
        model=GeoModel()
        native ? add_point!(model,0.,0,0;tag=1) :
            Tessella.Model.add_discrete_entity!(model,0,1)
        Tessella.Model.add_discrete_nodes!(model,0,1,[11,12],[0.,0,0,2.,0,0])
        if explicit===:first
            Tessella.Model.add_discrete_elements!(model,0,1,[15],[[101]],[[11]])
        elseif explicit===:both
            Tessella.Model.add_discrete_elements!(model,0,1,[15],[[101,102]],[[11,12]])
        end
        mesh,owners,cells=_mesh_identity_homology_parts(model;curve_parts=[])
        @test nnodes(mesh)==2
        @test length(cells[(0,1)])==(explicit===:both ? 2 : 1)
        cell_coords=[Tuple(mesh.coords[:,only(nodes)]) for (_,nodes) in cells[(0,1)]]
        expected=explicit===:none ? [(2.,0.,0.)] :
            explicit===:first ? [(0.,0.,0.)] : [(0.,0.,0.),(2.,0.,0.)]
        @test sort(cell_coords)==sort(expected)
        Tessella.Model.add_physical_group!(model,0,[1];tag=7)
        Tessella.Model.add_homology_request!(model;kind="Homology",domain_tags=[7],dims=[0])
        Tessella.Model.compute_homology!(model,cells)
        expected_names=explicit===:both ? ["H_0{7}1","H_0{7}2"] : ["H_0{7}1"]
        @test Set(values(model.physical_names))==Set(expected_names)
    end
end

@testset "geometric point mesh identity" begin
    # Independent Gmsh 4.15.2 fixtures: an orphan point is a separate
    # GVertex/MVertex, including at an existing edge/face mesh position.
    points=_mesh_identity_execute("""
        Point(1)={0,0,0,1};Point(2)={0,0,0,1};Mesh 1;
        """)
    @test nnodes(points.mesh)==2
    @test points.mesh.coords[:,1]==points.mesh.coords[:,2]
    line=_mesh_identity_execute("""
        Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={0.5,0,0,1};
        Line(1)={1,2};Transfinite Curve{1}=3;Mesh 1;
        """)
    @test nnodes(line.mesh)==4
    @test count(j->Tuple(line.mesh.coords[:,j])==(0.5,0.0,0.0),
                1:nnodes(line.mesh))==2
    @test length(Set(line.mesh.segs))==3
    square=_mesh_identity_execute("""
        Point(1)={0,0,0,1};Point(2)={1,0,0,1};
        Point(3)={1,1,0,1};Point(4)={0,1,0,1};Point(5)={0.5,0.5,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Transfinite Curve{:}=3;Transfinite Surface{1};
        Recombine Surface{1};Mesh 2;
        """)
    @test square.mesh isa MixedMesh
    @test nnodes(square.mesh)==10
    @test count(j->Tuple(square.mesh.coords[:,j])==(0.5,0.5,0.0),
                1:nnodes(square.mesh))==2
    quads=only(filter(b->b.msh==3,square.mesh.blocks))
    @test length(Set(quads.nodes))==9
    @test validate(square.mesh).ok
    # An interval and a distinct orphan vertex have two connected
    # components, even when the vertex lies at the interval's midpoint.
    homology=_mesh_identity_execute("""
        Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={0.5,0,0,1};
        Line(1)={1,2};Transfinite Curve{1}=3;
        Physical Point(7)={3};Physical Curve(8)={1};
        Homology(0){{7,8},{}};Mesh 1;
        """)
    @test Set(values(homology.model.physical_names))==
        Set(["H_0{7,8}1","H_0{7,8}2"])
    separate=_mesh_identity_execute("""
        Point(1)={0,0,0,1};Point(2)={1,0,0,1};
        Point(3)={0,0,0,1};Point(4)={-1,0,0,1};
        Line(1)={1,2};Line(2)={3,4};Transfinite Curve{:}=3;
        Physical Curve(8)={1,2};Homology(0){{8},{}};Mesh 1;
        """)
    @test nnodes(separate.mesh)==6
    @test isempty(intersect(Set(separate.mesh.segs[:,1:2]),Set(separate.mesh.segs[:,3:4])))
    @test Set(values(separate.model.physical_names))==Set(["H_0{8}1","H_0{8}2"])
    source=""
    for offset in (0,4)
        for (index,(x,y)) in enumerate(((0,0),(1,0),(1,1),(0,1)))
            source*="Point($(offset+index))={$x,$y,0,1};\n"
        end
        for (index,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
            source*="Line($(offset+index))={$(offset+a),$(offset+b)};\n"
        end
        tag=offset==0 ? 1 : 2
        source*="Curve Loop($tag)={$(offset+1),$(offset+2),$(offset+3),$(offset+4)};Plane Surface($tag)={$tag};\n"
    end
    coincident=_mesh_identity_execute(source*"""
        Transfinite Curve{:}=3;Transfinite Surface{:};Recombine Surface{:};
        Physical Surface(8)={1,2};Homology(0){{8},{}};Mesh 2;
        """)
    @test nnodes(coincident.mesh)==18
    @test sum(size(b.nodes,2) for b in coincident.mesh.blocks if b.msh==3)==8
    @test validate(coincident.mesh).ok
    @test Set(values(coincident.model.physical_names))==Set(["H_0{8}1","H_0{8}2"])
end
