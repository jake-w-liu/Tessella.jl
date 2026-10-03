# Gmsh retains orphan GVertex mesh nodes independently of coincident edge
# and face nodes. Compare multiplicity as well as coordinates.
using Tessella
using Tessella.MeshTypes: nnodes
binding=get(ENV,"GMSH_JULIA_API",Sys.iswindows() ?
    joinpath(homedir(),"AppData","Local","Programs","Python","Python312","Lib","gmsh.jl") :
    "/opt/homebrew/lib/gmsh.jl")
isfile(binding) || error("set GMSH_JULIA_API to pinned Gmsh 4.15.2")
include(binding)
gmsh.GMSH_API_VERSION=="4.15.2" || error("Gmsh 4.15.2 required")
sources=[
    "Point(1)={0,0,0,1};Point(2)={0,0,0,1};Mesh 1;",
    "Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={0.5,0,0,1};Line(1)={1,2};Transfinite Curve{1}=3;Mesh 1;",
    "Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={1,1,0,1};Point(4)={0,1,0,1};Point(5)={0.5,0.5,0,1};Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};Transfinite Curve{:}=3;Transfinite Surface{1};Recombine Surface{1};Mesh 2;"
]
push!(sources,"""
    Point(1)={0,0,0,1};Point(2)={1,0,0,1};
    Point(3)={0,0,0,1};Point(4)={-1,0,0,1};
    Line(1)={1,2};Line(2)={3,4};Transfinite Curve{:}=3;Mesh 1;
    """)
coincident_surfaces=""
for offset in (0,4)
    for (index,(x,y)) in enumerate(((0,0),(1,0),(1,1),(0,1)))
        global coincident_surfaces*="Point($(offset+index))={$x,$y,0,1};\n"
    end
    for (index,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        global coincident_surfaces*="Line($(offset+index))={$(offset+a),$(offset+b)};\n"
    end
    tag=offset==0 ? 1 : 2
    global coincident_surfaces*="Curve Loop($tag)={$(offset+1),$(offset+2),$(offset+3),$(offset+4)};Plane Surface($tag)={$tag};\n"
end
push!(sources,coincident_surfaces*"""
    Transfinite Curve{:}=3;Transfinite Surface{:};Recombine Surface{:};Mesh 2;
    """)
gmsh.initialize(String[],false,false)
try
    gmsh.option.getString("General.Version")=="4.15.2" || error("loaded Gmsh 4.15.2 required")
    gmsh.option.setNumber("General.Terminal",0)
    for (source,expected) in zip(sources,(2,4,10,6,18))
        mktempdir() do directory
            path=joinpath(directory,"identity.geo")
            write(path,source*"\n")
            native=execute_geo(path)
            gmsh.clear();gmsh.parser.clear();gmsh.merge(path)
            tags,xyz,_=gmsh.model.mesh.getNodes()
            length(tags)==expected==nnodes(native.mesh) || error("point identity count mismatch")
            oracle=sort([Tuple(xyz[3j-2:3j]) for j in eachindex(tags)])
            actual=sort([Tuple(native.mesh.coords[:,j]) for j in 1:nnodes(native.mesh)])
            # The existing transfinite numerical-integration spacing
            # difference is about 2e-12; identity multiplicities stay exact.
            maximum(abs(oracle[i][k]-actual[i][k]) for i in eachindex(oracle),k in 1:3)<=1e-7 ||
                error("point identity coordinate multiset mismatch")
        end
    end
    # MeshOnlyEmpty retains raw curve-owned vertices, including vertices
    # coincident with separately meshed model endpoints. Stored node order
    # does not change those classifications for native or discrete curves.
    for native_curve in (true,false)
        gmsh.clear();gmsh.model.add("raw_curve_identity")
        gmsh.option.setNumber("Mesh.MeshOnlyEmpty",1)
        model=GeoModel()
        for (point,x) in ((1,0.),(2,1.))
            add_point!(model,x,0,0;tag=point)
            gmsh.model.geo.addPoint(x,0,0,1,point)
        end
        if native_curve
            add_line!(model,1,2;tag=1)
            gmsh.model.geo.addLine(1,2,1)
        else
            Tessella.Model.add_discrete_entity!(model,1,1,[(0,1),(0,2)])
        end
        gmsh.model.geo.synchronize()
        native_curve || gmsh.model.addDiscreteEntity(1,1,[1,2])
        node_tags=[11,12,13];coordinates=[0.5,0,0,1,0,0,0,0,0]
        Tessella.Model.add_discrete_nodes!(model,1,1,node_tags,coordinates)
        Tessella.Model.add_discrete_elements!(model,1,1,[1],[[101,102]],[[13,11,11,12]])
        gmsh.model.mesh.addNodes(1,1,node_tags,coordinates)
        gmsh.model.mesh.addElementsByType(1,1,[101,102],[13,11,11,12])
        gmsh.model.mesh.generate(1)
        parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[]
        append!(parts,Tessella.Model._model_curve_mesh_parts(model,"oracle"))
        append!(parts,Tessella.Model._model_point_mesh_parts(model,"oracle"))
        part_owners=Tessella.GeoExec._geo_mesh_part_node_entities(model,parts,"oracle")
        mesh,node_owners=Tessella.GeoExec._geo_merge_entity_meshes(parts;
            node_entities=part_owners)
        tags,xyz,_=gmsh.model.mesh.getNodes()
        length(tags)==nnodes(mesh)==5 || error("raw curve node identity count mismatch")
        oracle=sort([Tuple(xyz[3j-2:3j]) for j in eachindex(tags)])
        actual=sort([Tuple(mesh.coords[:,j]) for j in 1:nnodes(mesh)])
        actual==oracle || error("raw curve coordinate multiplicity mismatch")
        for entity in ((0,1),(0,2),(1,1))
            count(==(entity),values(node_owners))==
                length(gmsh.model.mesh.getNodes(entity...)[1]) ||
                error("raw curve node classification differs for entity $entity")
        end
        _,oracle_connectivity=gmsh.model.mesh.getElementsByType(1)
        length(oracle_connectivity)==length(mesh.segs)==4 ||
            error("raw curve connectivity count mismatch")
    end
    # Distinct tags within one entity can describe disconnected coincident
    # chains. Gmsh preserves them; the current native coordinate merger must
    # reject this input explicitly rather than return a welded H_0 rank one.
    for native_curve in (true,false)
        gmsh.clear();gmsh.model.add("raw_coincident_tags")
        gmsh.model.mesh.clearHomologyRequests()
        gmsh.option.setNumber("Mesh.MeshOnlyEmpty",1)
        model=GeoModel()
        for (point,x) in ((1,0.),(2,1.))
            add_point!(model,x,0,0;tag=point)
            gmsh.model.geo.addPoint(x,0,0,1,point)
        end
        if native_curve
            add_line!(model,1,2;tag=1)
            gmsh.model.geo.addLine(1,2,1)
        else
            Tessella.Model.add_discrete_entity!(model,1,1,[(0,1),(0,2)])
        end
        gmsh.model.geo.synchronize()
        native_curve || gmsh.model.addDiscreteEntity(1,1,[1,2])
        node_tags=collect(11:16)
        xyz=[0.,0,0,0.5,0,0,1,0,0,0,0,0,0.5,0,0,1,0,0]
        connectivity=[11,12,12,13,14,15,15,16]
        Tessella.Model.add_discrete_nodes!(model,1,1,node_tags,xyz)
        Tessella.Model.add_discrete_elements!(model,1,1,[1],
            [[101,102,103,104]],[connectivity])
        gmsh.model.mesh.addNodes(1,1,node_tags,xyz)
        gmsh.model.mesh.addElementsByType(1,1,[101,102,103,104],connectivity)
        gmsh.model.addPhysicalGroup(1,[1],7)
        gmsh.model.mesh.generate(1)
        length(gmsh.model.mesh.getNodes()[1])==8 || error("raw coincident node tags were lost")
        length(gmsh.model.mesh.getNodes(1,1,false,false)[1])==6 ||
            error("raw coincident curve node classification differs")
        gmsh.model.mesh.addHomologyRequest("Homology",[7],[],[0])
        gmsh.model.mesh.computeHomology()
        names=[gmsh.model.getPhysicalName(dim,tag) for
               (dim,tag) in gmsh.model.getPhysicalGroups()]
        sort!(filter(name->startswith(name,"H_0"),names))==["H_0{7}1","H_0{7}2"] ||
            error("raw coincident chains must have two connected components")
        native_error=try
            Tessella.Model._model_curve_mesh_parts(model,"oracle")
            nothing
        catch err
            err
        end
        native_error isa ArgumentError && occursin(
            "preserving coincident nodes within one entity is not implemented",
            sprint(showerror,native_error)) || error("raw coincident tags require precise native blocker")
    end
    # Existing native point mesh vertices replace the synthetic geometry
    # vertex upstream, including a raw vertex located off the point geometry.
    gmsh.clear();gmsh.model.add("native_point_attachment")
    gmsh.model.mesh.clearHomologyRequests()
    gmsh.model.geo.addPoint(0,0,0,1,1);gmsh.model.geo.synchronize()
    gmsh.model.mesh.addNodes(0,1,[11],[2.,0,0])
    gmsh.model.addPhysicalGroup(0,[1],7)
    gmsh.model.mesh.generate(1)
    model=GeoModel();add_point!(model,0.,0,0;tag=1)
    Tessella.Model.add_discrete_nodes!(model,0,1,[11],[2.,0,0])
    parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[
        Tessella.Model._model_point_mesh_parts(model,"oracle")...]
    part_owners=Tessella.GeoExec._geo_mesh_part_node_entities(model,parts,"oracle")
    mesh,node_owners=Tessella.GeoExec._geo_merge_entity_meshes(parts;node_entities=part_owners)
    tags,xyz,_=gmsh.model.mesh.getNodes()
    length(tags)==nnodes(mesh)==1 && vec(mesh.coords)==xyz==[2.,0,0] ||
        error("native point attachment must replace its geometry vertex")
    cells=Tessella.GeoExec._geo_homology_cells(model,mesh,parts,"oracle";
        node_owner=node_owners,part_node_entities=part_owners)
    Tessella.Model.add_physical_group!(model,0,[1];tag=7)
    Tessella.Model.add_homology_request!(model;kind="Homology",domain_tags=[7],dims=[0])
    Tessella.Model.compute_homology!(model,cells)
    gmsh.model.mesh.addHomologyRequest("Homology",[7],[],[0])
    gmsh.model.mesh.computeHomology()
    names=[gmsh.model.getPhysicalName(dim,tag) for (dim,tag) in gmsh.model.getPhysicalGroups()]
    filter(name->startswith(name,"H_0"),names)==["H_0{7}1"] &&
        Set(values(model.physical_names))==Set(["H_0{7}1"]) ||
        error("native point attachment homology must not create a ghost vertex")
    # More stored vertices than MPoint elements is legal. Mesh0D creates a
    # single cell on the last vertex when no explicit cells exist; explicit
    # cell connectivity remains authoritative even when other vertices exist.
    for native_point in (true,false), explicit in (:none,:first,:both)
        gmsh.clear();gmsh.model.add("point_cells_and_orphan_vertices")
        gmsh.model.mesh.clearHomologyRequests()
        model=GeoModel()
        if native_point
            gmsh.model.geo.addPoint(0,0,0,1,1);gmsh.model.geo.synchronize()
            add_point!(model,0.,0,0;tag=1)
        else
            gmsh.model.addDiscreteEntity(0,1)
            Tessella.Model.add_discrete_entity!(model,0,1)
        end
        gmsh.model.mesh.addNodes(0,1,[11,12],[0.,0,0,2.,0,0])
        Tessella.Model.add_discrete_nodes!(model,0,1,[11,12],[0.,0,0,2.,0,0])
        if explicit!==:none
            elements=explicit===:first ? [101] : [101,102]
            connectivity=explicit===:first ? [11] : [11,12]
            gmsh.model.mesh.addElementsByType(1,15,elements,connectivity)
            Tessella.Model.add_discrete_elements!(model,0,1,[15],[elements],[connectivity])
        end
        gmsh.model.addPhysicalGroup(0,[1],7);gmsh.model.mesh.generate(1)
        parts=Tuple{Int,Int,Union{Mesh,MixedMesh}}[
            Tessella.Model._model_point_mesh_parts(model,"oracle")...]
        part_owners=Tessella.GeoExec._geo_mesh_part_node_entities(model,parts,"oracle")
        mesh,node_owners=Tessella.GeoExec._geo_merge_entity_meshes(parts;node_entities=part_owners)
        cells=Tessella.GeoExec._geo_homology_cells(model,mesh,parts,"oracle";
            node_owner=node_owners,part_node_entities=part_owners)
        _,oracle_nodes=gmsh.model.mesh.getElementsByType(15)
        actual=sort([Tuple(mesh.coords[:,only(nodes)]) for (_,nodes) in cells[(0,1)]])
        expected=sort([Tuple(gmsh.model.mesh.getNode(node)[1]) for node in oracle_nodes])
        nnodes(mesh)==length(gmsh.model.mesh.getNodes()[1])==2 && actual==expected ||
            error("point cell connectivity differs from stored vertex classification")
        Tessella.Model.add_physical_group!(model,0,[1];tag=7)
        Tessella.Model.add_homology_request!(model;kind="Homology",domain_tags=[7],dims=[0])
        Tessella.Model.compute_homology!(model,cells)
        gmsh.model.mesh.addHomologyRequest("Homology",[7],[],[0]);gmsh.model.mesh.computeHomology()
        names=[gmsh.model.getPhysicalName(dim,tag) for (dim,tag) in gmsh.model.getPhysicalGroups()]
        Set(filter(name->startswith(name,"H_0"),names))==Set(values(model.physical_names)) ||
            error("point cell homology includes unreferenced stored vertices")
    end
    # Only the canonical first Point vertex participates in carrier
    # embedding. Additional attached vertices remain stored but unused.
    # Transfinite curve samples drift by about 1e-12 upstream, so this
    # fixture checks exact Point incidence without assuming curve aliasing.
    gmsh.clear();gmsh.model.add("embedded_point_extra_vertex")
    gmsh.model.mesh.clearHomologyRequests()
    for (tag,x,y) in ((1,-1.,-1.),(2,2.,-1.),(3,2.,1.),(4,-1.,1.),
                     (5,0.,0.),(6,1.,0.),(7,0.5,0.))
        gmsh.model.geo.addPoint(x,y,0.,0.5,tag)
    end
    for (tag,a,b) in ((1,1,2),(2,2,3),(3,3,4),(4,4,1),(5,5,6))
        gmsh.model.geo.addLine(a,b,tag)
    end
    gmsh.model.geo.addCurveLoop([1,2,3,4],1)
    gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(0,[7],2,1)
    gmsh.model.mesh.embed(1,[5],2,1)
    gmsh.model.mesh.setTransfiniteCurve(5,5)
    gmsh.model.mesh.addNodes(0,7,[101,102],[0.5,0,0,0.25,0,0])
    gmsh.model.mesh.generate(2)
    point_tags,point_coords,_=gmsh.model.mesh.getNodes(0,7)
    point_nodes=Dict(Tuple(point_coords[3i-2:3i])=>point_tags[i] for i in eachindex(point_tags))
    Set(keys(point_nodes))==Set([(0.5,0.,0.),(0.25,0.,0.)]) ||
        error("embedded Point additional vertices were not retained")
    curve_used=Set(vcat(gmsh.model.mesh.getElements(1,5)[3]...))
    carrier_used=Set(vcat(gmsh.model.mesh.getElements(2,1)[3]...))
    first_point=point_nodes[(0.5,0.,0.)];extra_point=point_nodes[(0.25,0.,0.)]
    first_point in carrier_used && !(extra_point in carrier_used) && !(extra_point in curve_used) ||
        error("embedded Point must transmit only its first vertex to the carrier")
    # Point/Curve constraints recovered on a common embedded surface must
    # share one node through all three native element dimensions.
    source="""
        SetFactory("OpenCASCADE");Box(1)={0,0,0,1,1,1};
        Point(101)={0.2,0.2,0.5,0.5};Point(102)={0.8,0.2,0.5,0.5};
        Point(103)={0.5,0.8,0.5,0.5};
        Line(101)={101,102};Line(102)={102,103};Line(103)={103,101};
        Curve Loop(101)={101,102,103};Plane Surface(101)={101};
        Point(104)={0.3,0.35,0.5,0.5};Point(105)={0.7,0.35,0.5,0.5};
        Point(106)={0.5,0.35,0.5,0.5};Line(104)={104,105};
        Point{106} In Surface{101};Curve{104} In Surface{101};
        Surface{101} In Volume{1};
        """
    mktempdir() do directory
        path=joinpath(directory,"recovered_constraints.geo");write(path,source)
        native=execute_geo(path;mesh_dim=3)
        hits=findall(n->Tuple(native.mesh.coords[:,n])==(0.5,0.35,0.5),
                     1:nnodes(native.mesh))
        length(hits)==1 || error("recovered Point/Curve constraints have split identity")
        node=only(hits)
        node in native.mesh.segs && node in native.mesh.tris && node in native.mesh.tets ||
            error("recovered Point incidence is missing from a mesh dimension")
        gmsh.clear();gmsh.model.mesh.clearHomologyRequests();gmsh.merge(path)
        gmsh.model.mesh.generate(3)
        tags,xyz,_=gmsh.model.mesh.getNodes()
        oracle_hits=[tags[i] for i in eachindex(tags) if
            all(k->abs(xyz[3i-3+k]-(0.5,0.35,0.5)[k])<=1e-13,1:3)]
        length(oracle_hits)==1 && only(oracle_hits) in gmsh.model.mesh.getNodes(0,106)[1] ||
            error("Gmsh common-carrier incidence must be a single Point[106] node")
    end
    println("GEO_MESH_IDENTITY_DIFFERENTIAL_OK gmsh=4.15.2 cases=",length(sources)+13,
        " raw_coincident_tag_blockers=2")
finally
    gmsh.finalize()
end
