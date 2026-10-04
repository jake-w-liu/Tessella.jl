using Test, Tessella, TOML
using Tessella.Elements: MixedMesh, ElementBlock, msh_dimension
using Tessella.MeshTypes: Mesh, nnodes
using Tessella.Model: mesh_model_surface, mesh_model_volume, model_to_mixed
if !isdefined(@__MODULE__,:QuadTriNoNewRectGridCertificates)
    include("quadtri_nonew_rect_grid_certificates.jl")
end
const _RGSC=QuadTriNoNewRectGridCertificates

function _rgsc_saved_mesh(stage,dimension,tag)
    entity=only(e for e in stage["entities"] if (e["dim"],e["tag"])==(dimension,tag))
    ids=sort!(unique(Int(node) for cell in entity["cells"] for node in cell["nodes"]))
    position=Dict(node=>slot for (slot,node) in enumerate(ids))
    records=Dict(Int(node["tag"])=>node for node in stage["nodes"])
    coords=hcat((Float64.(records[node]["coordinates"]) for node in ids)...)
    families=sort!(unique(Int(cell["type"]) for cell in entity["cells"]))
    blocks=[ElementBlock(msh,hcat(([Int32(position[Int(node)]) for node in cell["nodes"]]
        for cell in entity["cells"] if cell["type"]==msh)...)) for msh in families]
    return MixedMesh(coords,blocks),ids
end

function _rgsc_saved_fixture(record)
    audit=record["audit"];raw=record["payload"]
    f=_RGSC.fixture(raw["name"];grid=Tuple(Int.(audit["grid"])),
        layers=audit["intervals"]==1 ? :one : :three,
        height=Float64(audit["direction"]),laterals=Bool(audit["recombine_laterals"]))
    if haskey(raw,"variant")
        variant=raw["variant"];axes=Tuple(Int.(variant["axes"]).+1)
        corners=Tuple(ntuple(k->k==axes[1] ? Float64(xy[1]) : k==axes[2] ? Float64(xy[2]) : 0.,3)
            for xy in variant["corners"])
        return merge(f,(;source=raw["input_geo"],axes,corners,
            height=Float64(variant["height"]),levels=Float64.(raw["requested_levels"])))
    end
    return merge(f,(;source=raw["input_geo"]))
end

function _rgsc_saved_supports(raw)
    primary=Dict(parse(Int,node)=>Int(mapped) for (node,mapped) in raw["primary_remap"])
    original=Dict(Int(node["tag"])=>node for node in raw["p1"]["nodes"])
    quadratic=Dict(Int(node["tag"])=>node for node in raw["p2"]["nodes"])
    @test Set(keys(primary))==Set(keys(original))
    for (old,new) in primary
        @test original[old]["coordinate_bits"]==quadratic[new]["coordinate_bits"]
        @test original[old]["owner"]==quadratic[new]["owner"]
    end
    owners=Dict{Tuple,Tuple{Int,Int}}()
    function add(nodes,owner)
        support=_RGSC.key(primary[Int(node)] for node in nodes)
        previous=get(owners,support,nothing)
        if previous===nothing || owner[1]<previous[1];owners[support]=owner
        elseif owner[1]==previous[1];@test owner==previous
        end
    end
    for entity in raw["p1"]["entities"]
        dim=Int(entity["dim"]);dim==0 && continue
        owner=(dim,Int(entity["tag"]))
        for cell in entity["cells"]
            nodes=Int.(cell["nodes"]);msh=Int(cell["type"])
            if dim==1;add(nodes,owner);continue;end
            faces=dim==2 ? (Tuple(eachindex(nodes)),) : _RGSC.FACES[msh]
            for pattern in faces
                face=Tuple(nodes[k] for k in pattern)
                for k in eachindex(face);add((face[k],face[mod1(k+1,length(face))]),owner);end
                add(face,owner)
            end
            msh==5 && add(nodes,owner)
        end
    end
    supports=Dict(parse(Int,node)=>Tuple(Int.(support)) for (node,support) in raw["supports"])
    @test Set(keys(supports))==setdiff(Set(keys(quadratic)),Set(values(primary)))
    @test length(Set(values(supports)))==length(supports)
    for (node,support) in supports
        @test Tuple(Int.(quadratic[node]["owner"]))==owners[_RGSC.key(support)]
        expected=[sum(quadratic[source]["coordinates"][k] for source in support)/length(support) for k in 1:3]
        @test maximum(abs.(Float64.(quadratic[node]["coordinates"]).-expected))<=2e-11
    end
    # Full typed P2 face supports follow actual P1 face identity. A common
    # coordinate is not used as evidence that two interpolation nodes agree.
    basis=Dict{Int,Tuple}(node=>(node,) for node in values(primary));merge!(basis,supports)
    faces=Dict{Tuple,Vector{Tuple}}()
    for entity in raw["p2"]["entities"]
        entity["dim"]==3 || continue
        for cell in entity["cells"]
            family=Int(cell["type"])-7;nodes=Int.(cell["nodes"])
            for pattern in _RGSC.FACES[family]
                support=_RGSC.key(nodes[k] for k in pattern)
                actual=_RGSC.key(node for node in nodes if Set(basis[node])⊆Set(support))
                @test length(actual)==(length(support)==3 ? 6 : 9)
                push!(get!(faces,support,Tuple[]),actual)
            end
        end
    end
    @test all(length(rows) in (1,2) for rows in values(faces))
    @test all(rows[1]==rows[2] for rows in values(faces) if length(rows)==2)
    boundary=Dict(support=>only(rows) for (support,rows) in faces if length(rows)==1)
    lower=Dict(_RGSC.key(cell["nodes"][1:(Int(cell["type"])==9 ? 3 : 4)])=>_RGSC.key(cell["nodes"])
        for entity in raw["p2"]["entities"] if entity["dim"]==2 for cell in entity["cells"])
    @test boundary==lower
end

function _rgsc_saved_oracles()
    artifact=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","quadtri_nonew_rect_grid_oracle.toml"))
    @test artifact["gmsh_version"]=="4.15.2" && artifact["fixture_count"]==16
    masks=Set(Tuple(Int.(row["states"])) for row in artifact["factory_masks"])
    @test length(artifact["factory_masks"])==315 && length(masks)==253
    variants=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","quadtri_nonew_rect_grid_variants_oracle.toml"))
    @test variants["fixture_count"]==8 && variants["gmsh_version"]=="4.15.2"
    for record in vcat(artifact["fixtures"],variants["fixtures"])
        raw=record["payload"];f=_rgsc_saved_fixture(record)
        @test raw["status"]=="certified" && isempty(raw["warnings_or_errors"])
        source,_=_rgsc_saved_mesh(raw["p1"],2,1)
        volume,_=_rgsc_saved_mesh(raw["p1"],3,1)
        cert=_RGSC.certify(volume,f,source;oracle=true)
        @test all(state in masks for state in values(cert.states))
        @test length(cert.domains)==prod(f.grid)*f.intervals
        @test length(raw["p2"]["nodes"])==(2f.grid[1]+1)*(2f.grid[2]+1)*(2f.intervals+1)
        _rgsc_saved_supports(raw)
        # Reconstruct the actual quadratic reference polynomial from the
        # captured Float64 nodes. The independent Python primary proof stores
        # its exact minimum Bernstein coefficient, not a quadrature sample.
        quadratic=Dict(Int(node["tag"])=>Tuple(Float64.(node["coordinates"])) for node in raw["p2"]["nodes"])
        recorded=Dict((Tuple(Int.(row["entity"])),Int(row["cell"]))=>row for row in raw["p2_whole_maps"])
        for entity in raw["p2"]["entities"]
            entity["dim"]==3 || continue
            for cell in entity["cells"]
                map=_RGSC.quadratic_map_certificate(Int(cell["type"]),Tuple(quadratic[Int(node)] for node in cell["nodes"]))
                reference=recorded[((3,Int(entity["tag"])),Int(cell["tag"]))]
                @test map.minimum==parse(_RGSC.Q,reference["minimum_exact"])
                @test map.coefficients==reference["coefficients"]
                @test collect(map.degrees)==reference["degrees"]
            end
        end
    end
end

function _rgsc_quadratic_reference()
    for msh in (11,12,13,14)
        vertices=[Tuple(_RGSC.Q(v)/2 for v in p) for p in _RGSC._P2_REFERENCE_TWICE[msh]]
        map=_RGSC.quadratic_map_certificate(msh,vertices)
        @test map.minimum==1 && map.maximum==1
        @test map.volume==(msh==11 ? 1//6 : msh==12 ? 8 : msh==13 ? 1 : 4//3)
    end
    # Positive primary corners alone do not certify an actual P2 map. These
    # independent nonprimary-node changes leave every primary corner intact.
    for (msh,node,axis,value) in ((11,5,1,2),(12,27,1,4),(13,7,1,2),(14,14,3,2))
        vertices=[Tuple(_RGSC.Q(v)/2 for v in p) for p in _RGSC._P2_REFERENCE_TWICE[msh]]
        vertices[node]=ntuple(k->k==axis ? _RGSC.Q(value) : vertices[node][k],3)
        @test_throws r"whole quadratic reference determinant" _RGSC.quadratic_map_certificate(msh,vertices)
    end
end

function _rgsc_product(f;entry=false)
    geometry=_RGSC.execute(f.source;dim=entry ? 3 : 0)
    out=Int.(geometry.lists["sweep"]);m=geometry.model
    source=mesh_model_surface(m,f.surface)
    volume=entry ? geo_entity_mesh(geometry,3,out[2]) : mesh_model_volume(m,out[2])
    cert=_RGSC.certify(volume,f,source)
    @test validate(volume).ok
    @test nnodes(volume)==prod(f.grid.+1)*(f.intervals+1)
    projected=model_to_mixed(m,volume,3,out[2])
    @test validate(projected).ok
    @test _RGSC.Quad.typed_signature(projected)==_RGSC.Quad.typed_signature(volume)
    a,b=f.grid;nb=2(a+b);ni=(a-1)*(b-1);n=f.intervals
    expected=[8,2(nb-4)+4(n-1),2ni+(nb-4)*(n-1),ni*(n-1)]
    @test [count(owner->owner[1]==dim,projected.entity_data.node_entities) for dim in 0:3]==expected
    parts=[(2,tag,mesh_model_surface(m,tag)) for tag in (f.surface,out[1],out[3:end]...)]
    @test _RGSC.Quad.surface_faces(parts,volume)==Set(keys(cert.boundary))
    @test _RGSC.Quad.surface_faces(_RGSC.execute(f.source;dim=2).mesh_parts,volume)==Set(keys(cert.boundary))
    return (;geometry,source,volume,cert,projected)
end

function _rgsc_match_primary_variant(record)
    f=_rgsc_saved_fixture(record);product=_rgsc_product(f;entry=true);raw=record["payload"]
    reference,ids=_rgsc_saved_mesh(raw["p1"],2,1)
    @test nnodes(reference)==nnodes(product.source)
    mapped=Dict{Int,Int}()
    for native in 1:nnodes(product.source)
        matches=findall(reference_node->maximum(abs.(product.source.coords[:,native].-reference.coords[:,reference_node]))<=2e-11,
            1:nnodes(reference))
        @test length(matches)==1
        length(matches)==1 || continue
        mapped[native]=only(matches)
    end
    @test length(Set(values(mapped)))==nnodes(reference)
    actual=Set(_RGSC.cycle(Tuple(mapped[Int(node)] for node in cell)) for block in _RGSC.blocks(product.source) for cell in eachcol(block.nodes))
    expected=Set(_RGSC.cycle(Tuple(Int.(cell))) for block in _RGSC.blocks(reference) for cell in eachcol(block.nodes))
    @test actual==expected
    # These are raw captured native Curve samples, not analytic fractions.
    for chain in raw["source_chains"]
        native=product.geometry.model.curve_params[Int(chain["curve"])]
        captured=Float64.(chain["stored_owned_parameters"])
        @test length(native)==length(captured)+2
        @test first(native)==0. && last(native)==1.
        @test maximum(abs.(sort(native[2:end-1]).-sort(captured)))<=2e-11
    end
    return product
end

function _rgsc_original_frames(f)
    M=Tessella.Model;g=_RGSC.execute(f.source;dim=0);m=g.model;out=Int.(g.lists["sweep"])
    t=out[2];caller="independent rectangular source frames"
    completed=M._extrude_nonew_plan(deepcopy(m),t,caller)
    original=completed.source_mesh;nodes=only(_RGSC.blocks(original)).nodes
    params=M._extrude_gate(m,3,t);spec=m.meshing.extrude_specs[(3,t)]
    baseline=_RGSC.certify(completed.sweep.volume,f,original)
    ncells=size(nodes,2);orders=[collect(1:ncells),collect(ncells:-1:1),circshift(collect(1:ncells),1)]
    for order in orders,rotation in 0:3
        rotated=hcat((circshift(nodes[:,order[column]],mod(rotation+column-1,4)) for column in 1:ncells)...)
        source=MixedMesh(original.coords,[ElementBlock(3,rotated)])
        independent=_RGSC.source_complex(source,f)
        data=M._extrude_nonew_rect_grid_source(m,f.surface,source,spec,caller)
        @test data.source_cells==[Tuple(Int32.(cell)) for cell in eachcol(rotated)]
        @test data.source_coordinates==[Tuple(original.coords[:,node]) for node in 1:nnodes(original)]
        @test Set(findall(data.boundary_vertices))==independent.border
        catalog=M._extrude_nonew_rect_grid_plan(data,completed.catalog.levels,
            completed.catalog.layer_refs,f.laterals,caller)
        M._extrude_nonew_rect_grid_product_certify(completed.sweep.cols,catalog,spec,caller)
        plan=M._extrude_nonew_rect_grid_finish(m,t,params,f.surface,source,catalog,
            completed.sweep.cols,out[1],Tuple(out[3:end]),caller)
        cert=_RGSC.certify(plan.sweep.volume,f,source)
        @test length(cert.domains)==ncells*f.intervals
        @test Set(keys(cert.boundary))==Set(keys(baseline.boundary))
    end
end

function _rgsc_malformed_source(f)
    M=Tessella.Model;g=_RGSC.execute(f.source;dim=0);m=g.model
    source=mesh_model_surface(m,f.surface);spec=m.meshing.extrude_specs[(3,Int(g.lists["sweep"][2]))]
    nodes=only(_RGSC.blocks(source)).nodes
    stable=repr(m);coordinates=copy(source.coords);connectivity=copy(nodes)
    bad=deepcopy(source);bad.coords[:,2]=bad.coords[:,1]
    @test_throws r"physically coincident" M._extrude_nonew_rect_grid_source(m,f.surface,bad,spec,"malformed source")
    @test repr(m)==stable && source.coords==coordinates && nodes==connectivity
    bad=deepcopy(source);only(_RGSC.blocks(bad)).nodes[:,2]=nodes[:,1]
    @test_throws ArgumentError M._extrude_nonew_rect_grid_source(m,f.surface,bad,spec,"malformed source")
    @test repr(m)==stable && source.coords==coordinates && nodes==connectivity
    bad=deepcopy(source);only(_RGSC.blocks(bad)).nodes[:,1]=reverse(nodes[:,1])
    @test_throws ArgumentError M._extrude_nonew_rect_grid_source(m,f.surface,bad,spec,"malformed source")
    @test repr(m)==stable && source.coords==coordinates && nodes==connectivity
    # A simple count/Euler test would accept the untouched incidence here.
    # The actual strict convexity/embedding certificate must reject this
    # displaced interior vertex without modifying native curve parameters.
    cert=_RGSC.source_complex(source,f);bad=deepcopy(source)
    interior=minimum(cert.interior);bad.coords[f.axes[1],interior]=maximum(source.coords[f.axes[1],:])+1.
    @test_throws ArgumentError M._extrude_nonew_rect_grid_source(m,f.surface,bad,spec,"malformed source")
    @test repr(m)==stable && source.coords==coordinates && nodes==connectivity
end

function _rgsc_atomic(f,diagnostic)
    g=_RGSC.execute(f.source;dim=0);m=g.model;out=Int.(g.lists["sweep"])
    before=repr(m);params=deepcopy(m.curve_params)
    @test_throws diagnostic mesh_model_volume(m,out[2])
    @test repr(m)==before && m.curve_params==params
    context=Tessella.IO._GeoNumericContext();sentinel=Mesh(zeros(3,1))
    push!(context.mesh_parts,(0,9001,sentinel));parts=context.mesh_parts;options=copy(context.option_numbers)
    @test_throws diagnostic Tessella.GeoExec._geo_mesh_model(m,3,context)
    @test repr(m)==before && m.curve_params==params
    @test context.mesh_parts===parts && only(parts)[3]===sentinel
    @test context.option_numbers==options
end

@testset "Dynamic rectangular Quad-grid independent certificates" begin
    @testset "Exact quadratic reference spaces and interior folds" begin
        _rgsc_quadratic_reference()
    end
    @testset "Sixteen unit and eight geometric primary products with actual P2 identities" begin
        _rgsc_saved_oracles()
    end
    @testset "Native B3/B2/B0 maps, complete shells and all public surface routes" begin
        for grid in ((2,3),(3,3)),plane in (:XY,:YZ,:ZX),height in (-1.,1.),laterals in (false,true),layers in (:one,:three,:graded)
            _rgsc_product(_RGSC.fixture("rectangular";grid,plane,height,laterals,layers))
        end
        for grid in ((2,4),(4,3),(4,4)),laterals in (false,true)
            _rgsc_product(_RGSC.fixture("dynamic_sizes";grid,laterals,layers=:graded))
        end
        for shape in (:skew,:rounded),laterals in (false,true)
            _rgsc_product(_RGSC.fixture("sampled_frames";grid=(2,3),shape,laterals,
                plane=:XZ,height=-.75,layers=:nonbinary,winding=-1,pins=:opposite,
                curve_reverse=5,tags=:sparse,offset=(.2,-.3,.4),progression=4.);entry=true)
        end
        variants=TOML.parsefile(joinpath(@__DIR__,"..","artifacts","quadtri_nonew_rect_grid_variants_oracle.toml"))
        for record in variants["fixtures"];_rgsc_match_primary_variant(record);end
    end
    @testset "Original source IDs survive cell orders and local cyclic frames" begin
        for grid in ((2,3),(3,3)),laterals in (false,true)
            _rgsc_original_frames(_RGSC.fixture("original_identity";grid,laterals,layers=:graded,height=-1.))
        end
    end
    @testset "Malformed actual complexes and unsupported transforms fail atomically" begin
        _rgsc_malformed_source(_RGSC.fixture("malformed";grid=(3,3),layers=:one))
        for laterals in (false,true)
            f=_RGSC.fixture("collapsed_layers";grid=(2,3),laterals,layers=:four,offset=(0.,0.,2.0^50),height=.25)
            _rgsc_atomic(f,r"strictly ordered|collapse")
            f=_RGSC.fixture("oblique";grid=(2,3),laterals,layers=:one)
            _rgsc_atomic(merge(f,(;source=replace(f.source,"Extrude{0.0,0.0,1.0}"=>"Extrude{0.25,0.0,1.0}"))),r"axis-normal|one finite nonzero")
        end
    end
end
