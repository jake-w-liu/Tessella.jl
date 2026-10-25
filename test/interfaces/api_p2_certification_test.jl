using Test,Tessella
using Tessella.Elements: lagrange_nodes,msh_dimension

const _P2C=Tessella.API
if !isdefined(@__MODULE__,:FOLDED_CAD_ROWS)
    include(joinpath(@__DIR__,"..","..","validation","native_cad_support",
                     "folded_cad_primary_fixtures_v1.jl"))
end

# These payloads were captured from Gmsh 4.15.2 before the native CAD
# admission change. They preserve actual supports, owners, parameters and IDs.
function _p2c_native_payload(state)
    @test Set(_P2C.mesh.get_nodes()[1])==Set(UInt64[node.tag for node in state.nodes])
    for node in state.nodes
        xyz,parameters,dimension,tag=_P2C.mesh.get_node(node.tag)
        @test reinterpret(UInt64,xyz)==reinterpret(UInt64,node.xyz)
        @test reinterpret(UInt64,parameters)==reinterpret(UInt64,Float64[node.parameters...])
        @test (dimension,tag)==Tuple(node.owner)
    end
    for entity in state.entities
        dimension,tag=entity.owner
        nodes,xyz,parameters=_P2C.mesh.get_nodes(dimension,tag,false,true)
        @test nodes==entity.owned_nodes
        @test reinterpret(UInt64,xyz)==reinterpret(UInt64,Float64[entity.owned_coordinates...])
        @test reinterpret(UInt64,parameters)==reinterpret(UInt64,Float64[entity.owned_parameters...])
        types,cells,connectivity=_P2C.mesh.get_elements(dimension,tag)
        @test types==entity.types
        @test cells==entity.cell_tags
        @test connectivity==entity.cell_nodes
    end
    @test (_P2C.mesh.get_max_node_tag(),_P2C.mesh.get_max_element_tag())==Tuple(state.maxima)
end
const _P2C_TYPES=(8,9,10,11,12,13,14)
const _P2C_REVERSED=((2,1),(2,1,3),(2,1,4,3),(2,1,3,4),
    (2,1,4,3,6,5,8,7),(2,1,3,5,4,6),(1,4,3,2,5))
# Actual public payload captured before certification was added.
const _P2C_CAPTURED_PYR14=[
    1e16 1e16 1e16 1e16 1e16+2 1e16 1e16 1e16 1e16 1e16 1e16 1e16 1e16 1e16;
    0. 4. 4. 0. 2. 2. 0. 1. 4. 3. 2. 3. 1. 2.;
    0. 0. 1. 1. .5 0. .5 .25 .5 .25 1. .75 .75 .5]

function _p2c_primary(msh,offset=0.)
    r=lagrange_nodes(msh)
    msh==7 && return [offset offset offset offset offset+2;
        0. 4. 4. 0. 2.;0. 0. 1. 1. .5]
    x=msh in (1,3,5) ? r[1,:].+1 : 2r[1,:]
    y=msh in (3,5) ? 2(r[2,:].+1) : 4r[2,:]
    z=msh in (5,6) ? (r[3,:].+1)./2 : r[3,:]
    return permutedims(hcat(offset.+x,y,z))
end

function _p2c_raw(msh;offset=0.,reverse=false,history=false,renumber=false)
    _P2C.initialize()
    _P2C.option("Mesh.Renumber",renumber ? 1 : 0)
    dim=msh_dimension(msh);points=_p2c_primary(msh,offset)
    tags=UInt64[500,100,900,300,700,200,1100,400][1:size(points,2)]
    order=reverse ? collect(_P2C_REVERSED[msh]) : collect(eachindex(tags))
    _P2C.model.add_discrete_entity(dim,7)
    _P2C.model.add_physical_group(dim,[7];tag=31,name="source")
    _P2C.mesh.add_nodes(dim,7,tags,vec(points))
    _P2C.mesh.add_elements_by_type(7,msh,[701],tags[order])
    if history
        _P2C.mesh.add_nodes(dim,7,[9000],[99.,99.,99.])
        _P2C.mesh.add_elements_by_type(7,msh,[9001],tags[order])
        _P2C.mesh.remove_elements(dim,7,[9001])
    end
    _P2C.mesh.generate(0)
    return (;dim,points,tags,order)
end

function _p2c_snapshot()
    payload=(_P2C.mesh.get_nodes(),_P2C.mesh.get_elements())
    slot=_P2C.MODEL_SLOTS[_P2C._current_slot_index_locked()]
    return (;model=_P2C.CURRENT[],model_value=repr(_P2C.CURRENT[]),
        cache=_P2C.LAST_MESH[],class=_P2C.LAST_MESH_CLASS[],payload,
        overlay=_P2C.LAST_MESH_HIGH_ORDER[],overlay_parent=_P2C.LAST_MESH_HIGH_ORDER_MESH[],
        overlay_mids=_P2C.LAST_MESH_HIGH_ORDER_MIDS[],
        mids_value=deepcopy(_P2C.LAST_MESH_HIGH_ORDER_MIDS[]),
        history=(_P2C.NODE_TAG_MAX[],_P2C.ELEMENT_TAG_MAX[]),
        public_max=(_P2C.mesh.get_max_node_tag(),_P2C.mesh.get_max_element_tag()),
        display=copy(_P2C.ELEMENT_VISIBILITY[]),
        visibility=deepcopy(slot.visibility),
        windows=deepcopy(slot.window_element_visibility))
end

function _p2c_unchanged(before)
    @test _P2C.CURRENT[]===before.model
    @test repr(_P2C.CURRENT[])==before.model_value
    @test _P2C.LAST_MESH[]===before.cache && _P2C.LAST_MESH_CLASS[]===before.class
    @test (_P2C.mesh.get_nodes(),_P2C.mesh.get_elements())==before.payload
    @test _P2C.LAST_MESH_HIGH_ORDER[]===before.overlay
    @test _P2C.LAST_MESH_HIGH_ORDER_MESH[]===before.overlay_parent
    @test _P2C.LAST_MESH_HIGH_ORDER_MIDS[]===before.overlay_mids
    @test _P2C.LAST_MESH_HIGH_ORDER_MIDS[]==before.mids_value
    @test (_P2C.NODE_TAG_MAX[],_P2C.ELEMENT_TAG_MAX[])==before.history
    @test (_P2C.mesh.get_max_node_tag(),_P2C.mesh.get_max_element_tag())==before.public_max
    @test _P2C.ELEMENT_VISIBILITY[]==before.display
    slot=_P2C.MODEL_SLOTS[_P2C._current_slot_index_locked()]
    @test slot.visibility==before.visibility
    @test slot.window_element_visibility==before.windows
end

# These exact maps interpolate the captured rounded P2 supports. They use
# independent polynomial formulas, not the production basis or certificate.
function _p2c_rounded_map(msh,u,v,w)
    x=msh in (1,3,5) ? u*u+u : msh==7 ? 4w*w-2w : 4u*u-2u
    y=msh in (3,5) ? 2(v+1) : msh==7 ? 2+2u : 4v
    z=msh in (5,6) ? (w+1)/2 : msh==7 ? (1+v)/2 : w
    return (big(10)^16+x,y,z)
end

function _p2c_bad_nodes(msh)
    msh==7 && return copy(_P2C_CAPTURED_PYR14)
    r=lagrange_nodes(_P2C_TYPES[msh])
    return [Float64(_p2c_rounded_map(msh,
        Rational{BigInt}(r[1,i]),Rational{BigInt}(r[2,i]),
        Rational{BigInt}(r[3,i]))[k]) for k in 1:3,i in axes(r,2)]
end

function _p2c_native_triangle(width=2.)
    _P2C.initialize()
    mktempdir() do directory
        path=joinpath(directory,"triangle.geo")
        write(path,"""
        Point(1)={0,0,0,1};Point(2)={$width,0,0,1};
        Point(3)={$width,4,0,1};Point(4)={0,4,0,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
        Transfinite Curve{1,2,3,4}=2;Transfinite Surface{1};
        """)
        _P2C.open_geo!(path;mesh_dim=0)
        return _P2C.mesh.generate(2)
    end
end

@testset "Legacy simplex elevation commits atomically" begin
    offset=(1.,0,0,1e16,0,1.,0,0,0,0,1.,0)
    source=_p2c_native_triangle()
    try
        @test source isa Tessella.MeshTypes.Mesh
        @test size(source.tris,2)==2
        # Native HighOrder constructs supports on the unchanged CAD even
        # when edited mesh coordinates make the resulting physical map fold.
        _P2C.mesh.affine_transform(offset)
        @test Set(_P2C.mesh.get_nodes()[2][1:3:end])==Set([1e16,1e16+2])
        cells=_P2C.mesh.get_elements(2,1)[2][1]
        _P2C.mesh.set_visibility(cells,5)
        _P2C.mesh.set_visibility_per_window(first(cells),0,2)
        @test _P2C.mesh.set_order(2)===nothing
        row=only(filter(row->row.name=="higher_p2_neighbors_primary_v4.json:2",FOLDED_CAD_ROWS))
        _p2c_native_payload(row.state)
        elevated=_P2C.mesh.get_elements(2,1)[2][1]
        @test _P2C.mesh.get_visibility(elevated)==fill(1,length(elevated))
    finally
        _P2C.finalize()
    end
end

@testset "Legacy manual edits transform actual P2 nodes without reinterpolation" begin
    offset=(1.,0,0,1e16,0,1.,0,0,0,0,1.,0)
    for route in (:affine_overlay,:primary_node)
        source=_p2c_native_triangle(route===:primary_node ? 4. : 2.)
        try
            route===:primary_node && _P2C.mesh.affine_transform(offset)
            # This fixture tests exact affine interpolation followed by manual
            # edits; CAD stays untranslated. SOL1 selects that explicit policy.
            route===:primary_node && _P2C.option("Mesh.SecondOrderLinear",1)
            _P2C.mesh.set_order(2)
            tags,flat,_=_P2C.mesh.get_nodes();xyz=reshape(flat,3,:)
            before=Dict(tag=>Tuple(xyz[:,i]) for (i,tag) in enumerate(tags))
            element_tags,connectivity=_P2C.mesh.get_elements_by_type(9)
            cells=reshape(connectivity,6,:)
            mids=Set(vec(cells[4:6,:]))
            expected=copy(before)
            if route===:affine_overlay
                # An explicit transform applies IEEE arithmetic to every
                # actual node, including caller-edited or singular geometry.
                # It does not request fresh midpoint interpolation.
                expected=Dict(tag=>(p[1]+1e16,p[2],p[3]) for (tag,p) in before)
                _P2C.mesh.affine_transform(offset)
            else
                target=only(tag for (tag,p) in before if p==(1e16+4,0.,0.))
                expected[target]=(1e16+2,0.,0.)
                _P2C.mesh.set_node(target,collect(expected[target]))
            end
            new_tags,new_flat,_=_P2C.mesh.get_nodes();new_xyz=reshape(new_flat,3,:)
            actual=Dict(tag=>Tuple(new_xyz[:,i]) for (i,tag) in enumerate(new_tags))
            @test actual==expected
            @test _P2C.mesh.get_elements_by_type(9)==(element_tags,connectivity)
            for tag in mids
                @test Tuple(_P2C.mesh.get_node(tag)[1])==expected[tag]
            end
            @test _P2C.mesh.get_element(first(element_tags))[1]==9
        finally
            _P2C.finalize()
        end
    end
end

function _p2c_weights(msh,u,v,w)
    msh==1 && return ((1-u)/2,(1+u)/2)
    msh==2 && return (1-u-v,u,v)
    msh==3 && return ((1-u)*(1-v)/4,(1+u)*(1-v)/4,
        (1+u)*(1+v)/4,(1-u)*(1+v)/4)
    msh==4 && return (1-u-v-w,u,v,w)
    if msh==5
        signs=((-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
            (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1))
        return Tuple((1+a*u)*(1+b*v)*(1+c*w)/8 for (a,b,c) in signs)
    end
    msh==6 && return ((1-u-v)*(1-w)/2,u*(1-w)/2,v*(1-w)/2,
        (1-u-v)*(1+w)/2,u*(1+w)/2,v*(1+w)/2)
    q=1-w
    iszero(q) && return (zero(w),zero(w),zero(w),zero(w),one(w))
    return (((q-u)*(q-v))/(4q),((q+u)*(q-v))/(4q),
        ((q+u)*(q+v))/(4q),((q-u)*(q+v))/(4q),w)
end

@testset "Rounded affine P2 supports can reverse the reference map" begin
    for msh in 1:7
        nodes=_p2c_bad_nodes(msh);reference=lagrange_nodes(_P2C_TYPES[msh])
        for i in axes(reference,2),axis in 1:3
            p=_p2c_rounded_map(msh,Rational{BigInt}(reference[1,i]),
                Rational{BigInt}(reference[2,i]),Rational{BigInt}(reference[3,i]))
            @test Rational{BigInt}(nodes[axis,i])==p[axis]
        end
        # At a strictly interior point, the signed speed/area/determinant is
        # negative, while the corresponding primary map has positive sign.
        derivative=msh in (1,3,5) ? 2*(big(-3)//4)+1 : 8*(big(1)//8)-2
        multiplier=msh in (2,4) ? 4 : msh in (3,6) ? 2 : 1
        @test multiplier*derivative<0
    end
    # Captured actual Pyr14 data gives X=1e16+4w²-2w, Y=2+2u,
    # Z=(1+v)/2. Its exact determinant at (0,0,1/8) is -1.
    @test 8*(big(1)//8)-2 == -1
    @test _p2c_bad_nodes(7)[1,5]==1e16+2
    @test all(==(1e16),_p2c_bad_nodes(7)[1,[8,10,12,13]])
end

@testset "Raw P2 elevation rejects precision folds atomically" begin
    for msh in 1:7,renumber in (false,true)
        fixture=_p2c_raw(msh;offset=1e16,history=true,renumber)
        try
            tags=_P2C.mesh.get_elements_by_type(msh,7)[1]
            primary_tags=_P2C.mesh.get_elements_by_type(msh,7)[2]
            _P2C.mesh.set_visibility(tags,5)
            _P2C.mesh.set_visibility_per_window(only(tags),0,2)
            before=_p2c_snapshot()
            @test_throws r"P2.*Jacobian|whole-domain P2" _P2C.mesh.set_order(2)
            _p2c_unchanged(before)
            # Failed construction must not consume hidden allocator history.
            _P2C.mesh.add_nodes(fixture.dim,7,Int[],[100.,100.,100.])
            @test _P2C.mesh.get_max_node_tag()==before.history[1]+1
            _P2C.mesh.add_elements_by_type(7,msh,Int[],primary_tags)
            @test _P2C.mesh.get_max_element_tag()==before.history[2]+1
        finally
            _P2C.finalize()
        end
    end
end

@testset "Well-conditioned P2 construction preserves both reference orientations" begin
    for msh in 1:7,reverse in (false,true)
        fixture=_p2c_raw(msh;reverse)
        try
            _P2C.mesh.set_order(2)
            cache=_P2C.mesh.get();block=only(cache.blocks)
            @test block.msh==_P2C_TYPES[msh]
            actual=cache.coords[:,block.nodes[:,1]]
            reference=lagrange_nodes(block.msh)
            primary=fixture.points[:,fixture.order]
            for slot in axes(reference,2),axis in 1:3
                weights=_p2c_weights(msh,Rational{BigInt}(reference[1,slot]),
                    Rational{BigInt}(reference[2,slot]),Rational{BigInt}(reference[3,slot]))
                expected=sum(weights[i]*Rational{BigInt}(primary[axis,i]) for i in eachindex(weights))
                @test Rational{BigInt}(actual[axis,slot])==expected
            end
            @test validate(cache).ok
            @test _P2C.mesh.get_nodes_for_physical_group(fixture.dim,31)[1]==
                _P2C.mesh.get_nodes(fixture.dim,7,true)[1]
        finally
            _P2C.finalize()
        end
    end
end

@testset "Native near-ULP line generation and elevation are atomic" begin
    for reverse in (false,true),route in (:generate0,:generate1,:setorder)
        _P2C.initialize()
        try
            _P2C.model.add_point(0.,0,0;tag=1)
            _P2C.model.add_point(2.,0,0;tag=2)
            _P2C.model.add_line(reverse ? 2 : 1,reverse ? 1 : 2;tag=1)
            _P2C.mesh.set_transfinite_curve(1,2)
            _P2C.model.set_coordinates(1,1e16,0,0)
            _P2C.model.set_coordinates(2,1e16+2,0,0)
            _P2C.mesh.generate(1)
            @test length(_P2C.mesh.get_elements_by_type(1,1)[1])==1
            _P2C.option("Mesh.ElementOrder",2)
            operation=route===:setorder ? ()->_P2C.mesh.set_order(2) :
                ()->_P2C.mesh.generate(route===:generate0 ? 0 : 1)
            operation()
            row=only(filter(row->row.kind=="curve" && row.reverse==reverse &&
                row.route==String(route),FOLDED_CAD_ROWS))
            _p2c_native_payload(row.state)
        finally
            _P2C.finalize()
        end
    end
end

@testset "Existing user P2 geometry remains queryable and retainable" begin
    for msh in (1,7),order in (1,2)
        _P2C.initialize()
        try
            _P2C.option("Mesh.Renumber",0)
            _P2C.option("Mesh.ElementOrder",order)
            dim=msh_dimension(msh);points=_p2c_bad_nodes(msh)
            tags=collect(101:100+size(points,2))
            _P2C.model.add_discrete_entity(dim,7)
            _P2C.mesh.add_nodes(dim,7,tags,vec(points))
            _P2C.mesh.add_elements_by_type(7,_P2C_TYPES[msh],[701],tags)
            @test _P2C.mesh.get_element(701)[1]==_P2C_TYPES[msh]
            @test _P2C.mesh.get_nodes()[2]==vec(points)
            before=(_P2C.mesh.get_nodes(),_P2C.mesh.get_elements())
            _P2C.mesh.generate(0)
            @test (_P2C.mesh.get_nodes(),_P2C.mesh.get_elements())==before
            @test _P2C.mesh.get_element(701)[2]==UInt64.(tags)
        finally
            _P2C.finalize()
        end
    end
end

_p2c_edge(a,b)=a<b ? (a,b) : (b,a)
_p2c_cell(nodes)=Tuple(sort(collect(nodes)))

# This compatibility fixture deliberately installs the legacy simplex cache.
# Its identities and owners come from the completed producer's explicit trace.
function _p2c_legacy_carrier!()
    actual=_P2C.LAST_MESH[]
    actual_class=_P2C.LAST_MESH_CLASS[]
    authority=_P2C.CURRENT[].meshing.mesh_authority
    actual isa Tessella.Elements.MixedMesh || error("fixture requires actual native mesh")
    authority!==nothing && authority.mesh===actual || error("fixture lacks actual authority")
    authority.carrier isa Tessella.MeshTypes.Mesh || error("fixture lacks simplex carrier")
    carrier=_P2C._copy_mesh(authority.carrier)
    rows=copy(authority.carrier_rows)
    length(rows)==Tessella.MeshTypes.nnodes(carrier) && length(unique(rows))==length(rows) ||
        error("fixture carrier rows are not an explicit injection")
    bindings=authority.carrier_cells
    blocks=((1,carrier.segs),(2,carrier.tris),(4,carrier.tets))
    length(bindings)==length(blocks) || error("fixture lacks cell birth bindings")
    node_tags=actual_class.public_tags.node_tags
    inverse=Dict(row=>Int32(column) for (column,row) in enumerate(rows))
    owners=Vector{Int32}[]
    for ((msh,cells),tags) in zip(blocks,bindings)
        length(tags)==size(cells,2) || error("fixture cell births are misaligned")
        block_owners=Int32[]
        for (column,tag) in enumerate(tags)
            tag!=0 || error("fixture simplex parent lacks actual identity")
            actual_msh,nodes,dimension,entity=_P2C.mesh.get_element(tag)
            actual_msh==msh || error("fixture parent is a proxy of another cell family")
            expected=node_tags[rows[cells[:,column]]]
            sort(nodes)==sort(expected) || error("fixture cell corners lost their birth identities")
            push!(block_owners,Int32(entity))
        end
        push!(owners,block_owners)
    end
    edges,faces,quads=_P2C._api_mesh_support_entities(carrier,owners...)
    # Include explicitly emitted Curve carriers even if the top geometry view
    # has no boundary Line block. The key is actual UID incidence, not XYZ.
    for (msh,dimension,offset,cells,block_owners) in _P2C._cache_catalog(actual,actual_class)
        msh==1 || continue
        for column in axes(cells,2)
            a=get(inverse,cells[1,column],Int32(0))
            b=get(inverse,cells[2,column],Int32(0))
            a!=0 && b!=0 || continue
            edges[minmax(a,b)]=(1,block_owners[column])
        end
    end
    parameters=Union{Nothing,Vector{Float64}}[
        value===nothing ? nothing : copy(value) for value in actual_class.node_parametric[rows]]
    legacy_class=_P2C._MeshClassification(carrier,actual_class.entity,
        copy(actual_class.entities),copy(actual_class.node_entities[rows]),
        deepcopy(actual_class.boundaries),owners...,
        Dict(1=>owners[1],2=>owners[2],4=>owners[3]),nothing,edges,faces,quads,parameters)
    _P2C.mesh.clear()
    lock(_P2C.STATE_LOCK) do
        _P2C._replace_mesh_cache_locked!(carrier,legacy_class;reconcile_records=false)
    end
    @test _P2C.LAST_MESH[]===carrier
    @test _P2C.LAST_MESH_CLASS[]===legacy_class
    @test _P2C.LAST_MESH_HIGH_ORDER[]===nothing
    return carrier
end

function _p2c_curved_triangles(;coincident=false,nTF=2)
    if !coincident
        _p2c_native_triangle(4.)
    else
        _P2C.initialize()
        mktempdir() do directory
            path=joinpath(directory,"independent_squares.geo")
            write(path,"""
            Point(1)={0,0,0};Point(2)={4,0,0};Point(3)={4,4,0};Point(4)={0,4,0};
            Point(101)={10,0,0};Point(102)={14,0,0};Point(103)={14,4,0};Point(104)={10,4,0};
            Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
            Line(101)={101,102};Line(102)={102,103};Line(103)={103,104};Line(104)={104,101};
            Curve Loop(1)={1,2,3,4};Curve Loop(2)={101,102,103,104};
            Plane Surface(1)={1};Plane Surface(2)={2};
            Transfinite Curve{1,2,3,4,101,102,103,104}=$nTF;
            Transfinite Surface{1,2};
            """)
            _P2C.open_geo!(path;mesh_dim=0)
            for (tag,p) in collect(_P2C.CURRENT[].points)
                tag>100 || continue
                _P2C.model.set_coordinates(tag,p[1]-10,p[2],p[3])
            end
            _P2C.mesh.generate(2)
        end
    end
    _p2c_legacy_carrier!()
    _P2C.mesh.set_order(2)
    for entity in (coincident ? (1,2) : (1,))
        tags,flat=_P2C.mesh.get_elements_by_type(9,entity)
        @test length(tags)==2(nTF-1)^2
        cells=reshape(flat,6,:)
        occurrences=Dict{Tuple{UInt64,UInt64},Vector{UInt64}}()
        for column in axes(cells,2),(slot,(a,b)) in enumerate(((1,2),(2,3),(3,1)))
            edge=_p2c_edge(cells[a,column],cells[b,column])
            push!(get!(occurrences,edge,UInt64[]),cells[slot+3,column])
        end
        location=2/(nTF-1)
        # The shared-edge midnode rides on first-order transfinite coordinates,
        # which carry Gmsh's own ~1e-12 residue (verified bitwise-identical
        # against the pinned 4.15.2 CLI); probe with a tight tolerance rather
        # than exact equality. The nearest other support is >=1 unit away, so
        # the match stays unique.
        shared=only(nodes for nodes in values(occurrences) if length(nodes)==2 &&
            isapprox(_P2C.mesh.get_node(first(nodes))[1][1:2],[location,location];
                     atol=1e-9))
        @test shared[1]==shared[2]
        p=_P2C.mesh.get_node(shared[1])[1]
        _P2C.mesh.set_node(shared[1],[p[1],p[2],entity==1 ? .05 : -.04])
    end
    return nothing
end

function _p2c_curved_payload(;certify_projection=true)
    tags,flat,_=_P2C.mesh.get_nodes();xyz=reshape(flat,3,:)
    coordinates=Dict(tag=>Tuple(xyz[:,i]) for (i,tag) in enumerate(tags))
    element_tags,flat=_P2C.mesh.get_elements_by_type(9)
    cells=reshape(flat,6,:)
    edges=Dict{Tuple{UInt64,UInt64},NTuple{3,Float64}}()
    primary=Dict{UInt64,NTuple{3,Float64}}()
    owners=Dict{NTuple{3,UInt64},Int}()
    element_keys=Dict{UInt64,NTuple{3,UInt64}}()
    for (column,tag) in enumerate(element_tags)
        key=_p2c_cell(cells[1:3,column])
        @test !haskey(owners,key)
        owners[key]=_P2C.mesh.get_element(tag)[4];element_keys[tag]=key
        points=Tuple(coordinates[cells[row,column]] for row in 1:3)
        for row in 1:3;primary[cells[row,column]]=points[row];end
        for (slot,(a,b)) in enumerate(((1,2),(2,3),(3,1)))
            edge=_p2c_edge(cells[a,column],cells[b,column])
            p=coordinates[cells[slot+3,column]]
            @test !haskey(edges,edge) || edges[edge]==p
            edges[edge]=p
            # All XY supports remain affine midpoints to parametric placement
            # precision; upstream transfinite division leaves a ~1e-12 residue
            # on the endpoint coordinates themselves, so an exact rational
            # midpoint is not representable. Consequently the XY minor of the
            # full quadratic derivative is constant to that bound, so det(J'J)
            # >= minor^2 > 0 even with the explicitly curved Z supports.
            if certify_projection
                for axis in 1:2
                    @test isapprox(p[axis],(points[a][axis]+points[b][axis])/2;atol=1e-9)
                end
            end
        end
        dx1=Rational{BigInt}(points[2][1])-Rational{BigInt}(points[1][1])
        dy1=Rational{BigInt}(points[2][2])-Rational{BigInt}(points[1][2])
        dx2=Rational{BigInt}(points[3][1])-Rational{BigInt}(points[1][1])
        dy2=Rational{BigInt}(points[3][2])-Rational{BigInt}(points[1][2])
        certify_projection && @test (dx1*dy2-dx2*dy1)^2>0
    end
    return (;edges,primary,owners,element_keys)
end

function _p2c_compare_curved(before,after;mapping=identity,removed=(),
                             coordinate_map=identity,primary_updates=Dict())
    selected_keys=Set(key for key in keys(before.owners) if !(key in removed))
    oldnodes=Set(node for key in selected_keys for node in key)
    expected_primary=Dict(mapping(node)=>get(primary_updates,node,coordinate_map(p))
        for (node,p) in before.primary if node in oldnodes)
    expected_owners=Dict(_p2c_cell(mapping.(key))=>before.owners[key] for key in selected_keys)
    wanted_edges=Set(_p2c_edge(key[a],key[b]) for key in selected_keys for (a,b) in ((1,2),(2,3),(3,1)))
    expected_edges=Dict(_p2c_edge(mapping(a),mapping(b))=>coordinate_map(p) for ((a,b),p) in before.edges
        if (a,b) in wanted_edges)
    @test after.primary==expected_primary
    @test after.owners==expected_owners
    @test after.edges==expected_edges
    @test any(p->!iszero(p[3]),values(after.edges))
end

@testset "Legacy permutations and selections preserve edited quadratic maps" begin
    for coincident in (false,true),route in (:reverse,:reverse_entity,:reverse_elements,
            :renumber_nodes,:renumber_elements,:reorder,:remove,:affine,:primary_node,:repeat_order)
        _p2c_curved_triangles(;coincident)
        try
            before=_p2c_curved_payload()
            @test length(before.owners)==(coincident ? 4 : 2)
            if coincident
                first_nodes=Set(node for (cell,owner) in before.owners if owner==1 for node in cell)
                second_nodes=Set(node for (cell,owner) in before.owners if owner==2 for node in cell)
                @test isempty(intersect(first_nodes,second_nodes))
                @test Set(before.primary[node] for node in first_nodes)==
                    Set(before.primary[node] for node in second_nodes)
                @test Set(p[3] for p in values(before.edges) if !iszero(p[3]))==Set([.05,-.04])
            end
            mapping=identity;removed=();coordinate_map=identity;primary_updates=Dict()
            tags=sort(collect(keys(before.element_keys)))
            if route===:reverse
                _P2C.mesh.reverse()
            elseif route===:reverse_entity
                _P2C.mesh.reverse([(2,1)])
            elseif route===:reverse_elements
                _P2C.mesh.reverse_elements([first(tags)])
            elseif route===:renumber_nodes
                n=Tessella.MeshTypes.nnodes(_P2C.mesh.get())
                old=collect(1:n);new=reverse(old)
                node_map=Dict(UInt64(a)=>UInt64(b) for (a,b) in zip(old,new))
                mapping=Base.Fix1(getindex,node_map)
                _P2C.mesh.renumber_nodes(old,new)
            elseif route===:renumber_elements
                _P2C.mesh.renumber_elements(tags[1:2],reverse(tags[1:2]))
            elseif route===:reorder
                _P2C.mesh.reorder_elements(9,1,[1,0])
            elseif route===:remove
                tag=first(_P2C.mesh.get_elements_by_type(9,1)[1])
                removed=(before.element_keys[tag],)
                _P2C.mesh.remove_elements(2,1,[tag])
            elseif route===:affine
                coordinate_map=p->(p[1],p[2],p[3]+.25)
                _P2C.mesh.affine_transform((1.,0,0,0,0,1.,0,0,0,0,1.,.25))
            elseif route===:repeat_order
                _P2C.mesh.set_order(2)
            else
                node=minimum(keys(before.primary));p=before.primary[node]
                updated=(p[1],p[2],p[3]+.01)
                primary_updates[node]=updated
                _P2C.mesh.set_node(node,collect(updated))
            end
            _p2c_compare_curved(before,_p2c_curved_payload();mapping,removed,coordinate_map,primary_updates)
        finally
            _P2C.finalize()
        end
    end
end

@testset "Explicit finite singular geometry remains inspectable" begin
    _p2c_curved_triangles()
    try
        before=_p2c_curved_payload()
        tags=sort(collect(keys(before.primary)))
        target,other=tags[1:2];requested=before.primary[other]
        _P2C.mesh.set_node(target,collect(requested))
        @test Tuple(_P2C.mesh.get_node(target)[1])==requested
        # Distinct primary identities remain distinct even when a manual
        # edit gives them equal coordinates. No positivity claim is made
        # for this intentionally singular caller-supplied geometry.
        _P2C.mesh.reverse()
        after=_p2c_curved_payload(;certify_projection=false)
        _p2c_compare_curved(before,after;primary_updates=Dict(target=>requested))
        @test target!=other && after.primary[target]==after.primary[other]
        @test _P2C.mesh.get_element(first(keys(after.element_keys)))[1]==9
    finally
        _P2C.finalize()
    end
end

if !isdefined(@__MODULE__,:NativeP2CertificationDefaultCurvedFixturesV1)
    include(joinpath(@__DIR__,"..","..","validation","native_geometry",
                     "native_p2_certification_default_curved_fixtures_v1.jl"))
end

function _p2c_native_default_bits(values)
    return [string(word;base=16,pad=16) for word in reinterpret(UInt64,Float64[values...])]
end

function _p2c_native_default_snapshot()
    tags,xyz,uv=_P2C.mesh.get_nodes(-1,-1,false,true)
    types,numbers,connectivity=_P2C.mesh.get_elements()
    nodes=Any[]
    for tag in sort(tags)
        point,parameters,dimension,entity=_P2C.mesh.get_node(tag)
        push!(nodes,Dict{String,Any}("tag"=>Int(tag),"xyz"=>point,
            "xyz_bits"=>_p2c_native_default_bits(point),"parameters"=>parameters,
            "parameter_bits"=>_p2c_native_default_bits(parameters),"owner"=>[dimension,entity]))
    end
    flat_cells=Int[tag for group in numbers for tag in group]
    cells=Any[]
    for tag in sort(flat_cells)
        kind,node_tags,dimension,entity=_P2C.mesh.get_element(tag)
        push!(cells,Dict{String,Any}("tag"=>tag,"msh"=>kind,"nodes"=>Int.(node_tags),
                                   "owner"=>[dimension,entity]))
    end
    entities=Any[]
    for (dimension,entity) in _P2C.model.get_entities()
        labels,points,parameters=_P2C.mesh.get_nodes(dimension,entity,false,true)
        kinds,cell_tags,cell_nodes=_P2C.mesh.get_elements(dimension,entity)
        push!(entities,Dict{String,Any}("owner"=>[dimension,entity],"node_tags"=>Int.(labels),
            "xyz"=>points,"xyz_bits"=>_p2c_native_default_bits(points),
            "parameters"=>parameters,"parameter_bits"=>_p2c_native_default_bits(parameters),
            "element_types"=>Int.(kinds),"element_tags"=>[Int.(row) for row in cell_tags],
            "element_nodes"=>[Int.(row) for row in cell_nodes],
            "visibility"=>_P2C.model.get_visibility(dimension,entity)))
    end
    return Dict{String,Any}("node_tags"=>Int.(tags),"xyz"=>xyz,"xyz_bits"=>_p2c_native_default_bits(xyz),
        "parameters"=>uv,"parameter_bits"=>_p2c_native_default_bits(uv),
        "element_types"=>Int.(types),"element_tags"=>[Int.(row) for row in numbers],
        "element_nodes"=>[Int.(row) for row in connectivity],"nodes_by_tag"=>nodes,
        "elements_by_tag"=>cells,"entities"=>entities,"element_visibility_tags"=>flat_cells,
        "element_visibility"=>Int.(_P2C.mesh.get_visibility(flat_cells)),
        "physical_groups"=>[collect(pair) for pair in _P2C.model.get_physical_groups()],
        "max_node_tag"=>Int(_P2C.mesh.get_max_node_tag()),
        "max_element_tag"=>Int(_P2C.mesh.get_max_element_tag()))
end

function _p2c_default_display_unchanged(before)
    @test _P2C.ELEMENT_VISIBILITY[]==before.display
    slot=_P2C.MODEL_SLOTS[_P2C._current_slot_index_locked()]
    @test slot.visibility==before.visibility
    # Gmsh has model/view window flags, but no public per-element window API.
    # Keep this existing Julia session contract separate from native evidence.
    @test slot.window_element_visibility==before.windows
end

function _p2c_prepare_default_optimizer_fixture()
    _p2c_curved_triangles()
    elements=_P2C.mesh.get_elements_by_type(9,1)[1]
    _P2C.mesh.set_visibility(elements,5)
    _P2C.mesh.set_visibility_per_window(first(elements),0,2)
    return _p2c_snapshot(),_p2c_curved_payload(),_p2c_native_default_snapshot()
end

@testset "Curved simplex optimizer dispatch matches upstream order gates and native default retains geometry" begin
    _p2c_curved_triangles()
    try
        elements=_P2C.mesh.get_elements_by_type(9,1)[1]
        _P2C.mesh.set_visibility(elements,5)
        _P2C.mesh.set_visibility_per_window(first(elements),0,2)
        before=_p2c_snapshot()
        @test_throws r"curved quadratic simplex refinement" _P2C.mesh.refine()
        _p2c_unchanged(before)
        # Upstream `laplaceSmoothing` refuses every order>1 face with its
        # first-order diagnostic — the check precedes its niter early-return,
        # so the error fires at niter=0 as well.  Nothing is published.
        for niter in (1,0)
            @test_throws r"only valid for first order mesh" (
                _P2C.mesh.optimize("Laplace2D",false,niter))
            _p2c_unchanged(before)
        end
        # `RelocateVertices` has no order gate: niter=0 returns without work;
        # at niter=1 the certificate-gated parametric arm relocates the bent
        # face-owned midpoint back onto its owning plane through its stored
        # (u,v) support.
        @test _P2C.mesh.optimize("Relocate2D",false,0)===nothing
        _p2c_unchanged(before)
        node_tags,_=_P2C.mesh.get_nodes()
        bent=only(tag for tag in node_tags if
            abs(_P2C.mesh.get_node(tag)[1][3])>0.01)
        @test _P2C.mesh.optimize("Relocate2D")===nothing
        @test abs(_P2C.mesh.get_node(bent)[1][3])<0.05
        @test _P2C.mesh.get_elements_by_type(9,1)[1]==elements
    finally
        _P2C.finalize()
    end

    fixtures=NativeP2CertificationDefaultCurvedFixturesV1
    @test fixtures.NATIVE_VERSION=="4.15.2" && length(fixtures.CASE_ROWS)==27
    @test length(fixtures.SEQUENCE_ROWS)==3
    for row in fixtures.CASE_ROWS
        before,curved_before,public_before=_p2c_prepare_default_optimizer_fixture()
        try
            # Private carrier iteration order has no native object-identity
            # contract. Known scalar identities, UV and incidence match the
            # exact public reconstruction bound to the 33 native captures.
            @test public_before["nodes_by_tag"]==fixtures.BEFORE_RECONSTRUCTION["nodes_by_tag"]
            @test public_before["elements_by_tag"]==fixtures.BEFORE_RECONSTRUCTION["elements_by_tag"]
            @test public_before["max_node_tag"]==fixtures.BEFORE_RECONSTRUCTION["max_node_tag"]
            @test public_before["max_element_tag"]==fixtures.BEFORE_RECONSTRUCTION["max_element_tag"]
            @test _P2C.mesh.optimize(row["method"],false,row["niter"])===nothing
            @test _p2c_native_default_snapshot()==fixtures.AFTER_STATE
            @test _P2C.logger.get_last_error()=="" && _P2C.NATIVE_LOOKUP_QUERY_ERRORS[]==0
            _p2c_default_display_unchanged(before)
            _p2c_compare_curved(curved_before,_p2c_curved_payload();
                mapping=node->UInt64((1,2,4,3)[Int(node)]))
        finally
            _P2C.finalize()
        end
    end
    for row in fixtures.SEQUENCE_ROWS
        before,curved_before,_=_p2c_prepare_default_optimizer_fixture()
        try
            for (method,niter) in row["operations"]
                @test _P2C.mesh.optimize(method,false,niter)===nothing
                @test _p2c_native_default_snapshot()==fixtures.AFTER_STATE
                @test _P2C.logger.get_last_error()=="" && _P2C.NATIVE_LOOKUP_QUERY_ERRORS[]==0
                _p2c_default_display_unchanged(before)
                _p2c_compare_curved(curved_before,_p2c_curved_payload();
                    mapping=node->UInt64((1,2,4,3)[Int(node)]))
            end
        finally
            _P2C.finalize()
        end
    end
end

@testset "Real primary compaction preserves coincident component supports" begin
    _p2c_curved_triangles(;coincident=true,nTF=3)
    try
        before=_p2c_curved_payload()
        old_count=Tessella.MeshTypes.nnodes(_P2C.mesh.get())
        @test old_count==18
        # The legacy cache retains orphan Point/Curve-owned primaries after
        # a component clear. Its actual Surface1-owned center is removed,
        # exercising a real shift of all later primary and midpoint IDs.
        removed_nodes=sort([node for node in keys(before.primary) if
            _P2C.mesh.get_node(node)[3:4]==(2,1)])
        @test length(removed_nodes)==1
        node_map=Dict(UInt64(old)=>UInt64(old-count(<=(old),removed_nodes))
            for old in 1:old_count if !(old in removed_nodes))
        removed_cells=Tuple(cell for (cell,owner) in before.owners if owner==1)
        _P2C.mesh.clear(vcat([(0,tag) for tag in 1:4],
            [(1,tag) for tag in 1:4],[(2,1)]))
        @test Tessella.MeshTypes.nnodes(_P2C.mesh.get())==old_count-1
        _p2c_compare_curved(before,_p2c_curved_payload();
            mapping=Base.Fix1(getindex,node_map),removed=removed_cells)
    finally
        _P2C.finalize()
    end
end

@testset "Conflicting merged quadratic supports reject atomically" begin
    _p2c_curved_triangles(;coincident=true)
    try
        _P2C.model.add_discrete_entity(0,30)
        _P2C.mesh.add_nodes(0,30,[9000],[0.,0.,0.])
        _P2C.mesh.add_elements_by_type(30,15,[9001],[9000])
        elements=_P2C.mesh.get_elements_by_type(9)[1]
        _P2C.mesh.set_visibility(elements,5)
        _P2C.mesh.set_visibility_per_window(first(elements),0,2)
        before=_p2c_snapshot()
        @test_throws r"conflicts with existing quadratic midpoint geometry" _P2C.mesh.remove_duplicate_nodes()
        _p2c_unchanged(before)
    finally
        _P2C.finalize()
    end
end
