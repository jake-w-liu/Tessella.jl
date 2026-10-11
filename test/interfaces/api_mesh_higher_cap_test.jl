using Test,Tessella
if !isdefined(@__MODULE__,:QuadTriNoNewB4FreeStripCertificates)
    include("../geometry/quadtri_nonew_b4_free_strip_certificates.jl")
end

function _hc_geometric_cells(mesh)
    values=mesh isa Tessella.MeshTypes.Mesh ?
        ((1,mesh.segs),(2,mesh.tris),(4,mesh.tets)) :
        ((Int(block.msh),block.nodes) for block in mesh.blocks)
    return [(msh,Tuple(Tuple(mesh.coords[:,node]) for node in cell))
        for (msh,cells) in values for cell in eachcol(cells)]
end

const _HC=Tessella.API
const _HCC=QuadTriNoNewB4FreeStripCertificates

function _hc_fixture(action;strip_length=5,layers=:one)
    fixture=_HCC.fixture("actual-higher-cap";strip_length,layers)
    _HC.initialize()
    try
        mktempdir() do directory
            path=joinpath(directory,"actual_higher_cap.geo");write(path,fixture.source)
            opened=_HC.open_geo!(path;mesh_dim=0)
            output=Int.(opened.lists["sweep"])
            _HC.option("Mesh.Renumber",0);_HC.option("Mesh.MeshOnlyEmpty",1)
            action(fixture,(;source=fixture.surface,cap=output[1],volume=output[2],laterals=output[3:end]))
        end
    finally
        _HC.finalize()
    end
end

function _hc_payload()
    return (;nodes=_HC.mesh.get_nodes(),cells=_HC.mesh.get_elements(),
        model=_HC.CURRENT[],model_repr=repr(_HC.CURRENT[]),cache=_HC.LAST_MESH[],
        class=_HC.LAST_MESH_CLASS[],history=(_HC.NODE_TAG_MAX[],_HC.ELEMENT_TAG_MAX[]))
end

function _hc_shell(volume,surfaces)
    types,_,connections=_HC.mesh.get_elements(3,volume)
    faces=Dict{Tuple,Int}()
    for (msh,values) in zip(types,connections)
        width=Tessella.Elements.msh_num_nodes(msh)
        for cell in eachcol(reshape(values,width,:)),pattern in _HCC.FACES[msh]
            face=Tuple(sort!(collect(cell[collect(pattern)])))
            faces[face]=get(faces,face,0)+1
        end
    end
    @test all(count->count in (1,2),Base.values(faces))
    boundary=Dict(face=>count for (face,count) in faces if count==1)
    shell=Dict{Tuple,Int}()
    for surface in surfaces
        types,_,connections=_HC.mesh.get_elements(2,surface)
        for (msh,values) in zip(types,connections)
            width=Tessella.Elements.msh_num_nodes(msh)
            for cell in eachcol(reshape(values,width,:))
                face=Tuple(sort!(collect(cell)));shell[face]=get(shell,face,0)+1
            end
        end
    end
    @test shell==boundary
    return nothing
end

function _hc_quadratic_shell(volume,surfaces)
    traces=Dict{Tuple,Int}()
    for msh in _HC.mesh.get_element_types(3,volume),width in (3,4)
        family=Int(msh)-7
        any(length(face)==width for face in _HCC.FACES[family]) || continue
        nodes=_HC.mesh.get_element_face_nodes(msh,width,volume,false)
        for face in eachcol(reshape(nodes,width==3 ? 6 : 9,:))
            key=Tuple(sort!(collect(face)))
            traces[key]=get(traces,key,0)+1
        end
    end
    @test all(count->count in (1,2),values(traces))
    shell=Dict{Tuple,Int}()
    for surface in surfaces
        types,_,connections=_HC.mesh.get_elements(2,surface)
        for (msh,nodes) in zip(types,connections)
            @test msh in (9,10)
            for face in eachcol(reshape(nodes,msh==9 ? 6 : 9,:))
                key=Tuple(sort!(collect(face)))
                shell[key]=get(shell,key,0)+1
            end
        end
    end
    @test shell==Dict(face=>count for (face,count) in traces if count==1)
end

@testset "Actual Source and cap retain primary identity through quadratic generation" begin
    for layers in (:one,:three)
        _hc_fixture(;layers) do fixture,owners
            _HC.mesh.generate(2)
            source=Tessella.Model._model_mesh_product(_HC.CURRENT[].meshing.mesh_authority,2,owners.source)
            original_source=_HC.mesh.get_elements_by_type(3,owners.source)[2]
            original_cap=_HC.mesh.get_elements_by_type(2,owners.cap)[2]
            primary=Dict(tag=>_HC.mesh.get_node(tag)[1] for tag in
                unique(vcat(original_source,original_cap)))
            _HC.option("Mesh.ElementOrder",2)
            _HC.mesh.generate(3)
            quadratic_source=reshape(_HC.mesh.get_elements_by_type(10,owners.source)[2],9,:)
            quadratic_cap=reshape(_HC.mesh.get_elements_by_type(9,owners.cap)[2],6,:)
            @test vec(quadratic_source[1:4,:])==original_source
            @test vec(quadratic_cap[1:3,:])==original_cap
            @test all(_HC.mesh.get_node(tag)[1]==point for (tag,point) in primary)
            volume=_HCC.public_volume(_HC,owners.volume)
            certificate=_HCC.certify_quadratic(volume,fixture,source;native_top=false)
            @test certificate.total==1 && certificate.nsupport==Tessella.MeshTypes.nnodes(volume)
            _hc_quadratic_shell(owners.volume,[owners.source,owners.cap,owners.laterals...])
            # Parameters come from the queried native carrier, including closure
            # vertices owned by a lower-dimensional entity.
            for (dim,tag) in _HC.model.get_entities()
                dim in (1,2) || continue
                nodes,xyz,parameters=_HC.mesh.get_nodes(dim,tag,true,true)
                @test length(parameters)==dim*length(nodes)
                @test isempty(nodes) || maximum(abs.(_HC.model.get_value(dim,tag,parameters).-xyz))<=2e-11
            end
            # The next pre-SetOrder1 pass must remove the previous P2 support
            # rows before the actual linear Source/cap admission runs again.
            _HC.mesh.generate(3)
            @test vec(reshape(_HC.mesh.get_elements_by_type(10,owners.source)[2],9,:)[1:4,:])==original_source
            @test vec(reshape(_HC.mesh.get_elements_by_type(9,owners.cap)[2],6,:)[1:3,:])==original_cap
            repeated=_HCC.certify_quadratic(_HCC.public_volume(_HC,owners.volume),fixture,source;native_top=false)
            @test repeated.total==certificate.total && repeated.ncell==certificate.ncell
            _hc_quadratic_shell(owners.volume,[owners.source,owners.cap,owners.laterals...])
        end
    end
end

@testset "Hidden native extrusion precedes OnlyEmpty and visible order elevation" begin
    # Independent Gmsh 4.15.2 rows in onlyempty_hidden_extrusion_primary_v1:
    # both visible and hidden supported native extrusions allocate fresh volume
    # cells. Hidden volumes stay linear during the visible-only P2 postpass.
    # Pointer-dependent initial subdivisions differ between separate primary
    # processes; the per-row before/after family counts remain identical.
    for order in (1,2),hidden in (false,true)
        _hc_fixture() do fixture,owners
            _HC.mesh.generate(2)
            source=Tessella.Model._model_mesh_product(_HC.CURRENT[].meshing.mesh_authority,2,owners.source)
            _HC.option("Mesh.ElementOrder",order);_HC.mesh.generate(3)
            types,tags,connections=_HC.mesh.get_elements(3,owners.volume)
            before=Dict(Int(kind)-(order==2 ? 7 : 0)=>length(values)
                for (kind,values) in zip(types,tags))
            primary=Dict{UInt64,Vector{Float64}}()
            for (kind,nodes) in zip(types,connections)
                family=Int(kind)-(order==2 ? 7 : 0)
                corners=Tessella.Elements.msh_num_nodes(family)
                for cell in eachcol(reshape(nodes,Tessella.Elements.msh_num_nodes(kind),:)),
                    node in cell[1:corners]
                    primary[node]=_HC.mesh.get_node(node)[1]
                end
            end
            maximum_cell=_HC.mesh.get_max_element_tag()
            _HC.model.set_visibility([(3,owners.volume)],hidden ? 0 : 1)
            _HC.option("Mesh.MeshOnlyVisible",1);_HC.mesh.generate(3)
            types,tags,_=_HC.mesh.get_elements(3,owners.volume)
            elevated=order==2 && !hidden
            @test Set(types)==Set(Int32.(elevated ? (11,14) : (4,7)))
            @test Dict(Int(kind)-(elevated ? 7 : 0)=>length(values)
                for (kind,values) in zip(types,tags))==before
            @test minimum(vcat(tags...))>maximum_cell
            @test all(_HC.mesh.get_node(node)[1]==point for (node,point) in primary)
            @test _HC.mesh.get_element_types(2,owners.source)==Int32[order==2 ? 10 : 3]
            volume=_HCC.public_volume(_HC,owners.volume)
            certificate=elevated ? _HCC.certify_quadratic(volume,fixture,source;native_top=false) :
                _HCC.certify(volume,fixture,source;native_top=false)
            @test certificate.total==1
            @test _HC.model.get_visibility(3,owners.volume)==(hidden ? 0 : 1)
        end
    end
end

@testset "Actual higher generation exposes retained lower products" begin
    # Independent Gmsh 4.15.2 M5/N1 proof: 28 Line + 34 Triangle + 5 Quad
    # + 8 Point + 1 foreign Tet = 76 cells. Virgin mesh CRC pins stay separate.
    _hc_fixture() do fixture,owners
        _HC.model.add_discrete_entity(3,90)
        _HC.mesh.add_nodes(3,90,6001:6004,[10.,0,0,11,0,0,10,1,0,10,0,1])
        _HC.mesh.add_elements_by_type(90,4,[9001],6001:6004)
        _HC.mesh.generate(2)
        types,tags,_=_HC.mesh.get_elements()
        @test Dict(kind=>length(elements) for (kind,elements) in zip(types,tags))==
              Dict(1=>28,2=>34,3=>5,4=>1,15=>8)
        @test sum(length,tags)==76
        @test _HC.mesh.get_element(9001)==(4,UInt64.(6001:6004),3,90)
        @test length(_HC.mesh.get_elements(2,owners.source)[2][1])==5
        @test length(_HC.mesh.get_elements(2,owners.cap)[2][1])==10
    end
end

@testset "Actual Source and cap vectors control native volume generation" begin
    _hc_fixture() do fixture,owners
        _HC.mesh.generate(2)
        node_tags=_HC.mesh.get_nodes()[1]
        _HC.mesh.renumber_nodes(node_tags,UInt64[10000+17i for i in eachindex(node_tags)])
        _HC.mesh.reorder_elements(3,owners.source,reverse(0:4))
        # A Source-only vector edit is inconsistent with the copied-cap pairs.
        before=_hc_payload()
        @test_throws ArgumentError _HC.mesh.generate(3)
        @test _HC.CURRENT[]===before.model && repr(_HC.CURRENT[])==before.model_repr
        @test _HC.LAST_MESH[]===before.cache && _HC.LAST_MESH_CLASS[]===before.class
        @test _HC.mesh.get_nodes()==before.nodes && _HC.mesh.get_elements()==before.cells
        @test (_HC.NODE_TAG_MAX[],_HC.ELEMENT_TAG_MAX[])==before.history
        cap_order=[2cell+offset for cell in reverse(0:4) for offset in (0,1)]
        _HC.mesh.reorder_elements(2,owners.cap,cap_order)
        source_connections=_HC.mesh.get_elements(2,owners.source)[3]
        cap_connections=_HC.mesh.get_elements(2,owners.cap)[3]
        boundary_nodes=Dict(entity=>_HC.mesh.get_nodes(entity...)
            for entity in _HC.model.get_entities() if entity[1]<=2)
        cell_highwater=_HC.mesh.get_max_element_tag()
        _HC.mesh.generate(3)
        @test _HC.mesh.get_elements(2,owners.source)[3]==source_connections
        @test _HC.mesh.get_elements(2,owners.cap)[3]==cap_connections
        @test all(entity->_HC.mesh.get_nodes(entity...)==boundary_nodes[entity],keys(boundary_nodes))
        @test all(tag->tag>cell_highwater,Iterators.flatten(_HC.mesh.get_elements(3,owners.volume)[2]))
        authority=_HC.CURRENT[].meshing.mesh_authority
        source=Tessella.Model._model_mesh_product(authority,2,owners.source)
        volume=_HCC.public_volume(_HC,owners.volume)
        certificate=_HCC.certify(volume,fixture,source;native_top=false)
        @test certificate.ncell>0 && certificate.total==1
        _hc_shell(owners.volume,[owners.source,owners.cap,owners.laterals...])
    end
end

@testset "Legacy volume cache retains original completed Source and cap products" begin
    _hc_fixture() do fixture,owners
        _HC.option("Mesh.MeshOnlyEmpty",0)
        # Admit the real Model frontend's legacy Volume-only carrier and its
        # completed boundary products. Public generation now publishes all
        # products, so it no longer creates this historical cache layout.
        model=_HC.CURRENT[]
        scope=Tessella.Model._extrude_nonew_scope(model,"legacy actual input";
            _working_model=true)
        legacy=Tessella.Model.mesh_model_volume(model,owners.volume;
            _extrude_scope=scope,_extrude_working_model=true)
        class=_HC._classify_cached_mesh(model,_HC._classification_skeleton(legacy,false),
            3,owners.volume,legacy;_extrude_scope=scope)
        previous=Tessella.Model._model_mesh_completed(model,3;native_inputs=scope.surfaces)
        model.meshing.mesh_authority=previous
        _HC._replace_mesh_cache_locked!(legacy,class;preserve_source_inputs=true)
        @test previous.mesh===nothing
        source=deepcopy(previous.native_inputs[(2,Int32(owners.source))])
        cap=deepcopy(previous.native_inputs[(2,Int32(owners.cap))])
        source_cells=_hc_geometric_cells(source);cap_cells=_hc_geometric_cells(cap)
        _HC.option("Mesh.MeshOnlyEmpty",1)
        _HC.mesh.generate(3)
        authority=_HC.CURRENT[].meshing.mesh_authority
        actual_source=Tessella.Model._model_mesh_product(authority,2,owners.source)
        actual_cap=Tessella.Model._model_mesh_product(authority,2,owners.cap)
        @test _hc_geometric_cells(actual_source)==source_cells
        @test _hc_geometric_cells(actual_cap)==cap_cells
        volume=_HCC.public_volume(_HC,owners.volume)
        certificate=_HCC.certify(volume,fixture,actual_source;native_top=false)
        @test certificate.total==1 && certificate.ncell>0
        _hc_shell(owners.volume,[owners.source,owners.cap,owners.laterals...])
    end
end

@testset "Actual copied-cap products remain authoritative across NoNew families" begin
    families=(
        ("two_tri",:QuadTriNoNewTwoTriCertificates),
        ("quad_patch",:QuadTriNoNewQuadPatchCertificates),
        ("quad_strip",:QuadTriNoNewQuadStripCertificates),
        ("three_quad_strip",:QuadTriNoNewThreeQuadStripCertificates),
        ("four_quad_strip",:QuadTriNoNewFourQuadStripCertificates),
        ("rect_grid",:QuadTriNoNewRectGridCertificates),
        ("b4_rec_strip",:QuadTriNoNewB4RecStripCertificates),
        ("b4_free_strip",:QuadTriNoNewB4FreeStripCertificates))
    for (name,symbol) in families
        @testset "$name actual Source and cap" begin
            if !isdefined(@__MODULE__,symbol)
                include("../geometry/quadtri_nonew_$(name)_certificates.jl")
            end
            certificate_module=getfield(@__MODULE__,symbol)
            fixture=certificate_module.fixture("actual-cap-$name";
                layers=:one,laterals=name=="b4_rec_strip")
            try
                _HC.initialize()
                mktempdir() do directory
                    path=joinpath(directory,"actual_cap.geo");write(path,fixture.source)
                    opened=_HC.open_geo!(path;mesh_dim=0)
                    sweep=Int.(opened.lists["sweep"])
                    source=hasproperty(fixture,:surface) ? fixture.surface : 1
                    _HC.option("Mesh.Renumber",0);_HC.option("Mesh.MeshOnlyEmpty",1)
                    _HC.mesh.generate(2)
                    source_connections=_HC.mesh.get_elements(2,source)[3]
                    cap_connections=_HC.mesh.get_elements(2,sweep[1])[3]
                    boundary_nodes=Dict(entity=>_HC.mesh.get_nodes(entity...)
                        for entity in _HC.model.get_entities() if entity[1]<=2)
                    _HC.mesh.generate(3)
                    @test _HC.mesh.get_elements(2,source)[3]==source_connections
                    @test _HC.mesh.get_elements(2,sweep[1])[3]==cap_connections
                    @test all(entity->_HC.mesh.get_nodes(entity...)==boundary_nodes[entity],keys(boundary_nodes))
                    authority=_HC.CURRENT[].meshing.mesh_authority
                    source_mesh=Tessella.Model._model_mesh_product(authority,2,source)
                    volume=_HCC.public_volume(_HC,sweep[2])
                    # The cap is the retained public Surface product. A native
                    # regeneration formula based on current source row numbers
                    # cannot predict its original pointer identity after merging.
                    # The original independent TwoTri certificate accepts a
                    # Mesh carrier. Expose its exact actual Tri3 matrix and
                    # coordinates without changing its proof or source rows.
                    certificate_source=name=="two_tri" ? Tessella.MeshTypes.Mesh(
                        source_mesh.coords;tris=only(source_mesh.blocks).nodes,
                        tri_tag=only(source_mesh.blocks).tags) : source_mesh
                    certificate=name in ("quad_strip","three_quad_strip","four_quad_strip",
                        "rect_grid","b4_rec_strip","b4_free_strip") ?
                        certificate_module.certify(volume,fixture,certificate_source;native_top=false) :
                        certificate_module.certify(volume,fixture,certificate_source)
                    @test certificate.ncell>0
                    _hc_shell(sweep[2],[source,sweep[1],sweep[3:end]...])
                end
            finally
                _HC.finalize()
            end
        end
    end
end
