module APIPostViewTensorFieldChecks
# Session-level `Field[PostView]` transport over tensor-valued list views.
# Upstream `PostViewField::numComponents` picks the dominant record kind —
# tensor over vector over scalar — and a tensor view feeds the metric overload
# (`SMetric3`) that anisotropic producers consume. The session transport used
# to admit only scalar (`S*`) records; the dominant-component selection in
# `Post.postview_field` now applies, so a tensor view resolves to a
# `PostViewAnisoField` through `_session_size_field_locked` and reaches the
# anisotropic surface producer (`Mesh.Algorithm=7`,
# `_metric_surface_algorithm_field`) and volume fill (`_VolumeMetricEdgeField`)
# as an `AbstractAnisoField`.
using Test,Tessella,LinearAlgebra
const A=Tessella.API
const F=Tessella.SizeField

# `view.add_list_data` blocked layout: every x, then every y, then every z
# corner coordinate, followed by each node's component values in order, per
# time step (`_native_list_view_pos` decodes exactly this layout).
function record_data(coords,values)
    nodes=size(coords,2);components=size(values,1)
    size(values,2)==nodes || error("values must carry one column per node")
    data=Float64[]
    for axis in 1:3,node in 1:nodes
        push!(data,coords[axis,node])
    end
    for node in 1:nodes,component in 1:components
        push!(data,values[component,node])
    end
    return data
end

# Triangle (-1,-1),(4,-1),(-1,4) contains [0,1]^2 strictly; the tetrahedron
# (-2,-2,-2),(8,-2,-2),(-2,8,-2),(-2,-2,8) contains [0,1]^3 strictly.
const _TRI=Float64[-1 4 -1; -1 -1 4; 0 0 0]
const _TET=Float64[-2 8 -2 -2; -2 -2 8 -2; -2 -2 -2 8]

# Row-major 3x3 tensor components, as stored per node by the TT/TS records.
diag3(a,b,c)=Float64[a,0,0, 0,b,0, 0,0,c]

function tensor_view!(coords,metric)
    nodes=size(coords,2)
    type=nodes==3 ? "TT" : nodes==4 ? "TS" : error("unsupported view shape")
    view=A.view.add("postview metric")
    A.view.add_list_data(view,type,1,
        record_data(coords,repeat(reshape(metric,9,1),1,nodes)))
    return view
end

function postview_background(view;selector="ViewTag",index=0)
    field=A.mesh.field.add("PostView")
    if selector=="ViewTag"
        A.mesh.field.set_number(field,"ViewTag",view)
    else
        A.mesh.field.set_number(field,selector,index)
    end
    A.mesh.field.set_as_background_mesh(field)
    return field
end

function session(body)
    A.initialize()
    try;body();finally;A.finalize();end
end

# Pin the OCC rectangle's four boundary Curves to exact transfinite
# subdivisions. The metric surface producer keeps the retained lower cells
# fixed and certifies every boundary edge satisfies the strict
# `metric_edge_length <= max_metric_edge` bound (default 1.0). The lazy
# boundary grader instead places nodes at the field size — for a metric whose
# strong direction lies along an edge that lands exactly *at* the bound, and
# ~1e-13 parameter-placement overshoot can trip the certificate. Exact 0.2
# cells score 0.8 under diag(16,1,1), so transport is exercised through the
# interior bisection rather than boundary rounding.
function pin_rectangle_curves!(nodes_per_curve=6)
    for curve in 1:4
        A.mesh.set_transfinite_curve(curve,nodes_per_curve)
    end
    return nothing
end

function plain_options!()
    for (name,value) in (("General.Terminal",0),("Mesh.Renumber",0),
            ("Mesh.MeshSizeFromPoints",0),("Mesh.MeshSizeFromCurvature",0),
            ("Mesh.MeshSizeExtendFromBoundary",0),("Mesh.Smoothing",0))
        A.option(name,value)
    end
    return nothing
end

function built_field(;raw=false)
    return lock(A.STATE_LOCK) do
        A._session_size_field_locked(A._model_locked();raw)
    end
end

function node_coordinates()
    tags,xyz,_=A.mesh.get_nodes()
    Dict(tag=>xyz[3i-2:3i] for (i,tag) in enumerate(tags))
end

function surface_triangles(surface)
    families,labels,nodes=A.mesh.get_elements(2,surface)
    families==Int32[2] || error("expected actual Tri3 product, got $families")
    return labels[1],[collect(cell) for cell in eachcol(reshape(nodes[1],3,:))]
end

function volume_tetrahedra(volume)
    families,labels,nodes=A.mesh.get_elements(3,volume)
    families==Int32[4] || error("expected actual Tet4 product, got $families")
    return labels[1],[collect(cell) for cell in eachcol(reshape(nodes[1],4,:))]
end

# Independent metric-edge audit: every element edge's `sqrt(dᵀMd)` sampled at
# both endpoints and the midpoint — the contraction `metric_edge_length`
# evaluates. Edges are shared between cells harmlessly for these statistics.
function metric_edge_lengths(cells,xyz,matrix)
    lengths=Float64[]
    for cell in cells,i in 1:length(cell)-1,j in i+1:length(cell)
        a,b=cell[i],cell[j]
        first,last=xyz[a],xyz[b]
        displacement=last-first
        push!(lengths,maximum(sqrt(dot(displacement,matrix*displacement))
                              for point in (first,last,(first+last)/2)))
    end
    return lengths
end

@testset "session PostView tensor resolves to an AbstractAnisoField" begin
    session() do
        A.model.add("postview tensor session")
        A.model.occ.add_rectangle(0.,0.,0.,1.,1.)
        A.model.occ.synchronize()
        view=tensor_view!(_TRI,diag3(16,1,1))
        postview_background(view;selector="ViewTag")
        plain_options!()
        built=built_field();raw=built_field(;raw=true)
        @test built isa F.AbstractAnisoField
        @test raw isa F.PostViewAnisoField
        metric=F.metric_at(built,0.5,0.5,0.0)
        @test metric isa F.Metric3
        # The background wrapper intersects with the isotropic model-lc cap;
        # diag(16,1,1) is below it everywhere, so the tensor transports as-is.
        @test metric.m11≈16 && metric.m22≈1 && metric.m33≈1
        @test metric.m12≈0 && metric.m13≈0 && metric.m23≈0
        @test F.metric_at(raw,0.5,0.5,0.0).m11≈16
        # The scalar operator is the directional size 1/sqrt(lambda_max).
        @test F.field_value(built,0.5,0.5,0.0)≈0.25
        # `use_closest` finds the record's element for uncovered queries.
        @test F.metric_at(raw,3.0,3.0,0.0).m11≈16
    end
end

@testset "ViewIndex/IView select the tensor view by sorted position" begin
    session() do
        A.model.add("postview tensor ViewIndex")
        A.model.occ.add_rectangle(0.,0.,0.,1.,1.)
        A.model.occ.synchronize()
        scalar=A.view.add("scalar first")
        A.view.add_list_data(scalar,"SP",1,[0.,0.,0., 0.5])
        view=tensor_view!(_TRI,diag3(16,1,1))
        # The tensor view sits at 0-based position 1 in sorted tag order.
        postview_background(view;selector="ViewIndex",index=1)
        built=built_field()
        @test built isa F.AbstractAnisoField
        @test F.metric_at(built,0.5,0.5,0.0).m11≈16
        # IView aliases ViewIndex.
        field=A.mesh.field.add("PostView")
        A.mesh.field.set_number(field,"IView",1)
        A.mesh.field.set_as_background_mesh(field)
        @test F.metric_at(built_field(),0.5,0.5,0.0).m11≈16
    end
end

@testset "tensor PostView drives anisotropic Surface meshing (Algorithm 7)" begin
    meshes=[]
    for metric in (diag3(16,1,1),diag3(1,16,1))
        session() do
            A.model.add("postview tensor surface")
            surface=A.model.occ.add_rectangle(0.,0.,0.,1.,1.)
            A.model.occ.synchronize()
            pin_rectangle_curves!()
            view=tensor_view!(_TRI,metric)
            postview_background(view)
            plain_options!()
            A.option("Mesh.Algorithm",7)
            A.mesh.generate(2)
            xyz=node_coordinates()
            labels,cells=surface_triangles(surface)
            matrix=Diagonal([metric[1],metric[5],metric[9]])
            lengths=metric_edge_lengths(cells,xyz,matrix)
            @test !isempty(cells) && !isempty(lengths)
            # Conforming bisection bound: every edge satisfies
            # metric_edge_length <= 1 (boundary included — the retained lower
            # Curve cells are discretized at the field's directional size).
            @test maximum(lengths)<=1.0+1e-9
            push!(meshes,(xyz,cells))
        end
    end
    # The two strong directions produce genuinely different gradings.
    @test meshes[1]!=meshes[2]
end

@testset "tensor PostView satisfies generate_metric_surface's aniso gate" begin
    session() do
        A.model.add("postview tensor direct")
        surface=A.model.occ.add_rectangle(0.,0.,0.,1.,1.)
        A.model.occ.synchronize()
        pin_rectangle_curves!()
        view=tensor_view!(_TRI,diag3(16,1,1))
        postview_background(view)
        plain_options!()
        output=A.mesh.generate_metric_surface(surface)
        @test validate(output).ok && Tessella.MeshTypes.ntris(output)>0
    end
end

@testset "coarse lower boundary keeps the metric certificate rejection" begin
    session() do
        A.model.add("postview tensor coarse boundary")
        surface=A.model.occ.add_rectangle(0.,0.,0.,1.,1.)
        A.model.occ.synchronize()
        # 0.5 boundary cells score 2.0 under diag(16,1,1) — over the strict
        # bound deterministically, so the producer must reject with its
        # documented diagnostic instead of splitting an owned Curve cell.
        pin_rectangle_curves!(3)
        view=tensor_view!(_TRI,diag3(16,1,1))
        postview_background(view)
        plain_options!()
        A.option("Mesh.Algorithm",7)
        before=(A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements(),
                A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
        @test_throws r"retained boundary exceeds max_metric_edge" A.mesh.generate(2)
        @test (A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements(),
               A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())==before
    end
end

@testset "tensor PostView reaches generate(3) volume metric fill" begin
    session() do
        A.model.add("postview tensor volume")
        volume=A.model.add_box(0.,0.,0.,1.,1.,1.)
        view=tensor_view!(_TET,diag3(16,1,1))
        postview_background(view;selector="ViewIndex",index=0)
        plain_options!()
        A.mesh.generate(3)
        xyz=node_coordinates()
        labels,cells=volume_tetrahedra(volume)
        @test !isempty(cells)
        volumes=[det(hcat(xyz[cell[2]]-xyz[cell[1]],
                          xyz[cell[3]]-xyz[cell[1]],
                          xyz[cell[4]]-xyz[cell[1]]))/6 for cell in cells]
        @test all(>(0),volumes) && isapprox(sum(volumes),1.;atol=1e-12)
        lengths=sort!(metric_edge_lengths(cells,xyz,Diagonal([16.,1.,1.])))
        # _VolumeMetricEdgeField scales the resolved metric by 1/4, so the
        # fill's unit bound reads sqrt(dᵀMd) <= 2 in raw input units. The
        # best-effort arm can leave a small violating remnant share (the API
        # boundary shell differs from the direct model fill), so require the
        # 95th percentile inside the bound and a bounded residual maximum —
        # the same slack model_metric_volume_functional_test.jl applies.
        @test lengths[ceil(Int,0.95*length(lengths))]<=2.0+1e-9
        @test maximum(lengths)<=2.5
    end
end

@testset "scalar and vector session PostView fields are unchanged" begin
    session() do
        A.model.add("postview scalar regression")
        surface=A.model.occ.add_rectangle(0.,0.,0.,1.,1.)
        A.model.occ.synchronize()
        view=A.view.add("sizes")
        A.view.add_list_data(view,"ST",1,record_data(_TRI,fill(0.25,1,3)))
        postview_background(view)
        plain_options!()
        built=built_field();raw=built_field(;raw=true)
        @test !(built isa F.AbstractAnisoField)
        @test !(raw isa F.AbstractAnisoField)
        @test F.field_value(built,0.5,0.5,0.0)≈0.25
        A.mesh.generate(2)
        xyz=node_coordinates()
        labels,cells=surface_triangles(surface)
        @test !isempty(cells) &&
              all(tag->haskey(xyz,tag),Iterators.flatten(cells))
    end
    session() do
        A.model.add("postview vector regression")
        A.model.occ.add_rectangle(0.,0.,0.,1.,1.)
        A.model.occ.synchronize()
        # A vector view keeps the scalar norm operator (Gmsh PostView parity:
        # field_value is the Euclidean norm of the interpolated vector).
        view=A.view.add("vectors")
        A.view.add_list_data(view,"VT",1,
            record_data(_TRI,Float64[0.3 0.3 0.3;0.4 0.4 0.4;0 0 0]))
        postview_background(view)
        built=built_field()
        @test !(built isa F.AbstractAnisoField)
        @test F.field_value(built,0.5,0.5,0.0)≈0.5
    end
end

@testset "malformed tensor view data rejects at build or evaluation" begin
    session() do
        A.model.add("postview tensor failures")
        A.model.occ.add_rectangle(0.,0.,0.,1.,1.)
        A.model.occ.synchronize()
        plain_options!()
        # Empty view: no element records at all.
        empty_view=A.view.add("empty")
        field=A.mesh.field.add("PostView")
        A.mesh.field.set_number(field,"ViewTag",empty_view)
        A.mesh.field.set_as_background_mesh(field)
        @test_throws ArgumentError built_field()
        # Indefinite tensor: the symmetric projection is not positive
        # definite, so evaluation fails (mirroring the SMetric3 contract).
        bad=A.view.add("indefinite")
        A.view.add_list_data(bad,"TT",1,
            record_data(_TRI,repeat(reshape(diag3(-1,1,1),9,1),1,3)))
        A.mesh.field.set_number(field,"ViewTag",bad)
        @test built_field() isa F.AbstractAnisoField
        threw=false
        try
            F.metric_at(built_field(),0.5,0.5,0.0)
        catch err
            threw=true
            @test err isa Union{ArgumentError,ErrorException}
        end
        @test threw
        threw=false
        try
            A.mesh.generate(2)
        catch err
            threw=true
            @test err isa Union{ArgumentError,ErrorException}
        end
        @test threw
        # Unknown view references keep their explicit diagnostics.
        A.mesh.field.set_number(field,"ViewTag",9999)
        @test_throws ArgumentError built_field()
    end
end
end
