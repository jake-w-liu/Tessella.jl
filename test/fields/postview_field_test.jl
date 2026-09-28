using Test
using Tessella
using Tessella.Post: PosElement, PosView, read_pos, postview_field
using Tessella.SizeField: metric_raw, PostViewAnisoField, geo_postview_context,
                          GMSH_MAX_SIZE
using Tessella.IO: GeoParams, GeoFieldSpec

# Element records are duck-typed: anything with .kind/.coords/.values works.
_line(kind,coords,values)=(kind=kind,coords=coords,values=values)

const _LINE=[0.0 2.0;0.0 0.0;0.0 0.0]

@testset "PostViewField(records) construction and evaluation" begin
    records=[_line(:line,_LINE,[0.2 0.8]),
             _line(:point,reshape([3.0,0,0],3,1),[1.4;;])]
    field=PostViewField(records;use_closest=false)
    @test field isa PostViewField
    @test field.num_components==1
    @test field_value(field,0.5,0.0,0.0)≈0.35
    @test field_value(field,3.0,0.0,0.0)==1.4
    @test field_value(field,5.0,5.0,5.0)==GMSH_MAX_SIZE

    # PosElement records from the .pos reader behave identically.
    elements=[PosElement(:line,_LINE,[0.2 0.8])]
    parsed=PostViewField(elements;use_closest=false)
    @test field_value(parsed,0.5,0.0,0.0)≈0.35

    # Mixed element families: Gmsh searches 3-D records first.
    tri=_line(:triangle,[0.0 1.0 0.0;0.0 0.0 1.0;0.0 0.0 0.0],[0.2 0.5 0.8])
    line=_line(:line,[0.0 0.0;0.0 0.0;0.0 1.0],[9.0 9.5])
    mixed=PostViewField([tri,line];use_closest=false)
    @test field_value(mixed,0.25,0.25,0.0)≈0.425
    # (0,0,0.5) lies on the line; the co-leaf flat triangle's own thickened
    # bounding box does not claim it (Gmsh checks each element's box).
    @test field_value(mixed,0.0,0.0,0.5)≈9.25
    @test field_value(PostViewField([tri];use_closest=false),
                      0.0,0.0,0.5)==GMSH_MAX_SIZE

    vector_records=[_line(:line,_LINE,[3.0 0.0;0.0 4.0;0.0 0.0])]
    vector_field=PostViewField(vector_records;use_closest=false)
    @test vector_field.num_components==3
    @test field_value(vector_field,1.0,0.0,0.0)==2.5

    # Validation: consistent component/step counts, kinds, shapes, limits.
    @test_throws ArgumentError PostViewField(
        [_line(:line,_LINE,[0.2 0.8]),
         _line(:line,_LINE,[1.0 1.0;1.0 1.0;1.0 1.0])])
    @test_throws ArgumentError PostViewField(
        [_line(:line,_LINE,reshape([0.2,0.8,0.4,1.6],1,2,2)),
         _line(:line,_LINE,[0.2 0.8])])
    @test_throws ArgumentError PostViewField([_line(:bogus,_LINE,[0.2 0.8])])
    @test_throws ArgumentError PostViewField(
        [_line(:line,_LINE,zeros(4,2,1))])
    @test_throws ArgumentError PostViewField(
        [_line(:line,[0.0 1.0;0.0 0.0],[0.2 0.8])])
    @test_throws ArgumentError PostViewField(
        [_line(:line,_LINE,[NaN 0.8])])
    @test_throws ArgumentError PostViewField(PosElement[];time=2)
    @test_throws ArgumentError PostViewField(
        [_line(:line,_LINE,[0.2 0.8]),
         _line(:line,_LINE,[0.4 0.6]),
         _line(:point,reshape([3.0,0,0],3,1),[1.4;;])];max_elements=2)
    @test_throws ArgumentError PostViewField(
        [_line(:line,_LINE,[0.2 0.8]),
         _line(:point,reshape([3.0,0,0],3,1),[1.4;;]),
         _line(:point,reshape([4.0,0,0],3,1),[1.6;;]),
         _line(:point,reshape([5.0,0,0],3,1),[1.8;;]),
         _line(:point,reshape([6.0,0,0],3,1),[2.0;;])];max_nodes=4)

    # Time-step selection.
    steps=reshape([0.2,0.8, 0.4,1.6, 0.6,3.2],1,2,3)
    stepped=PostViewField([_line(:line,_LINE,steps)];time=3,use_closest=false)
    @test field_value(stepped,0.5,0.0,0.0)≈1.25
    @test_throws ArgumentError PostViewField(
        [_line(:line,_LINE,steps)];time=4)
end

@testset "tensor PostViewField: metric_raw, metric_at, PostViewAnisoField" begin
    tensor_records=[_line(:line,_LINE,
        reshape([1.0,0,0,0,1,0,0,0,1, 4.0,0,0,0,4,0,0,0,4],9,2))]
    tensor_field=PostViewField(tensor_records;use_closest=false)
    @test tensor_field.num_components==9
    @test field_value(tensor_field,0.5,0.0,0.0)==GMSH_MAX_SIZE

    # metric_raw interpolates the 9 row-major components.
    raw=metric_raw(tensor_field,1.0,0.0,0.0)
    @test collect(raw)≈[2.5,0,0, 0,2.5,0, 0,0,2.5]
    @test metric_raw(tensor_field,5.0,5.0,5.0)===nothing
    @test_throws ArgumentError metric_raw(
        PostViewField([_line(:line,_LINE,[0.2 0.8])]),0.0,0.0,0.0)

    # metric_at symmetrizes and enforces SPD.
    m=metric_at(tensor_field,1.0,0.0,0.0)
    @test m isa Metric3
    @test m.m11≈2.5 && m.m22≈2.5 && m.m33≈2.5
    @test metric_at(tensor_field,5.0,5.0,5.0)==
          Metric3(1.0,1.0,1.0,0.0,0.0,0.0)
    skew_records=[_line(:point,reshape([0.0,0,0],3,1),
        reshape([0.0,1.0,0, -1.0,0.0,0, 0,0,1.0],9,1,1))]
    skew=PostViewField(skew_records)
    @test_throws ArgumentError metric_at(skew,0.0,0.0,0.0)
    indefinite=PostViewField([_line(:point,reshape([0.0,0,0],3,1),
        reshape([-1.0,0,0, 0,1.0,0, 0,0,1.0],9,1,1))])
    @test_throws ArgumentError metric_at(indefinite,0.0,0.0,0.0)

    # use_closest falls back to the closest node's containing element.
    closest=PostViewField(tensor_records)
    @test metric_at(closest,5.0,0.0,0.0).m11≈4.0
    @test metric_raw(closest,5.0,0.0,0.0)[1]≈4.0

    # PostViewAnisoField composes as a metric field.
    aniso=PostViewAnisoField(tensor_field)
    @test field_value(aniso,1.0,0.0,0.0)≈inv(sqrt(2.5))
    @test field_value(aniso,1.0,0.0,0.0,(2,7))≈inv(sqrt(2.5))
    @test metric_at(aniso,1.0,0.0,0.0).m11≈2.5
    @test_throws ArgumentError PostViewAnisoField(
        PostViewField([_line(:line,_LINE,[0.2 0.8])]))

    # CropNegativeValues does not apply to tensor metrics (Gmsh parity).
    nocrop=PostViewField([_line(:point,reshape([0.0,0,0],3,1),
        reshape([-4.0,0,0, 0,1.0,0, 0,0,1.0],9,1,1))];crop_negative=true)
    @test_throws ArgumentError metric_at(nocrop,0.0,0.0,0.0)

    # Parsed tensor view drives the whole pipeline.
    view=read_pos(IOBuffer("""
        View "m" {
        TT(0,0,0, 1,0,0, 0,1,0){
        4,0,0,0,4,0,0,0,4, 4,0,0,0,4,0,0,0,4, 4,0,0,0,4,0,0,0,4};
        };
        """))[1]
    parsed_aniso=postview_field(view)
    @test parsed_aniso isa PostViewAnisoField
    @test field_value(parsed_aniso,0.25,0.25,0.0)≈0.5
end

@testset "geo_postview_context resolver" begin
    views=[PosView("sizes";elements=PosElement[
               PosElement(:line,_LINE,[0.2 0.8])]),
           PosView("metrics";elements=PosElement[
               PosElement(:line,_LINE,
                   reshape([1.0,0,0,0,1,0,0,0,1,
                            4.0,0,0,0,4,0,0,0,4],9,2))])]
    resolver=geo_postview_context(views)
    scalar_params=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
        Dict(1=>GeoFieldSpec(1,"PostView",
            Dict("ViewIndex"=>"0","CropNegativeValues"=>"0",
                 "UseClosest"=>"0"))),1)
    scalar_field=build_geo_size_field(scalar_params,Dict();
                                      context_fields=resolver)
    @test scalar_field isa Tessella.SizeField.AbstractSizeField
    @test field_value(scalar_field,0.5,0.0,0.0)≈0.35
    # A view miss returns MAX_LC, which the background wrapper clamps to the
    # model characteristic length (1.0 here, no model_bbox).
    @test field_value(scalar_field,3.0,0.0,0.0)==1.0

    # IView aliases ViewIndex; the tensor view resolves as aniso.
    iview_params=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
        Dict(1=>GeoFieldSpec(1,"PostView",Dict("IView"=>"1"))),1)
    iview_field=build_geo_size_field(iview_params,Dict();
                                     context_fields=resolver)
    @test iview_field isa Tessella.SizeField.AbstractAnisoField
    @test metric_at(iview_field,1.0,0.0,0.0).m11≈2.5

    # Tensor view resolves through the metric overload path.
    tensor_params=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
        Dict(1=>GeoFieldSpec(1,"PostView",Dict("ViewIndex"=>"1"))),1)
    tensor_field=build_geo_size_field(tensor_params,Dict();
                                      context_fields=resolver)
    @test tensor_field isa Tessella.SizeField.AbstractAnisoField
    @test metric_at(tensor_field,1.0,0.0,0.0).m11≈2.5
    @test field_value(tensor_field,1.0,0.0,0.0)≈inv(sqrt(2.5))

    # ViewTag selects by tag and wins over ViewIndex.
    tagged=[PosView("a";tag=5,elements=views[1].elements),
            PosView("b";tag=9,elements=views[2].elements)]
    tag_resolver=geo_postview_context(tagged)
    tag_params=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
        Dict(1=>GeoFieldSpec(1,"PostView",
            Dict("ViewIndex"=>"0","ViewTag"=>"9"))),1)
    tag_field=build_geo_size_field(tag_params,Dict();
                                   context_fields=tag_resolver)
    @test tag_field isa Tessella.SizeField.AbstractAnisoField
    @test metric_at(tag_field,1.0,0.0,0.0).m11≈2.5

    # Errors are strict.
    @test_throws ArgumentError geo_postview_context(PosView[])
    @test_throws ArgumentError geo_postview_context(views;time=0)
    @test_throws ArgumentError build_geo_size_field(
        GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"PostView",Dict{String,String}())),1),
        Dict();context_fields=geo_postview_context(
            [(tag=1,elements=Any[])]))
    @test_throws ArgumentError build_geo_size_field(
        GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"PostView",Dict("ViewIndex"=>"9"))),1),
        Dict();context_fields=resolver)
    @test_throws ArgumentError build_geo_size_field(
        GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
            Dict(1=>GeoFieldSpec(1,"PostView",
                Dict("ViewTag"=>"42"))),1),
        Dict();context_fields=resolver)

    # A resolver rejects non-PostView specs.
    nonpost=GeoParams(NaN,NaN,1.0,0,Dict{Tuple{Int,Int},String}(),
        Dict(1=>GeoFieldSpec(1,"Distance",Dict("PointsList"=>"{1}"))),1)
    @test_throws ArgumentError build_geo_size_field(nonpost,Dict();
                                                    context_fields=resolver)
end
