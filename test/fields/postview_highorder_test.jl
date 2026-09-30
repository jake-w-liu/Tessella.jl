using Test
using Tessella
using Tessella.Post: PosElement, PosScheme, PosView, read_pos, write_pos,
                     postview_field
using Tessella.Elements: lagrange_nodes
using Tessella.SizeField: metric_raw, GMSH_MAX_SIZE, field_value

_parse(s)=read_pos(IOBuffer(s))[1]
function _roundtrip(view::PosView)
    io=IOBuffer()
    write_pos(io,[view])
    return read_pos(IOBuffer(String(take!(io))))[1]
end
_flat(v)=join((string(x) for x in v),",")
_nodes(mshtype::Int)=join((_flat(lagrange_nodes(mshtype)[:,j])
                          for j in axes(lagrange_nodes(mshtype),2)),",")
_values(v)=join((string(x) for x in v),",")

@testset "high-order .pos records" begin
    @testset "order-2 suffixed records parse and bind the shared scheme" begin
        for (record,mshtype) in (("SL2",8),("ST2",9),("SQ2",10),("SS2",11),
                                 ("SH2",12),("SI2",13),("SY2",14))
            n=size(lagrange_nodes(mshtype),2)
            view=_parse("""View "v" {
              $record($(_nodes(mshtype))){$(_values(1.0:n))};
            };""")
            element=view.elements[1]
            @test element.suffix=="2"
            @test element.scheme!==nothing
            @test size(element.coords,2)==n
            @test size(element.values,2)==n
        end
        # The order-2 scheme is one shared instance per family.
        view=_parse("""View "v" {
          SL2(0,0,0, 1,0,0, 0.5,0,0){1,2,3};
          SL2(0,1,0, 1,1,0, 0.5,1,0){4,5,6};
        };""")
        @test view.elements[1].scheme===view.elements[2].scheme
        # Order-2 suffixes on vector/tensor records work the same way.
        view=_parse("""View "v" {
          VL2(0,0,0, 1,0,0, 0.5,0,0){$(join(fill("1.0",9),","))};
        };""")
        @test size(view.elements[1].values)==(3,3,1)
        view=_parse("""View "v" {
          TL2(0,0,0, 1,0,0, 0.5,0,0){$(join(fill("1.0",27),","))};
        };""")
        @test size(view.elements[1].values)==(9,3,1)
        # Element widths check against the bound scheme, not the name arity.
        @test_throws ArgumentError _parse("""View "v" {
          SL2(0,0,0, 1,0,0){1,2};
        };""")
        @test_throws ArgumentError _parse("""View "v" {
          SL2(0,0,0, 1,0,0, 0.5,0,0){1,2};
        };""")
        # Unknown record names still reject; the point family has no order-2
        # scheme, so a suffixed point parses verbatim but binds nothing.
        @test_throws ArgumentError _parse("""View "v" {
          XX2(0,0,0){1};
        };""")
        view=_parse("""View "v" {
          SP2(0,0,0){1};
        };""")
        @test view.elements[1].suffix=="2"
        @test view.elements[1].scheme===nothing
    end

    @testset "INTERPOLATION_SCHEME binding follows Gmsh precedence+position" begin
        # A scheme binds the first-in-precedence family present at its file
        # position — here line (before triangle in precedence) even though
        # the triangle record came first in the file.
        view=_parse("""View "v" {
          ST(0,0,0, 1,0,0, 0,1,0){1,2,3};
          SL(0,0,0, 1,0,0){4,5};
          INTERPOLATION_SCHEME{{1,0,-1},{0,0,1}}{{0},{2},{1}};
        };""")
        @test view.scheme_positions==[2]
        @test length(view.interpolations)==1
        @test view.elements[1].scheme===nothing  # triangle stays first-order
        @test view.elements[2].scheme!==nothing  # line bound retroactively
        # Four-matrix schemes cannot bind pyramids/prisms: with a pyramid and
        # a hexahedron present, the hexahedron binds.
        view=_parse("""View "v" {
          SY(-1,-1,0, 1,-1,0, 1,1,0, -1,1,0, 0,0,1){1,1,1,1,1};
          SH(-1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
             -1,-1,1, 1,-1,1, 1,1,1, -1,1,1){1,2};
          INTERPOLATION_SCHEME{{0.5,-0.5},{0.5,0.5}}{{0},{1}}
            {{1,0,0},{0,1,0},{0,0,1},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}}
            {{0},{0},{0}};
        };""")
        @test view.elements[1].scheme===nothing  # pyramid: outside 4-matrix precedence
        @test view.elements[2].scheme!==nothing  # hexahedron bound
        # First binding wins: a second scheme record is dead for an
        # already-bound family.
        view=_parse("""View "v" {
          SL(0,0,0, 1,0,0){1,2};
          INTERPOLATION_SCHEME{{1,0,-1},{0,0,1}}{{0},{2},{1}};
          SL(0,1,0, 1,1,0){3,4};
          INTERPOLATION_SCHEME{{9,9},{9,9}}{{0},{9}};
        };""")
        scheme=view.elements[1].scheme
        @test view.elements[2].scheme===scheme
        @test size(scheme.coefval)==(2,3)  # scheme 1, not scheme 2
        # A suffixed record binds order-2 at its own position; a later
        # INTERPOLATION_SCHEME can no longer rebind that family.
        view=_parse("""View "v" {
          SL2(0,0,0, 1,0,0, 0.5,0,0){1,2,3};
          INTERPOLATION_SCHEME{{1,0,-1},{0,0,1}}{{0},{2},{1}};
          SL2(0,1,0, 1,1,0, 0.5,1,0){4,5,6};
        };""")
        @test view.elements[1].scheme===view.elements[2].scheme
        @test size(view.elements[1].scheme.coefval)==(3,3)
    end

    @testset "round-trip preserves schemes, positions, and elements" begin
        source="""View "v" {
          SL(0,0,0, 1,0,0){1,3};
          INTERPOLATION_SCHEME{{1,0,-1},{0,0,1}}{{0},{2},{1}};
          SL(0,1,0, 1,1,0){0.5,2.5};
          ST2(0,0,0, 1,0,0, 0,1,0, 0.5,0,0, 0.5,0.5,0, 0,0.5,0){1,2,3,4,5,6};
        };"""
        view=_parse(source)
        @test view.scheme_positions==[1]
        reparsed=_roundtrip(view)
        @test reparsed.scheme_positions==view.scheme_positions
        @test length(reparsed.elements)==length(view.elements)
        @test reparsed.elements[1].scheme==view.elements[1].scheme
        @test reparsed.elements[2].scheme==view.elements[2].scheme
        @test reparsed.elements[3].scheme==view.elements[3].scheme
        @test reparsed.elements[3].suffix=="2"
        @test reparsed.elements[3].coords==view.elements[3].coords
        @test reparsed.elements[3].values==view.elements[3].values
        @test reparsed.time==view.time
        # Second round-trip is byte-identical (deterministic serialization).
        io1=IOBuffer();write_pos(io1,[view])
        io2=IOBuffer();write_pos(io2,[reparsed])
        @test String(take!(io1))==String(take!(io2))
    end

    @testset "PosView constructor resolution rules" begin
        coords=[0.0 1.0;0.0 0.0;0.0 0.0]
        plain=PosElement(:line,coords,reshape([1.0,2.0],1,2))
        scheme=PosScheme([1.0 0 -1;0 0 1],[0.0;2.0;1.0;;])
        # An element-bound custom scheme requires an interpolations entry
        # binding the family at a position after one of its elements.
        custom=PosElement(:line,coords,reshape([1.0,2.0],1,2);scheme=scheme)
        @test_throws ArgumentError PosView("v";elements=[custom])
        @test PosView("v";elements=[custom],
                      interpolations=[[scheme.coefval,scheme.expval]],
                      scheme_positions=[1]) isa PosView
        # Contradictory explicit schemes reject.
        other=PosScheme([0.0 1;0.0 0],[0.0;1.0;;])
        clash=PosElement(:line,coords,reshape([1.0,2.0],1,2);scheme=other)
        @test_throws ArgumentError PosView(
            "v";elements=[custom,clash],
            interpolations=[[scheme.coefval,scheme.expval]],
            scheme_positions=[1])
        # Point elements cannot carry a scheme at all.
        @test_throws ArgumentError PosElement(
            :point,reshape([0.0,0.0,0.0],3,1),reshape([1.0],1,1);scheme=scheme)
        # Suffix validity and the implicit order-2 binding.
        @test_throws ArgumentError PosElement(:line,coords,[1.0 2.0];
                                             suffix="bad-char!")
        @test PosElement(:line,[0.0 1.0 0.5;0.0 0.0 0.0;0.0 0.0 0.0],
                         reshape([1.0,2.0,3.0],1,3);
                         suffix="2").scheme!==nothing
        # scheme_positions must be nondecreasing and inside 0:nelements.
        @test_throws ArgumentError PosView(
            "v";elements=[plain],
            interpolations=[[scheme.coefval,scheme.expval]],
            scheme_positions=[2])
        @test_throws ArgumentError PosView(
            "v";elements=[plain],
            interpolations=[[scheme.coefval,scheme.expval],
                            [scheme.coefval,scheme.expval]],
            scheme_positions=[1,0])
    end

    @testset "malformed records reject with useful errors" begin
        # Matrix count must be 2 or 4.
        @test_throws ArgumentError _parse("""View "v" {
          SL(0,0,0, 1,0,0){1,2};
          INTERPOLATION_SCHEME{{1,0,-1},{0,0,1}};
        };""")
        @test_throws ArgumentError _parse("""View "v" {
          SL(0,0,0, 1,0,0){1,2};
          INTERPOLATION_SCHEME{{1,0},{0,1}}{{0},{1}}{{1,0},{0,1}};
        };""")
        # Ragged matrix rows reject.
        @test_throws ArgumentError _parse("""View "v" {
          SL(0,0,0, 1,0,0){1,2};
          INTERPOLATION_SCHEME{{1,0,-1},{0,0}}{{0},{2},{1}};
        };""")
        # Scheme shape mismatches reject (coefficients vs exponents).
        @test_throws ArgumentError PosScheme([1.0 0;0 1],[0.0;1.0;2.0;;])
        @test_throws ArgumentError PosScheme([1.0;2.0;;],[0.0;1.0;;],
                                           [1.0;;])
        @test_throws ArgumentError PosScheme([1.0;2.0;;],[0.0;1.0;;],
                                           [1.0;;],[0.0 0.0 0.0 0.0])
        # A bound family validates every element's widths.
        @test_throws ArgumentError _parse("""View "v" {
          SL(0,0,0, 1,0,0){1,2};
          INTERPOLATION_SCHEME{{1,0,-1},{0,0,1}}{{0},{2},{1}};
          SL(0,1,0, 1,1,0){3};
        };""")
    end
end

@testset "PostViewField high-order evaluation" begin
    @testset "built-in order-2 elements interpolate exactly" begin
        # Straight P2 line: exact quadratic interpolation.
        view=_parse("""View "v" {
          SL2(0,0,0, 1,0,0, 0.5,0,0){1.0,2.0,3.0};
        };""")
        field=postview_field(view;crop_negative=false)
        @test length(field.scheme_cells)==1
        @test field_value(field,0.5,0.0,0.0)==3.0
        @test field_value(field,0.25,0.0,0.0)≈2.375
        @test field_value(field,0.75,0.0,0.0)≈2.875
        @test field_value(field,0.0,0.0,0.0)==1.0
        @test field_value(field,1.0,0.0,0.0)==2.0
        # Quadratic tet: values = x² sampled on the P2 nodes, exact everywhere.
        nodes=lagrange_nodes(11)
        vals=[n[1]^2 for n in eachcol(nodes)]
        view=_parse("""View "v" {
          SS2($(_nodes(11))){$(_values(vals))};
        };""")
        field=postview_field(view;crop_negative=false)
        for (x,y,z) in ((0.5,0.25,0.25),(0.1,0.1,0.1),(0.25,0.25,0.0))
            @test field_value(field,x,y,z)≈x^2 atol=1e-12
        end
        # Quadratic hexahedron (tensor-product basis): value x·y.
        nodes=lagrange_nodes(12)
        vals=[n[1]*n[2] for n in eachcol(nodes)]
        view=_parse("""View "v" {
          SH2($(_nodes(12))){$(_values(vals))};
        };""")
        field=postview_field(view;crop_negative=false)
        for (x,y,z) in ((0.0,0.0,0.0),(0.4,0.5,0.2),(-0.6,0.8,-0.4))
            @test field_value(field,x,y,z)≈x*y atol=1e-11
        end
        # Quadratic prism: value x²+z on the reference prism.
        nodes=lagrange_nodes(13)
        vals=[n[1]^2+n[3] for n in eachcol(nodes)]
        view=_parse("""View "v" {
          SI2($(_nodes(13))){$(_values(vals))};
        };""")
        field=postview_field(view;crop_negative=false)
        for (x,y,z) in ((0.2,0.3,0.5),(0.0,0.0,-0.5),(0.5,0.25,1.0))
            @test field_value(field,x,y,z)≈x^2+z atol=1e-11
        end
        # Quadratic Bergot pyramid: value x²+w.
        nodes=lagrange_nodes(14)
        vals=[n[1]^2+n[3] for n in eachcol(nodes)]
        view=_parse("""View "v" {
          SY2($(_nodes(14))){$(_values(vals))};
        };""")
        field=postview_field(view;crop_negative=false)
        for (x,y,z) in ((0.0,0.0,0.5),(0.25,0.0,0.25),(0.3,-0.2,0.4))
            @test field_value(field,x,y,z)≈x^2+z atol=1e-11
        end
        # Quadratic quadrangle in the plane: value x²+x·y/2.
        nodes=lagrange_nodes(10)
        vals=[n[1]^2+0.5n[1]*n[2] for n in eachcol(nodes)]
        view=_parse("""View "v" {
          SQ2($(_nodes(10))){$(_values(vals))};
        };""")
        field=postview_field(view;crop_negative=false)
        for (x,y) in ((0.0,0.0),(0.3,0.4),(-0.5,0.5))
            @test field_value(field,x,y,0.0)≈x^2+0.5x*y atol=1e-11
        end
        # Quadratic triangle: value x²+y.
        nodes=lagrange_nodes(9)
        vals=[n[1]^2+n[2] for n in eachcol(nodes)]
        view=_parse("""View "v" {
          ST2($(_nodes(9))){$(_values(vals))};
        };""")
        field=postview_field(view;crop_negative=false)
        for (x,y) in ((0.25,0.25),(0.5,0.25),(0.1,0.8))
            @test field_value(field,x,y,0.0)≈x^2+y atol=1e-11
        end
    end

    @testset "curved geometry inverts by Newton iteration" begin
        # Parabolic P2 arc: nodes (0,0),(1,0),(0.5,0.3).
        view=_parse("""View "v" {
          SL2(0,0,0, 1,0,0, 0.5,0.3,0){1.0,2.0,3.0};
        };""")
        field=postview_field(view;crop_negative=false,use_closest=false)
        @test field_value(field,0.5,0.3,0.0)==3.0        # mid node
        @test field_value(field,0.25,0.225,0.0)≈2.375   # u=-0.5 on the arc
        @test field_value(field,0.75,0.225,0.0)≈2.875   # u=+0.5
        @test field_value(field,0.5,0.5,0.0)==GMSH_MAX_SIZE  # off the arc
        # The same view with use_closest recovers the nearest node's value.
        field=postview_field(view;crop_negative=false,use_closest=true)
        @test field_value(field,0.5,0.5,0.0)==3.0
    end

    @testset "custom two- and four-matrix schemes" begin
        # Two-matrix: linear value basis {u⁰,u²,u¹} on u∈[-1,1] with
        # N1=(1-u)/2, N2=(1+u)/2.
        view=_parse("""View "v" {
          SL(0,0,0, 1,0,0){1.0,3.0};
          INTERPOLATION_SCHEME{{0.5,0,-0.5},{0.5,0,0.5}}{{0},{2},{1}};
        };""")
        field=postview_field(view;crop_negative=false)
        @test field_value(field,0.0,0.0,0.0)≈1.0
        @test field_value(field,1.0,0.0,0.0)≈3.0
        @test field_value(field,0.5,0.0,0.0)≈2.0
        # Four-matrix: order-2 value basis on a curved geometry map.
        order2=PosElement(:line,[0.0 1.0 0.5;0.0 0.0 0.0;0.0 0.0 0.0],
                          reshape([1.0,2.0,3.0],1,3);suffix="2")
        coef=order2.scheme.coefval
        exps=order2.scheme.expval
        view=PosView("v";elements=[
                PosElement(:line,[0.0 1.0 0.5;0.0 0.0 0.3;0.0 0.0 0.0],
                           reshape([1.0,2.0,3.0],1,3);
                           scheme=order2.scheme)],
                     interpolations=[[coef,exps,coef,exps]],
                     scheme_positions=[1])
        field=postview_field(view;crop_negative=false,use_closest=false)
        @test field_value(field,0.5,0.3,0.0)==3.0
        @test field_value(field,0.25,0.225,0.0)≈2.375
        @test field_value(field,0.5,0.5,0.0)==GMSH_MAX_SIZE
    end

    @testset "vector, tensor, and time-step handling" begin
        # Vector order-2 record: interpolated vector norm.
        nodes=lagrange_nodes(8)
        vals=Matrix{Float64}(undef,3,3)
        for j in 1:3
            vals[1,j]=nodes[1,j];vals[2,j]=1.0;vals[3,j]=0.0
        end
        view=_parse("""View "v" {
          VL2($(_nodes(8))){$(_values(vec(vals)))};
        };""")
        field=postview_field(view;crop_negative=false)
        @test field.num_components==3
        # lagrange_nodes coords span u∈[-1,1], so the arc maps x∈[-1,1]:
        # node 1 (u=-1 ↔ x=-1) has (-1,1,0), node 3 (u=0 ↔ x=0) has (0,1,0).
        @test field_value(field,-1.0,0.0,0.0)≈sqrt(2.0)
        @test field_value(field,0.0,0.0,0.0)≈1.0
        # Tensor order-2 record: metric_raw interpolates the 9 components.
        tvals=Matrix{Float64}(undef,9,3)
        for j in 1:3,c in 1:9
            tvals[c,j]=c+j
        end
        view=_parse("""View "v" {
          TL2($(_nodes(8))){$(_values(vec(tvals)))};
        };""")
        aniso=postview_field(view)
        @test aniso.field.num_components==9
        raw=metric_raw(aniso.field,0.5,0.0,0.0)
        @test raw!==nothing
        # Nodal basis: at node 1 (u=-1 ↔ x=-1) components equal its values.
        @test collect(metric_raw(aniso.field,-1.0,0.0,0.0))≈tvals[:,1]
        # Multiple time steps select through `time`.
        view=_parse("""View "v" {
          TIME{0.0,1.0};
          SL2(0,0,0, 1,0,0, 0.5,0,0){1,2,3, 4,5,6};
        };""")
        field=postview_field(view;time=2)
        @test field_value(field,0.5,0.0,0.0)==6.0
        field=postview_field(view;time=1)
        @test field_value(field,0.5,0.0,0.0)==3.0
        @test_throws ArgumentError postview_field(view;time=3)
    end

    @testset "mixed scheme + first-order families" begin
        view=_parse("""View "v" {
          SL2(0,0,0, 1,0,0, 0.5,0,0){1.0,2.0,3.0};
          ST(0,0,1, 1,0,1, 0,1,1){10.0,20.0,30.0};
        };""")
        field=postview_field(view;crop_negative=false,use_closest=false)
        @test field_value(field,0.25,0.0,0.0)≈2.375     # scheme line
        # Barycentric interp on the canonical triangle: .5·10+.25·20+.25·30
        @test field_value(field,0.25,0.25,1.0)≈17.5
        @test length(field.scheme_cells)==1
        @test size(field.triangles,2)==1
    end

    @testset "validation and determinism" begin
        # Degenerate scheme elements reject at field construction.
        @test_throws ArgumentError postview_field(_parse("""View "v" {
          SL2(0,0,0, 0,0,0, 0,0,0){1.0,2.0,3.0};
        };"""))
        # A planar P2 tet has a singular geometry map everywhere.
        @test_throws ArgumentError postview_field(_parse("""View "v" {
          SS2(0,0,0, 1,0,0, 0,1,0, 1,1,0,
              0.5,0,0, 0.5,0.5,0, 0,0.5,0, 0.5,0,0, 1,0.5,0, 0.5,1,0){
              $(_values(fill(1.0,10)))};
        };"""))
        # Deterministic: repeated queries agree bit-for-bit.
        view=_parse("""View "v" {
          SL2(0,0,0, 1,0,0, 0.5,0.3,0){1.0,2.0,3.0};
        };""")
        field=postview_field(view;crop_negative=false)
        @test field_value(field,0.25,0.225,0.0)==
              field_value(field,0.25,0.225,0.0)
        # The scheme-cell path is allocation-free after warm-up.
        measure(f)=(f(); @allocated f())
        @test measure(()->field_value(field,0.25,0.225,0.0))==0
        @test measure(()->field_value(field,0.5,0.5,0.0))==0
    end
end
