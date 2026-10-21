# Gmsh parsed .pos: T2/T3 text records carry `FExpr` numeric arguments
# (option reads such as `General.GraphicsWidth/4`, arithmetic, calls) —
# `read_pos` evaluates them through live-override + DefaultOptions tables,
# and session import stores the records as T2/T3 string-list entries that
# re-export losslessly through `view.write`.
module PosTextExprChecks
using Test,Tessella
const A=Tessella.API
const P=Tessella.Post

# The verbatim t8 tutorial appendix view that the 2026-10-19 tutorial audit
# hit: a scalar graph view plus a text-annotation view embedding
# `General.GraphicsWidth/4`.
const VIEW4=joinpath(@__DIR__,"..","..","validation","gmsh_parity","view4.pos")

@testset "view4.pos T2/T3 expression records parse" begin
    views=P.read_pos(VIEW4)
    @test length(views)==2
    graph,annotations=views
    @test length(graph.elements)==101 && isempty(graph.texts)
    @test [e.kind for e in graph.elements]==fill(:point,101)
    @test length(annotations.elements)==1 && length(annotations.texts)==4
    t1,t2,t3,t4=annotations.texts
    @test t1.dim==2 && t1.coords==(10.0,-12.0,0.0) && t1.style==0.0 &&
        t1.strings==["File created on Fri Oct 18 23:50:20 2002"]
    # `General.GraphicsWidth` defaults to 800: `General.GraphicsWidth/4` = 200.
    @test t2.dim==2 && t2.coords==(200.0,-28.0,0.0) && t2.style==0.0 &&
        length(t2.strings)==5 && t2.strings[1]=="First time step" &&
        t2.strings[5]=="Last time step!"
    @test t3.dim==3 && t3.coords==(0.1,0.0,0.0) &&
        t3.strings==["This is a 3D string, defined in model coordinates"]
    @test t4.dim==3 && t4.coords==(0.1,0.3,0.0) && t4.strings==["Test $i" for i in 1:7]
    # The lone SP record in the text view declares five value steps.
    @test size(annotations.elements[1].values,3)==5 && annotations.time==[0.,1,2,3,4]
end

@testset "read_pos evaluates FExpr in all record lists" begin
    text="""
    View "exprs" {
    TIME{0,2*0.5};
    T2(1+2*3,-28,0){"sum"};
    T2(2^3^2,0,0){"pow-assoc"};
    T2(-2^2,0,0){"unary-binds-below-pow"};
    T2(1>0 ? 5 : 6,0,0){"ternary"};
    T2(Sqrt(9)+Min(3,1,2),0,0){"calls"};
    T2(Max(1,4)/2+Floor(0.9),0,0){"calls2"};
    T2(Pi,0,0){"constant"};
    T2(General.GraphicsHeight/3,0,0){"option"};
    T2(General.Verbosity>=5 && 1==1 || 0,0,0){"logic"};
    T2(!0,0,0){"not"};
    T2(Fmod(7,4)+Modulo(7,4),0,0){"fmod"};
    T3(ATan2(0,1),Hypot(3,4),0,0){"3d"};
    };
    """
    view=only(P.read_pos(IOBuffer(text)))
    @test view.time==[0.0,1.0]
    texts=view.texts
    @test [t.coords[1] for t in texts[1:11]]≈
        [7.0,512.0,-4.0,5.0,4.0,2.0,pi,200.0,1.0,1.0,6.0]
    @test texts[12].dim==3 && texts[12].coords==(0.0,5.0,0.0)
    # Element coordinate and value lists are expressions too.
    el=only(P.read_pos(IOBuffer(
        "View \"v\" {SL(0,0,0,1+1,2*0,3/1){Pi,2};};")))
    @test el.elements[1].coords==[0.0 2.0;0.0 0.0;0.0 3.0] &&
        el.elements[1].values==reshape([pi,2.0],1,2,1)
    # Live option overrides win over the DefaultOptions table.
    over=P.read_pos(IOBuffer(
        "View \"v\" {T2(General.GraphicsWidth/4,0,0){\"x\"};};");
        options=Dict{String,Float64}("General.GraphicsWidth"=>400.0))
    @test over[1].texts[1].coords[1]==100.0
end

@testset "read_pos expression failures stay ArgumentError" begin
    io(text;kw...)=P.read_pos(IOBuffer(text);kw...)
    @test_throws ArgumentError io("View \"v\" {T2(Foo,0,0){\"x\"};};")
    @test_throws ArgumentError io("View \"v\" {T2(General.Nope,0,0){\"x\"};};")
    @test_throws ArgumentError io("View \"v\" {T2(BogusCat.X,0,0){\"x\"};};")
    @test_throws ArgumentError io("View \"v\" {T2(Nosuch(1),0,0){\"x\"};};")
    @test_throws ArgumentError io("View \"v\" {T2(Sqrt(1,2),0,0){\"x\"};};")
    @test_throws ArgumentError io("View \"v\" {T2(1+,0,0){\"x\"};};")
    @test_throws ArgumentError io("View \"v\" {T2(1?2,0,0){\"x\"};};")
    @test_throws ArgumentError io("View \"v\" {T2(1 = 1,0,0){\"x\"};};")
    @test_throws ArgumentError io("View \"v\" {SP(0,0,/0){1};};")
    @test_throws ArgumentError io("View \"v\" {SP(0,0,0){nan};};")
    @test_throws ArgumentError io(
        "View \"v\" {T2(1,2,3){\"x\"};};";options=Dict{String,Float64}("k"=>NaN))
    @test_throws ArgumentError io("View \"v\" {};";options=[1.0])
    @test_throws ArgumentError io("View \"v\" {};";options=Dict(1=>2.0))
end

@testset "write_pos re-emits evaluated T2/T3 and round-trips" begin
    dir=mktempdir()
    source="""View "annotated" {
SL(0,0,0,1,0,0){0.2,0.8};
T2(General.GraphicsWidth/4,-28,0){"first","second"};
T3(0.1,0.3,0,8){"world","again"};
};
"""
    in_path=joinpath(dir,"in.pos");Base.write(in_path,source)
    parsed=only(P.read_pos(in_path))
    out_path=joinpath(dir,"out.pos")
    P.write_pos(out_path,[parsed])
    reread=only(P.read_pos(out_path))
    @test [(t.dim,t.coords,t.style,t.strings) for t in reread.texts]==
        [(t.dim,t.coords,t.style,t.strings) for t in parsed.texts] &&
        [(e.kind,e.coords,e.values) for e in reread.elements]==
        [(e.kind,e.coords,e.values) for e in parsed.elements]
    text_out=read(out_path,String)
    @test occursin("T2(200,-28,0){\"first\",\"second\"};",text_out) &&
        occursin("T3(0.1,0.3,0,8){\"world\",\"again\"};",text_out)
    # write∘read is a fixed point on evaluated output.
    out2=joinpath(dir,"out2.pos");P.write_pos(out2,[reread])
    @test read(out2,String)==text_out
end

@testset "session merge retains T2/T3 records and exports them" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        dir=mktempdir()
        path=joinpath(dir,"annotated.pos")
        Base.write(path,"""View "annotated" {
SL(0,0,0,1,0,0){0.2,0.8};
T2(General.GraphicsWidth/4,-28,0){"first","second"};
T2(10,-12,0){"note"};
T3(0.1,0.3,0,8){"world","again"};
};
""")
        A.model.add("carrier")
        A.merge(path)
        tag=only(A.view.get_tags())
        types,data=A.view.get_list_data(tag)
        # List data keeps the axis-major corners then node values (SL has 2).
        @test types==["SL"] && data==[[0.,1,0,0,0,0,0.2,0.8]]
        # Two T2 records: counts are records, strings are '\0'-joined
        # per-record payloads (upstream T2C), coords hold x/y/z/style.
        keys,strings,coords=A.view.get_list_data_strings(tag,2,-1,true)
        @test keys==["T2"] && strings==["first\0second","note"] &&
            coords==[200.0,-28,0,0, 10.0,-12,0,0]
        keys3,strings3,coords3=A.view.get_list_data_strings(tag,1,-1,true)
        @test keys3==["T3"] && strings3==["world\0again"] &&
            coords3==[0.1,0.3,0.0,8.0]
        # .pos export restores the records verbatim (evaluated arguments).
        out=joinpath(dir,"exported.pos");A.view.write(tag,out)
        view=only(P.read_pos(out))
        @test length(view.elements)==1 && length(view.texts)==3
        @test view.texts[1].coords==(200.0,-28.0,0.0) &&
            view.texts[1].strings==["first","second"] &&
            view.texts[2].strings==["note"] &&
            view.texts[3].dim==3 && view.texts[3].style==8.0 &&
            view.texts[3].strings==["world","again"]
        # The tutorial audit file itself merges as two session views.
        A.merge(VIEW4)
        @test length(A.view.get_tags())==3
    finally;A.finalize();end
end

@testset "session merge rejects bad expressions atomically" begin
    A.initialize()
    try
        A.option("General.Terminal",0);A.model.add("kept")
        prior=A.view.get_tags()
        dir=mktempdir();bad=joinpath(dir,"bad.pos")
        for body in ("View \"v\" {T2(Unknown.Width,0,0){\"a\"};};",
                     "View \"v\" {T2(1,2){\"a\"};};",
                     "View \"v\" {INTERPOLATION_SCHEME{{1}}{{0}};};")
            Base.write(bad,"View \"ok\" {SP(0,0,0){1};};\n"*body)
            @test_throws ArgumentError A.merge(bad)
            @test A.view.get_tags()==prior
        end
    finally;A.finalize();end
end
end
