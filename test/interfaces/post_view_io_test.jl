using Test
using Tessella
using Tessella.Post: View, PosElement, PosText, PosView, read_pos, write_pos,
                     postview_field
using Tessella.SizeField: PostViewField, PostViewAnisoField, GMSH_MAX_SIZE

const _POS_SCALAR_LINE="""
View "line" {
SL(0,0,0, 1,0,0){0.2,0.8};
};
"""

@testset "PosElement/PosText/PosView validation" begin
    element=PosElement(:line,[0.0 1.0;0.0 0.0;0.0 0.0],[0.2 0.8])
    @test element.kind==:line && size(element.coords)==(3,2)
    @test size(element.values)==(1,2,1)
    tensor=PosElement(:point,reshape([0.0,0,0],3,1),reshape(collect(1.0:9.0),9,1,1))
    @test size(tensor.values)==(9,1,1)
    @test_throws ArgumentError PosElement(:bogus,[0.0;0;0;;],[1.0;;])
    @test_throws ArgumentError PosElement(:line,[0.0 1.0;0.0 0.0],[0.2 0.8])
    @test_throws ArgumentError PosElement(:line,[0.0 1.0 2.0;zeros(2,3)],[0.2 0.8 0.4])
    @test_throws ArgumentError PosElement(:line,[0.0 1.0;0.0 0.0;0.0 0.0],[0.2 0.8;0.3 0.4])
    @test_throws ArgumentError PosElement(:line,[0.0 NaN;0.0 0.0;0.0 0.0],[0.2 0.8])
    @test_throws ArgumentError PosElement(:line,[0.0 1.0;0.0 0.0;0.0 0.0],[0.2 Inf])
    @test_throws ArgumentError PosElement(:line,[0.0 1.0;0.0 0.0;0.0 0.0],
                                          zeros(4,2,1))
    @test_throws ArgumentError PosElement(:line,[0.0 1.0;0.0 0.0;0.0 0.0],
                                          zeros(1,2,0))
    @test_throws ArgumentError PosElement(:line,fill(true,3,2),[0.2 0.8])
    @test_throws ArgumentError PosElement(:line,
        [0.0 1.0;0.0 0.0;0.0 0.0],fill(true,1,2))

    text=PosText(3,(1.0,2.0,3.0),10.0,["a","b"])
    @test text.dim==3 && text.coords==(1.0,2.0,3.0) && text.style==10.0
    @test text.strings==["a","b"]
    @test PosText(2,[1,2,0],10,["a"]).coords==(1.0,2.0,0.0)
    @test_throws ArgumentError PosText(4,(1,2,3),10,["a"])
    @test_throws ArgumentError PosText(2,(1,2,1),10,["a"])
    @test_throws ArgumentError PosText(3,(1,2,3),10,String[])
    @test_throws ArgumentError PosText(3,(1,2,3),10,["a\"b"])
    @test_throws ArgumentError PosText(3,(1,2,NaN),10,["a"])
    @test_throws ArgumentError PosText(3,(1,2,3),Inf,["a"])

    view=PosView("v";elements=PosElement[element])
    @test view.name=="v" && view.tag==0 && length(view.elements)==1
    @test view.time==[0.0]
    @test_throws ArgumentError PosView("a\"b")
    @test_throws ArgumentError PosView("v";tag=-1)
    @test_throws ArgumentError PosView("v";tag=typemax(Int64))
    @test_throws ArgumentError PosView("v";elements=Any["bogus"])
    @test_throws ArgumentError PosView("v";texts=Any["bogus"])
    @test_throws ArgumentError PosView("v";elements=PosElement[element],
                                       time=[0.0,1.0])
    @test_throws ArgumentError PosView("v";time=[NaN])
    @test_throws ArgumentError PosView("v";interpolations=[[ones(2,2)]])
    @test_throws ArgumentError PosView("v";interpolations=[[ones(0,2),ones(0,2)]])

    steps=Array{Float64,3}(undef,1,2,2)
    steps[1,:,1].=[0.2,0.8];steps[1,:,2].=[0.4,1.6]
    stepped=PosElement(:line,[0.0 1.0;0.0 0.0;0.0 0.0],steps)
    view2=PosView("v2";tag=7,elements=PosElement[stepped])
    @test view2.time==[0.0,1.0]
    @test PosView("v2";tag=7,elements=PosElement[stepped],time=[5.0,6.0]).time==
          [5.0,6.0]
end

@testset "read_pos parses all record families and comments" begin
    text="""
    // leading comment
    # hash comment
    /* block
       comment */
    View "everything" {
    SP(1,2,3){0.5};
    VP(1,2,3){1,2,3};
    TP(1,2,3){1,0,0, 0,1,0, 0,0,1};
    SL(0,0,0, 1,0,0){0.1,0.9};
    VL(0,0,0, 1,0,0){1,0,0, 0,1,0};
    TL(0,0,0, 1,0,0){1,0,0,0,1,0,0,0,1, 2,0,0,0,2,0,0,0,2};
    ST(0,0,0, 1,0,0, 0,1,0){1,2,3};
    VT(0,0,0, 1,0,0, 0,1,0){1,0,0, 0,1,0, 0,0,1};
    TT(0,0,0, 1,0,0, 0,1,0){1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1,
                           1,0,0,0,1,0,0,0,1};
    SQ(0,0,0, 1,0,0, 1,1,0, 0,1,0){1,2,3,4};
    VQ(0,0,0, 1,0,0, 1,1,0, 0,1,0){1,0,0, 0,1,0, 0,0,1, 1,1,1};
    TQ(0,0,0, 1,0,0, 1,1,0, 0,1,0){1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1,
                                  1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1};
    SS(0,0,0, 1,0,0, 0,1,0, 0,0,1){1,2,3,4};
    VS(0,0,0, 1,0,0, 0,1,0, 0,0,1){1,0,0, 0,1,0, 0,0,1, 1,1,1};
    TS(0,0,0, 1,0,0, 0,1,0, 0,0,1){1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1,
                                  1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1};
    SH(0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1){1,2,3,4,5,6,7,8};
    VH(0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1){
      1,0,0, 0,1,0, 0,0,1, 1,1,1, 1,0,0, 0,1,0, 0,0,1, 1,1,1};
    TH(0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1){
      1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1,
      1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1,
      1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1};
    SI(0,0,0, 1,0,0, 0,1,0, 0,0,1, 1,0,1, 0,1,1){1,2,3,4,5,6};
    VI(0,0,0, 1,0,0, 0,1,0, 0,0,1, 1,0,1, 0,1,1){
      1,0,0, 0,1,0, 0,0,1, 1,1,1, 1,0,0, 0,1,0};
    TI(0,0,0, 1,0,0, 0,1,0, 0,0,1, 1,0,1, 0,1,1){
      1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1,
      1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1};
    SY(0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1){1,2,3,4,5};
    VY(0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1){1,0,0, 0,1,0, 0,0,1, 1,1,1, 2,2,2};
    TY(0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1){
      1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1,
      1,0,0,0,1,0,0,0,1, 1,0,0,0,1,0,0,0,1};
    T2(1,2,10){"corner","label"};
    T3(1,2,3,8){"world"};
    INTERPOLATION_SCHEME{{0,0},{1,0}}{0.5,0.5};
    };
    """
    views=read_pos(IOBuffer(text))
    @test length(views)==1
    view=views[1]
    @test view.name=="everything" && view.tag==1
    @test length(view.elements)==24
    @test [element.kind for element in view.elements[[1,4,7,10,13,16,19,22]]]==
          [:point,:line,:triangle,:quadrangle,:tetrahedron,:hexahedron,:prism,
           :pyramid]
    @test [size(element.values,1) for element in view.elements[1:3]]==[1,3,9]
    @test size(view.elements[13].values)==(1,4,1)   # SS tetrahedron
    @test view.elements[4].values[1,:,1]==[0.1,0.9]
    @test view.elements[5].values[:,1,1]==[1.0,0.0,0.0]
    @test view.elements[6].values[:,2,1]==[2.0,0.0,0.0,0.0,2.0,0.0,0.0,0.0,2.0]
    @test length(view.texts)==2
    @test view.texts[1].dim==2 && view.texts[1].strings==["corner","label"]
    @test view.texts[2].coords==(1.0,2.0,3.0) && view.texts[2].style==8.0
    @test length(view.interpolations)==1 &&
          view.interpolations[1][1]==[0.0 0.0;1.0 0.0] &&
          view.interpolations[1][2]==[0.5;0.5;;]
    @test view.time==[0.0]
end

@testset "read_pos multi-view, TIME, and single-quoted names" begin
    text="""
    View "first" { SP(0,0,0){1}; };
    View 'second' {
    TIME{0,5};
    SL(0,0,0, 1,0,0){1,2, 3,4};
    };
    View "third" {};
    """
    views=read_pos(IOBuffer(text))
    @test length(views)==3
    @test [view.tag for view in views]==[1,2,3]
    @test views[1].name=="first" && views[2].name=="second"
    @test views[2].time==[0.0,5.0]
    @test size(views[2].elements[1].values)==(1,2,2)
    @test views[2].elements[1].values[1,2,2]==4.0
    @test isempty(views[3].elements) && isempty(views[3].time)
    paren_time="""
    View "v" {
    TIME(0,5);
    SL(0,0,0, 1,0,0){1,2, 3,4};
    };
    """
    @test read_pos(IOBuffer(paren_time))[1].time==[0.0,5.0]
end

@testset "read_pos rejects malformed and unsupported input" begin
    io(text)=read_pos(IOBuffer(text))
    @test_throws ArgumentError io("View \"v\" { SL(0,0,0, 1,0,0){1,2}; }")
    @test_throws ArgumentError io("View \"v\" { SL(0,0,0, 1,0,0){1,2}; }; extra")
    @test_throws ArgumentError io("Vieu \"v\" { };")
    # SP has one value per node, so {1,2} is a legal 2-step record, not a
    # malformed count.
    @test size(io("View \"v\" { SP(0,0,0){1,2}; };")[1].elements[1].values,3)==2
    @test_throws ArgumentError io("View \"v\" { SL(0,0,0){1,2}; };")
    @test_throws ArgumentError io("View \"v\" { XX(0,0,0){1}; };")
    # SL2 is a supported order-2 record (3 geometry + 3 value nodes).
    @test io("View \"v\" { SL2(0,0,0, 1,0,0, 0.5,0,0){1,2,3}; };"
        )[1].elements[1].suffix=="2"
    @test_throws ArgumentError io(
        "View \"v\" { SL2(0,0,0, 1,0,0){1,2}; };")
    @test_throws ArgumentError io("View \"v\" { TIME{0,5}; SP(0,0,0){1}; };")
    @test_throws ArgumentError io("View \"v\" { TIME{}; SP(0,0,0){1}; };")
    @test_throws ArgumentError io("View \"v\" { TIME{0}; TIME{0}; };")
    @test_throws ArgumentError io("View \"v\" { SL(0,0,0,1,0,0){1,2, 3}; };")
    @test_throws ArgumentError io("View \"v\" { SL(0,0,0,1,0,0){1,2, 3,4}; " *
        "SL(0,0,0,1,0,0){1,2}; };")
    @test_throws ArgumentError io("View \"v\" { SP(0,0,0){1e400}; };")
    @test_throws ArgumentError io("View \"v\" { SP(0,0,0){nan}; };")
    @test_throws ArgumentError io("View \"v\" { SP(0,0,0){1};; };")
    @test_throws ArgumentError io("View \"v\" { T2(1,2){\"a\"}; };")
    @test_throws ArgumentError io("View \"v\" { T2(1,2,10){}; };")
    @test_throws ArgumentError io("View \"v\" { T2(1,2,10){\"a\"; };")
    @test_throws ArgumentError io("View \"v\" { INTERPOLATION_SCHEME{{1}}; };")
    @test_throws ArgumentError io("View \"unterminated { };")
    @test_throws ArgumentError io("View \"v\" { SP(0,0,0){1}; }; garbage")
    @test_throws ArgumentError io("/* never closed")
    @test_throws ArgumentError io("View \"v\" { SP(0,0,/0){1}; };")
    @test read_pos(IOBuffer(""))==PosView[]
    @test read_pos(IOBuffer(" \t\r\n# only comments\n"))==PosView[]

    bounded="View \"v\" { SP(0,0,0){1}; };"
    @test_throws ArgumentError read_pos(IOBuffer(bounded);max_views=0)
    @test_throws ArgumentError read_pos(IOBuffer(bounded);max_elements=0)
    @test_throws ArgumentError read_pos(IOBuffer(bounded);max_values=2)
    @test_throws ArgumentError read_pos(IOBuffer(bounded);max_name_bytes=0)
    timed="View \"v\" { TIME{0,1,2}; };"
    @test_throws ArgumentError read_pos(IOBuffer(timed);max_steps=2)
    @test_throws ArgumentError read_pos(IOBuffer(timed);max_file_bytes=3)
    @test_throws ArgumentError read_pos(
        IOBuffer("View \"v\" { T2(1,2,10){\"abcd\"}; };");max_string_bytes=3)
    @test_throws ArgumentError read_pos(IOBuffer(bounded);max_views=-1)
    @test_throws ArgumentError read_pos(IOBuffer(bounded);max_views=1.5)
    @test_throws ArgumentError read_pos(IOBuffer(bounded);max_views=true)
end

@testset "write_pos deterministic output and round-trip" begin
    elements=PosElement[
        PosElement(:point,[1.0;2.0;3.0;;],[0.5;;]),
        PosElement(:line,[0.0 1.0;0.0 0.0;0.0 0.0],[0.2 0.8]),
        PosElement(:triangle,[0.0 1.0 0.0;0.0 0.0 1.0;0.0 0.0 0.0],
                   reshape(collect(1.0:9.0),3,3)),
        PosElement(:tetrahedron,
                   [0.0 1.0 0.0 0.0;0.0 0.0 1.0 0.0;0.0 0.0 0.0 1.0],
                   [0.1 0.2 0.3 0.4])]
    view=PosView("mixed";tag=4,elements=elements,
                 texts=PosText[PosText(2,(5.0,6.0,0.0),12.0,["tag"])])
    io=IOBuffer()
    write_pos(io,view)
    expected="View \"mixed\" {\n" *
        "SP(1,2,3){0.5};\n" *
        "SL(0,0,0,1,0,0){0.2,0.8};\n" *
        "VT(0,0,0,1,0,0,0,1,0){1,2,3,4,5,6,7,8,9};\n" *
        "SS(0,0,0,1,0,0,0,1,0,0,0,1){0.1,0.2,0.3,0.4};\n" *
        "T2(5,6,12){\"tag\"};\n" *
        "};\n"
    @test String(take!(io))==expected

    # write → read → write is a fixed point.
    again=read_pos(IOBuffer(expected))
    io2=IOBuffer();write_pos(io2,again)
    @test String(take!(io2))==expected

    timed=PosView("steps";elements=PosElement[
        PosElement(:line,[0.0 1.0;0.0 0.0;0.0 0.0],
                   reshape([0.2,0.8,0.4,1.6],1,2,2))],
        time=[0.0,5.0])
    io3=IOBuffer();write_pos(io3,timed)
    text3=String(take!(io3))
    @test occursin("TIME{0,5};",text3)
    @test occursin("SL(0,0,0,1,0,0){0.2,0.8,0.4,1.6};",text3)
    rt=read_pos(IOBuffer(text3))
    @test rt[1].time==[0.0,5.0]
    @test rt[1].elements[1].values==timed.elements[1].values
    io4=IOBuffer();write_pos(io4,rt)
    @test String(take!(io4))==text3

    # Legacy scalar View converts to SP records; Vector of views writes too.
    legacy=View("legacy",[0.0 2.0;0.0 0.0;0.0 0.0],[0.25,0.75])
    io5=IOBuffer();write_pos(io5,legacy)
    @test String(take!(io5))==
        "View \"legacy\" {\nSP(0,0,0){0.25};\nSP(2,0,0){0.75};\n};\n"
    io6=IOBuffer();write_pos(io6,[view,timed])
    @test String(take!(io6))==expected*text3
    @test_throws ArgumentError write_pos(IOBuffer(),Any["bogus"])

    # Path-based I/O.
    path=joinpath(mktempdir(),"out.pos")
    write_pos(path,rt)
    @test read(path,String)==text3
    reread=read_pos(path)
    @test reread[1].name=="steps" && reread[1].time==[0.0,5.0]
    @test_throws ArgumentError read_pos(path;max_file_bytes=filesize(path)-1)
    @test_throws ArgumentError read_pos(joinpath(mktempdir(),"missing.pos"))
end

@testset "postview_field builds field views" begin
    view=read_pos(IOBuffer("""
        View "sizes" { SL(0,0,0, 2,0,0){0.2,0.8}; };
        """))[1]
    field=postview_field(view;use_closest=false)
    @test field isa PostViewField
    @test field_value(field,0.5,0.0,0.0)≈0.35
    @test field_value(field,3.0,0.0,0.0)==GMSH_MAX_SIZE

    tensor_view=read_pos(IOBuffer("""
        View "metrics" {
        TL(0,0,0, 2,0,0){
        1,0,0,0,1,0,0,0,1, 4,0,0,0,4,0,0,0,4};
        };
        """))[1]
    tensor_field=postview_field(tensor_view)
    @test tensor_field isa PostViewAnisoField
    @test_throws ArgumentError postview_field(PosView("empty"))
    # Mixed component kinds follow Gmsh precedence: vector records win over
    # scalar records, tensor records win over both.
    mixed=PosView("mixed";elements=PosElement[
        PosElement(:line,[0.0 2.0;0.0 0.0;0.0 0.0],[0.2 0.8]),
        PosElement(:line,[0.0 2.0;0.0 0.0;0.0 0.0],[3.0 0.0;0.0 4.0;0.0 0.0])])
    mixed_field=postview_field(mixed;use_closest=false)
    @test field_value(mixed_field,1.0,0.0,0.0)==2.5
end
