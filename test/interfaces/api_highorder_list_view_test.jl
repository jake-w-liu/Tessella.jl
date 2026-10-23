module HighOrderListViewChecks
using Test,Tessella,LinearAlgebra
const A=Tessella.API
const E=Tessella.Elements

# Minimal MSH 2.2 section splitter for the list-view mesh payload. Repeated
# $ElementNodeData sections keep their order.
function msh22(path)
    lines=readlines(path)
    @test lines[1]=="\$MeshFormat" && lines[2]=="2.2 0 8" &&
        lines[3]=="\$EndMeshFormat"
    sections=Pair{String,Vector{String}}[];index=4
    while index<=length(lines)
        name=lines[index];@test startswith(name,"\$")
        stop=findfirst(==(replace(name,"\$"=>"\$End")),lines[index+1:end])
        @test stop!==nothing
        push!(sections,name=>lines[index+1:index+stop-1]);index+=stop+1
    end
    return sections
end
sec1(sections,name)=only(s[2] for s in sections if s[1]==name)
secall(sections,name)=[s[2] for s in sections if s[1]==name]

function exercise(directory)
    @testset "SL2 import, queries and mesh/pos/txt export" begin
        A.initialize()
        try
            A.option("General.Terminal",0)
            path=joinpath(directory,"sl2.pos")
            Base.write(path,"View \"quad line\" {TIME{0.25,1.5};"*
                "SL2(0,0,0,1,0,0,.5,.05,0){1,2,3,4,5,6};SP(9,9,9){7,8};};\n")
            A.merge(path)
            tag=only(A.view.get_tags())
            # Suffixed record names stay verbatim; the flat payload keeps the
            # x/y/z coord blocks then node-major component values per step.
            types,data=A.view.get_list_data(tag)
            @test types==["SP","SL2"]
            @test data[1]==[9.,9,9,7,8]
            @test data[2]==[0.,1,.5, 0,0,.05, 0,0,0, 1,2,3,4,5,6]
            @test A.view.get_number(tag,"NbTimeStep")==2
            # .pos refuses like upstream rather than discarding the bound
            # scheme; the error is booked and the file is untouched.
            refused=joinpath(directory,"refused.pos")
            Base.write(refused,"sentinel")
            @test_throws ArgumentError A.view.write(tag,refused)
            @test occursin("interpolation matrices in old POS format",
                           A.logger.get_last_error())
            @test read(refused,String)=="sentinel"
            # .msh carries MSH 2.2 nodes/elements plus the bound scheme.
            msh=joinpath(directory,"sl2.msh")
            A.view.write(tag,msh)
            sections=msh22(msh)
            @test sec1(sections,"\$Nodes")==
                ["4","1 9 9 9","2 0 0 0","3 1 0 0","4 0.5 0.05 0"]
            elements=sec1(sections,"\$Elements")
            @test elements==["2","1 15 2 1 1 1","2 8 2 1 1 2 3 4"]
            scheme=sec1(sections,"\$InterpolationScheme")
            @test scheme[1:4]==["\"INTERPOLATION_SCHEME\"","1","2","2"]
            @test parse.(Int,split(scheme[5]))==[3,3]
            @test parse.(Int,split(scheme[9]))==[3,1]
            @test rstrip.(scheme[10:12])==["0","2","1"]
            datas=secall(sections,"\$ElementNodeData")
            @test length(datas)==2
            @test datas[1][1:3]==["2","\"quad line\"","\"INTERPOLATION_SCHEME\""]
            @test datas[1][5]=="0.25" && datas[2][5]=="1.5"
            @test datas[1][6:9]==["3","0","1","2"]
            @test datas[1][10:11]==["1 1 7","2 3 1 2 3"]
            mesh=E.read_mixed_msh(msh)
            @test sort!([Int(block.msh) for block in mesh.blocks])==[8,15]
            # Raw text keeps every interpolation node; CSV mirrors it.
            txt=joinpath(directory,"sl2.txt");csv=joinpath(directory,"sl2.csv")
            A.view.write(tag,txt);A.view.write(tag,csv)
            rows=readlines(txt)
            @test length(rows)==4 && occursin("0.05",rows[2])
            @test length(readlines(csv))==4
            # .vtk has no high-order corner layout: refused before open.
            vtk=joinpath(directory,"sl2.vtk")
            @test_throws ArgumentError A.view.write(tag,vtk)
            @test !isfile(vtk)
            # The decoded scheme feeds the postview field evaluator.
            parsed=only(Tessella.Post.read_pos(path))
            field=Tessella.Post.postview_field(parsed)
            @test Tessella.SizeField.field_value(field,.5,0,0)≈3
        finally;A.finalize();end
    end

    @testset "INTERPOLATION_SCHEME binds a family and reaches .msh" begin
        A.initialize()
        try
            A.option("General.Terminal",0)
            path=joinpath(directory,"scheme.pos")
            # 2-matrix order-2 value scheme bound to the triangle family;
            # geometry stays linear (ngeo=3), values live on 3 DOFs.
            Base.write(path,"View \"bound\" {"*
                "ST(0,0,0,1,0,0,0,1,0){1,2,3};"*
                "INTERPOLATION_SCHEME{{0,.5,-.5},{0,.5,.5},{1,-1,0}}{{0,0},{2,0},{1,1}};"*
                "ST(0,0,0,1,0,0,0,1,0){4,5,6};};\n")
            A.merge(path)
            tag=only(A.view.get_tags())
            types,data=A.view.get_list_data(tag)
            @test types==["ST"] && length(data[1])==24
            msh=joinpath(directory,"scheme.msh")
            A.view.write(tag,msh)
            sections=msh22(msh)
            @test sec1(sections,"\$Elements")==["2","1 2 2 1 1 1 2 3","2 2 2 1 1 1 2 3"]
            scheme=sec1(sections,"\$InterpolationScheme")
            @test scheme[3]=="3" && parse.(Int,split(scheme[5]))==[3,3] &&
                parse.(Int,split(scheme[9]))==[3,2]
            @test sec1(sections,"\$ElementNodeData")[3]=="\"INTERPOLATION_SCHEME\""
            refused=joinpath(directory,"s.pos")
            @test_throws ArgumentError A.view.write(tag,refused)
            @test !isfile(refused)
        finally;A.finalize();end
    end

    @testset "complete pyramid SY2 and Pyramid13 (MSH 19) serialization" begin
        A.initialize()
        try
            A.option("General.Terminal",0)
            path=joinpath(directory,"pyr.pos")
            # SY2 binds the complete 14-node Bergot order-2 pyramid:
            # 5 corners, 8 edge mids, one quadrilateral-face interior.
            verts=[(0.,0,0),(1,0,0),(1,1,0),(0,1,0),(.5,.5,1)]
            edges=[0.5.*(verts[a].+verts[b]) for (a,b) in
                ((1,2),(1,4),(1,5),(2,3),(2,5),(3,4),(3,5),(4,5))]
            nodes=vcat(verts,edges,[(.5,.5,0)])
            @assert length(nodes)==14
            coords=join((join(node,",") for node in nodes),",")
            Base.write(path,"View \"pyr\" {SY2($coords)"*
                "{1,2,3,4,5,6,7,8,9,10,11,12,13,14};};\n")
            A.merge(path)
            tag=only(A.view.get_tags())
            @test first(A.view.get_list_data(tag))==["SY2"]
            msh=joinpath(directory,"pyr.msh")
            A.view.write(tag,msh)
            sections=msh22(msh)
            @test sec1(sections,"\$Nodes")[1]=="14"
            @test startswith(sec1(sections,"\$Elements")[2],"1 14 2 1 1")
            mesh=E.read_mixed_msh(msh)
            @test only(mesh.blocks).msh==14 && size(only(mesh.blocks).nodes,2)==1
            # The incomplete quadratic pyramid (MSH 19) serializes through
            # API-set interpolation matrices plus canonical SY list data.
            tag13=A.view.add("pyr13")
            basis=A._native_view_writer_pyramid_basis(19)
            coef=basis.coefficients;expo=Float64.(basis.exponents)
            A.view.set_interpolation_matrices(tag13,"Pyramid",13,
                vec(coef'),vec(expo'),13,vec(coef'),vec(expo'))
            corners=[0. 1 1 0 .5;0 0 1 1 .5;0 0 0 0 1]
            mids=[0.5.*(corners[:,a].+corners[:,b]) for (a,b) in
                ((1,2),(1,4),(1,5),(2,3),(2,5),(3,4),(3,5),(4,5))]
            c13=hcat(corners,mids...)
            flat=vcat(vec(c13[1,:]),vec(c13[2,:]),vec(c13[3,:]),collect(1.:13.))
            A.view.add_list_data(tag13,"SY",1,flat)
            msh13=joinpath(directory,"pyr13.msh")
            A.view.write(tag13,msh13)
            sections13=msh22(msh13)
            @test sec1(sections13,"\$Elements")[2]=="1 19 2 1 1 "*join(1:13,' ')
            mesh13=E.read_mixed_msh(msh13)
            @test only(mesh13.blocks).msh==19
            # An order->2 incomplete pyramid has no Gmsh-consumable type:
            # refused before any file is touched.
            tag21=A.view.add("pyr21")
            A.view.set_interpolation_matrices(tag21,"Pyramid",21,
                vec(Matrix(1.0I,21,21)),zeros(63),21,
                vec(Matrix(1.0I,21,21)),zeros(63))
            A.view.add_list_data(tag21,"SY",1,zeros(84))
            bad=joinpath(directory,"bad.msh");Base.write(bad,"sentinel")
            @test_throws ArgumentError A.view.write(tag21,bad)
            @test read(bad,String)=="sentinel"
        finally;A.finalize();end
    end

    @testset "endpoint validation, visibility and aliases" begin
        A.initialize()
        try
            A.option("General.Terminal",0)
            tag=A.view.add("plain")
            A.view.add_list_data(tag,"SL",1,[0.,1,0,0,0,0,1,2])
            # Canonical list names only — upstream's "Unknown data type" error.
            @test_throws ArgumentError A.view.add_list_data(tag,"SL2",1,zeros(12))
            @test_throws ArgumentError A.view.add_list_data(tag,"NOPE",1,zeros(4))
            # Per-window visibility is retained display state.
            A.view.set_visibility_per_window(tag,0,0)
            @test A.view.get_number(tag,"VisibilityPerWindow:0")==0
            A.view.set_visibility_per_window(tag,1,0)
            @test A.view.get_number(tag,"VisibilityPerWindow:0")==1
            # "Visible" drives the combine-visible filter.
            A.view.set_number(tag,"Visible",0)
            second=A.view.add("other")
            A.view.add_list_data(second,"SL",1,[1.,2,0,0,0,0,3,4])
            A.view.combine("elements","visible")
            @test length(A.view.get_tags())==2
            A.view.set_number(tag,"Visible",1)
            # set/delete interpolation matrices toggles .pos refusal.
            A.view.set_interpolation_matrices(tag,"Line",2,
                [1.,0,0,1],[0,0,0,1,0,0])
            refused=joinpath(directory,"r.pos")
            @test_throws ArgumentError A.view.write(tag,refused)
            @test !isfile(refused)
            A.view.set_interpolation_matrices(tag,"Line",0,Float64[],Float64[])
            A.view.write(tag,refused)
            @test isfile(refused)
            # An alias exports the owner's data under its own name.
            alias=A.view.add_alias(tag)
            out=joinpath(directory,"alias.msh")
            A.view.write(alias,out)
            sections=msh22(out)
            @test sec1(sections,"\$Elements")==["1","1 1 2 1 1 1 2"]
            @test last(sec1(sections,"\$ElementNodeData"))=="1 2 1 2"
        finally;A.finalize();end
    end
end

directory=get(ENV,"TESSELLA_HIGHORDER_VIEW_DIR","")
if isempty(directory);mktempdir(exercise);else;mkpath(directory);exercise(directory);end
end
