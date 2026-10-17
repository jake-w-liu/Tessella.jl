# `model.view` alias/combine/interpolation endpoints — shared view data
# (`addAlias`), grouped view merging (`combine`) and per-family
# interpolation matrix storage (`setInterpolationMatrices`).
module ViewAliasChecks
using Test, Tessella
const A=Tessella.API

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

# Scalar triangle row: 3 nodes -> 9 coordinates then 3 values per step.
st(values...) = vcat([0.,1,0, 0,0,1, 0,0,0], collect(Float64, values))

@testset "view.add_alias shares data, keeps options per-view" begin
    A.initialize()
    try
        A.option("General.Terminal",0)

        src=A.view.add("source")
        A.view.add_list_data(src,"ST",2,
            vcat(st(1.,2,3),[0.,0,1, 0,1,1, 0,0,0, 4.,5,6]))
        types,data=A.view.get_list_data(src)

        alias=A.view.add_alias(src)
        @test alias!=src && alias in A.view.get_tags()
        @test A.view.get_index(alias)>=0
        @test A.view.get_list_data(alias)==(types,data)
        @test A.view.get_number(alias,"NbTimeStep")==1.0

        # Writes through the alias land on the shared data of the source.
        A.view.add_list_data(alias,"SQ",1,
            [0.,0,1,1, 0,0,0,1, 0,0,0,0, 7.,8,9,10])
        types2,data2=A.view.get_list_data(src)
        @test types2==["ST","SQ"]
        @test data2[2]==[0.,0,1,1, 0,0,0,1, 0,0,0,0, 7.,8,9,10]

        # The reverse direction: writes through the source reach the alias.
        A.view.add_list_data(src,"SL",1,[0.,1, 0,0, 0,0, 11.,12])
        types3,data3=A.view.get_list_data(alias)
        @test types3==["ST","SQ","SL"]
        @test data3[3]==[0.,1, 0,0, 0,0, 11.,12]

        # String list blocks share as well.
        A.view.add_list_data_string(alias,"ST",["e0","e1"],
            [0.,0,0, 1,1,0])
        got_types,got_strings,got_coords=
            A.view.get_list_data_strings(src,2,-1,true)
        @test got_types==["ST"] && got_strings==["e0","e1"]
        @test got_coords==[0.,0,0, 1,1,0]

        # Options and the display name stay per-view.
        @test A.view.get_string(alias,"Name")=="source"
        A.view.set_string(alias,"Name","renamed")
        @test A.view.get_string(src,"Name")=="source"
        A.view.set_number(alias,"Scale",2.5)
        @test A.view.get_number(src,"Scale")==0.0

        # copyOptions=true copies options at creation; false gives defaults.
        A.view.set_number(src,"InterpolationScheme",1.0)
        copied=A.view.add_alias(src,true)
        fresh=A.view.add_alias(src,false)
        @test A.view.get_number(copied,"InterpolationScheme")==1.0
        @test A.view.get_number(fresh,"InterpolationScheme")==0.0

        # An explicit tag is honored; an occupied tag is replaced.
        forced=A.view.add_alias(src,false,777)
        @test forced==777
        occupied=A.view.add("dup")
        @test A.view.get_list_data(occupied)==(String[],Vector{Float64}[])
        @test A.view.add_alias(src,false,occupied)==occupied
        @test A.view.get_list_data(occupied)[1]==["ST","SQ","SL"]

        # Alias-of-alias resolves to the shared canonical data.
        deep=A.view.add_alias(alias)
        A.view.add_list_data(deep,"SP",1,[0.,0,0, 42.])
        @test A.view.get_list_data(src)[1]==["ST","SQ","SL","SP"]

        # Unknown reference tag and bad argument types reject.
        @test occursin("unknown view tag",
            _errtext(()->A.view.add_alias(99999)))
        @test occursin("tag must be",
            _errtext(()->A.view.add_alias(src,false,"x")))
        @test occursin("copyOptions must be Bool",
            _errtext(()->A.view.add_alias(src,1)))

        # Removing the source leaves the shared data alive for the aliases.
        A.view.remove(src)
        @test !(src in A.view.get_tags())
        survived=A.view.get_list_data(alias)
        @test survived[1]==["ST","SQ","SL","SP"]
        @test A.view.get_list_data(deep)==survived
        A.view.add_list_data(alias,"SP",2,[0.,0,0,1., 1,0,0,2.])
        @test A.view.get_list_data(deep)[1][end]=="SP"
    finally
        A.finalize()
    end
end

@testset "view alias shares model data; step display stays per-view" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.occ.add_box(0,0,0,1,1,1)
        A.model.occ.synchronize()
        mname=A.model.get_current()

        src=A.view.add("model source")
        A.view.add_model_data(src,0,mname,"NodeData",[1,2],[5.0,6.0])
        alias=A.view.add_alias(src)
        mt,mtags,mdata=A.view.get_model_data(alias,0)
        @test mt==["NodeData"] && mtags==[1,2] && mdata==[5.0,6.0]

        # Writes through the alias extend the shared step range.
        A.view.add_model_data(alias,2,mname,"NodeData",[1],[9.0])
        @test A.view.get_number(src,"NbTimeStep")==3.0
        @test A.view.get_model_data(src,2)[3]==[9.0]
        @test A.view.get_number(alias,"NbTimeStep")==3.0

        # Displayed step selection stays per-view even over shared data.
        A.view.set_number(alias,"TimeStep",2)
        @test A.view.get_number(alias,"TimeStep")==2.0
        @test A.view.get_number(src,"TimeStep")==0.0
    finally
        A.finalize()
    end
end

@testset "view.combine elements merges payloads in tag order" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        v1=A.view.add("first")
        st1=st(1.,2,3)
        A.view.add_list_data(v1,"ST",1,st1)
        v2=A.view.add("second")
        st2=[0.,1,1, 0,0,1, 0,0,0, 4.,5,6]
        sq=[0.,0,1,1, 0,0,0,1, 0,0,0,0, 7.,8,9,10]
        A.view.add_list_data(v2,"ST",1,st2)
        A.view.add_list_data(v2,"SQ",1,sq)

        A.view.combine("elements","all")
        @test !(v1 in A.view.get_tags()) && !(v2 in A.view.get_tags())
        merged=only(A.view.get_tags())
        @test A.view.get_string(merged,"Name")=="all_Combine"
        types,data=A.view.get_list_data(merged)
        @test types==["ST","SQ"]
        @test data[1]==vcat(st1,st2)
        @test data[2]==sq
        @test A.view.get_number(merged,"NbTimeStep")==1.0
    finally
        A.finalize()
    end
end

@testset "view.combine groups by name and visibility" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        n1=A.view.add("grp");A.view.add_list_data(n1,"ST",1,st(1.))
        n2=A.view.add("grp");A.view.add_list_data(n2,"ST",1,st(2.))
        solo=A.view.add("other");A.view.add_list_data(solo,"ST",1,st(3.))

        A.view.combine("elements","name")
        tags=A.view.get_tags()
        @test !(n1 in tags) && !(n2 in tags) && solo in tags
        @test length(tags)==2
        named=only(setdiff(tags,[solo]))
        @test A.view.get_string(named,"Name")=="grp_Combine"
        @test A.view.get_list_data(named)[2][1]==vcat(st(1.),st(2.))
    finally
        A.finalize()
    end

    A.initialize()
    try
        A.option("General.Terminal",0)
        keep=A.view.add("visA");A.view.add_list_data(keep,"ST",1,st(4.))
        other=A.view.add("visB");A.view.add_list_data(other,"ST",1,st(5.))
        hidden=A.view.add("hid");A.view.add_list_data(hidden,"ST",1,st(6.))
        A.view.set_number(hidden,"Visible",0)

        A.view.combine("elements","visible")
        tags=A.view.get_tags()
        @test hidden in tags && !(keep in tags) && !(other in tags)
        @test length(tags)==2
        combined=only(setdiff(tags,[hidden]))
        @test A.view.get_string(combined,"Name")=="visible_Combine"
        @test A.view.get_list_data(combined)[2][1]==vcat(st(4.),st(5.))
    finally
        A.finalize()
    end

    A.initialize()
    try
        A.option("General.Terminal",0)
        o1=A.view.add("opt");A.view.add_list_data(o1,"ST",1,st(6.))
        A.view.set_number(o1,"Scale",3.25)
        o2=A.view.add("opt2");A.view.add_list_data(o2,"ST",1,st(7.))

        A.view.combine("elements","all",false,true)
        tags=A.view.get_tags()
        @test o1 in tags && o2 in tags && length(tags)==3
        combined=only(setdiff(tags,[o1,o2]))
        @test A.view.get_number(combined,"Scale")==3.25

        # A lone view is a group of one: never combined, never removed.
        A.view.remove(o2)
        before=A.view.get_tags()
        A.view.combine("elements","name")
        @test A.view.get_tags()==before

        @test occursin("what must be",
            _errtext(()->A.view.combine("bogus","all")))
        @test occursin("how must be",
            _errtext(()->A.view.combine("elements","bogus")))
        @test occursin("remove must be Bool",
            _errtext(()->A.view.combine("elements","all",1)))
    finally
        A.finalize()
    end
end

@testset "view.combine steps concatenates value payload" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        s1=A.view.add("s1");A.view.add_list_data(s1,"ST",1,st(1.,2,3))
        s2=A.view.add("s2");A.view.add_list_data(s2,"ST",1,st(4.,5,6, 7,8,9))
        A.view.combine("steps","all")
        merged=only(A.view.get_tags())
        types,data=A.view.get_list_data(merged)
        @test types==["ST"] && A.view.get_number(merged,"NbTimeStep")==3.0
        @test only(data)==vcat([0.,1,0, 0,0,1, 0,0,0],
                               [1.,2,3],[4.,5,6],[7.,8,9])

        # Reference view (lowest tag) sets the element count; mismatched
        # views contribute nothing to the type but are still removed.
        A.view.remove(merged)
        b=A.view.add("b");A.view.add_list_data(b,"ST",1,st(5.,5,5))
        c=A.view.add("c");A.view.add_list_data(c,"ST",1,st(9.,9,9))
        a=A.view.add("a");A.view.add_list_data(a,"ST",2,vcat(st(0.,0,0),st(1,1,1)))
        A.view.combine("steps","all")
        merged2=only(A.view.get_tags())
        types2,data2=A.view.get_list_data(merged2)
        @test types2==["ST"]
        @test only(data2)==vcat([0.,1,0, 0,0,1, 0,0,0],[5.,5,5],[9.,9,9])
        @test A.view.get_number(merged2,"NbTimeStep")==2.0
    finally
        A.finalize()
    end
end

@testset "view.combine reports per-group failures without removing views" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        st1=st(1.,2,3)
        st2=vcat([0.,1,0, 0,0,1, 0,0,0],[1.,2,3, 4,5,6])
        one=A.view.add("one step");A.view.add_list_data(one,"ST",1,st1)
        two=A.view.add("two steps");A.view.add_list_data(two,"ST",1,st2)
        A.view.combine("elements","all")
        @test occursin("different number of time steps",
            A.NATIVE_MSG_LAST_ERROR[])
        @test one in A.view.get_tags() && two in A.view.get_tags()
    finally
        A.finalize()
    end

    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.occ.add_box(0,0,0,1,1,1);A.model.occ.synchronize()
        mname=A.model.get_current()
        mv=A.view.add("modelled")
        A.view.add_model_data(mv,0,mname,"ElementData",[1],[1.5])
        listed=A.view.add("listed")
        A.view.add_list_data(listed,"ST",1,st(1.,2,3))
        before=A.view.get_tags()
        A.view.combine("elements","all")
        @test occursin("hybrid list/mesh-based datasets",
            A.NATIVE_MSG_LAST_ERROR[])
        @test A.view.get_tags()==before
    finally
        A.finalize()
    end

    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.occ.add_box(0,0,0,1,1,1);A.model.occ.synchronize()
        mname=A.model.get_current()
        mv=A.view.add("modelled")
        A.view.add_model_data(mv,0,mname,"ElementData",[1],[1.5])
        mv2=A.view.add("modelled2")
        A.view.add_model_data(mv2,3,mname,"ElementData",[1],[2.5])

        A.view.combine("elements","all")
        @test occursin("not implemented for this type of data",
            A.NATIVE_MSG_LAST_ERROR[])
        @test mv in A.view.get_tags() && mv2 in A.view.get_tags()

        # Step combination renumbers all populated steps contiguously.
        A.view.combine("steps","all")
        combined=only(A.view.get_tags())
        @test !(mv in A.view.get_tags())
        @test A.view.get_number(combined,"NbTimeStep")==2.0
        @test A.view.get_model_data(combined,0)[3]==[1.5]
        @test A.view.get_model_data(combined,1)[3]==[2.5]
        @test A.view.get_string(combined,"Name")=="all_Combine"
    finally
        A.finalize()
    end
end

@testset "view.set_interpolation_matrices stores and validates" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        v=A.view.add("interp")
        coef=Float64[i+j for i in 0:2 for j in 0:2]
        exp=Float64[0,0,0, 1,0,0, 0,1,0]
        A.view.set_interpolation_matrices(v,"Triangle",3,coef,exp)
        stored=Tessella.API.SESSION_VIEWS[v].interp
        @test length(stored[Int32(3)])==2
        # Row-major storage: F[i,j]=coef[3*(i-1)+j], P[i,j]=exp[3*(i-1)+j].
        @test stored[Int32(3)][1][2,3]==coef[3*(2-1)+3]
        @test stored[Int32(3)][2][2,2]==exp[3*(2-1)+2]
        @test size(stored[Int32(3)][1])==(3,3)
        @test size(stored[Int32(3)][2])==(3,3)

        # Geometry interpolation appends the two extra matrices.
        A.view.set_interpolation_matrices(v,"Line",2,ones(4),zeros(6),
                                          2,ones(4),zeros(6))
        @test length(Tessella.API.SESSION_VIEWS[v].interp[Int32(2)])==4

        # Re-setting replaces; non-positive d removes.
        A.view.set_interpolation_matrices(v,"Triangle",2,ones(4),zeros(6))
        @test size(Tessella.API.SESSION_VIEWS[v].interp[Int32(3)][1])==(2,2)
        A.view.set_interpolation_matrices(v,"Triangle",-1,Float64[],Float64[])
        @test !haskey(Tessella.API.SESSION_VIEWS[v].interp,Int32(3))

        # Stored on the shared data when called through an alias.
        alias=A.view.add_alias(v)
        A.view.set_interpolation_matrices(alias,"Prism",1,ones(1),zeros(3))
        @test haskey(Tessella.API.SESSION_VIEWS[v].interp,Int32(7))
        @test !haskey(Tessella.API.SESSION_VIEWS[alias].interp,Int32(7))

        # Rejections.
        @test occursin("unknown view tag",
            _errtext(()->A.view.set_interpolation_matrices(
                99999,"Line",2,ones(4),zeros(6))))
        @test occursin("Unknown element family type",
            _errtext(()->A.view.set_interpolation_matrices(
                v,"Circle",2,ones(4),zeros(6))))
        @test occursin("Wrong number of coefficients",
            _errtext(()->A.view.set_interpolation_matrices(
                v,"Line",2,ones(5),zeros(6))))
        @test occursin("Wrong number of exponents",
            _errtext(()->A.view.set_interpolation_matrices(
                v,"Line",2,ones(4),zeros(5))))
        @test occursin("Wrong number of coefficients",
            _errtext(()->A.view.set_interpolation_matrices(
                v,"Line",2,ones(4),zeros(6),2,ones(5),zeros(6))))
        @test occursin("dGeo must be",
            _errtext(()->A.view.set_interpolation_matrices(
                v,"Line",2,ones(4),zeros(6),1.5,ones(1),zeros(3))))
    finally
        A.finalize()
    end
end
end
