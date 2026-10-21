# `plugin` module UI01: Isosurface (`scalar value == Value` level sets of the
# view's own data) on small list-view fixtures with analytically known
# sections — per-step views, ExtractVolume side selection, verbatim
# on-levelset elements — plus a t9-style .pos tet view run.
module PluginIsosurfaceChecks
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

# Records of a scalar list view as (type,node coordinates,node values).
function _records(types,data)
    arity=Dict("SP"=>1,"SL"=>2,"ST"=>3,"SQ"=>4,"SS"=>4,"SH"=>8,"SI"=>6,"SY"=>5)
    out=Tuple{String,Vector{NTuple{3,Float64}},Vector{Float64}}[]
    for (i,type) in enumerate(types)
        a=arity[type];stride=3a+a;vec=data[i]
        for element in 0:length(vec)÷stride-1
            base=element*stride
            pts=NTuple{3,Float64}[(vec[base+node],vec[base+a+node],
                                   vec[base+2a+node]) for node in 1:a]
            vals=Float64[vec[base+3a+node] for node in 1:a]
            push!(out,(type,pts,vals))
        end
    end
    return out
end

@testset "Isosurface plugin" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        latest()=last(A.view.get_tags())

        # Single tet, u = x+2y+3z on (0,0,0),(1,0,0),(0,1,0),(0,0,1) with
        # corner values (0,1,2,3): the u = 1.5 level set cuts edges
        # (1,3),(1,4),(2,3),(2,4) — a four-crossing quadrangle.
        tet=A.view.add("itet")
        A.view.add_list_data(tet,"SS",1,
            [0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1., 0.,1.,2.,3.])
        A.plugin.set_number("Isosurface","Value",1.5)
        A.plugin.set_number("Isosurface","View",tet)
        A.plugin.run("Isosurface")
        types,data=A.view.get_list_data(latest())
        @test types==["SQ"]
        recs=_records(types,data)
        @test length(recs)==1
        @test Set(recs[1][2])==Set([(0.,.75,0.),(0.,0.,.5),
                                    (.75,0.,.25),(.5,.5,0.)])
        @test recs[1][3]≈fill(1.5,4)
        # Every crossing satisfies the interpolated field exactly.
        for p in recs[1][2]
            @test p[1]+2p[2]+3p[3]≈1.5 atol=1e-12
        end

        # The same tet with a single corner below: the u = 0.5 level set is
        # the mid-edge triangle on the three edges leaving node 1.
        tet2=A.view.add("itet2")
        A.view.add_list_data(tet2,"SS",1,
            [0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1., 0.,1.,1.,1.])
        A.plugin.set_number("Isosurface","Value",0.5)
        A.plugin.set_number("Isosurface","View",tet2)
        A.plugin.run("Isosurface")
        types,data=A.view.get_list_data(latest())
        @test types==["ST"]
        recs=_records(types,data)
        @test length(recs)==1
        @test Set(recs[1][2])==Set([(0.5,0.,0.),(0.,0.5,0.),(0.,0.,0.5)])
        @test recs[1][3]≈fill(0.5,3)

        # Unit hexahedron with u = z (bottom 0, top 1): the z = 0.5 level
        # set through the six sub-tetrahedra tiles the mid-square once.
        hex=A.view.add("ihex")
        A.view.add_list_data(hex,"SH",1,
            [0.,1.,1.,0., 0.,1.,1.,0.,   # x1..x8
             0.,0.,1.,1., 0.,0.,1.,1.,   # y
             0.,0.,0.,0., 1.,1.,1.,1.,   # z
             zeros(4)..., ones(4)...])
        A.plugin.set_number("Isosurface","View",hex)
        A.plugin.run("Isosurface")
        types,data=A.view.get_list_data(latest())
        recs=_records(types,data)
        @test !isempty(recs)
        area=0.0
        for (_,pts,vals) in recs
            for p in pts
                @test p[3]≈0.5 atol=1e-12
            end
            @test vals≈fill(0.5,length(pts)) atol=1e-12
            for k in 3:length(pts)
                e1=pts[k-1].-pts[1];e2=pts[k].-pts[1]
                area+=0.5hypot(e1[2]*e2[3]-e1[3]*e2[2],
                               e1[3]*e2[1]-e1[1]*e2[3],
                               e1[1]*e2[2]-e1[2]*e2[1])
            end
        end
        @test area≈1.0 atol=1e-12
        corners=Set{NTuple{3,Float64}}()
        foreach(r->union!(corners,r[2]),recs)
        @test issubset(Set([(0.,0.,.5),(1.,0.,.5),(1.,1.,.5),(0.,1.,.5)]),
                       corners)

        # A surface triangle (u = x+y) cut at 0.5: a level *line* between the
        # crossing midpoints.
        tri=A.view.add("itri")
        A.view.add_list_data(tri,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 0.,1.,0.])
        A.plugin.set_number("Isosurface","View",tri)
        A.plugin.run("Isosurface")
        types,data=A.view.get_list_data(latest())
        @test types==["SL"]
        recs=_records(types,data)
        @test length(recs)==1
        @test Set(recs[1][2])==Set([(0.5,0.,0.),(0.5,0.5,0.)])
        @test recs[1][3]≈[0.5,0.5]

        # A unit quadrangle with u = x+y = 0.5: each sub-triangle emits a
        # segment; their union is the straight x+y = 0.5 line from (0.5,0)
        # to (0,0.5) through the shared diagonal crossing.
        qua=A.view.add("iqua")
        A.view.add_list_data(qua,"SQ",1,
            [0.,1.,1.,0., 0.,0.,1.,1., 0.,0.,0.,0., 0.,1.,2.,1.])
        A.plugin.set_number("Isosurface","View",qua)
        A.plugin.run("Isosurface")
        types,data=A.view.get_list_data(latest())
        @test types==["SL"]
        recs=_records(types,data)
        @test length(recs)==2
        allpts=NTuple{3,Float64}[]
        for (_,pts,vals) in recs
            for p in pts
                @test p[1]+p[2]≈0.5 atol=1e-12
                push!(allpts,p)
            end
            @test vals≈[0.5,0.5]
        end
        @test Set(allpts)==Set([(0.5,0.,0.),(0.,0.5,0.),(0.25,0.25,0.)])

        # A 1-D line element emits a level *point*.
        lin=A.view.add("ilin")
        A.view.add_list_data(lin,"SL",1,[0.,2., 0.,0., 0.,0., 0.,1.])
        A.plugin.set_number("Isosurface","View",lin)
        A.plugin.run("Isosurface")
        types,data=A.view.get_list_data(latest())
        @test types==["SP"]
        @test data[1]≈[1.,0.,0.,0.5]

        # Off-levelset data emits nothing; an element lying entirely on the
        # level set re-emits verbatim.
        flat=A.view.add("iflat")
        A.view.add_list_data(flat,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 4.,5.,6.])
        A.plugin.set_number("Isosurface","Value",0.0)
        A.plugin.set_number("Isosurface","View",flat)
        A.plugin.run("Isosurface")
        types,_=A.view.get_list_data(latest())
        @test isempty(types)
        flat2=A.view.add("iflat2")
        A.view.add_list_data(flat2,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 0.,0.,0.])
        A.plugin.set_number("Isosurface","View",flat2)
        A.plugin.run("Isosurface")
        types,data=A.view.get_list_data(latest())
        @test types==["ST"]
        @test data[1]≈[0.,1.,0., 0.,0.,1., 0.,0.,0., 0.,0.,0.]

        # ExtractVolume keeps whole elements on the selected side: the
        # straddling tet goes to neither, the all-above one to +1, the
        # all-below one to -1 — verbatim.
        pair=A.view.add("ipair")
        A.view.add_list_data(pair,"SS",3,
            [0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1., 0.,1.,1.,1.,   # straddle
             0.,1.,0.,0., 0.,0.,1.,0., 1.,1.,1.,2., 1.,1.,1.,1.,   # above
             0.,1.,0.,0., 0.,0.,1.,0., -2.,-2.,-2.,-1., 0.,0.,0.,-1.]) # below
        A.plugin.set_number("Isosurface","Value",0.5)
        A.plugin.set_number("Isosurface","ExtractVolume",1.0)
        A.plugin.set_number("Isosurface","View",pair)
        A.plugin.run("Isosurface")
        types,data=A.view.get_list_data(latest())
        @test types==["SS"]
        @test data[1]≈[0.,1.,0.,0., 0.,0.,1.,0., 1.,1.,1.,2., 1.,1.,1.,1.]
        A.plugin.set_number("Isosurface","ExtractVolume",-1.0)
        A.plugin.run("Isosurface")
        types,data=A.view.get_list_data(latest())
        @test types==["SS"]
        @test data[1]≈[0.,1.,0.,0., 0.,0.,1.,0., -2.,-2.,-2.,-1., 0.,0.,0.,-1.]
        A.plugin.set_number("Isosurface","ExtractVolume",0.0)

        # Two-step data emits one `*_Isosurface_<step>` view per step: the
        # level set follows each step's own values.
        two=A.view.add("itwo")
        A.view.add_list_data(two,"ST",1,
            [0.,1.,0., 0.,0.,1., 0.,0.,0., 0.,1.,0., 0.,0.,1.])
        A.plugin.set_number("Isosurface","Value",0.5)
        A.plugin.set_number("Isosurface","View",two)
        A.plugin.run("Isosurface")
        tags=A.view.get_tags()
        types,data=A.view.get_list_data(tags[end-1])
        @test types==["SL"]
        recs=_records(types,data)
        @test Set(recs[1][2])==Set([(0.5,0.,0.),(0.5,0.5,0.)])
        @test recs[1][3]≈[0.5,0.5]
        types,data=A.view.get_list_data(tags[end])
        @test types==["SL"]
        recs=_records(types,data)
        @test Set(recs[1][2])==Set([(0.,0.5,0.),(0.5,0.5,0.)])
        @test recs[1][3]≈[0.5,0.5]

        # Registry plumbing: upstream options accepted, unknown ones and
        # missing views reject loudly.
        for (opt,v) in (("RecurLevel",2.0),("TargetError",1e-3),
                        ("Visible",0.0))
            A.plugin.set_number("Isosurface",opt,v)
        end
        @test_throws ArgumentError A.plugin.set_number("Isosurface","Bogus",1.0)
        @test_throws ArgumentError A.plugin.set_string("Isosurface","View","x")
        @test occursin("unknown view tag",
            _errtext(()->(A.plugin.set_number("Isosurface","View",77);
                          A.plugin.run("Isosurface"))))
    finally
        A.finalize()
    end
end

@testset "Isosurface plugin on a .pos tet view (t9-style)" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        # A small parsed scalar-tetrahedra view like the tutorial's
        # view3.pos: two tets carrying u = x+y+z — one crossing into a
        # triangle, one into a quadrangle.
        pos=joinpath(mktempdir(),"iso.pos")
        open(pos,"w") do io
            print(io,"""View "iso_source" {
SS(0,1,0,0, 0,0,1,0, 0,0,0,1){0,1,1,1};
SS(0,0,0,-1, 0,1,0,0, 0,0,1,0){0,1,1,-1};
};
""")
        end
        A.merge(pos)
        src=last(A.view.get_tags())
        A.plugin.set_number("Isosurface","Value",0.67)
        A.plugin.set_number("Isosurface","View",src)
        A.plugin.run("Isosurface")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        recs=_records(types,data)
        @test !isempty(recs)
        @test sort!(unique!(first.(recs)))==["SQ","ST"]
        for (_,pts,vals) in recs
            for p in pts
                @test all(-1-1e-12.<p.<1+1e-12)
            end
            @test vals≈fill(0.67,length(pts)) atol=1e-12
        end
    finally
        A.finalize()
    end
end

end # module PluginIsosurfaceChecks
