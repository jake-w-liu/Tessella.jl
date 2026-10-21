# Tessella-side dump mirroring `oracle_gmsh.py`: runs the same plugin cases
# against the same fixtures and writes `tessella_out.txt`; `diff_compare.py`
# then diffs it against the Gmsh 4.15.2 `oracle_out.json`.
using Tessella
const A=Tessella.API
FX=joinpath(@__DIR__)

function view_name(tag)
    s=A.view.get_string(tag,"Name")
    return isempty(s) ? A.view.get_string(tag,"Name0") : s
end

function run_case(io,fixture,plugin,numbers,strings=Dict{String,String}())
    A.initialize()
    try
        A.option("General.Terminal",0)
        fixtures=fixture isa AbstractVector ? fixture : Any[fixture]
        srctags=Int[]
        for f in fixtures
            prev=Set(A.view.get_tags())
            A.merge(joinpath(FX,f))
            # A fixture file may declare several views; keep them in tag
            # order so oracle "view index i" maps to the i-th view.
            append!(srctags,sort!(collect(setdiff(A.view.get_tags(),prev))))
        end
        src=last(srctags)
        for (k,v) in numbers
            # Gmsh View options carry 0-based view indices; Tessella uses
            # view tags. Index i refers to the i-th merged fixture view.
            vv=v
            if (startswith(k,"View") || k=="OtherView") && isfinite(v)
                idx=Int(v)+1
                idx<=length(srctags) && (vv=Float64(srctags[idx]))
            end
            A.plugin.set_number(plugin,k,vv)
        end
        for (k,v) in strings
            A.plugin.set_string(plugin,k,v)
        end
        before=Set(A.view.get_tags())
        err=""
        try
            A.plugin.run(plugin)
        catch e
            err=sprint(showerror,e)
        end
        println(io,"CASE ",plugin," ",join(fixtures,"+"))
        isempty(err) || println(io,"ERROR ",err)
        for t in A.view.get_tags()
            t in before && !(plugin in ("NearestNeighbor","LongitudeLatitude")
                             && t==src) && continue
            types,datas=A.view.get_list_data(t)
            println(io,"VIEW ",t," name=\"",view_name(t),"\"")
            for (ty,d) in zip(types,datas)
                cnt=A.SESSION_VIEWS[t].list_counts[ty]
                println(io,"REC ",ty," ",cnt," ",join(string.(d),","))
            end
        end
    finally
        A.finalize()
    end
end

cases=[
    ("scalar_tri.pos","CutGrid",
     Dict("X0"=>0.,"Y0"=>0.,"Z0"=>0.,"X1"=>1.,"Y1"=>0.,"Z1"=>0.,
          "X2"=>0.,"Y2"=>1.,"Z2"=>0.,"NumPointsU"=>3.,"NumPointsV"=>2.,
          "ConnectPoints"=>1.,"View"=>0.)),
    ("scalar_tri.pos","CutGrid",
     Dict("X0"=>0.,"Y0"=>0.,"Z0"=>0.,"X1"=>1.,"Y1"=>0.,"Z1"=>0.,
          "X2"=>0.,"Y2"=>1.,"Z2"=>0.,"NumPointsU"=>3.,"NumPointsV"=>2.,
          "ConnectPoints"=>0.,"View"=>0.)),
    ("scalar_tri.pos","CutBox",
     Dict("X0"=>0.,"Y0"=>0.,"Z0"=>0.,"X1"=>1.,"Y1"=>0.,"Z1"=>0.,
          "X2"=>0.,"Y2"=>1.,"Z2"=>0.,"X3"=>0.,"Y3"=>0.,"Z3"=>1.,
          "NumPointsU"=>3.,"NumPointsV"=>2.,"NumPointsW"=>2.,
          "ConnectPoints"=>1.,"Boundary"=>0.,"View"=>0.)),
    ("scalar_tri.pos","CutBox",
     Dict("X0"=>0.,"Y0"=>0.,"Z0"=>0.,"X1"=>1.,"Y1"=>0.,"Z1"=>0.,
          "X2"=>0.,"Y2"=>1.,"Z2"=>0.,"X3"=>0.,"Y3"=>0.,"Z3"=>1.,
          "NumPointsU"=>3.,"NumPointsV"=>2.,"NumPointsW"=>2.,
          "ConnectPoints"=>0.,"Boundary"=>1.,"View"=>0.)),
    ("scalar_tri.pos","CutBox",
     Dict("X0"=>0.,"Y0"=>0.,"Z0"=>0.,"X1"=>1.,"Y1"=>0.,"Z1"=>0.,
          "X2"=>0.,"Y2"=>1.,"Z2"=>0.,"X3"=>0.,"Y3"=>0.,"Z3"=>1.,
          "NumPointsU"=>3.,"NumPointsV"=>2.,"NumPointsW"=>2.,
          "ConnectPoints"=>1.,"Boundary"=>1.,"View"=>0.)),
    ("scalar_tri.pos","CutGrid",
     Dict("X0"=>0.,"Y0"=>0.,"Z0"=>0.,"X1"=>1.,"Y1"=>0.,"Z1"=>0.,
          "X2"=>0.,"Y2"=>1.,"Z2"=>0.,"NumPointsU"=>3.,"NumPointsV"=>1.,
          "ConnectPoints"=>1.,"View"=>0.)),
    ("scalar_tri.pos","CutParametric",
     Dict("MinU"=>0.,"MaxU"=>1.,"NumPointsU"=>4.,"MinV"=>0.,"MaxV"=>1.,
          "NumPointsV"=>3.,"ConnectPoints"=>1.,"View"=>0.),
     Dict("X"=>"u","Y"=>"v","Z"=>"0")),
    ("vector_tet.pos","StreamLines",
     Dict("X0"=>0.25,"Y0"=>0.25,"Z0"=>0.1,"X1"=>0.5,"Y1"=>0.25,"Z1"=>0.1,
          "X2"=>0.25,"Y2"=>0.5,"Z2"=>0.1,"NumPointsU"=>2.,"NumPointsV"=>1.,
          "DT"=>0.1,"MaxIter"=>5.,"TimeStep"=>0.,"View"=>0.)),
    ("vector_scalar.pos","StreamLines",
     Dict("X0"=>0.25,"Y0"=>0.25,"Z0"=>0.1,"X1"=>0.5,"Y1"=>0.25,"Z1"=>0.1,
          "X2"=>0.25,"Y2"=>0.5,"Z2"=>0.1,"NumPointsU"=>2.,"NumPointsV"=>1.,
          "DT"=>0.1,"MaxIter"=>5.,"TimeStep"=>0.,"View"=>0.,"OtherView"=>1.)),
    ("vector_tet.pos","Particles",
     Dict("X0"=>0.25,"Y0"=>0.25,"Z0"=>0.1,"NumPointsU"=>2.,"NumPointsV"=>1.,
          "A2"=>1.,"A1"=>0.,"A0"=>0.,"DT"=>0.1,"MaxIter"=>5.,"TimeStep"=>0.,
          "View"=>0.)),
    ("points5.pos","Triangulate",Dict("View"=>0.)),
    ("points5.pos","Tetrahedralize",Dict("View"=>0.)),
    ("points4.pos","NearestNeighbor",Dict("View"=>0.)),
    ("scalar_tri.pos","Scal2Tens",
     Dict("NumberOfComponents"=>3.,"View0"=>0.,"View1"=>0.,"View2"=>0.)),
    ("tensor_tet.pos","Eigenvectors",
     Dict("ScaleByEigenvalues"=>1.,"View"=>0.)),
    ("lonlat.pos","LongitudeLatitude",Dict("View"=>0.)),
    (["scalar_tri.pos","scalar_b.pos"],"Summation",
     Dict("View 0"=>0.,"View 1"=>1.),
     Dict("Resuling View Name"=>"sum")),
]

open(joinpath(FX,"tessella_out.txt"),"w") do io
    for case in cases
        fixture,plugin,numbers=case[1],case[2],case[3]
        strings=length(case)>3 ? case[4] : Dict{String,String}()
        run_case(io,fixture,plugin,numbers,strings)
    end
end
println("tessella_out.txt written: ",length(cases)," cases")
