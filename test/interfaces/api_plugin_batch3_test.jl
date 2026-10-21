# `plugin` module UI01 batch 3: CutBox/CutGrid/CutParametric grid probes,
# StreamLines/Particles integrators, Triangulate/Tetrahedralize point
# remeshing, NearestNeighbor, Scal2Tens, Eigenvectors, Summation,
# MeshVolume, LongitudeLatitude, ShowNeighborElements, GaussPoints,
# AnalyseMeshQuality, MeshSizeFieldView and SimplePartition — verified
# against Gmsh 4.15.2 outputs on the fixtures in `fixtures_plugin3/`.
module PluginBatch3Checks
using Test, Tessella
const A=Tessella.API
const FX=joinpath(@__DIR__,"fixtures_plugin3")

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

@testset "CutGrid/CutBox/CutParametric grid probes" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.merge(joinpath(FX,"scalar_tri.pos"))
        src=last(A.view.get_tags())

        A.plugin.set_number("CutGrid","NumPointsU",3)
        A.plugin.set_number("CutGrid","NumPointsV",2)
        A.plugin.set_number("CutGrid","View",src)
        A.plugin.run("CutGrid")
        out=last(A.view.get_tags())
        types,data=A.view.get_list_data(out)
        @test types==["SQ"]
        @test length(data[1])==2*(3*4+4)          # 2 quads, 4 scalar values each
        @test data[1][1:4]≈[0.0,0.5,0.5,0.0]      # first quad x-block

        A.plugin.set_number("CutGrid","ConnectPoints",0)
        A.plugin.set_number("CutGrid","View",src)
        A.plugin.run("CutGrid")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["SP"]
        @test length(data[1])==6*4                # 6 points, 1 value each

        A.plugin.set_number("CutBox","NumPointsU",3)
        A.plugin.set_number("CutBox","NumPointsV",2)
        A.plugin.set_number("CutBox","NumPointsW",2)
        A.plugin.set_number("CutBox","Boundary",0)
        A.plugin.set_number("CutBox","View",src)
        A.plugin.run("CutBox")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["SH"]
        @test length(data[1])==2*(3*8+8)

        A.plugin.set_number("CutBox","ConnectPoints",0)
        A.plugin.set_number("CutBox","Boundary",1)
        A.plugin.set_number("CutBox","View",src)
        A.plugin.run("CutBox")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["SP"]
        @test length(data[1])==32*4               # boundary points only

        A.plugin.set_number("CutParametric","MinU",0.0)
        A.plugin.set_number("CutParametric","MaxU",1.0)
        A.plugin.set_number("CutParametric","NumPointsU",4)
        A.plugin.set_number("CutParametric","MinV",0.0)
        A.plugin.set_number("CutParametric","MaxV",1.0)
        A.plugin.set_number("CutParametric","NumPointsV",3)
        A.plugin.set_number("CutParametric","ConnectPoints",1)
        A.plugin.set_string("CutParametric","X","u")
        A.plugin.set_string("CutParametric","Y","v")
        A.plugin.set_string("CutParametric","Z","0")
        A.plugin.set_number("CutParametric","View",src)
        A.plugin.run("CutParametric")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["SQ"]
        @test length(data[1])==6*(3*4+4)          # (4-1)*(3-1) quads
    finally
        A.finalize()
    end
end

@testset "StreamLines and Particles integrators" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.merge(joinpath(FX,"vector_tet.pos"))
        src=last(A.view.get_tags())
        for opt in ("X0","Y0","Z0","X1","Y1","Z1","X2","Y2","Z2")
            A.plugin.set_number("StreamLines",opt,
                opt=="X0" ? 0.25 : opt=="Y0" ? 0.25 : opt=="Z0" ? 0.1 :
                opt=="X1" ? 0.5  : opt=="Y1" ? 0.25 : opt=="Z1" ? 0.1 :
                opt=="X2" ? 0.25 : opt=="Y2" ? 0.5  : 0.1)
        end
        A.plugin.set_number("StreamLines","NumPointsU",2)
        A.plugin.set_number("StreamLines","NumPointsV",1)
        A.plugin.set_number("StreamLines","MaxIter",5)
        A.plugin.set_number("StreamLines","View",src)
        A.plugin.run("StreamLines")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["VP"]
        @test length(data[1])==2*(3+15)           # 2 seeds × (xyz + 5×3)
        # Constant (1,0,0) field: displacement grows until the seed exits.
        # VP layout: 3 coord slots per record, then 5 iters × 3 comps.
        dx=data[1][4:3:16]
        @test issorted(dx) && dx[1]>0

        A.plugin.set_number("Particles","X0",0.25)
        A.plugin.set_number("Particles","Y0",0.25)
        A.plugin.set_number("Particles","Z0",0.1)
        A.plugin.set_number("Particles","NumPointsU",2)
        A.plugin.set_number("Particles","NumPointsV",1)
        A.plugin.set_number("Particles","A2",1)
        A.plugin.set_number("Particles","A1",0)
        A.plugin.set_number("Particles","A0",0)
        A.plugin.set_number("Particles","MaxIter",5)
        A.plugin.set_number("Particles","View",src)
        A.plugin.run("Particles")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["VP"]
        @test length(data[1])==2*(3+15)
        # Upstream 3-term Newmark recurrence on constant F=(1,0,0):
        # Xⁿ⁺¹ = 2Xⁿ - Xⁿ⁻¹ + DT²·F → displacements 0.01, 0.03, 0.06, 0.10, 0.15.
        @test data[1][4:3:16]≈[0.01,0.03,0.06,0.10,0.15] atol=1e-12
    finally
        A.finalize()
    end
end

@testset "Triangulate/Tetrahedralize/NearestNeighbor point ops" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.merge(joinpath(FX,"points5.pos"))
        src=last(A.view.get_tags())

        A.plugin.set_number("Triangulate","View",src)
        A.plugin.run("Triangulate")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["ST"]
        # Upstream jitters input points by 1e-12·lc before its Delaunay, so
        # the exact triangle count is not deterministic on degenerate sets;
        # Tessella's deterministic kernel may legitimately pick a different
        # valid triangulation. Invariants: arity-3 records whose vertices are
        # input points carrying the source scalar.
        @test length(data[1])%12==0
        srcdata=A.view.get_list_data(src)[2][1]
        ptval=Dict(Tuple(srcdata[4i+1:4i+3])=>srcdata[4i+4] for i in 0:6)
        for r in 0:(length(data[1])÷12)-1
            rec=data[1][12r+1:12r+12]
            for k in 0:2
                p=(rec[k+1],rec[3+k+1],rec[6+k+1])
                @test haskey(ptval,p) && ptval[p]==rec[9+k+1]
            end
        end

        A.plugin.set_number("Tetrahedralize","View",src)
        A.plugin.run("Tetrahedralize")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["SS"]
        @test length(data[1])%16==0
        vol=0.0
        for r in 0:(length(data[1])÷16)-1
            rec=data[1][16r+1:16r+16]
            p=[(rec[k+1],rec[4+k+1],rec[8+k+1]) for k in 0:3]
            for k in 0:3
                @test haskey(ptval,p[k+1]) && ptval[p[k+1]]==rec[12+k+1]
            end
            e1=collect(p[2]).-collect(p[1]);e2=collect(p[3]).-collect(p[1])
            e3=collect(p[4]).-collect(p[1])
            vol+=abs(e1[1]*(e2[2]*e3[3]-e2[3]*e3[2])-
                     e1[2]*(e2[1]*e3[3]-e2[3]*e3[1])+
                     e1[3]*(e2[1]*e3[2]-e2[2]*e3[1]))/6
        end
        # The emitted tets partition the convex hull (the corner tet: 1/6).
        @test vol≈1/6 atol=1e-12

        A.merge(joinpath(FX,"points4.pos"))
        srcnn=last(A.view.get_tags())
        A.plugin.set_number("NearestNeighbor","View",srcnn)
        A.plugin.run("NearestNeighbor")
        types,data=A.view.get_list_data(srcnn)
        @test types==["SP"]
        # In-place rewrite: every first component becomes sqrt(3)/2 here.
        vals=data[1][4:4:end]
        @test vals≈fill(sqrt(3)/2,6) atol=1e-12
    finally
        A.finalize()
    end
end

@testset "Scal2Tens/Summation/Eigenvectors view algebra" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.merge(joinpath(FX,"scalar_tri.pos"))
        v0=last(A.view.get_tags())
        A.merge(joinpath(FX,"scalar_b.pos"))
        v1=last(A.view.get_tags())

        A.plugin.set_number("Scal2Tens","NumberOfComponents",3)
        A.plugin.set_number("Scal2Tens","View0",v0)
        A.plugin.set_number("Scal2Tens","View1",v0)
        A.plugin.set_number("Scal2Tens","View2",v0)
        A.plugin.run("Scal2Tens")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["VT"]
        @test length(data[1])==4*(3*3+9)          # 4 tris × 3 comps

        A.plugin.set_number("Summation","View 0",v0)
        A.plugin.set_number("Summation","View 1",v1)
        A.plugin.set_string("Summation","Resuling View Name","sum")
        A.plugin.run("Summation")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["ST"]
        @test data[1][10:12]≈[5.0,7.0,8.0]        # node-wise sums, first tri

        A.merge(joinpath(FX,"tensor_tet.pos"))
        vt=last(A.view.get_tags())
        A.plugin.set_number("Eigenvectors","View",vt)
        A.plugin.run("Eigenvectors")
        # Three new vector views: Min/Mid/Max eigenvector fields.
        tail=last(A.view.get_tags(),3)
        for t in tail
            types,data=A.view.get_list_data(t)
            @test types==["VS"]
            @test length(data[1])==1*(3*4+12)
        end
    finally
        A.finalize()
    end
end

@testset "LongitudeLatitude in-place projection" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.merge(joinpath(FX,"lonlat.pos"))
        src=last(A.view.get_tags())
        A.plugin.set_number("LongitudeLatitude","View",src)
        A.plugin.run("LongitudeLatitude")
        types,data=A.view.get_list_data(src)
        @test types==["ST","VT"]
        st=data[1]
        # First triangle: nodes at lon 0/90/0 deg, lat 0/0/90 deg.
        @test st[1:9]≈[0.0,90.0,0.0, 0.0,0.0,90.0, 0.0,0.0,0.0]
        # Antimeridian unwrap on the second record: the -90° node moves to
        # +270° (upstream `xmax*180/pi - x > 180 → x += 360`).
        @test st[13]≈177.4208518 atol=1e-6
        @test st[14]≈2.5791482 atol=1e-6
        @test st[15]≈270.0 atol=1e-12
        vt=data[2]
        # Vector components rotated into lon/lat directions (oracle-checked).
        @test vt[10:12]≈[0.0,1.0,0.0] atol=1e-12
        @test vt[13:15]≈[-1.0,0.0,0.0] atol=1e-12
        @test vt[16:18]≈[1.0,0.0,0.0] atol=1e-12
    finally
        A.finalize()
    end
end

@testset "Model plugins on a unit box mesh" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.model.add("m")
        A.model.add_box(0,0,0,1,1,1)
        A.option("Mesh.MeshSizeMin",1.0)
        A.option("Mesh.MeshSizeMax",1.0)
        A.mesh.generate(3)

        # MeshVolume: unit box volume = 1, SP record at the model bbox center.
        A.plugin.set_number("MeshVolume","Dimension",3)
        A.plugin.set_number("MeshVolume","PhysicalGroup",-1)
        A.plugin.run("MeshVolume")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["SP"]
        @test data[1][1:4]≈[0.5,0.5,0.5,1.0]

        # GaussPoints on surfaces: every surface triangle contributes its
        # Order-0 Gauss point(s) valued by the element tag.
        A.plugin.set_number("GaussPoints","Order",0)
        A.plugin.set_number("GaussPoints","Dimension",2)
        A.plugin.run("GaussPoints")
        types,data=A.view.get_list_data(last(A.view.get_tags()))
        @test types==["SP"]
        nrec=length(data[1])÷4
        @test nrec>0
        for i in 0:nrec-1
            x,y,z=data[1][4i+1:4i+3]
            @test 0.0-1e-12<=x<=1.0+1e-12 && 0.0-1e-12<=y<=1.0+1e-12 &&
                  0.0-1e-12<=z<=1.0+1e-12
            @test min(abs(x),abs(x-1),abs(y),abs(y-1),abs(z),abs(z-1))<1e-8
        end

        # AnalyseMeshQuality with JacobianDeterminant prints stats; a good
        # box mesh keeps all tets visible and reports minJ>0.
        A.plugin.set_number("AnalyseMeshQuality","JacobianDeterminant",1)
        A.plugin.set_number("AnalyseMeshQuality","CreateView",1)
        A.plugin.set_number("AnalyseMeshQuality","DimensionOfElements",3)
        A.plugin.run("AnalyseMeshQuality")
        tags=A.view.get_tags()
        # The published view is an ElementData minJ/maxJ quality view.
        dtype,etags,vals=A.view.get_model_data(tags[end],0)
        @test dtype==["ElementData"]
        @test length(etags)==length(vals)>0
        @test all(v->0.0<v[1]<=1.0+1e-12,vals)

        # ShowNeighborElements: seed the first 3-D element; neighbors sharing
        # >2 vertices become visible, everything else hidden.
        A.plugin.set_number("ShowNeighborElements","NumLayers",1)
        _,etags3,_=A.mesh.get_elements(3,1)
        seed=etags3[1][1]
        A.plugin.set_number("ShowNeighborElements","Element1",seed)
        A.plugin.run("ShowNeighborElements")
        vis=Tessella.API.ELEMENT_VISIBILITY[]
        @test !isempty(vis)
        # NumLayers=1: only the seed and elements sharing >2 of its vertices
        # stay visible; everything else is hidden.
        @test get(vis,seed,Int32(0))==Int32(1)
        @test any(==(Int32(0)),values(vis))

        # SimplePartition: NumSlicesX=2 splits the box into 2 partitions.
        A.plugin.set_number("SimplePartition","NumSlicesX",2)
        A.plugin.run("SimplePartition")
        @test A.model.get_partitions(3,1)==[1,2]

        # MeshSizeFieldView on a NodeData view: every component is rewritten
        # with field(x,y,z) (upstream `putOnView` ignores `Component`).
        nd=A.view.add("sizes")
        node_ids,coords,_=A.mesh.get_nodes()
        A.view.add_model_data(nd,0,"m","NodeData",
            node_ids,[[0.0] for _ in node_ids])
        field=A.mesh.field.add("MathEval")
        A.mesh.field.set_string(field,"F","x + y + z + 0.25")
        A.plugin.set_number("MeshSizeFieldView","MeshSizeField",field)
        A.plugin.set_number("MeshSizeFieldView","View",nd)
        A.plugin.run("MeshSizeFieldView")
        _,ntags,vals=A.view.get_model_data(nd,0)
        xyz=Dict(node_ids[i]=>coords[3i-2:3i] for i in eachindex(node_ids))
        @test all(eachindex(ntags)) do i
            x,y,z=xyz[ntags[i]]
            vals[i][1]≈x+y+z+0.25
        end
        # Unknown field tags are rejected, mirroring upstream's
        # `Msg::Error("Unknown mesh size field %d")`.
        A.plugin.set_number("MeshSizeFieldView","MeshSizeField",0)
        @test occursin("field",lowercase(_errtext(
            ()->A.plugin.run("MeshSizeFieldView"))))
    finally
        A.finalize()
    end
end
end
