module OCCNaturalSeamChecks
using Test,Tessella
using LinearAlgebra: norm
const A=Tessella.API
const M=Tessella.Model
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function session(f,name)
    A.initialize()
    try
        @test A.model.occ.import_shapes(joinpath(FIXTURES,name))==[(2,1)]
        A.model.occ.synchronize()
        return f(A.CURRENT[].surface_geometry[1])
    finally
        A.finalize()
    end
end

@testset "Natural tensor closure is independent of world translation" begin
    for shift in (1e9,-1e12)
        session("brep_native_natural_bspline.brep") do geometry
            @test !M._occ_natural_surface_closed(geometry)
            A.model.occ.translate([(2,1)],shift,0.,0.)
            A.model.occ.synchronize()
            translated=A.CURRENT[].surface_geometry[1]
            @test M._brep_uperiod(translated)==0
            @test M._brep_vperiod(translated)==0
            @test !M._occ_natural_surface_closed(translated)
            @test_throws ArgumentError M._validate_surface_loop(
                A.CURRENT[],[1],"translated open natural patch")
        end
    end
end

@testset "A seam certificate covers the complete polynomial domain" begin
    # The degree-nine displacement is 2^20*prod(v-j/8,j=0:8).
    # All nine formerly used samples coincide while the actual seam is open.
    session("brep_native_sampled_seam_gap.brep") do geometry
        for v in range(0.,1.;length=9)
            @test norm(A.model.get_value(2,1,[0.,v]).-
                A.model.get_value(2,1,[1.,v]))<1e-8
        end
        @test norm(A.model.get_value(2,1,[0.,1/16]).-
            A.model.get_value(2,1,[1.,1/16]))≈30.929946899414062 atol=1e-10
        @test M._brep_uperiod(geometry)==0
        @test !M._occ_natural_surface_closed(geometry)
        @test_throws ArgumentError M._validate_surface_loop(
            A.CURRENT[],[1],"polynomial seam gap")
    end
end

function patch(controls,weights,U,V,pu,pv)
    surface=Tessella.NURBS.NURBSSurface(pu,pv,U,V,controls,weights)
    return (occ=:tensor_bspline,patch=M._ModelNURBSSurface(:bspline,surface))
end

@testset "Exact seams support rational and nonclamped representations" begin
    # The endpoint homogeneous curves have different factors (1+v) and
    # (2-v), but both dehomogenize to the same Cartesian curve (v,0,0).
    controls=[(0.,0.,0.) (1/3,0.,0.) (1.,0.,0.);
              (0.,1.,0.) (.5,1.,0.) (1.,1.,0.);
              (0.,0.,0.) (2/3,0.,0.) (1.,0.,0.)]
    weights=[1. 1.5 2.;1. 1. 1.;2. 1.5 1.]
    rational=patch(controls,weights,[0.,0.,0.,1.,1.,1.],
        [0.,0.,0.,1.,1.,1.],2,2)
    @test M._model_nurbs_patch_seam_period(rational,1)==1.
    @test M._model_nurbs_patch_seam_period(rational,2)==0.
    # The active U endpoints average different neighboring control rows.
    # The first and last control rows themselves are deliberately different.
    ring=((1.,0.),(0.,1.),(-1.,0.),(2.,1.))
    controls2=[j==1 ? (0.,0.,-1.) : j==3 ? (0.,0.,1.) :
        (ring[i]...,0.) for i in 1:4,j in 1:3]
    nonclamped=patch(controls2,ones(4,3),[-2.,-1.,0.,1.,2.,3.,4.],
        [0.,0.,0.,1.,1.,1.],2,2)
    @test M._model_nurbs_patch_seam_period(nonclamped,1)==2.
    @test M._model_nurbs_patch_seam_period(nonclamped,2)==0.
    for v in (0.,.13,.63,1.)
        @test collect(first(M._model_nurbs_surface_d2(nonclamped,0.,v)))≈
            collect(first(M._model_nurbs_surface_d2(nonclamped,2.,v))) atol=1e-14
    end
end

@testset "A conic-reparametrized natural revolution retains both seams" begin
    session("brep_native_natural_revolution.brep") do geometry
        @test geometry.reparam.u!==nothing && geometry.reparam.v!==nothing
        @test M._brep_uperiod(geometry)==2π
        @test M._brep_vperiod(geometry)==2π
        @test M._occ_natural_surface_closed(geometry)
        @test M._validate_surface_loop(A.CURRENT[],[1],"natural revolution")===nothing
        for (u,v) in ((.3,.7),(2.,3.),(5.,4.))
            @test A.model.get_value(2,1,[u,v])≈
                [(1+.25cos(v))*cos(u),(1+.25cos(v))*sin(u),.25sin(v)] atol=2e-14
        end
        A.model.occ.convert_to_nurbs([(2,1)])
        A.model.occ.synchronize()
        @test M._occ_natural_surface_closed(A.CURRENT[].surface_geometry[1])
    end
end

@testset "A closed weighted natural patch retains a closed published mesh" begin
    for name in ("brep_native_weighted_closed.brep","brep_native_nonclamped_closed.brep")
    session(name) do geometry
        @test M._brep_uperiod(geometry)==1.
        @test M._brep_vperiod(geometry)==0.
        @test M._occ_natural_surface_closed(geometry)
        @test M._validate_surface_loop(A.CURRENT[],[1],"weighted closed patch")===nothing
        A.option("Mesh.MeshSizeMin",1.)
        A.option("Mesh.MeshSizeMax",1.)
        A.option("Mesh.Algorithm",1)
        A.mesh.generate(2)
        tags,xyz,parameters=A.mesh.get_nodes(2,1)
        @test A.model.get_value(2,1,parameters)≈xyz atol=2e-12
        families,_,groups=A.mesh.get_elements(2,1)
        @test families==Int32[2]
        incidence=Dict{Tuple{UInt64,UInt64},Int}()
        for cells in groups,i in 1:3:length(cells),pair in ((1,2),(2,3),(3,1))
            edge=minmax(cells[i+pair[1]-1],cells[i+pair[2]-1])
            incidence[edge]=get(incidence,edge,0)+1
        end
        @test !isempty(tags) && all(==(2),values(incidence))
        @test length(tags)-length(incidence)+sum(length,groups)÷3==2
        coordinates=reshape(xyz,3,:)
        rows=Dict(uid=>row for (row,uid) in enumerate(tags))
        cells=reduce(vcat,groups;init=UInt64[])
        mesh=Tessella.Mesh(coordinates;
            tris=reshape(Int32[rows[uid] for uid in cells],3,:))
        @test Tessella.validate(mesh).ok
        @test M._plc_self_intersection3(mesh)===nothing
        @test maximum(norm(coordinates[:,rows[a]].-coordinates[:,rows[b]])
            for (a,b) in keys(incidence))<=1+1e-12
        data=A.CURRENT[].meshing.mesh_authority.mesh.entity_data
        @test all(==((2,Int32(1))),data.node_entities)
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(0,0,1,0)
    end
    end
end
end
