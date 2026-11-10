module HighOrderElasticTests
using Test, LinearAlgebra, Tessella
const H=Tessella.HighOrder
const E=Tessella.Elements
const A=Tessella.API

function assemble(msh,maps)
    nodes=NTuple{3,Float64}[]
    index=Dict{NTuple{3,Float64},Int32}()
    reference=E.lagrange_nodes(msh)
    cells=Matrix{Int32}(undef,size(reference,2),length(maps))
    for (t,map) in enumerate(maps), k in axes(reference,2)
        p=map(Tuple(reference[:,k]))
        key=ntuple(d->p[d]+0.0,3)
        cells[k,t]=get!(index,key) do
            push!(nodes,key)
            Int32(length(nodes))
        end
    end
    return reduce(hcat,collect.(nodes)),cells,index
end

function tetra_fan_maps()
    corners=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.))
    center=Float64[.25,.25,.25]
    maps=Function[]
    for face in ((2,3,4),(1,3,4),(1,2,4),(1,2,3))
        a,b,c=(collect(corners[i]) for i in face)
        det(hcat(b-a,c-a,center-a))<0 && ((b,c)=(c,b))
        push!(maps,let a=a,b=b,c=c
            p->Tuple(a+(b-a)*p[1]+(c-a)*p[2]+(center-a)*p[3])
        end)
    end
    return maps
end

@testset "Elastic FEM restores straight references in quadratic families" begin
    cases=((9,[identity,p->(1-p[1],1-p[2],p[3])],(.5,.5,0.)),
           (10,[identity],(0.,0.,0.)),
           (11,tetra_fan_maps(),(.125,.125,.125)),
           (12,[identity],(0.,0.,0.)),
           (13,[identity,p->(-p[2],p[1],p[3])],(0.,.5,0.)),
           (14,[identity,p->(-p[1],p[2],-p[3])],(0.,0.,0.)))
    for (msh,maps,position) in cases
        coords,cells,index=assemble(msh,maps)
        original=copy(coords)
        v=Int(index[position])
        coords[:,v].+=(.03,-.02,.01)
        # Surface displacement remains in its plane.
        msh in (9,10) && (coords[3,v]=0.)
        stats,_=H.elastic_p2!(coords,((Int32(msh),cells),))
        @test stats[1]>=1
        @test coords≈original atol=2e-12
        @test all(H._p2opt_certify(H._p2opt_family(Int32(msh)).kind,coords,cells,t)
                  for t in axes(cells,2))
    end
end

@testset "Elastic FEM validates before mutation and respects ownership" begin
    coords,cells,index=assemble(10,[identity])
    center=Int(index[(0.,0.,0.)])
    coords[1,center]=.1
    before=copy(coords)
    for nu in (-1.,.5,Inf,NaN)
        @test_throws ArgumentError H.elastic_p2!(coords,((Int32(10),cells),);poisson=nu)
        @test isequal(coords,before)
    end
    stats,_=H.elastic_p2!(coords,((Int32(10),cells),);ownership=ones(Int,size(coords,2)))
    @test stats[1]==0
    @test isequal(coords,before)
    @test_throws ErrorException H.elastic_p2!(coords,((Int32(10),cells),);
        project=(v,p)->error("projection failure"))
    @test isequal(coords,before)
end

@testset "Legacy elastic overlays retain embedded curve mids" begin
    coords,cells,index=assemble(11,tetra_fan_maps())
    primaries=sort!(unique(vec(cells[1:4,:])))
    permutation=vcat(primaries,setdiff(collect(axes(coords,2)),primaries))
    inverse=invperm(permutation)
    coords=coords[:,permutation]
    cells=Int32.(inverse[cells])
    v=inverse[Int(index[(.125,.125,.125)])]
    coords[1,v]+=.01
    original=copy(coords)
    overlay=H.P2Mesh(coords,cells;_check_geometry=false)
    mids=[(3,Int32(7)) for _ in length(primaries)+1:size(coords,2)]
    mids[v-length(primaries)]=(1,Int32(9))
    stats=A._high_order_optimize_overlay_stages!(overlay,mids;
        elastic=true,banded=false)
    # The prescribed curve displacement drives the surrounding free DOFs.
    @test stats[1]>0
    @test isequal(overlay.coords[:,v],original[:,v])
    boundary=H._elastic_boundary_nodes(((Int32(11),cells),))
    @test all(node->isequal(overlay.coords[:,node],original[:,node]),boundary)
end

@testset "Elastic repair certifies the complete element" begin
    coords,cells,index=assemble(10,[identity])
    center=Int(index[(0.,0.,0.)])
    coords[1,center]=1.5
    @test !H._p2opt_certify(Val(10),coords,cells,1)
    stats,_=H.elastic_p2!(coords,((Int32(10),cells),))
    @test stats[1]>=1
    @test H._p2opt_certify(Val(10),coords,cells,1)
    @test coords[:,center]≈zeros(3) atol=1e-12
end

@testset "Elastic reference normalization and rigid-frame covariance" begin
    for scale in (1e-150,1.,1e150)
        map=p->(scale*(p[1]+.2p[2]),scale*(.8p[2]),scale*(.6p[2]))
        coords,cells,index=assemble(10,[map])
        original=copy(coords)
        center=Int(index[(0.,0.,0.)])
        coords[:,center].+=scale.*(.03,.016,.012)
        stats,_=H.elastic_p2!(coords,((Int32(10),cells),);poisson=-.5)
        @test stats[1]>=1
        @test maximum(abs.(coords-original))/scale<1e-12
        @test all(isfinite,coords)
    end
end

@testset "Thin elastic references retain finite weighted gradients" begin
    for aspect in (1e-100,1e-200,1e-300)
        coords=E.lagrange_nodes(12)
        coords[2:3,:].*=aspect
        cells=reshape(Int32.(1:27),27,1)
        center=findfirst(k->all(iszero,coords[:,k]),1:27)
        coords[:,center].+=(.03,.01aspect,-.02aspect)
        stats,_=H.elastic_p2!(coords,((Int32(12),cells),))
        @test stats[1]>0
        @test abs(coords[1,center])<1e-12
        @test maximum(abs.(coords[2:3,center]))/aspect<1e-12
    end
end

@testset "Disconnected elastic domains retain their own finite scale" begin
    reference=E.lagrange_nodes(12)
    big=copy(reference)
    big[1,:].+=10
    for smallscale in (1e-20,1e-150,1e-300)
        original=hcat(big,reference.*smallscale)
        coords=copy(original)
        cells=hcat(Int32.(1:27),Int32.(28:54))
        center=findfirst(k->all(iszero,reference[:,k]),1:27)
        coords[:,center].+=(.03,.01,-.02)
        coords[:,center+27].+=smallscale.*(.03,.01,-.02)
        stats,_=H.elastic_p2!(coords,((Int32(12),cells),))
        @test stats[1]>=2
        @test maximum(abs.(coords[:,center]-original[:,center]))<1e-12
        @test maximum(abs.((coords[:,center+27]-original[:,center+27])./smallscale))<1e-12
    end
end

@testset "Disconnected elastic Dirichlet loads remain scale covariant" begin
    reference=E.lagrange_nodes(10)
    center=findfirst(k->all(iszero,reference[:,k]),1:9)
    bottom=findfirst(k->reference[:,k]==[0.,-1.,0.],1:9)
    reference[2,bottom]-=.1
    isolated=copy(reference)
    H.elastic_p2!(isolated,((Int32(10),reshape(Int32.(1:9),9,1)),))
    @test abs(isolated[2,center])>.001
    for smallscale in (1e-150,1e-300)
        big=copy(reference)
        big[1,:].+=10
        coords=hcat(big,reference.*smallscale)
        cells=hcat(Int32.(1:9),Int32.(10:18))
        fixed=copy(coords[:,bottom+9])
        H.elastic_p2!(coords,((Int32(10),cells),))
        @test coords[:,center]-[10.,0.,0.]≈isolated[:,center] atol=2e-12
        @test coords[:,center+9]./smallscale≈isolated[:,center] atol=2e-12
        @test isequal(coords[:,bottom+9],fixed)
    end
end

@testset "Public elastic Poisson controls and boundary preservation" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(joinpath(@__DIR__,"..","fixtures","msh_high_order_elastic_disk.msh"))
        A.option("Mesh.HighOrderPoissonRatio",0.0)
        A.option("Mesh.HighOrderThresholdMin",.99)
        before=A.mesh.get_nodes()
        boundary=A.mesh.get_nodes(1,1,true)[1]
        A.mesh.optimize("HighOrderElastic",true,1)
        tags,coords,_=A.mesh.get_nodes()
        @test tags==before[1]
        @test any(coords .!= before[2])
        for tag in boundary
            k=findfirst(==(tag),tags)
            @test coords[3k-2:3k]==before[2][3k-2:3k]
        end
        @test coords[79:81]≈Float64[0.20175685171910943, -0.5320220614384539, 0.0] atol=1e-8 rtol=0
        @test coords[103:105]≈Float64[0.14053908224064962, -0.04310945992554645, 0.0] atol=1e-8 rtol=0
        @test coords[148:150]≈Float64[0.2829488895169657, 0.1670132678276144, 0.0] atol=1e-8 rtol=0
        @test coords[223:225]≈Float64[-0.2787775764037904, 0.26382097950255884, 0.0] atol=1e-8 rtol=0
        @test coords[286:288]≈Float64[-0.25347445425899795, -0.6682642207031829, 0.0] atol=1e-8 rtol=0
    finally
        A.finalize()
    end
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(joinpath(@__DIR__,"..","fixtures","msh_high_order_elastic_disk.msh"))
        A.option("Mesh.HighOrderPoissonRatio",0.2)
        A.option("Mesh.HighOrderThresholdMin",.99)
        before=A.mesh.get_nodes()
        boundary=A.mesh.get_nodes(1,1,true)[1]
        A.mesh.optimize("HighOrderElastic",true,1)
        tags,coords,_=A.mesh.get_nodes()
        @test tags==before[1]
        @test any(coords .!= before[2])
        for tag in boundary
            k=findfirst(==(tag),tags)
            @test coords[3k-2:3k]==before[2][3k-2:3k]
        end
        @test coords[79:81]≈Float64[0.20173595854019555, -0.5319348385167016, 0.0] atol=1e-8 rtol=0
        @test coords[103:105]≈Float64[0.1404843591192392, -0.043062207180750094, 0.0] atol=1e-8 rtol=0
        @test coords[148:150]≈Float64[0.2828303447155906, 0.16699882634023475, 0.0] atol=1e-8 rtol=0
        @test coords[223:225]≈Float64[-0.2787887160623528, 0.2637955620822512, 0.0] atol=1e-8 rtol=0
        @test coords[286:288]≈Float64[-0.2534039407533041, -0.6680901902845942, 0.0] atol=1e-8 rtol=0
    finally
        A.finalize()
    end
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(joinpath(@__DIR__,"..","fixtures","msh_high_order_elastic_disk.msh"))
        A.option("Mesh.HighOrderPoissonRatio",0.33)
        A.option("Mesh.HighOrderThresholdMin",.99)
        before=A.mesh.get_nodes()
        boundary=A.mesh.get_nodes(1,1,true)[1]
        A.mesh.optimize("HighOrderElastic",true,1)
        tags,coords,_=A.mesh.get_nodes()
        @test tags==before[1]
        @test any(coords .!= before[2])
        for tag in boundary
            k=findfirst(==(tag),tags)
            @test coords[3k-2:3k]==before[2][3k-2:3k]
        end
        @test coords[79:81]≈Float64[0.20171200258373592, -0.53180403073948, 0.0] atol=1e-8 rtol=0
        @test coords[103:105]≈Float64[0.1404272155128989, -0.04298153194539482, 0.0] atol=1e-8 rtol=0
        @test coords[148:150]≈Float64[0.28268691854538697, 0.16698082048202523, 0.0] atol=1e-8 rtol=0
        @test coords[223:225]≈Float64[-0.2788906213464496, 0.26380289232468535, 0.0] atol=1e-8 rtol=0
        @test coords[286:288]≈Float64[-0.2533495779728241, -0.667994074251216, 0.0] atol=1e-8 rtol=0
    finally
        A.finalize()
    end
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.open(joinpath(@__DIR__,"..","fixtures","msh_high_order_elastic_disk.msh"))
        A.option("Mesh.HighOrderPoissonRatio",0.49)
        A.option("Mesh.HighOrderThresholdMin",.99)
        before=A.mesh.get_nodes()
        boundary=A.mesh.get_nodes(1,1,true)[1]
        A.mesh.optimize("HighOrderElastic",true,1)
        tags,coords,_=A.mesh.get_nodes()
        @test tags==before[1]
        @test any(coords .!= before[2])
        for tag in boundary
            k=findfirst(==(tag),tags)
            @test coords[3k-2:3k]==before[2][3k-2:3k]
        end
        @test coords[79:81]≈Float64[0.2018250247328992, -0.5303706162729214, 0.0] atol=1e-8 rtol=0
        @test coords[103:105]≈Float64[0.14040321505746742, -0.04284194324930699, 0.0] atol=1e-8 rtol=0
        @test coords[148:150]≈Float64[0.2820644564622929, 0.16677276140885366, 0.0] atol=1e-8 rtol=0
        @test coords[223:225]≈Float64[-0.28081428179996115, 0.26415257398328157, 0.0] atol=1e-8 rtol=0
        @test coords[286:288]≈Float64[-0.25354843539618843, -0.6696500399334403, 0.0] atol=1e-8 rtol=0
    finally
        A.finalize()
    end
end
end
