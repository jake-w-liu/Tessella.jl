module HighOrderElasticKernelTests
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

end
