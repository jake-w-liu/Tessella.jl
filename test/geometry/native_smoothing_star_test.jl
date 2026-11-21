module NativeSmoothingStarChecks
using Test, Tessella
const M=Tessella.Model
const atan2=Tessella.GmshLibm._gm87_atan2

# The Laplacian smoothing of the native planar producer measures the star of
# every node six times per pass (parameter area and worst `qmTriangle::eta`).
# The measure now evaluates one emulated x87 `atan2` where a single corner
# decides it, and reuses an accepted trial's measure for the next trial. Both
# are only allowed to save work: the numbers must be the ones the plain
# formula gives, bit for bit, because the admitted moves decide the mesh.
@testset "native smoothing star measure equals the plain formula" begin
    # `qmTriangle::eta` as upstream writes it: all three angles.
    function plain_eta(world,t)
        p=ntuple(i->(world[1,t[i]],world[2,t[i]],world[3,t[i]]),3)
        angles=ntuple(3) do i
            e=p[i].-p[mod1(i+1,3)];f=p[mod1(i+2,3)].-p[mod1(i+1,3)]
            cross=(e[2]*f[3]-e[3]*f[2],e[3]*f[1]-e[1]*f[3],e[1]*f[2]-e[2]*f[1])
            180.0*atan2(sqrt(cross[1]^2+cross[2]^2+cross[3]^2),
                        e[1]*f[1]+e[2]*f[2]+e[3]*f[3])/pi
        end
        return 1.0-abs(60.0-min(min(angles[1],angles[2]),angles[3]))*(1.0/60.0)
    end
    function plain_star(state,incident)
        area=0.0;minimum=1.0
        for (_,cell) in incident
            a,b,c=state.triangles[cell]
            area+=0.5*abs((state.uv[1,b]-state.uv[1,a])*(state.uv[2,c]-state.uv[2,a])-
                          (state.uv[1,c]-state.uv[1,a])*(state.uv[2,b]-state.uv[2,a]))
            minimum=min(plain_eta(state.world,state.triangles[cell]),minimum)
        end
        return area,minimum
    end
    # A hexagonal star around node 1; `ring(k)` places the six neighbours.
    function star(centre,ring)
        uv=zeros(2,7);uv[:,1].=centre
        for k in 1:6
            uv[:,k+1].=ring(k)
        end
        world=vcat(uv,zeros(1,7))
        triangles=NTuple{3,Int32}[(1,k+1,mod1(k+1,6)+1) for k in 1:6]
        state=M._NativeSurfaceRecombine(world,copy(uv),fill((2,Int32(1)),7),
            UInt64.(1:7),triangles,UInt64.(1:6),NTuple{4,Int32}[],UInt64[],
            trues(7),UInt64(7),(u,v)->(u,v,0.0),"test")
        return state,[(false,cell) for cell in 1:6]
    end
    regular(k)=(cospi((k-1)/3),sinpi((k-1)/3))
    cases=(
        # One corner decides the star: the single-`atan2` path.
        star((0.13,-0.07),k->regular(k).*(1.0+0.11*k)),
        # An exactly regular star: every corner ties, so all are evaluated.
        star((0.0,0.0),regular),
        # Two mirror-image cells tie for the smallest angle.
        star((0.0,0.2),k->(regular(k)[1],regular(k)[2]*0.8)),
        # A needle next to ordinary cells.
        star((0.9,0.02),regular))
    for (state,incident) in cases
        @test M._native_recombine_star(state,incident)===plain_star(state,incident)
        for cell in 1:6
            @test M._native_recombine_triangle_eta(state,state.triangles[cell])===
                plain_eta(state.world,state.triangles[cell])
        end
    end
    # The fast path is taken where it may be and declined on a tie.
    @test M._native_recombine_triangle_star(cases[1]...)!==nothing
    @test M._native_recombine_triangle_star(cases[2]...)===nothing
    # A trial returns the measure at the trial position for the next one.
    state,incident=cases[1]
    actual=(state.uv[1,1],state.uv[2,1]);trial=(0.05,0.01)
    admitted,measure=M._native_recombine_laplace_admit!(state,Int32(1),actual,
        trial,incident)
    @test (state.uv[1,1],state.uv[2,1])==actual
    state.uv[:,1].=trial;state.world[1:2,1].=trial
    @test measure===plain_star(state,incident)
    @test admitted isa Bool

    # The same measure through the node-to-triangle table of the
    # triangles-only loop.
    for (state,incident) in cases
        cells=Int32[cell for (_,cell) in incident]
        @test M._native_recombine_star(state,cells,1,length(cells))===
            plain_star(state,incident)
    end
end

# A mesh of triangles is smoothed by a loop that keeps its node-to-triangle
# table in two arrays and allocates nothing per node (the general loop made
# 28 heap objects per triangle). It must move every node exactly as the
# general loop does.
@testset "triangles-only smoothing equals the general loop" begin
    # A 6x6 lattice over the unit square, interior nodes displaced, each
    # square cut along alternating diagonals; only interior nodes are owned
    # by the face.
    function lattice()
        n=6;uv=zeros(2,n*n);owners=Vector{Tuple{Int,Int32}}(undef,n*n)
        row(i,j)=Int32((j-1)*n+i)
        for j in 1:n,i in 1:n
            inner=1<i<n && 1<j<n
            shift=inner ? 0.07*sinpi(0.37*i*j) : 0.0
            uv[:,row(i,j)].=((i-1)/(n-1)+shift,(j-1)/(n-1)-0.6*shift)
            owners[row(i,j)]=inner ? (2,Int32(1)) : (1,Int32(1))
        end
        triangles=NTuple{3,Int32}[]
        for j in 1:n-1,i in 1:n-1
            a,b,c,d=row(i,j),row(i+1,j),row(i+1,j+1),row(i,j+1)
            if iseven(i+j)
                push!(triangles,(a,b,c),(a,c,d))
            else
                push!(triangles,(a,b,d),(b,c,d))
            end
        end
        # Node labels in an order unrelated to the rows: the loop visits
        # nodes and sums neighbours by label.
        labels=UInt64[UInt64(mod(17*k,n*n+1)) for k in 1:n*n]
        return M._NativeSurfaceRecombine(vcat(uv,zeros(1,n*n)),copy(uv),owners,
            labels,triangles,UInt64.(1:length(triangles)),NTuple{4,Int32}[],
            UInt64[],trues(n*n),UInt64(n*n),(u,v)->(u,v,0.0),"test")
    end
    frame=M._GmshNativePlaneFrame((0.0,0.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0),
                                  (0.0,0.0,1.0))
    for passes in (1,3)
        general=lattice();fast=lattice()
        M._native_recombine_laplace_general!(general,passes,frame)
        M._native_recombine_laplace!(fast,passes,frame)
        @test fast.uv==general.uv && fast.world==general.world
        @test fast.uv!=lattice().uv
    end
    # The suite runs interpreted, where every call allocates; the bound holds
    # for compiled code. The loop owns five arrays, whatever the node count.
    if Base.JLOptions().compile_enabled==1
        state=lattice()
        M._native_recombine_laplace!(lattice(),1,frame)
        @test (@allocated M._native_recombine_laplace!(state,3,frame))<=4096
    end
end
end
