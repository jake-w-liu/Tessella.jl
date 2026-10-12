using Test
using Tessella
const A=Tessella.API
const M=Tessella.Model

function _occ_setup()
    A.initialize()
    return A.model.occ
end

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

function _occ_curve(tag)
    m=A._with_model(identity)
    id=m.occ_kernel.bindings[(1,tag)]
    g=m.occ_kernel.shapes[id].geometry.geometry
    return g.curve,g.t0,g.t1
end

function _occ_uses(tag)
    m=A._with_model(identity)
    id=m.occ_kernel.bindings[(1,tag)]
    return m.occ_kernel.shapes[id].children
end

# Independent oracle: the unique clamped C2 cubic spline solved by the classic
# tridiagonal second-derivative system (different linear system, same unique
# interpolant as the B-spline collocation solve).
function _complete_spline(points,taus,T1,Tn,t)
    n=length(points)
    h=[taus[i+1]-taus[i] for i in 1:n-1]
    A_=zeros(n,n)
    b=zeros(n,3)
    A_[1,1]=2h[1];A_[1,2]=h[1]
    b[1,:].=6 .*([(points[2][k]-points[1][k])/h[1]-T1[k] for k in 1:3])
    for i in 2:n-1
        A_[i,i-1]=h[i-1];A_[i,i]=2(h[i-1]+h[i]);A_[i,i+1]=h[i]
        for k in 1:3
            b[i,k]=6((points[i+1][k]-points[i][k])/h[i]-
                     (points[i][k]-points[i-1][k])/h[i-1])
        end
    end
    A_[n,n-1]=h[n-1];A_[n,n]=2h[n-1]
    b[n,:].=6 .*([Tn[k]-(points[n][k]-points[n-1][k])/h[n-1] for k in 1:3])
    Mv=A_\b
    i=min(searchsortedlast(taus,clamp(t,taus[1],taus[n])),n-1)
    dt=t-taus[i]
    return ntuple(k->(Mv[i,k]*(taus[i+1]-t)^3/(6h[i])+Mv[i+1,k]*dt^3/(6h[i])+
        (points[i][k]/h[i]-Mv[i,k]*h[i]/6)*(taus[i+1]-t)+
        (points[i+1][k]/h[i]-Mv[i+1,k]*h[i]/6)*dt),3)
end

function _lagrange_tangent(points,taus,t)
    n=length(points)
    result=(0.0,0.0,0.0)
    for j in 1:n
        numerator=1.0
        for k in 1:n
            k==j || (numerator*=taus[j]-taus[k])
        end
        derivative=0.0
        for m in 1:n
            m==j && continue
            term=1.0
            for k in 1:n
                (k==j || k==m) && continue
                term*=t-taus[k]
            end
            derivative+=term/numerator
        end
        result=result.+derivative.*points[j]
    end
    return result
end

function _chord_params(points;closed=false)
    n=length(points)
    taus=[0.0]
    for i in 1:n-1
        push!(taus,taus[end]+sqrt(sum((points[i+1][k]-points[i][k])^2 for k in 1:3)))
    end
    closed && push!(taus,taus[end]+sqrt(sum((points[1][k]-points[n][k])^2 for k in 1:3)))
    return taus
end

@testset "occ.add_spline — non-periodic" begin
    occ=_occ_setup()
    points=[(0.0,0.0,0.0),(1.0,1.0,0.0),(2.0,0.5,0.0),(3.0,0.0,0.0)]
    pids=[occ.add_point(p...) for p in points]
    # Two points: degree-1 segment through the given vertices.
    line=occ.add_spline([pids[1],pids[4]])
    curve,t0,t1=_occ_curve(line)
    @test curve.curve.degree==1
    @test curve.curve.controls==[points[1],points[4]]
    m=A._with_model(identity)
    @test collect(_occ_uses(line))==
        [m.occ_kernel.bindings[(0,pids[1])],m.occ_kernel.bindings[(0,pids[4])]]
    a,=M._model_occ_curve_d2(curve,t0)
    b,=M._model_occ_curve_d2(curve,t1)
    @test collect(a)≈collect(points[1])
    @test collect(b)≈collect(points[4])
    @test occ.get_mass(1,line)≈3.0
    # Three points: degree-2 interpolant.
    quad=occ.add_spline([pids[1],pids[2],pids[3]])
    qcurve,qt0,qt1=_occ_curve(quad)
    @test qcurve.curve.degree==2
    qtaus=_chord_params(points[1:3])
    @test qt0≈qtaus[1] && qt1≈qtaus[3]
    for i in 1:3
        p,=M._model_occ_curve_d2(qcurve,qtaus[i])
        @test collect(p)≈collect(points[i]) atol=1e-12
    end
    # Four points: degree-3 clamped spline; unique interpolant cross-checked
    # against the independent tridiagonal system at 200 samples.
    spline=occ.add_spline(pids)
    scurve,st0,st1=_occ_curve(spline)
    @test scurve.curve.degree==3
    @test length(scurve.curve.controls)==6
    taus=_chord_params(points)
    @test st0≈taus[1] && st1≈taus[4]
    for i in 1:4
        p,=M._model_occ_curve_d2(scurve,taus[i])
        @test collect(p)≈collect(points[i]) atol=1e-12
    end
    T1=_lagrange_tangent(points,taus,taus[1])
    Tn=_lagrange_tangent(points,taus,taus[4])
    _,d1,_=M._model_occ_curve_d2(scurve,taus[1])
    @test collect(d1)≈collect(T1) atol=1e-10
    _,dn,_=M._model_occ_curve_d2(scurve,taus[4])
    @test collect(dn)≈collect(Tn) atol=1e-10
    worst=0.0
    for u in range(st0,st1,length=200)
        p,=M._model_occ_curve_d2(scurve,u)
        reference=_complete_spline(points,taus,T1,Tn,u)
        worst=max(worst,sqrt(sum((p[k]-reference[k])^2 for k in 1:3)))
    end
    @test worst<1e-10
    # C2 across interior knots.
    for i in 2:3
        _,_,dl=M._model_occ_curve_d2(scurve,taus[i]-1e-7)
        _,_,dr=M._model_occ_curve_d2(scurve,taus[i]+1e-7)
        @test maximum(abs.(collect(dr).-collect(dl)))<1e-4
    end
end

@testset "occ.add_spline — periodic" begin
    occ=_occ_setup()
    points=[(0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,1.0,0.0),(0.0,1.0,0.0)]
    pids=[occ.add_point(p...) for p in points]
    spline=occ.add_spline([pids;pids[1]])
    curve,t0,t1=_occ_curve(spline)
    @test curve.curve.degree==3
    taus=_chord_params(points;closed=true)
    @test t0≈taus[1] && t1≈taus[5]
    for i in 1:4
        p,=M._model_occ_curve_d2(curve,taus[i])
        @test collect(p)≈collect(points[i]) atol=1e-12
    end
    start,=M._model_occ_curve_d2(curve,t0)
    stop,=M._model_occ_curve_d2(curve,t1)
    @test collect(start)≈collect(stop) atol=1e-12
    _,ds,_=M._model_occ_curve_d2(curve,t0)
    _,de,_=M._model_occ_curve_d2(curve,t1)
    @test collect(ds)≈collect(de) atol=1e-10
    m=A._with_model(identity)
    id=m.occ_kernel.bindings[(0,pids[1])]
    @test collect(_occ_uses(spline))==[id,id]
    # Closed loop meshes as a wire.
    wire=occ.add_wire([spline])
    A.model.occ.synchronize()
    A.mesh.generate(1)
    @test length(A.mesh.get_nodes()[1])>0
    # Degree-1 closed spline: [a,b,a].
    q1=occ.add_point(0.0,0.0,1.0)
    q2=occ.add_point(1.0,0.0,1.0)
    closed=occ.add_spline([q1,q2,q1])
    cc,ct0,ct1=_occ_curve(closed)
    @test cc.curve.degree==1
    @test length(cc.curve.controls)==3
    @test collect(cc.curve.controls[3])≈collect(cc.curve.controls[1])
end

@testset "occ.add_spline — errors, atomicity, queries, transforms" begin
    occ=_occ_setup()
    m=A._with_model(identity)
    p1=occ.add_point(0.0,0.0,0.0)
    p2=occ.add_point(1.0,0.0,0.0)
    p3=occ.add_point(2.0,0.0,0.0)
    near=occ.add_point(1.0+5e-9,0.0,0.0)
    spline=occ.add_spline([p1,p2,p3])
    maximums=copy(m.occ_kernel.maxima)
    @test occursin("OCC.add_spline: OCC spline requires at least two points",
                   _errtext(()->occ.add_spline([p1])))
    @test occursin("OCC.add_spline: unknown OCC Point[99]",
                   _errtext(()->occ.add_spline([p1,99])))
    @test occursin("OCC.add_spline: OCC Curve[$spline] already exists",
                   _errtext(()->occ.add_spline([p1,p2];tag=spline)))
    @test occursin("OCC.add_spline: OCC spline points must be at least 1.0e-8 apart",
                   _errtext(()->occ.add_spline([p1,p2,near])))
    @test occursin("OCC.add_spline: could not interpolate spline",
                   _errtext(()->occ.add_spline([p1,p1])))
    @test m.occ_kernel.maxima==maximums
    # Queries + transforms on the new curve kind.
    @test occ.get_mass(1,spline)≈2.0
    xmin,ymin,zmin,xmax,ymax,zmax=occ.get_bounding_box(1,spline)
    @test xmin≈0.0 atol=2e-7
    @test xmax≈2.0 atol=2e-7
    @test zmin≈0.0 atol=2e-7
    @test zmax≈0.0 atol=2e-7
    occ.translate([(1,spline)],0.0,0.0,5.0)
    xmin,ymin,zmin,xmax,ymax,zmax=occ.get_bounding_box(1,spline)
    @test zmin≈5.0 atol=2e-7
    @test zmax≈5.0 atol=2e-7
    A.model.occ.synchronize()
    @test A.model.get_type(1,spline)=="BSpline"
    xmin,ymin,zmin,xmax,ymax,zmax=A.model.get_bounding_box(1,spline)
    @test zmin≈5.0 atol=1e-12
    @test zmax≈5.0 atol=1e-12
end
