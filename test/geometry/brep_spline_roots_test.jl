module BRepSplineRootsChecks
using Test,Tessella
const M=Tessella.Model

# A root search over a fitted section trace sampled nine points in every
# knot span of the trace: 12,850 evaluations a search for a circle fitted on
# a sphere chart, 421,450 for one sphere/box Boolean (730 s interpreted).
# `_brep_spl_linear_roots` gives the same roots, bit for bit, and samples
# only the spans whose poles do not keep the linear form clear of zero.

# a closed C1 cubic through `count` nodes of an off-centre ellipse
function ring(count)
    ts=Float64[];points=NTuple{2,Float64}[];slopes=NTuple{2,Float64}[]
    for i in 0:count
        angle=2π*i/count
        push!(ts,angle)
        push!(points,(0.3+1.2cos(angle),-0.2+0.7sin(angle)))
        push!(slopes,(-1.2sin(angle),0.7cos(angle)))
    end
    knots,degree,poles=M._hermite_bspline_2d(ts,points,slopes)
    return (spl2d=(knots,degree,poles),)
end

# a cubic with full-multiplicity inner knots: its derivative jumps there
const KINKED=(spl2d=(Float64[0,0,0,0,1,1,1,2,2,2,3,3,3,3],3,
    NTuple{2,Float64}[(0,0),(0.3,0.5),(0.7,0.5),(1,0),(1.2,-0.6),
                      (1.8,-0.6),(2,0),(2.5,0.9),(2.6,-0.4),(3,0.2)]),)

function both(pc,order,a,b,c,lo,hi,tol;tangents=true)
    calls=Ref(0)
    form=order==0 ?
        (t->(calls[]+=1;uv=M._brep_pc_eval(pc,t);a*uv[1]+b*uv[2]+c)) :
        (t->(calls[]+=1;d=M._brep_pc_tangent(pc,t);a*d[1]+b*d[2]+c))
    full=M._brep_scalar_roots(form,M._brep_pc_grid(pc,lo,hi),tol;tangents=tangents)
    full_calls=calls[];calls[]=0
    pruned=M._brep_spl_linear_roots(form,pc,order,Float64(a),Float64(b),
                                    Float64(c),lo,hi,tol;tangents=tangents)
    return full,pruned,full_calls,calls[]
end

@testset "roots of a linear form on a spline pcurve" begin
    pc=ring(48)
    # (order, a, b, c, lo, hi, tangents): crossings of a line, of a
    # horizontal and of a seam value, an extremum (tangent contact), a
    # level through a knot, no root at all, and a range inside the knots
    cases=((0,0.0,1.0,0.1,0.0,2π,false),
           (0,1.0,0.0,-0.3,0.0,2π,false),
           (0,0.6,-0.8,0.05,0.0,2π,true),
           (0,0.0,1.0,-0.5,0.0,2π,true),      # v = 0.5 touches the top
           (0,1.0,0.0,-1.5,0.0,2π,true),      # u = 1.5 is the node at angle 0
           (0,0.0,1.0,5.0,0.0,2π,true),
           (1,1.0,0.0,0.0,0.0,2π,false),      # du/dt = 0
           (1,0.0,1.0,0.0,0.4,5.1,false),
           (0,0.3,0.4,-0.2,1.3,4.9,true))
    for (order,a,b,c,lo,hi,tangents) in cases
        full,pruned,full_calls,pruned_calls=both(pc,order,a,b,c,lo,hi,1e-9;
                                                 tangents=tangents)
        @test pruned==full
        @test pruned_calls<=full_calls
        # 48 spans, 441 grid nodes; a root costs its span and a bisection
        isempty(full) && @test pruned_calls==0
        length(full)<=2 && @test pruned_calls<=250
    end
    # where the derivative jumps, a span's last node belongs to the next one
    for (order,a,b,c) in ((1,0.0,1.0,0.0),(1,1.0,0.0,-0.9),(0,0.0,1.0,0.0),
                          (0,0.0,1.0,-0.5))
        full,pruned,_,_=both(KINKED,order,a,b,c,0.0,3.0,1e-9;tangents=false)
        @test pruned==full
    end
    # a record that is no plain spline takes the full grid
    line=(lin2d=((0.0,0.0),(1.0,0.5)),)
    @test M._brep_spl_linear_roots(t->M._brep_pc_eval(line,t)[2]-0.25,line,0,
        0.0,1.0,-0.25,0.0,1.0,1e-9)==[0.5]
end

@testset "a bisection stops on neighbouring numbers" begin
    calls=Ref(0)
    roots=M._brep_scalar_roots(t->(calls[]+=1;t-1/3),[0.0,1.0],1e-9)
    @test roots==[1/3]
    @test calls[]<=60          # 2 grid values and some 53 halvings; was 82
end

@testset "periodic copies of a line across a fitted trace" begin
    # `_brep_ix_sec_lin` asks the trace's uv extent only which copies of the
    # line can meet it: the pole hull and one sample at each extreme pole
    # bracket the extent of the sample grid
    pc=ring(48)
    outer,inner=M._brep_spl_extent_bounds(pc,0.0,2π)
    umin=vmin=Inf;umax=vmax=-Inf
    for t in M._brep_pc_grid(pc,0.0,2π)
        uv=M._brep_pc_eval(pc,t)
        umin=min(umin,uv[1]);umax=max(umax,uv[1])
        vmin=min(vmin,uv[2]);vmax=max(vmax,uv[2])
    end
    @test outer[1]<=umin<=inner[1] && inner[2]<=umax<=outer[2]
    @test outer[3]<=vmin<=inner[3] && inner[4]<=vmax<=outer[4]
    # the crossings of the line v = 0.1 + 0.2u, with and without a u period
    edge=((0.0,0.1),(1.0,0.2))
    plain=M._brep_ix_sec_lin(pc,edge,-2.0,2.0,1e-9)
    @test length(plain)==2
    for t in plain
        uv=M._brep_pc_eval(pc,t)
        @test uv[2]≈0.1+0.2uv[1] atol=1e-9
    end
    wrapped=M._brep_ix_sec_lin(pc,edge,-2.0,2.0,1e-9;pu=4.0)
    @test length(wrapped)>=2 && issubset(plain,wrapped)
end
end
