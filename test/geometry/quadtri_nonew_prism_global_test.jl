module QuadTriNoNewPrismGlobalTests

using Test,Tessella
using LinearAlgebra: norm

const M=Tessella.Model
const CELL=(Int32(1),Int32(2),Int32(3))
const TRIANGLE=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.))
const PRISM=(TRIANGLE...,(0.,0.,1.),(1.,0.,1.),(0.,1.,1.))

function columns(points=TRIANGLE;levels=(0.,1.,2.),delta=(0.,0.,1.))
    result=Matrix{NTuple{3,Float64}}(undef,length(levels),length(points))
    for (layer,level) in pairs(levels),corner in eachindex(points)
        result[layer,corner]=ntuple(axis->points[corner][axis]+level*delta[axis],3)
    end
    return result
end

# Determinants over the represented coordinates, independent of the production
# adaptive predicate. These provide exact strict halfspace/separation evidence.
function exact_side(p,q,r,s)
    v=ntuple(i->ntuple(k->Rational{BigInt}((p,q,r)[i][k])-
        Rational{BigInt}(s[k]),3),3)
    a,b,c=v
    return sign(a[1]*(b[2]*c[3]-b[3]*c[2])-
        a[2]*(b[1]*c[3]-b[3]*c[1])+a[3]*(b[1]*c[2]-b[2]*c[1]))
end

function helical_columns(pitch;intervals=52)
    result=Matrix{NTuple{3,Float64}}(undef,intervals+1,3)
    for layer in 0:intervals
        level=layer/intervals
        sine,cosine=sincos(13pi/6*level)
        for corner in 1:3
            x,y,_=TRIANGLE[corner]
            result[layer+1,corner]=((x+1)*cosine-1,
                y+pitch*level,-(x+1)*sine)
        end
    end
    return result
end

# Literal Pri6 interpolation and differentiation, independently of both the
# global hull algorithm and the native reference-coordinate query machinery.
function prism_parameters(vertices,point)
    parameters=[.4,.25,.5]
    for _ in 1:30
        u,v,t=parameters
        barycentric=(1-u-v,u,v)
        weights=(barycentric[1]*(1-t),barycentric[2]*(1-t),
            barycentric[3]*(1-t),barycentric[1]*t,barycentric[2]*t,barycentric[3]*t)
        residual=[sum(weights[i]*vertices[i][axis] for i in 1:6)-point[axis]
            for axis in 1:3]
        norm(residual,Inf)<=2e-13 && return parameters
        gradients=(-(1-t),1-t,0.,-t,t,0.),
            (-(1-t),0.,1-t,-t,0.,t),
            (-barycentric[1],-barycentric[2],-barycentric[3],
             barycentric[1],barycentric[2],barycentric[3])
        jacobian=[sum(gradients[axis][i]*vertices[i][row] for i in 1:6)
            for row in 1:3,axis in 1:3]
        parameters-=jacobian\residual
    end
    error("independent Pri6 witness inversion did not converge")
end

@testset "True prism hull bounds and exact separating planes" begin
    @test M._extrude_nonew_hull_bounds(PRISM)==(0.,0.,0.,1.,1.,1.)
    separated=ntuple(i->(PRISM[i][1]+3,PRISM[i][2],PRISM[i][3]),6)
    plane=((1.,0.,0.),(0.,1.,0.),(1.,0.,1.))
    signs=[exact_side(plane...,point) for point in PRISM]
    @test length(unique(filter(!iszero,signs)))==1
    @test all(point->exact_side(plane...,point)==-only(unique(filter(!iszero,signs))),separated)
    @test M._extrude_nonew_hull_support_separates(PRISM,separated)
    @test M._extrude_nonew_hull_support_separates(separated,PRISM)
    touching=ntuple(i->(PRISM[i][1],PRISM[i][2],PRISM[i][3]+1),6)
    @test !M._extrude_nonew_hull_support_separates(PRISM,touching)
    @test !M._extrude_nonew_hull_support_separates(PRISM,PRISM)
    contained=ntuple(i->ntuple(k->PRISM[i][k]/4+.25,3),6)
    # The whole smaller prism is contained in x>=.25,y>=.25,x+y<=.75,
    # .25<=z<=.5, strictly inside the unit prism.
    @test all(p->p[1]>0 && p[2]>0 && p[1]+p[2]<1 && 0<p[3]<1,contained)
    @test !M._extrude_nonew_hull_support_separates(PRISM,contained)
    @test !M._extrude_nonew_hull_support_separates(contained,PRISM)
end

@testset "Three-point shared cap requires every strict nonshared corner" begin
    for scale in (1.,1e-150,1e150),direction in (1.,-1.),cell in (CELL,reverse(CELL))
        points=ntuple(i->ntuple(k->TRIANGLE[i][k]*scale,3),3)
        cols=columns(points;delta=(0.,0.,direction*scale))
        p,q,r=Tuple(cols[2,cell[k]] for k in 1:3)
        low=[exact_side(p,q,r,cols[1,cell[k]]) for k in 1:3]
        high=[exact_side(p,q,r,cols[3,cell[k]]) for k in 1:3]
        @test all(==(low[1]),low) && low[1]!=0
        @test high==.-low
        @test isnothing(M._extrude_nonew_hull_adjacent(cols,cell,1,"test"))
    end
    tilted=((0.,0.,0.),(1.,0.,.25),(0.,1.,.125))
    @test isnothing(M._extrude_nonew_hull_adjacent(columns(tilted),CELL,1,"test"))
    backward=columns(;levels=(0.,1.,.5))
    @test_throws r"opposite strict cap halfspaces" M._extrude_nonew_hull_adjacent(backward,CELL,1,"test")
    touching=columns();touching[1,3]=touching[2,3]
    @test_throws r"opposite strict cap halfspaces" M._extrude_nonew_hull_adjacent(touching,CELL,1,"test")
    degenerate=columns();degenerate[2,3]=degenerate[2,2]
    @test_throws r"strict shared-cap halfspace" M._extrude_nonew_hull_adjacent(degenerate,CELL,1,"test")
    square=(TRIANGLE[1],TRIANGLE[2],(1.,1.,0.),TRIANGLE[3])
    nonplanar=columns(square);nonplanar[2,4]=(0.,1.,1.125)
    @test_throws r"not exactly planar" M._extrude_nonew_hull_adjacent(nonplanar,
        (Int32(1),Int32(2),Int32(3),Int32(4)),1,"test")
end

@testset "Prism global certificate rejects overlap and accepts separated helices" begin
    for count in (1,2,52),delta in ((.25,-.125,1.),(.25,-.125,-1.))
        cols=columns(;levels=collect(range(0.,1.;length=count+1)),delta)
        @test length(M._extrude_nonew_corners(cols,CELL,1))==6
        @test isnothing(M._extrude_nonew_global_certify(cols,CELL,nothing,"test"))
    end
    for pitch in (-.1,.1)
        cols=helical_columns(pitch)
        witness=(.4935836460303578,.25+pitch/104,-.09789464416503871)
        for interval in (1,49)
            vertices=Tuple(cols[layer,corner] for layer in (interval,interval+1) for corner in 1:3)
            parameters=prism_parameters(vertices,witness)
            u,v,t=parameters
            @test min(u,v,1-u-v,t,1-t)>1e-9
            @test isnothing(M._extrude_nonew_prism_jacobian_certify(vertices,"test"))
        end
        @test_throws r"global P1 cell-hull disjointness is not certified" M._extrude_nonew_global_certify(cols,CELL,nothing,"test")
    end
    for pitch in (-3.,3.)
        @test isnothing(M._extrude_nonew_global_certify(helical_columns(pitch),CELL,nothing,"test"))
    end
    @test_throws r"triangle or quadrangle" M._extrude_nonew_global_certify(columns(),(Int32(1),Int32(2)),nothing,"test")
    @test_throws r"at least one interval" M._extrude_nonew_global_certify(columns(;levels=(0.,)),CELL,nothing,"test")
end

@testset "Prism hull work budgets fail precisely" begin
    cols=columns()
    @test_throws r"point-test budget is exceeded" M._extrude_nonew_hull_range_side(
        cols,CELL,1,1,3,Ref(0),"test")
    nodes=M._ExtrudeNoNewHullNode[]
    root=M._extrude_nonew_hull_tree!(nodes,cols,CELL,1,2,Ref(1000),"test")
    @test length(nodes)==3
    @test_throws r"hull-pair traversal budget is exceeded" M._extrude_nonew_hull_pair!(
        nodes,cols,CELL,root,root,Ref(0),Ref(1000),"test")
end

end
