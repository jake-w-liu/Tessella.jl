module NativeCurvedRecombineChecks
using Test,Tessella,Serialization,SHA
const A=Tessella.API
const XY=[(-1.,0.),(-.7,-.7),(0.,-1.),(.7,-.7),(1.,0.),(.7,.7),(0.,1.),(-.7,.7),(-.3,0.),(0.,-.3),(.3,0.),(0.,.3)]
const CELLS=[9 9 10 10 11 12 12;10 1 2 4 5 6 8;11 2 3 5 6 7 1;12 10 4 11 12 8 9]
fingerprint(value)=begin io=IOBuffer();serialize(io,value);sha256(take!(io));end
function setup(kind;parameters=true,smoothing=0,boundary=false)
    A.initialize();A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    A.option("Mesh.Smoothing",smoothing);A.option("Mesh.RecombineOptimizeTopology",5)
    occ=A.model.occ;n=kind===:saddle ? 2 : 3
    control=Int[]
    for j in 1:n,i in 1:n
        x=-1+2*(i-1)/(n-1);y=-1+2*(j-1)/(n-1)
        z=kind===:saddle ? .2*x*y : .2*(1,-1,1)[i]+.1*(1,-1,1)[j]
        push!(control,occ.add_point(x,y,z))
    end
    face=occ.add_bspline_surface(control,n;degree_u=n-1,degree_v=n-1)
    occ.synchronize()
    xy=boundary ? [(-1.,-1.),(0.,-1.),(1.,-1.),(1.,0.),(1.,1.),(0.,1.),(-1.,1.),(-1.,0.),XY[9:12]...] : XY
    uv=Float64[a for (x,y) in xy for a in ((x+1)/2,(y+1)/2)]
    xyz=A.model.get_value(2,face,uv)
    if boundary
        edges=A.model.get_boundary([(2,face)],false,false)
        vertices=A.model.get_boundary([(2,face)],false,false,true)
        for i in 1:8
            point=xyz[3*i-2:3*i]
            if isodd(i)
                errors=[sum(abs2,A.model.get_value(d,t,Float64[]).-point) for (d,t) in vertices]
                d,t=vertices[argmin(errors)]
                A.mesh.add_nodes(d,t,[i],point)
            else
                candidates=[(A.model.get_value(d,t,[.5]),[.5]) for (d,t) in edges]
                errors=[sum(abs2,c[1].-point) for c in candidates]
                selected=argmin(errors);d,t=edges[selected]
                A.mesh.add_nodes(d,t,[i],point,candidates[selected][2])
            end
        end
    else
        outer=[occ.add_point(xyz[3*i-2:3*i]...) for i in 1:8]
        occ.synchronize()
        for i in 1:8;A.mesh.add_nodes(0,outer[i],[i],xyz[3*i-2:3*i]);end
    end
    A.mesh.add_nodes(2,face,collect(9:12),xyz[25:36],parameters ? uv[17:24] : Float64[])
    A.mesh.add_elements_by_type(face,3,collect(1:7),vec(CELLS))
    return face
end
# Public Gmsh 4.15.2 controls retain cells 2:7, consume node 11 and lift
# the surviving node at the average surface parameter, including on a
# quadratic paraboloid whose world-space chord midpoint is off the surface.
@testset "Native diamonds evaluate the actual nonlinear CAD surface" begin
    for kind in (:saddle,:paraboloid),parameters in (false,true)
        try
            face=setup(kind;parameters)
            before=[A.mesh.get_node(i) for i in 1:12]
            held=A.CURRENT[];digest=fingerprint(held)
            A.mesh.recombine()
            types,tags,cells=A.mesh.get_elements(2,face)
            expected=copy(CELLS[:,2:end]);expected[expected.==11].=9
            @test types==[3]
            @test tags==[UInt64.(2:7)]
            @test cells==[UInt64.(vec(expected))]
            @test sort(first(A.mesh.get_nodes()))==UInt64[i for i in 1:12 if i!=11]
            @test A.mesh.get_node(9)[1]≈A.model.get_value(2,face,[.5,.5]) atol=1e-9
            @test [A.mesh.get_node(i) for i in 1:8]==before[1:8]
            @test [A.mesh.get_node(i) for i in (10,12)]==before[[10,12]]
            parameter=A.mesh.get_node(9)[2]
            @test parameters ? isapprox(parameter,[.5,.5];atol=1e-9) : isempty(parameter)
            @test fingerprint(held)==digest
        finally
            A.finalize()
        end
    end
end
const BOUNDARY_CONTROL=Dict(
    :saddle=>[
        ([-1.6667334925779542e-05, -0.0015708059524401907, 5.236229782563144e-09],[0.4999916663325371, 0.4992145970237799]),
        ([0.28250187434767016, -0.31353320807710716, -0.017714743790404167],[0.6412509371738351, 0.3432333959614464]),
        ([-0.28567507660559954, 0.30630616754357964, -0.017500807575555947],[0.3571624616972002, 0.6531530837717898]),
    ],
    :paraboloid=>[
        ([0.0004017815550130255, -0.001206750126336864, 1.779102703263345e-07],[0.5002008907775065, 0.4993966249368315]),
        ([0.3090992303389566, -0.3241954118062478, 0.029618733342849334],[0.6545496151694783, 0.33790229409687617]),
        ([-0.31090888993330257, 0.3204154472780707, 0.029599473453352305],[0.34454555503334877, 0.6602077236390353]),
    ],
)
# Native Curve-owned vertices carry their Curve U. Smoothing uses that
# face's pcurve, including the actual nonlinear iso-Curves on this patch.
@testset "Native nonlinear boundary Curve relocation" begin
    for kind in (:saddle,:paraboloid)
        try
            face=setup(kind;parameters=true,smoothing=1,boundary=true)
            boundary=[A.mesh.get_node(i) for i in 1:8]
            held=A.CURRENT[];digest=fingerprint(held)
            A.mesh.recombine()
            types,tags,cells=A.mesh.get_elements(2,face)
            expected=copy(CELLS[:,2:end]);expected[expected.==11].=9
            @test types==[3]
            @test tags==[UInt64.(2:7)]
            @test cells==[UInt64.(vec(expected))]
            @test [A.mesh.get_node(i) for i in 1:8]==boundary
            @test sort(first(A.mesh.get_nodes()))==UInt64[i for i in 1:12 if i!=11]
            @test fingerprint(held)==digest
            for (node,(xyz,uv)) in zip((9,10,12),BOUNDARY_CONTROL[kind])
                actual,parameter,dimension,entity=A.mesh.get_node(node)
                @test actual≈xyz atol=2e-12
                @test parameter≈uv atol=2e-12
                @test (dimension,entity)==(2,face)
                @test actual≈A.model.get_value(2,face,parameter) atol=2e-12
            end
        finally
            A.finalize()
        end
    end
end
function setup_circle(stored)
  face=A.model.occ.add_disk(0,0,0,1,1);A.model.occ.synchronize()
  curve=A.model.get_boundary([(2,face)],false,false)[1][2]
  vertex=A.model.get_boundary([(1,curve)],false,false)[1][2]
  xy=[(cos(t),sin(t)) for t in (pi,5pi/4,3pi/2,7pi/4,0,pi/4,pi/2,3pi/4)]
  append!(xy,[(-.3,0.),(0.,-.3),(.3,0.),(0.,.3)])
  xyz=Float64[a for (x,y) in xy for a in (x,y,0.)]
  for i in 1:8
   d,t=i==5 ? (0,vertex) : (1,curve)
   point=xyz[3*i-2:3*i]
   par=d==0 || !stored ? Float64[] : A.model.get_parametrization(d,t,point)
   A.mesh.add_nodes(d,t,[i],point,par)
  end
  A.mesh.add_nodes(2,face,collect(9:12),xyz[25:end],A.model.get_parametrization(2,face,xyz[25:end]))
  A.mesh.add_elements_by_type(face,3,collect(1:7),[9,10,11,12,9,1,2,10,10,2,3,4,10,4,5,11,11,5,6,12,12,6,7,8,12,8,1,9])
  return face
end
const CIRCLE_CONTROL=Dict(
    true=>[
        ([-1.7575125118805831e-16, 3.0755417194061306e-17, 0.0],[-1.1165249497527285e-16, -5.66218834374304e-17]),
        ([-4.255547009734441e-16, -0.4142135729975656, 0.0],[-3.614559447606586e-16, -0.4142135729975657]),
        ([-8.253104968768911e-17, 0.41421357299756567, 0.0],[-1.8432293474903632e-17, 0.41421357299756556]),
    ],
    false=>[
        ([0.18313912161055637, 2.226160334420496e-05, 0.0],[0.18313912161055643, 2.2261603344117582e-05]),
        ([0.189351618074891, -0.2180360975566466, 0.0],[0.18935161807489106, -0.21803609755664669]),
        ([0.18949738649103004, 0.21874093267138872, 0.0],[0.1894973864910301, 0.21874093267138864]),
    ],
)
@testset "Circular boundary reparametrization retains stored and absent U" begin
    for stored in (true,false)
        A.initialize();A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        A.option("Mesh.Smoothing",1);A.option("Mesh.RecombineOptimizeTopology",5)
        try
            face=setup_circle(stored)
            boundary=[A.mesh.get_node(i) for i in 1:8]
            held=A.CURRENT[];digest=fingerprint(held)
            A.mesh.recombine()
            types,tags,cells=A.mesh.get_elements(2,face)
            expected=copy(CELLS[:,2:end]);expected[expected.==11].=9
            @test types==[3]
            @test tags==[UInt64.(2:7)]
            @test cells==[UInt64.(vec(expected))]
            @test [A.mesh.get_node(i) for i in 1:8]==boundary
            @test sort(first(A.mesh.get_nodes()))==UInt64[i for i in 1:12 if i!=11]
            @test fingerprint(held)==digest
            for (node,(xyz,uv)) in zip((9,10,12),CIRCLE_CONTROL[stored])
                actual,parameter,dimension,entity=A.mesh.get_node(node)
                @test actual≈xyz atol=2e-12
                @test parameter≈uv atol=2e-12
                @test (dimension,entity)==(2,face)
                @test actual≈A.model.get_value(2,face,parameter) atol=2e-12
            end
        finally
            A.finalize()
        end
    end
end
end
