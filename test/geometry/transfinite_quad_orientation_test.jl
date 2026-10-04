module TransfiniteQuadOrientationTests

using Test,Tessella
using Tessella.Elements: ElementBlock,MixedMesh

const CORNERS=((0.,0.,0.),(2.,0.,0.),(2.5,1.5,0.),(.5,1.5,0.))
const ARRANGEMENTS=(:left,:right,:alternate_left,:alternate_right)
const WORDS=Dict(:left=>"Left",:right=>"Right",:alternate_left=>"AlternateLeft",
                 :alternate_right=>"AlternateRight")

function source(;sign=1,axis=:xy,arrangement=:left,pins=(),recombine=false,
                count=2,reverse=false,directions=0)
    points=axis===:xy ? CORNERS : Tuple((p[1],p[3],p[2]) for p in CORNERS)
    point_text=join(("Point($i)={$(join(p,',')),1};" for (i,p) in pairs(points)))
    storage_directions=ntuple(k->iszero(directions&(1<<(k-1))) ? 1 : -1,4)
    curve_text=join("Line($k)={$(storage_directions[k]==1 ? k : mod1(k+1,4)),$(storage_directions[k]==1 ? mod1(k+1,4) : k)};" for k in 1:4)
    loop=sign==1 ? collect(1:4) : [-4,-3,-2,-1]
    signed=[s*storage_directions[abs(s)] for s in loop]
    pin_text=isempty(pins) ? "" : "={$(join(pins,','))}"
    return """
    $point_text
    $curve_text
    Curve Loop(1)={$(join(signed,','))};Plane Surface(1)={1};
    Transfinite Curve{:}=$count;Transfinite Surface{1}$pin_text $(WORDS[arrangement]);
    $(recombine ? "Recombine Surface{1};" : "")
    $(reverse ? "Reverse Surface{1};" : "")
    """
end

function geometry(text)
    mktempdir() do directory
        path=joinpath(directory,"quad.geo");write(path,text)
        return execute_geo(path;mesh_dim=0)
    end
end

blocks(mesh::Mesh)=(ElementBlock(2,mesh.tris),)
blocks(mesh::MixedMesh)=mesh.blocks
cycle(c)=minimum(ntuple(k->c[mod1(i+k-1,length(c))],length(c)) for i in 1:length(c))

function expected(spec)
    pins=get(spec,:pins,())
    sign=spec.sign;directions=get(spec,:directions,0)
    # Independent literal Gmsh4.15.2 GEdgeLoop frames for the four possible
    # first-curve storage directions; explicit pins replace that frame.
    frame=isempty(pins) ? sign==1 ?
        (iszero(directions&1) ? (1,2,3,4) : (2,1,4,3)) :
        (iszero(directions&8) ? (4,1,2,3) : (1,4,3,2)) : pins
    coords=[get(spec,:axis,:xy)===:xy ? CORNERS[i] :
            (CORNERS[i][1],CORNERS[i][3],CORNERS[i][2]) for i in frame]
    count=get(spec,:count,2);steps=count-1
    # Dyadic count2/3 affine grids are exact represented-coordinate fixtures.
    points=[ntuple(3) do d
        coords[1][d]+(coords[2][d]-coords[1][d])*i/steps+
                      (coords[4][d]-coords[1][d])*j/steps
    end for j in 0:steps for i in 0:steps]
    index(i,j)=i+1+j*count
    frame_sign=frame[2]==mod1(frame[1]+1,4) ? 1 : -1
    flip=(frame_sign!=sign) ⊻ get(spec,:reverse,false)
    cells=Tuple[]
    for i in 0:steps-1,j in 0:steps-1
        a=index(i,j);b=index(i+1,j);c=index(i+1,j+1);d=index(i,j+1)
        if get(spec,:recombine,false)
            push!(cells,(a,b,c,d))
        else
            right=spec.arrangement===:right ||
                (spec.arrangement===:alternate_right && isodd(i+j)) ||
                (spec.arrangement===:alternate_left && iseven(i+j))
            append!(cells,right ? [(a,b,c),(c,d,a)] : [(a,b,d),(d,b,c)])
        end
    end
    return points,[cycle(flip ? (c[1],reverse(c[2:end])...) : c) for c in cells]
end

function cases()
    result=NamedTuple[]
    for axis in (:xy,:xz),sign in (-1,1),arrangement in ARRANGEMENTS,
        pins in ((),(1,2,3,4),(1,4,3,2),(2,3,4,1)),recombine in (false,true)
        push!(result,(;axis,sign,arrangement,pins,recombine))
    end
    for sign in (-1,1),arrangement in ARRANGEMENTS,recombine in (false,true)
        push!(result,(;sign,arrangement,recombine,count=3))
        push!(result,(;sign,arrangement,recombine,reverse=true))
    end
    for directions in 0:15,sign in (-1,1),recombine in (false,true)
        push!(result,(;sign,arrangement=:alternate_left,recombine,directions))
    end
    return result
end

@testset "Public four-sided transfinite frames and CAD winding" begin
    # All128 TF2 arrangement/pin/type/plane cases were independently compared
    # with two fresh pinned Gmsh4.15.2 models; negative loops are deterministic.
    fixtures=cases()
    @test length(fixtures)==224
    for spec in fixtures
        @testset "$(spec)" begin
            execution=geometry(source(;spec...))
            mesh=Tessella.Model.mesh_model_surface(execution.model,1)
            coords,cells=expected(spec)
            @test validate(mesh).ok
            @test [Tuple(mesh.coords[:,i]) for i in axes(mesh.coords,2)]==coords
            actual=[cycle(Tuple(c)) for block in blocks(mesh) for c in eachcol(block.nodes)]
            @test actual==cells
            @test Set(block.msh for block in blocks(mesh))==Set([get(spec,:recombine,false) ? 3 : 2])
        end
    end
end

@testset "Unsigned sphere meshing frames preserve the CAD radius" begin
    # TransfiniteSph's radius uses the first CAD generatrix corner. Automatic
    # GEdgeLoop corners can start at a different corner when the loop reverses;
    # raising either affected corner makes an accidental radius change visible.
    for raised_corner in (1,4),sign in (-1,1),recombine in (false,true)
        points=ntuple(4) do k
            xy=((0.,0.),(1.,0.),(1.,1.),(0.,1.))[k]
            (xy...,k==raised_corner ? .125 : 0.)
        end
        point_text=join("Point($k)={$(join(p,',')),1};" for (k,p) in pairs(points))
        loop=sign==1 ? "1,2,3,4" : "-4,-3,-2,-1"
        text="""
        $point_text
        Point(5)={.5,.5,-2,1};
        Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
        Curve Loop(1)={$loop};Surface(1)={1} In Sphere{5};
        Transfinite Curve{1,2,3,4}=3;Transfinite Surface{1};
        $(recombine ? "Recombine Surface{1};" : "")
        """
        execution=geometry(text)
        mesh=Tessella.Model.mesh_model_surface(execution.model,1)
        center=only(filter(p->p[1]==.5 && p[2]==.5,
                           [Tuple(mesh.coords[:,i]) for i in axes(mesh.coords,2)]))
        radius=sqrt(.5^2+.5^2+(points[1][3]+2)^2)
        @test size(mesh.coords,2)==9
        @test collect(center)≈[.5,.5,radius-2] atol=4eps(Float64) rtol=0
        @test all(p->p in Set(Tuple(mesh.coords[:,i]) for i in axes(mesh.coords,2)),points)
        @test validate(mesh).ok
    end
end

end
