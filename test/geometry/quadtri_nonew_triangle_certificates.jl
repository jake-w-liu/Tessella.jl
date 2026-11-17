if !isdefined(@__MODULE__,:QuadTriNoNewCertificates)
    include("quadtri_nonew_certificates.jl")
end

module QuadTriNoNewTriangleCertificates

using Tessella
using Tessella.Elements: ElementBlock, MixedMesh, msh_dimension
using Tessella.MeshTypes: nnodes
using ..QuadTriNoNewCertificates
const Quad=QuadTriNoNewCertificates
const CORNERS=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.))
require(ok,message)=ok || throw(ArgumentError("Triangle NoNew certificate: "*message))
crc(mesh)=mixed_crc(mesh isa MixedMesh ? mesh : Tessella.Elements.simplex_to_mixed(mesh))

function root_source(winding=1)
    loop=winding==1 ? "1,2,3" : "-3,-2,-1"
    return """
    Mesh.TransfiniteTri=1;
    Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={0,1,0,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,1};
    Curve Loop(1)={$loop};Plane Surface(1)={1};
    Transfinite Curve{:}=2;Transfinite Surface{1};
    """
end

function fixture(name,layers,levels,laterals;winding=1,motion=:translation,
                 delta=(0.,0.,1.),angle=pi/6)
    literal=angle<0 ? "-Pi/6" : "Pi/6"
    transform=motion===:translation ? "Extrude{$(join(delta,','))}" :
        motion===:rotation ? "Extrude{{0,1,0},{-1,0,0},$literal}" :
                            "Extrude{{0,0,2},{0,0,1},{0,0,0},$literal}"
    source=root_source(winding)*"sweep[]=$transform{Surface{1};$layers;Recombine;"*
        "QuadTriNoNewVerts$(laterals ? " RecombLaterals" : "");};\n"
    return (;name,source,levels=Float64.(levels),laterals,winding,motion,delta,angle)
end

function fixtures()
    return [
        fixture("triangle_one_Z1_W1_Rfalse","Layers{1}",[0.,1.],false),
        fixture("triangle_one_Z-1_W-1_Rtrue","Layers{1}",[0.,1.],true;
            winding=-1,delta=(0.,0.,-1.)),
        fixture("triangle_graded_Z1_W-1_Rfalse","Layers{{2,1},{0.2,1.0}}",
            [0.,.1,.2,1.],false;winding=-1),
        fixture("triangle_graded_Z-1_W1_Rtrue","Layers{{2,1},{0.2,1.0}}",
            [0.,.1,.2,1.],true;delta=(0.,0.,-1.)),
        fixture("triangle_uniform3_W1_Rtrue","Layers{3}",[0.,1/3,2/3,1.],true),
        fixture("triangle_tilted_W1_Rfalse","Layers{3}",[0.,1/3,2/3,1.],false;
            delta=(.25,-.125,2.)),
        fixture("triangle_rotation_W1_Rfalse","Layers{3}",[0.,1/3,2/3,1.],false;
            motion=:rotation),
        fixture("triangle_rotation_W-1_Rtrue","Layers{3}",[0.,1/3,2/3,1.],true;
            winding=-1,motion=:rotation),
        fixture("triangle_rotation_negative_Rtrue","Layers{3}",[0.,1/3,2/3,1.],true;
            motion=:rotation,angle=-pi/6),
        fixture("triangle_twist_W1_Rtrue","Layers{3}",[0.,1/3,2/3,1.],true;
            motion=:twist),
        fixture("triangle_twist_W-1_Rfalse","Layers{3}",[0.,1/3,2/3,1.],false;
            winding=-1,motion=:twist),
        fixture("triangle_twist_negative_Rfalse","Layers{3}",[0.,1/3,2/3,1.],false;
            motion=:twist,angle=-pi/6),
    ]
end

function transformed(f,p,u)
    f.motion===:translation && return ntuple(k->p[k]+u*f.delta[k],3)
    c,s=cos(u*f.angle),sin(u*f.angle)
    f.motion===:rotation && return (-1+c*(p[1]+1)+s*p[3],p[2],-s*(p[1]+1)+c*p[3])
    return (c*p[1]-s*p[2],s*p[1]+c*p[2],p[3]+2u)
end

determinant(a,b,c)=a[1]*(b[2]*c[3]-b[3]*c[2])-
    a[2]*(b[1]*c[3]-b[3]*c[1])+a[3]*(b[1]*c[2]-b[2]*c[1])
exact_point(p)=Tuple(Rational{BigInt}(x) for x in p)

# Direct differentiation of the six reference shape functions lambda_i*(1-t)
# and lambda_i*t. No production interval bounds or template table is used.
function prism_jacobian(vertices,u,v,t)
    points=map(exact_point,vertices)
    u,v,t=Rational{BigInt}.((u,v,t))
    barycentric=(1-u-v,u,v)
    gradients=((-1,-1),(1,0),(0,1))
    derivative=ntuple(3) do axis
        ntuple(3) do coordinate
            sum(points[i][coordinate]*(axis==3 ?
                (i<=3 ? -barycentric[i] : barycentric[i-3]) :
                gradients[mod1(i,3)][axis]*(i<=3 ? 1-t : t)) for i in 1:6)
        end
    end
    return determinant(derivative...)
end

function positive_prism(vertices)
    # det(J) is affine over the reference triangle and quadratic in t.
    # Its minimum therefore occurs on one of these three triangle vertices;
    # an exact stationary-point calculation covers the entire axial interval.
    for (u,v) in ((0,0),(1,0),(0,1))
        low=prism_jacobian(vertices,u,v,0)
        middle=prism_jacobian(vertices,u,v,1//2)
        high=prism_jacobian(vertices,u,v,1)
        low>0 && high>0 || return false
        a=2(low+high-2middle);b=high-low-a
        if a>0
            stationary=-b/(2a)
            if 0<stationary<1
                a*stationary^2+b*stationary+low>0 || return false
            end
        end
    end
    return true
end

function full_jacobians(mesh)
    count=0
    for block in Quad.blocks(mesh)
        block isa ElementBlock && msh_dimension(block.msh)==3 || continue
        require(block.msh in (4,6),"unexpected isolated triangle cell family")
        for cell in eachcol(block.nodes)
            count+=1
            vertices=Tuple(Quad.point(mesh,node) for node in cell)
            if block.msh==4
                p=map(exact_point,vertices)
                require(determinant(p[2].-p[1],p[3].-p[1],p[4].-p[1])>0,
                    "nonpositive constant Tet4 Jacobian")
            else
                require(positive_prism(vertices),"nonpositive full Pri6 Jacobian")
            end
        end
    end
    return count
end

function certify(mesh,f)
    certificate=Quad.certify_complex(mesh)
    require(full_jacobians(mesh)==certificate.ncell,"full-map certificate count")
    expected=[transformed(f,p,u) for u in f.levels,p in CORNERS]
    grid_ids=Matrix{Int}(undef,size(expected))
    for i in eachindex(expected)
        matches=[j for j in 1:nnodes(mesh) if maximum(abs.(Quad.point(mesh,j).-expected[i]))<=2e-11]
        require(length(matches)==1,"swept triangle corner missing or duplicated")
        grid_ids[i]=only(matches)
    end
    require(length(unique(grid_ids))==length(grid_ids),"collapsed columns")
    require(Set(vec(grid_ids))==Set(1:nnodes(mesh)),"added or orphan volume vertices")
    remaining=Set(keys(certificate.boundary))
    for row in (1,length(f.levels))
        face=Quad.key(grid_ids[row,:])
        require(face in remaining,"triangular cap changed")
        delete!(remaining,face)
    end
    nlayer=length(f.levels)-1
    for layer in 1:nlayer,side in 1:3
        next=mod1(side+1,3)
        corners=Set((grid_ids[layer,side],grid_ids[layer,next],
                     grid_ids[layer+1,side],grid_ids[layer+1,next]))
        faces=[face for face in remaining if issubset(Set(face),corners)]
        if f.laterals
            require(length(faces)==1 && length(only(faces))==4,"recombined lateral changed")
        else
            require(length(faces)==2 && all(face->length(face)==3,faces),"free lateral missing split")
            common=Quad.key(intersect(collect(faces[1]),collect(faces[2])))
            require(common in (Quad.key((grid_ids[layer,side],grid_ids[layer+1,next])),
                               Quad.key((grid_ids[layer,next],grid_ids[layer+1,side]))),
                "free lateral has an invalid diagonal")
        end
        setdiff!(remaining,faces)
    end
    require(isempty(remaining),"unexpected exterior face")
    family_counts=Quad.counts(mesh)
    require(family_counts==(f.laterals ? Dict(6=>nlayer) : Dict(4=>3nlayer)),
        "unexpected triangle sweep families")
    level_of=Dict(grid_ids[row,column]=>row for row in axes(grid_ids,1),column in 1:3)
    per_layer=zeros(Int,nlayer)
    for block in Quad.blocks(mesh),cell in eachcol(block.nodes)
        lower=minimum(level_of[node] for node in cell)
        require(maximum(level_of[node] for node in cell)==lower+1,"cell spans unrelated intervals")
        per_layer[lower]+=1
    end
    require(all(==(f.laterals ? 1 : 3),per_layer),"interval cell counts differ")
    if f.motion===:translation
        require(isapprox(certificate.total,abs(f.delta[3])/2;atol=2e-12,rtol=2e-12),
            "analytic translation volume differs")
    end
    return merge(certificate,(;grid_ids,family_counts))
end

function certify_cap_winding(parts,volume,f,certificate)
    source=only(part for (dim,tag,part) in parts if dim==2 && tag==1)
    source_block=source isa MixedMesh ? only(b for b in source.blocks if b.msh==2) :
        ElementBlock(2,source.tris)
    source_points=Tuple(Quad.point(source,node) for node in source_block.nodes[:,1])
    lookup=Dict(Quad.point(volume,node)=>node for node in 1:nnodes(volume))
    lower=Tuple(lookup[p] for p in source_points)
    upper=Tuple(only(node for node in 1:nnodes(volume) if
        maximum(abs.(Quad.point(volume,node).-transformed(f,p,1.)))<=2e-11) for p in source_points)
    exact=map(exact_point,source_points)
    step=exact_point(Quad.point(volume,upper[1])).-exact[1]
    sigma=sign(determinant(exact[2].-exact[1],exact[3].-exact[1],step))
    require(sigma!=0,"source/top winding has zero sweep orientation")
    lower_same=Quad.canonical_cycle(lower)==Quad.canonical_cycle(certificate.boundary[Quad.key(lower)])
    upper_same=Quad.canonical_cycle(upper)==Quad.canonical_cycle(certificate.boundary[Quad.key(upper)])
    require(lower_same==(sigma<0) && upper_same==(sigma>0),"cap outward direction is inconsistent")
    # This compares the actual material-source order, independent of whether
    # a surface mesher normalizes an input loop's winding convention.
    source_winding=sign(determinant(exact[2].-exact[1],exact[3].-exact[1],(0,0,1)))
    require(source_winding==f.winding,"source winding differs from the explicit XY loop")
    return (;sigma,lower_same,upper_same,source_winding)
end

function paired_source()
    return root_source()*"""
    Point(101)={3,0,0};Point(102)={4,0,0};Point(103)={3,1,0};
    Line(101)={101,102};Line(102)={102,103};Line(103)={103,101};
    Curve Loop(101)={101,102,103};Plane Surface(101)={101};
    Transfinite Curve{101,102,103}=2;Transfinite Surface{101};
    left[]=Extrude{0,0,1}{Surface{1};Layers{3};Recombine;QuadTriNoNewVerts;};
    right[]=Extrude{0,0,1}{Surface{101};Layers{3};Recombine;QuadTriNoNewVerts RecombLaterals;};
    """
end


end
