using Test,Tessella,Random
using Tessella.Model:mesh_model_surface,mesh_model_volume,set_point_coordinates!
using Tessella.MeshTypes:nnodes
if !isdefined(@__MODULE__,:QuadTriNoNewQuadStripCertificates)
    include("quadtri_nonew_quad_strip_certificates.jl")
end

module QuadTriCentroidOverflowTests

using Test,Tessella,Random
using Tessella.Model:mesh_model_surface,mesh_model_volume,set_point_coordinates!
using Tessella.MeshTypes:nnodes
using ..QuadTriNoNewQuadStripCertificates
const Certificate=QuadTriNoNewQuadStripCertificates
const Model=Tessella.Model
const Q=Rational{BigInt}

# Keep the established binary addition order as an independent compatibility
# reference; exact means below establish the overflow geometry separately.
function original_mean(v)
    sums=zeros(3);count=0;half=length(v)÷2
    for k in eachindex(v)
        k>half && v[k]==v[k-half] && continue
        for axis in 1:3;sums[axis]+=v[k][axis];end
        count+=1
    end
    return Tuple(sums./count)
end
function exact_mean(v)
    kept=[k for k in eachindex(v) if k<=length(v)÷2 || v[k]!=v[k-length(v)÷2]]
    return ntuple(axis->sum(Q(v[k][axis]) for k in kept)/length(kept),3)
end
bits(v)=reinterpret.(UInt64,v)
function regular_batch(v,count)
    total=0.0
    for _ in 1:count
        p=Model._extrude_quadtri_centroid(v);total+=p[1]+p[2]+p[3]
    end
    return total
end
function boxes(x)
    x===Core.Box && return 1
    x isa GlobalRef && x.mod===Core && x.name===:Box && return 1
    x isa Core.CodeInfo && return sum(boxes,x.code;init=0)
    x isa Expr && return sum(boxes,x.args;init=0)
    return 0
end
function boxed_control(n)
    value=0;f=()->value
    for _ in 1:n;value+=1;end
    return f()
end

function macro_source(count,normal_sign,recombined)
    # Construct topology at ordinary scale. Moving existing CAD points and
    # their retained translation describes a finite exact axis product,
    # without relying on huge-coordinate CAD construction/coherence.
    return """
    Geometry.AutoCoherence=0;
    Point(1)={0,0,0};Point(2)={0,1,0};Point(3)={0,1,1};Point(4)={0,0,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
    Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
    Transfinite Curve{1,3}=$count;Transfinite Curve{2,4}=2;
    Transfinite Surface{1}={1,2,3,4};Recombine Surface{1};
    Extrude{$normal_sign,0,0}{Surface{1};Layers{1};Recombine;
        QuadTriNoNewVerts$(recombined ? " RecombLaterals" : "");}
    """
end

function public_product(count,base_sign,normal_sign,recombined)
    base=base_sign*1e308
    cap=normal_sign>0 ? nextfloat(base,1000) : prevfloat(base,1000)
    delta=cap-base
    mktempdir() do directory
        path=joinpath(directory,"finite-centroid.geo")
        write(path,macro_source(count,normal_sign,recombined))
        model=Tessella.execute_geo(path;mesh_dim=0).model
        for tag in keys(model.points)
            p=model.points[tag]
            set_point_coordinates!(model,tag,p[1]==0.0 ? base : cap,p[2],p[3])
        end
        for key in keys(model.meshing.extrude_specs)
            spec=model.meshing.extrude_specs[key]
            @test spec.type===:translate
            model.meshing.extrude_specs[key]=merge(spec,(T=(delta,0.,0.),))
        end
        source=mesh_model_surface(model,1)
        @test nnodes(source)==2count
        @test all(isfinite,source.coords) && all(source.coords[1,:].==base)
        source_cells=[Tuple(Int.(cell)) for b in Certificate.blocks(source) for cell in eachcol(b.nodes)]
        @test length(source_cells)==count-1 && all(length(cell)==4 for cell in source_cells)
        volume=mesh_model_volume(model,1)
        primary=2nnodes(source)
        @test nnodes(volume)==primary+(recombined ? count-1 : 0)
        @test all(isfinite,volume.coords)
        column_points=Set((x,source.coords[2,k],source.coords[3,k]) for x in (base,cap),k in 1:nnodes(source))
        actual_points=Set(Certificate.point(volume,k) for k in 1:nnodes(volume))
        @test issubset(column_points,actual_points)
        centers=setdiff(actual_points,column_points)
        @test length(centers)==(recombined ? count-1 : 0)
        expected_x=Float64((Q(base)+Q(cap))/2)
        @test all(p->min(base,cap)<p[1]<max(base,cap) && p[1]==expected_x &&
                    0<p[2]<1 && 0<p[3]<1,centers)
        if recombined
            @test sum(size(b.nodes,2) for b in Certificate.blocks(volume) if b.msh==4)==2(count-1)
            @test sum(size(b.nodes,2) for b in Certificate.blocks(volume) if b.msh==7)==5(count-1)
        end
        faces=Dict{Tuple,Vector{Tuple}}();used=Set{Int}();total=Q(0)
        for b in Certificate.blocks(volume),cell in eachcol(b.nodes)
            nodes=Tuple(Int.(cell));p=Tuple(Certificate.exact(Certificate.point(volume,k)) for k in nodes)
            @test Certificate.map_positive(p,Int(b.msh))
            total+=Certificate.cell_volume(p,Int(b.msh))
            union!(used,nodes)
            for pattern in Certificate.FACES[Int(b.msh)]
                face=Tuple(nodes[k] for k in pattern)
                push!(get!(faces,Certificate.key(face),Tuple[]),face)
            end
        end
        @test used==Set(1:nnodes(volume))
        @test total==abs(Q(delta))
        for incidence in values(faces)
            @test length(incidence) in (1,2)
            length(incidence)==2 && (@test Certificate.cycle(incidence[1])==Certificate.cycle(reverse(incidence[2])))
        end
        # Repeating the actual product retains both the selected geometry
        # and the stable native source parameter rows.
        parameters=deepcopy(model.curve_params)
        repeated=mesh_model_volume(model,1)
        @test repeated.coords==volume.coords
        @test [(b.msh,b.nodes) for b in Certificate.blocks(repeated)]==[(b.msh,b.nodes) for b in Certificate.blocks(volume)]
        @test model.curve_params==parameters
    end
end

@testset "QuadTri finite centroid overflow" begin
    @testset "Healthy binary order and fixed columns" begin
        rng=MersenneTwister(721)
        for width in (6,8),index in 1:250
            v=ntuple(_->Tuple(randn(rng,3)),width)
            index%3==0 && (v=ntuple(k->k>width÷2 ? v[k-width÷2] : v[k],width))
            @test bits(Model._extrude_quadtri_centroid(v))==bits(original_mean(v))
        end
        # A healthy coordinate keeps its original rounding bits even when
        # another coordinate's sum needs the exact fallback.
        v=ntuple(k->(1e308, k==1 ? 1.0 : k==2 ? 1e-16 : 0.0, Float64(k)),8)
        result=Model._extrude_quadtri_centroid(v)
        @test bits(result)[2:3]==bits(original_mean(v))[2:3]
    end
    @testset "Exact finite means, cancellation and input rejection" begin
        for width in (6,8),sign in (-1.,1.)
            v=ntuple(k->(sign*(k<=width÷2 ? 1e308 : nextfloat(1e308,1000)),Float64(k),Float64(k%2)),width)
            expected=Float64.(exact_mean(v));actual=Model._extrude_quadtri_centroid(v)
            @test all(isfinite,actual) && actual==expected
            @test min(v[1][1],v[end][1])<actual[1]<max(v[1][1],v[end][1])
        end
        for sign in (-1.,1.)
            v=ntuple(k->(sign*(k<=4 ? 1e308 : -1e308),0.,Float64(k)),8)
            @test !isfinite(original_mean(v)[1])
            @test Model._extrude_quadtri_centroid(v)[1]==0.0
        end
        ratio=(Q(1e308)+Q(nextfloat(0.0)))/8
        @test ndigits(numerator(ratio);base=2)>1024 && ndigits(denominator(ratio);base=2)>1024
        @test Float64(ratio)==1e308/8
        for nonfinite in (Inf,-Inf,NaN)
            v=ntuple(k->(k==1 ? nonfinite : 1.,0.,0.),8)
            @test isequal(Model._extrude_quadtri_centroid(v),original_mean(v))
        end
    end
    @testset "Real single and three Quad products" begin
        for count in (2,4),base_sign in (-1.,1.),normal_sign in (-1,1),recombined in (false,true)
            public_product(count,base_sign,normal_sign,recombined)
        end
    end
    @testset "Ordinary resource behavior" begin
        v=ntuple(k->(Float64(k),Float64(k%3),Float64(k%2)),8)
        for count in (1000,2000,4000)
            regular_batch(v,count);regular_batch(v,count)
            @test minimum(@allocated(regular_batch(v,count)) for _ in 1:3)==0
        end
        @test sum(boxes(Base.uncompressed_ast(m)) for m in methods(boxed_control))>0
        for helper in (Model._extrude_quadtri_centroid,Model._extrude_quadtri_centroid_exact),method in methods(helper)
            @test boxes(Base.uncompressed_ast(method))==0
        end
    end
end

end
