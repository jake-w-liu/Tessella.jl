using Test
using SHA
using Random
using Tessella

const _REFERENCE_GEOMETRY=Tessella.MeshReferenceGeometry

function _reference_geometry_fixture()
    coordinates=Float64[
        1 4  -3 -1  -2 1 0.5   2 3 -0.5   0 2 0 0   1 0 0 0;
       -2 2   1  4   0 0.2 3  -1 2  1     0 0 3 0   0 0 1 0;
        0.5 -1 2 3.5 0 1 -0.4  0.5 1 4    0 0 0 4   0 1 0 4]
    return Mesh(
        coordinates;
        segs=Int32[1 3;2 4],
        tris=Int32[5 8;6 9;7 10],
        tets=Int32[11 15;12 16;13 17;14 18])
end

@testset "fixed-seed affine Jacobian oracle" begin
    rng=Xoshiro(0x4ac9e13f72b8605d)
    for _ in 1:32
        a=ntuple(_->randn(rng),3)
        b=ntuple(i->a[i]+randn(rng),3)
        while hypot((b.-a)...)<0.25
            b=ntuple(i->a[i]+randn(rng),3)
        end
        u=3rand(rng)-1.5
        mesh=Mesh(hcat(collect(a),collect(b));
                  segs=reshape(Int32[1,2],2,1))
        jacobian,determinant,physical=
            _REFERENCE_GEOMETRY.mesh_jacobian(mesh,1,[u,0,0])
        derivative=ntuple(i->(b[i]-a[i])/2,3)
        @test isapprox(
            jacobian[1:3],collect(derivative);
            atol=2e-14,rtol=2e-14)
        @test isapprox(
            determinant[1],hypot(derivative...);
            atol=2e-14,rtol=2e-14)
        @test isapprox(
            physical,
            [((1-u)/2)*a[i]+((1+u)/2)*b[i] for i in 1:3];
            atol=2e-14,rtol=2e-14)
        @test isapprox(
            sum(jacobian[1:3].*jacobian[4:6]),0.0;
            atol=2e-14,rtol=0.0)
        @test isapprox(hypot(jacobian[4:6]...),1.0;
                       atol=2e-14,rtol=2e-14)
        @test isapprox(hypot(jacobian[7:9]...),1.0;
                       atol=2e-14,rtol=2e-14)
    end

    for _ in 1:32
        a=ntuple(_->randn(rng),3)
        b=ntuple(i->a[i]+randn(rng),3)
        c=ntuple(i->a[i]+randn(rng),3)
        edge1=b.-a;edge2=c.-a
        cross=(edge1[2]*edge2[3]-edge1[3]*edge2[2],
               edge1[3]*edge2[1]-edge1[1]*edge2[3],
               edge1[1]*edge2[2]-edge1[2]*edge2[1])
        while hypot(cross...)<0.25
            c=ntuple(i->a[i]+randn(rng),3)
            edge2=c.-a
            cross=(edge1[2]*edge2[3]-edge1[3]*edge2[2],
                   edge1[3]*edge2[1]-edge1[1]*edge2[3],
                   edge1[1]*edge2[2]-edge1[2]*edge2[1])
        end
        u=randn(rng);v=randn(rng)
        mesh=Mesh(hcat(collect(a),collect(b),collect(c));
                  tris=reshape(Int32[1,2,3],3,1))
        jacobian,determinant,physical=
            _REFERENCE_GEOMETRY.mesh_jacobian(mesh,1,[u,v,0])
        magnitude=hypot(cross...)
        @test isapprox(
            jacobian[1:6],collect((edge1...,edge2...));
            atol=2e-14,rtol=2e-14)
        @test isapprox(
            jacobian[7:9],collect(cross)./magnitude;
            atol=2e-14,rtol=2e-14)
        @test isapprox(determinant[1],magnitude;
                       atol=2e-14,rtol=2e-14)
        @test isapprox(
            physical,
            [a[i]+u*edge1[i]+v*edge2[i] for i in 1:3];
            atol=4e-14,rtol=4e-14)
    end

    for _ in 1:32
        a=ntuple(_->randn(rng),3)
        b=ntuple(i->a[i]+randn(rng),3)
        c=ntuple(i->a[i]+randn(rng),3)
        d=ntuple(i->a[i]+randn(rng),3)
        edge1=b.-a;edge2=c.-a;edge3=d.-a
        determinant_expected=
            edge1[1]*(edge2[2]*edge3[3]-edge2[3]*edge3[2])-
            edge1[2]*(edge2[1]*edge3[3]-edge2[3]*edge3[1])+
            edge1[3]*(edge2[1]*edge3[2]-edge2[2]*edge3[1])
        while abs(determinant_expected)<0.25
            d=ntuple(i->a[i]+randn(rng),3)
            edge3=d.-a
            determinant_expected=
                edge1[1]*(edge2[2]*edge3[3]-edge2[3]*edge3[2])-
                edge1[2]*(edge2[1]*edge3[3]-edge2[3]*edge3[1])+
                edge1[3]*(edge2[1]*edge3[2]-edge2[2]*edge3[1])
        end
        u=randn(rng);v=randn(rng);w=randn(rng)
        mesh=Mesh(hcat(collect(a),collect(b),collect(c),collect(d));
                  tets=reshape(Int32[1,2,3,4],4,1))
        jacobian,determinant,physical=
            _REFERENCE_GEOMETRY.mesh_jacobian(mesh,1,[u,v,w])
        @test isapprox(
            jacobian,collect((edge1...,edge2...,edge3...));
            atol=2e-14,rtol=2e-14)
        @test isapprox(
            determinant[1],determinant_expected;
            atol=4e-14,rtol=4e-14)
        @test isapprox(
            physical,
            [a[i]+u*edge1[i]+v*edge2[i]+w*edge3[i] for i in 1:3];
            atol=6e-14,rtol=6e-14)
    end
end

function _reference_geometry_digest(mesh)
    stream=IOBuffer()
    local_coordinates=Float64[-1,0,0, 0.2,0.3,0.1, 1,0,0]
    for element_type in (1,2,4)
        for values in _REFERENCE_GEOMETRY.mesh_jacobians(
            mesh,element_type,local_coordinates)
            write(stream,htol(UInt64(length(values))))
            foreach(value->write(
                stream,htol(reinterpret(UInt64,value))),values)
        end
    end
    for tag in (1,3,6)
        for values in _REFERENCE_GEOMETRY.mesh_jacobian(
            mesh,tag,local_coordinates)
            write(stream,htol(UInt64(length(values))))
            foreach(value->write(
                stream,htol(reinterpret(UInt64,value))),values)
        end
    end
    return bytes2hex(SHA.sha256(take!(stream)))
end

@noinline function _reference_geometry_batch_allocation(count::Int)
    coordinates=Float64[0 1 0 0;0 0 1 0;0 0 0 1]
    tetrahedra=repeat(reshape(Int32[1,2,3,4],4,1),1,count)
    mesh=Mesh(coordinates;tets=tetrahedra)
    local_coordinates=Float64[0.2,0.3,0.1]
    _REFERENCE_GEOMETRY.mesh_jacobians(mesh,4,local_coordinates)
    GC.gc()
    return @allocated _REFERENCE_GEOMETRY.mesh_jacobians(
        mesh,4,local_coordinates)
end

@testset "robust linear-simplex reference geometry" begin
    fixture=_reference_geometry_fixture()
    segment_local=Float64[-1,0,0, 0,9,-7, 1,0,0]
    segment_jacobians,segment_determinants,segment_coordinates=
        _REFERENCE_GEOMETRY.mesh_jacobians(fixture,1,segment_local)
    expected_segment_jacobian=[
        1.5,2.0,-0.75,0.8,-0.6,0.0,
        -0.17240873133980728,-0.22987830845307636,
        -0.9578262852211514]
    @test length(segment_jacobians)==54
    @test isapprox(
        segment_jacobians[1:9],expected_segment_jacobian;
        atol=4eps(Float64),rtol=4eps(Float64))
    @test segment_jacobians[1:9]==segment_jacobians[10:18]
    @test segment_jacobians[1:9]==segment_jacobians[19:27]
    @test segment_determinants[1:3]==fill(hypot(1.5,2.0,-0.75),3)
    @test segment_coordinates[1:9]==
          [1.0,-2.0,0.5, 2.5,0.0,-0.25, 4.0,2.0,-1.0]

    adjacent=Mesh(
        Float64[1.0 nextfloat(1.0);0 0;0 0];
        segs=reshape(Int32[1,2],2,1))
    adjacent_coordinates=_REFERENCE_GEOMETRY.mesh_jacobian(
        adjacent,1,[-1,0,0, 1,0,0])[3]
    @test adjacent_coordinates==
          [1.0,0.0,0.0,nextfloat(1.0),0.0,0.0]

    triangle_local=Float64[0,0,0, 0.2,0.3,99, 1,0,0]
    triangle_jacobians,triangle_determinants,triangle_coordinates=
        _REFERENCE_GEOMETRY.mesh_jacobians(fixture,2,triangle_local)
    expected_triangle_jacobian=[
        3.0,0.2,1.0,2.5,3.0,-0.4,
        -0.31529453733321783,0.3787629182249695,
        0.8701310283546596]
    @test isapprox(
        triangle_jacobians[1:9],expected_triangle_jacobian;
        atol=4eps(Float64),rtol=4eps(Float64))
    @test triangle_jacobians[1:9]==triangle_jacobians[10:18]
    @test triangle_determinants[1:3]==fill(9.768643713433303,3)
    @test isapprox(
        triangle_coordinates[1:9],
        [-2.0,0.0,0.0, -0.65,0.94,0.08, 1.0,0.2,1.0];
        atol=4eps(Float64),rtol=4eps(Float64))

    tetrahedron_local=Float64[0,0,0, 0.2,0.3,0.1]
    tetrahedron_jacobians,tetrahedron_determinants,tetrahedron_coordinates=
        _REFERENCE_GEOMETRY.mesh_jacobians(fixture,4,tetrahedron_local)
    @test tetrahedron_jacobians==[
        2.0,0,0,0,3.0,0,0,0,4.0,
        2.0,0,0,0,3.0,0,0,0,4.0,
        -1.0,0,1.0,-1.0,1.0,0,-1.0,0,4.0,
        -1.0,0,1.0,-1.0,1.0,0,-1.0,0,4.0]
    @test tetrahedron_determinants==[24.0,24.0,-3.0,-3.0]
    @test isapprox(
        tetrahedron_coordinates,
        [0.0,0.0,0.0, 0.4,0.9,0.4,
         1.0,0.0,0.0, 0.4,0.3,0.6];
        atol=4eps(Float64),rtol=4eps(Float64))

    single=_REFERENCE_GEOMETRY.mesh_jacobian(
        fixture,6,Float64[0.2,0.3,0.1])
    @test single[1]==
          [-1.0,0,1.0,-1.0,1.0,0,-1.0,0,4.0]
    @test single[2]==[-3.0]
    @test isapprox(
        single[3],[0.4,0.3,0.6];
        atol=4eps(Float64),rtol=4eps(Float64))

    detached=_REFERENCE_GEOMETRY.mesh_jacobians(
        fixture,4,Float64[0.2,0.3,0.1])
    detached[1][1]=99.0;detached[2][1]=99.0;detached[3][1]=99.0
    repeated=_REFERENCE_GEOMETRY.mesh_jacobian(
        fixture,5,Float64[0.2,0.3,0.1])
    @test repeated[1]==[2.0,0,0,0,3.0,0,0,0,4.0]
    @test repeated[2]==[24.0]
    @test isapprox(
        repeated[3],[0.4,0.9,0.4];
        atol=4eps(Float64),rtol=4eps(Float64))

    @test _REFERENCE_GEOMETRY.mesh_jacobians(
        fixture,3,Float64[0,0,0])==(Float64[],Float64[],Float64[])
    @test _REFERENCE_GEOMETRY.mesh_jacobians(
        fixture,4,Float64[])==(Float64[],Float64[],Float64[])
    @test _REFERENCE_GEOMETRY.mesh_jacobian(
        fixture,1,())==(Float64[],Float64[],Float64[])

    for call in (
        ()->_REFERENCE_GEOMETRY.mesh_jacobians(fixture,true,[0,0,0]),
        ()->_REFERENCE_GEOMETRY.mesh_jacobians(fixture,999,[0,0,0]),
        ()->_REFERENCE_GEOMETRY.mesh_jacobians(fixture,1,0),
        ()->_REFERENCE_GEOMETRY.mesh_jacobians(fixture,1,[0]),
        ()->_REFERENCE_GEOMETRY.mesh_jacobians(fixture,1,[0,0]),
        ()->_REFERENCE_GEOMETRY.mesh_jacobians(fixture,1,[0,0,0,0]),
        ()->_REFERENCE_GEOMETRY.mesh_jacobians(
            fixture,1,Any[true,0,0]),
        ()->_REFERENCE_GEOMETRY.mesh_jacobians(fixture,1,[missing,0,0]),
        ()->_REFERENCE_GEOMETRY.mesh_jacobians(fixture,1,[NaN,0,0]),
        ()->_REFERENCE_GEOMETRY.mesh_jacobians(fixture,1,[Inf,0,0]),
        ()->_REFERENCE_GEOMETRY.mesh_jacobian(fixture,true,[0,0,0]),
        ()->_REFERENCE_GEOMETRY.mesh_jacobian(fixture,0,[0,0,0]),
        ()->_REFERENCE_GEOMETRY.mesh_jacobian(fixture,7,[0,0,0]),
        ()->_REFERENCE_GEOMETRY.mesh_jacobian(
            fixture,big(typemax(Int))+1,[0,0,0]),
    )
        @test_throws ArgumentError call()
    end

    degenerate_segment=Mesh(
        zeros(3,2);segs=reshape(Int32[1,2],2,1))
    degenerate_triangle=Mesh(
        Float64[0 1 2;0 0 0;0 0 0];
        tris=reshape(Int32[1,2,3],3,1))
    degenerate_tetrahedron=Mesh(
        Float64[0 1 0 1;0 0 1 1;0 0 0 0];
        tets=reshape(Int32[1,2,3,4],4,1))
    for (element_type,bad) in (
        (1,degenerate_segment),(2,degenerate_triangle),
        (4,degenerate_tetrahedron))
        @test_throws ArgumentError _REFERENCE_GEOMETRY.mesh_jacobians(
            bad,element_type,[0,0,0])
        @test_throws ArgumentError _REFERENCE_GEOMETRY.mesh_jacobian(
            bad,1,[0,0,0])
        @test _REFERENCE_GEOMETRY.mesh_jacobians(
            bad,element_type,[])==(Float64[],Float64[],Float64[])
    end

    maximum=floatmax(Float64)
    extreme_segment=Mesh(
        Float64[-maximum maximum;0 0;0 0];
        segs=reshape(Int32[1,2],2,1))
    extreme_result=_REFERENCE_GEOMETRY.mesh_jacobian(
        extreme_segment,1,[-1,0,0, 0,0,0, 1,0,0])
    @test extreme_result[1][1]==maximum
    @test extreme_result[2]==fill(maximum,3)
    @test extreme_result[3]==[-maximum,0,0, 0,0,0, maximum,0,0]
    @test_throws ArgumentError _REFERENCE_GEOMETRY.mesh_jacobian(
        extreme_segment,1,[2,0,0])

    tiny=nextfloat(0.0)
    unrepresentable_segment=Mesh(
        Float64[0 tiny;0 0;0 0];
        segs=reshape(Int32[1,2],2,1))
    @test_throws ArgumentError _REFERENCE_GEOMETRY.mesh_jacobian(
        unrepresentable_segment,1,[0,0,0])

    for factor in (2.0^-300,2.0^300)
        scaled=Mesh(
            Float64[0 factor 0 0;0 0 factor 0;0 0 0 factor];
            tris=reshape(Int32[1,2,3],3,1),
            tets=reshape(Int32[1,2,3,4],4,1))
        triangle=_REFERENCE_GEOMETRY.mesh_jacobian(
            scaled,1,[0.2,0.3,0])
        tetrahedron=_REFERENCE_GEOMETRY.mesh_jacobian(
            scaled,2,[0.2,0.3,0.1])
        @test triangle[2]==[factor^2]
        @test tetrahedron[2]==[factor^3]
        @test isapprox(
            triangle[3],[0.2factor,0.3factor,0.0];
            atol=0.0,rtol=8eps(Float64))
        @test isapprox(
            tetrahedron[3],[0.2factor,0.3factor,0.1factor];
            atol=0.0,rtol=8eps(Float64))
    end

    overflowing_triangle=Mesh(
        Float64[0 2.0^600 0;0 0 2.0^600;0 0 0];
        tris=reshape(Int32[1,2,3],3,1))
    @test_throws ArgumentError _REFERENCE_GEOMETRY.mesh_jacobian(
        overflowing_triangle,1,[0,0,0])

    underflowing_triangle=Mesh(
        Float64[0 2.0^-600 0;0 0 2.0^-600;0 0 0];
        tris=reshape(Int32[1,2,3],3,1))
    underflowing_tetrahedron=Mesh(
        Float64[0 2.0^-400 0 0;0 0 2.0^-400 0;0 0 0 2.0^-400];
        tets=reshape(Int32[1,2,3,4],4,1))
    @test_throws ArgumentError _REFERENCE_GEOMETRY.mesh_jacobian(
        underflowing_triangle,1,[0,0,0])
    @test_throws ArgumentError _REFERENCE_GEOMETRY.mesh_jacobian(
        underflowing_tetrahedron,1,[0,0,0])

    underflowing_coordinate=Mesh(
        Float64[0 1 0;0 tiny 0;0 0 1];
        tris=reshape(Int32[1,2,3],3,1))
    @test_throws ArgumentError _REFERENCE_GEOMETRY.mesh_jacobian(
        underflowing_coordinate,1,[0.5,0,0])

    allocation_small=_reference_geometry_batch_allocation(10_000)
    allocation_large=_reference_geometry_batch_allocation(20_000)
    @test allocation_small<=1_150_000
    @test allocation_large<=2allocation_small+8_192

    @test _reference_geometry_digest(fixture)==
          "10188d3068aeabe3ac38e858d24673aea7fe1f9a6d001b5ed4030e601322d58f"
    @test isempty(Docs.undocumented_names(
        Tessella.MeshReferenceGeometry;private=false))
    @test isempty(Test.detect_ambiguities(
        Tessella.MeshReferenceGeometry;recursive=true))
end
