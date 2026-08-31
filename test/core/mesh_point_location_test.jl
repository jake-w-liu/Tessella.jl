using Test
using SHA
using Random
using Tessella

const _POINT_LOCATION=Tessella.MeshPointLocation

function _point_location_fixture()
    coordinates=Float64[0 1 0 0;
                        0 0 1 0;
                        0 0 0 1]
    return Mesh(
        coordinates;
        segs=reshape(Int32[1,2],2,1),
        tris=reshape(Int32[1,2,3],3,1),
        tets=reshape(Int32[1,2,3,4],4,1))
end

function _point_location_digest(mesh,locator)
    stream=IOBuffer()
    for (x,y,z,dim,strict) in (
        (0.25,0.0,0.0,-1,true),
        (0.2,0.3,0.0,-1,true),
        (0.1,0.2,0.3,-1,true),
        (-5.0e-6,0.2,0.2,3,false))
        tags=_POINT_LOCATION.locate_elements(
            locator,x,y,z;dim=dim,strict=strict)
        write(stream,htol(UInt64(length(tags))))
        foreach(tag->write(stream,htol(tag)),tags)
    end
    for (tag,x,y,z) in (
        (1,0.25,2.0,3.0),(2,0.2,0.3,5.0),(3,0.1,0.2,0.3),
        (3,-0.25,0.5,1.25))
        for value in _POINT_LOCATION.mesh_local_coordinates(mesh,tag,x,y,z)
            write(stream,htol(reinterpret(UInt64,value)))
        end
    end
    return bytes2hex(SHA.sha256(take!(stream)))
end

function _point_location_segment_fixture(count::Int)
    coordinates=zeros(3,count+1)
    segments=Matrix{Int32}(undef,2,count)
    for segment in 1:count
        coordinates[1,segment+1]=segment
        segments[1,segment]=Int32(segment)
        segments[2,segment]=Int32(segment+1)
    end
    return Mesh(coordinates;segs=segments)
end

@noinline function _point_location_build_allocation(mesh)
    GC.gc()
    return @allocated _POINT_LOCATION.SimplexLocator(mesh)
end

@testset "robust finalized-simplex point location" begin
    fixture=_point_location_fixture()
    @test _POINT_LOCATION.mesh_element_offsets(fixture)==(1,2,3)
    @test _POINT_LOCATION.mesh_element_block(fixture,1)==(0,fixture.segs)
    @test _POINT_LOCATION.mesh_element_block(fixture,2)==(1,fixture.tris)
    @test _POINT_LOCATION.mesh_element_block(fixture,4)==(2,fixture.tets)
    @test _POINT_LOCATION.mesh_element_block(fixture,Int32(4))==
          (2,fixture.tets)
    @test _POINT_LOCATION.mesh_element_block(fixture,3)===nothing
    @test_throws ArgumentError _POINT_LOCATION.mesh_element_block(fixture,true)

    segment=_POINT_LOCATION.mesh_element_record(fixture,1)
    triangle=_POINT_LOCATION.mesh_element_record(fixture,2)
    tetrahedron=_POINT_LOCATION.mesh_element_record(fixture,3)
    @test segment==(element_type=Int32(1),dimension=1,
                    node_tags=UInt64[1,2])
    @test triangle==(element_type=Int32(2),dimension=2,
                     node_tags=UInt64[1,2,3])
    @test tetrahedron==(element_type=Int32(4),dimension=3,
                        node_tags=UInt64[1,2,3,4])
    segment.node_tags[1]=99
    @test _POINT_LOCATION.mesh_element_record(fixture,1).node_tags==
          UInt64[1,2]

    locator=_POINT_LOCATION.SimplexLocator(fixture)
    @test locator.mesh===fixture
    @test _POINT_LOCATION.locate_elements(
        locator,0.25,0.0,0.0;strict=true)==UInt64[3,2,1]
    @test _POINT_LOCATION.locate_elements(
        locator,0.2,0.3,0.0;strict=true)==UInt64[3,2]
    @test _POINT_LOCATION.locate_elements(
        locator,0.1,0.2,0.3;strict=true)==UInt64[3]
    @test _POINT_LOCATION.locate_elements(
        locator,0.25,0.0,0.0;dim=1,strict=true)==UInt64[1]
    @test _POINT_LOCATION.locate_elements(
        locator,0.25,0.0,0.0;dim=2,strict=true)==UInt64[2]
    @test _POINT_LOCATION.locate_elements(
        locator,0.25,0.0,0.0;dim=3,strict=true)==UInt64[3]
    @test isempty(_POINT_LOCATION.locate_elements(
        locator,0.25,0.0,0.0;dim=0,strict=true))
    @test isempty(_POINT_LOCATION.locate_elements(
        locator,2.0,2.0,2.0;strict=true))
    @test isempty(_POINT_LOCATION.locate_elements(
        locator,2.0,2.0,2.0;strict=false))

    @test _POINT_LOCATION.mesh_local_coordinates(
        fixture,1,0.25,2.0,3.0)==(-0.5,0.0,0.0)
    @test _POINT_LOCATION.mesh_local_coordinates(
        fixture,2,0.2,0.3,5.0)==(0.2,0.3,0.0)
    @test _POINT_LOCATION.mesh_local_coordinates(
        fixture,3,0.1,0.2,0.3)==(0.1,0.2,0.3)
    @test _POINT_LOCATION.mesh_local_coordinates(
        fixture,3,-0.25,0.5,1.25)==(-0.25,0.5,1.25)

    oblique_segment=Mesh(
        Float64[1 3;2 1;3 5];segs=reshape(Int32[1,2],2,1))
    @test isapprox(collect(_POINT_LOCATION.mesh_local_coordinates(
        oblique_segment,1,2.6,1.7,2.6)),[-0.4,0.0,0.0];
        atol=16eps(Float64),rtol=16eps(Float64))

    oblique_triangle=Mesh(
        Float64[1 2 1;2 3 3;3 3 4];
        tris=reshape(Int32[1,2,3],3,1))
    triangle_point=(6.2,-2.5,8.3)
    @test isapprox(collect(_POINT_LOCATION.mesh_local_coordinates(
        oblique_triangle,1,triangle_point...)),[0.2,0.3,0.0];
        atol=32eps(Float64),rtol=32eps(Float64))

    oblique_tetrahedron=Mesh(
        Float64[1 2 1 2;2 3 3 2;3 3 4 4];
        tets=reshape(Int32[1,2,3,4],4,1))
    tetrahedron_point=(1.3,2.5,3.4)
    @test isapprox(collect(_POINT_LOCATION.mesh_local_coordinates(
        oblique_tetrahedron,1,tetrahedron_point...)),[0.2,0.3,0.1];
        atol=32eps(Float64),rtol=32eps(Float64))
    @test _POINT_LOCATION.locate_elements(
        _POINT_LOCATION.SimplexLocator(oblique_tetrahedron),
        tetrahedron_point...;strict=true)==UInt64[1]

    @test isempty(_POINT_LOCATION.locate_elements(
        locator,-5.0e-6,0.2,0.2;dim=3,strict=true))
    @test _POINT_LOCATION.locate_elements(
        locator,-5.0e-6,0.2,0.2;dim=3,strict=false)==UInt64[3]
    @test isempty(_POINT_LOCATION.locate_elements(
        locator,-1.1,0.1,0.1;dim=3,strict=false))

    relaxed_levels=Mesh(
        Float64[0 1 0 0 4.5e-5 1.000045 4.5e-5 4.5e-5;
                0 0 1 0 0 0 1 0;
                0 0 0 1 0 0 0 1];
        tets=Int32[1 5;2 6;3 7;4 8])
    relaxed_locator=_POINT_LOCATION.SimplexLocator(relaxed_levels)
    @test isempty(_POINT_LOCATION.locate_elements(
        relaxed_locator,-5.0e-6,0.2,0.2;strict=true))
    @test _POINT_LOCATION.locate_elements(
        relaxed_locator,-5.0e-6,0.2,0.2;strict=false)==UInt64[1]

    shared=Mesh(
        Float64[0 1 0 1;0 0 1 1;0 0 0 0];
        tris=Int32[1 2;2 4;3 3])
    @test _POINT_LOCATION.locate_elements(
        _POINT_LOCATION.SimplexLocator(shared),0.5,0.5,0.0;
        strict=true)==UInt64[1,2]

    epsilon=eps(Float64)
    skinny_triangle=Mesh(
        Float64[0 1 1;0 0 epsilon;0 0 0];
        tris=reshape(Int32[1,2,3],3,1))
    @test _POINT_LOCATION.mesh_local_coordinates(
        skinny_triangle,1,0.75,epsilon/4,0.0)==(0.5,0.25,0.0)
    @test _POINT_LOCATION.locate_elements(
        _POINT_LOCATION.SimplexLocator(skinny_triangle),
        0.75,epsilon/4,0.0;strict=true)==UInt64[1]

    tolerance_sensitive_triangle=Mesh(
        Float64[0 1 1;
                0 0 1.0e-5;
                0 0 3.6977196714458935e-6];
        tris=reshape(Int32[1,2,3],3,1))
    tolerance_sensitive_point=(
        0.8720278628389854,3.934735911487432e-6,
        2.469703846133699e-6)
    @test _POINT_LOCATION.locate_elements(
        _POINT_LOCATION.SimplexLocator(tolerance_sensitive_triangle),
        tolerance_sensitive_point...;strict=true)==UInt64[1]
    @test isapprox(collect(_POINT_LOCATION.mesh_local_coordinates(
        tolerance_sensitive_triangle,1,tolerance_sensitive_point...)),
        [0.44554509350225024,0.4264827693367351,0.0];
        atol=4eps(Float64),rtol=4eps(Float64))

    skinny_tetrahedron=Mesh(
        Float64[0 1 0 0;0 0 1 0;0 0 0 epsilon];
        tets=reshape(Int32[1,2,3,4],4,1))
    @test _POINT_LOCATION.mesh_local_coordinates(
        skinny_tetrahedron,1,0.25,0.25,epsilon/4)==(0.25,0.25,0.25)
    @test _POINT_LOCATION.locate_elements(
        _POINT_LOCATION.SimplexLocator(skinny_tetrahedron),
        0.25,0.25,epsilon/4;strict=true)==UInt64[1]

    tetrahedron_delta=2.0^-20
    cancellation_tetrahedron=Mesh(
        Float64[0 1 1 1;
                0 1 1+tetrahedron_delta 1;
                0 1 1 1+tetrahedron_delta];
        tets=reshape(Int32[1,2,3,4],4,1))
    cancellation_point=(
        0.75,0.75+tetrahedron_delta/4,0.75+tetrahedron_delta/4)
    @test _POINT_LOCATION.mesh_local_coordinates(
        cancellation_tetrahedron,1,cancellation_point...)==(0.25,0.25,0.25)
    @test _POINT_LOCATION.locate_elements(
        _POINT_LOCATION.SimplexLocator(cancellation_tetrahedron),
        cancellation_point...;strict=true)==UInt64[1]

    maximum=floatmax(Float64)
    extreme_segment=Mesh(
        Float64[-maximum maximum;0 0;0 0];
        segs=reshape(Int32[1,2],2,1))
    @test _POINT_LOCATION.mesh_local_coordinates(
        extreme_segment,1,0.0,0.0,0.0)==(0.0,0.0,0.0)
    @test _POINT_LOCATION.locate_elements(
        _POINT_LOCATION.SimplexLocator(extreme_segment),
        0.0,0.0,0.0;strict=true)==UInt64[1]

    tiny_segment=Mesh(
        Float64[0 nextfloat(0.0);0 0;0 0];
        segs=reshape(Int32[1,2],2,1))
    @test_throws ArgumentError _POINT_LOCATION.mesh_local_coordinates(
        tiny_segment,1,1.0,0.0,0.0)
    @test _POINT_LOCATION.mesh_local_coordinates(
        tiny_segment,1,0.0,floatmax(Float64),0.0)==(-1.0,0.0,0.0)
    @test isempty(_POINT_LOCATION.locate_elements(
        _POINT_LOCATION.SimplexLocator(tiny_segment),1.0,0.0,0.0;
        strict=false))

    relaxed_with_tiny_segment=Mesh(
        Float64[0 nextfloat(0.0) 10 11;
                0 0 0 0;
                0 0 0 0];
        segs=Int32[1 3;2 4])
    @test _POINT_LOCATION.locate_elements(
        _POINT_LOCATION.SimplexLocator(relaxed_with_tiny_segment),
        10.0-5.0e-6,0.0,0.0;strict=false)==UInt64[2]

    empty_mesh=Mesh(zeros(3,0))
    empty_locator=_POINT_LOCATION.SimplexLocator(empty_mesh)
    @test isempty(empty_locator.order)
    @test isempty(_POINT_LOCATION.locate_elements(
        empty_locator,0.0,0.0,0.0))

    for call in (
        ()->_POINT_LOCATION.mesh_element_record(fixture,0),
        ()->_POINT_LOCATION.mesh_element_record(fixture,4),
        ()->_POINT_LOCATION.mesh_element_record(fixture,true),
        ()->_POINT_LOCATION.mesh_element_record(
            fixture,big(typemax(Int))+1),
        ()->_POINT_LOCATION.mesh_local_coordinates(fixture,1,NaN,0,0),
        ()->_POINT_LOCATION.mesh_local_coordinates(fixture,1,0,Inf,0),
        ()->_POINT_LOCATION.mesh_local_coordinates(fixture,1,0,0,false),
        ()->_POINT_LOCATION.locate_elements(locator,missing,0,0),
        ()->_POINT_LOCATION.locate_elements(locator,0,0,0;dim=false),
        ()->_POINT_LOCATION.locate_elements(locator,0,0,0;dim=4),
        ()->_POINT_LOCATION.locate_elements(locator,0,0,0;strict=0),
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
    for bad in (degenerate_segment,degenerate_triangle,degenerate_tetrahedron)
        @test_throws ArgumentError _POINT_LOCATION.SimplexLocator(bad)
        @test_throws ArgumentError _POINT_LOCATION.mesh_local_coordinates(
            bad,1,0.0,0.0,0.0)
    end

    allocation_small=_point_location_segment_fixture(5_000)
    allocation_large=_point_location_segment_fixture(10_000)
    _POINT_LOCATION.SimplexLocator(
        _point_location_segment_fixture(10))
    small_build=_point_location_build_allocation(allocation_small)
    large_build=_point_location_build_allocation(allocation_large)
    @test small_build>0
    @test large_build>small_build
    @test large_build<=2.2small_build+65_536
    @test Base.summarysize(_POINT_LOCATION.SimplexLocator(allocation_large))<
          2_000_000

    @test _point_location_digest(fixture,locator)==
          "5c15e992466524cdc496c5f00d9e497bd7375eb88d7438a0573c5dc767f82d07"
    @test isempty(Docs.undocumented_names(
        Tessella.MeshPointLocation;private=false))
    @test isempty(Test.detect_ambiguities(
        Tessella.MeshPointLocation;recursive=true))
end

@testset "fixed-seed affine-coordinate oracle" begin
    rng=Xoshiro(0x91d4c72a6be50318)
    for _ in 1:64
        origin=ntuple(_->randn(rng),3)
        edge=ntuple(_->randn(rng),3)
        while hypot(edge...)<0.25
            edge=ntuple(_->randn(rng),3)
        end
        parameter=rand(rng)
        point=ntuple(axis->origin[axis]+parameter*edge[axis],3)
        coordinates=hcat(collect(origin),collect(origin .+ edge))
        mesh=Mesh(coordinates;segs=reshape(Int32[1,2],2,1))
        expected=(2parameter-1,0.0,0.0)
        @test isapprox(collect(_POINT_LOCATION.mesh_local_coordinates(
            mesh,1,point...)),collect(expected);atol=2e-12,rtol=2e-12)
        @test _POINT_LOCATION.locate_elements(
            _POINT_LOCATION.SimplexLocator(mesh),point...;strict=true)==UInt64[1]
    end

    for _ in 1:64
        origin=ntuple(_->randn(rng),3)
        edge1=ntuple(_->randn(rng),3)
        edge2=ntuple(_->randn(rng),3)
        cross=(edge1[2]*edge2[3]-edge1[3]*edge2[2],
               edge1[3]*edge2[1]-edge1[1]*edge2[3],
               edge1[1]*edge2[2]-edge1[2]*edge2[1])
        while hypot(cross...)<0.25
            edge2=ntuple(_->randn(rng),3)
            cross=(edge1[2]*edge2[3]-edge1[3]*edge2[2],
                   edge1[3]*edge2[1]-edge1[1]*edge2[3],
                   edge1[1]*edge2[2]-edge1[2]*edge2[1])
        end
        u=0.05+0.4rand(rng);v=0.05+(0.85-u)*rand(rng)
        point=ntuple(
            axis->origin[axis]+u*edge1[axis]+v*edge2[axis],3)
        coordinates=hcat(collect(origin),collect(origin .+ edge1),
                          collect(origin .+ edge2))
        mesh=Mesh(coordinates;tris=reshape(Int32[1,2,3],3,1))
        @test isapprox(collect(_POINT_LOCATION.mesh_local_coordinates(
            mesh,1,point...)),[u,v,0.0];atol=2e-12,rtol=2e-12)
        @test _POINT_LOCATION.locate_elements(
            _POINT_LOCATION.SimplexLocator(mesh),point...;strict=true)==UInt64[1]
    end

    for _ in 1:64
        origin=ntuple(_->randn(rng),3)
        edge1=ntuple(_->randn(rng),3)
        edge2=ntuple(_->randn(rng),3)
        edge3=ntuple(_->randn(rng),3)
        determinant=(
            edge1[1]*(edge2[2]*edge3[3]-edge2[3]*edge3[2])-
            edge1[2]*(edge2[1]*edge3[3]-edge2[3]*edge3[1])+
            edge1[3]*(edge2[1]*edge3[2]-edge2[2]*edge3[1]))
        while abs(determinant)<0.25
            edge3=ntuple(_->randn(rng),3)
            determinant=(
                edge1[1]*(edge2[2]*edge3[3]-edge2[3]*edge3[2])-
                edge1[2]*(edge2[1]*edge3[3]-edge2[3]*edge3[1])+
                edge1[3]*(edge2[1]*edge3[2]-edge2[2]*edge3[1]))
        end
        u=0.05+0.2rand(rng);v=0.05+0.2rand(rng)
        w=0.05+0.2rand(rng)
        point=ntuple(axis->origin[axis]+u*edge1[axis]+
                           v*edge2[axis]+w*edge3[axis],3)
        coordinates=hcat(
            collect(origin),collect(origin .+ edge1),
            collect(origin .+ edge2),collect(origin .+ edge3))
        mesh=Mesh(coordinates;tets=reshape(Int32[1,2,3,4],4,1))
        @test isapprox(collect(_POINT_LOCATION.mesh_local_coordinates(
            mesh,1,point...)),[u,v,w];atol=2e-12,rtol=2e-12)
        @test _POINT_LOCATION.locate_elements(
            _POINT_LOCATION.SimplexLocator(mesh),point...;strict=true)==UInt64[1]
    end
end
