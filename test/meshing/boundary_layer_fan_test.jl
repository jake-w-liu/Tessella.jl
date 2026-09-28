using Test
using Tessella
using Tessella.MeshTypes: Mesh, nnodes, ntris, tet_volume
using Tessella.Elements: validate, MixedMesh
using Tessella.Geometry: box_surface
using Tessella.BoundaryLayer: mesh_boundary_layer_fan, _fan_prism_volume

# ── independent oracle helpers ────────────────────────────────────────────────

@inline function _srt(a,b,c)
    a>b && ((a,b)=(b,a))
    b>c && ((b,c)=(c,b))
    a>b && ((a,b)=(b,a))
    return (a,b,c)
end

@inline _pt(c,i) = (c[1,i],c[2,i],c[3,i])

# Sorted-vertex union boundary faces of a MixedMesh, recomputed independently:
# prisms contribute bottom/top tris plus the canonically split side quads,
# tets their four faces. Returns (count-1 face set, max-incidence<=2 flag).
function _fan_boundary(m::MixedMesh)
    inc=Dict{NTuple{3,Int32},Int}()
    for b in m.blocks, c in 1:size(b.nodes,2)
        v=b.nodes[:,c]
        if b.msh==6
            faces=NTuple{3,Int32}[_srt(v[1],v[2],v[3]),_srt(v[4],v[5],v[6])]
            for (i,j) in ((1,2),(2,3),(3,1))
                a,b2,c2,d=v[i],v[j],v[j+3],v[i+3]
                if min(a,c2)<min(b2,d)
                    push!(faces,_srt(a,b2,c2)); push!(faces,_srt(a,c2,d))
                else
                    push!(faces,_srt(b2,c2,d)); push!(faces,_srt(b2,d,a))
                end
            end
            for f in faces
                inc[f]=get(inc,f,0)+1
            end
        elseif b.msh==4
            for f in (_srt(v[2],v[3],v[4]),_srt(v[1],v[3],v[4]),
                      _srt(v[1],v[2],v[4]),_srt(v[1],v[2],v[3]))
                inc[f]=get(inc,f,0)+1
            end
        end
    end
    bnd=Set{NTuple{3,Int32}}()
    ok=true
    for (f,cnt) in inc
        cnt>2 && (ok=false)
        cnt==1 && push!(bnd,f)
    end
    return bnd,ok
end

function _fan_input_boundary(surface::Mesh)
    Set{NTuple{3,Int32}}(_srt(surface.tris[1,t],surface.tris[2,t],surface.tris[3,t])
                         for t in 1:ntris(surface))
end

function _fan_volume(m::MixedMesh)
    v=0.0
    for b in m.blocks, c in 1:size(b.nodes,2)
        if b.msh==6
            v+=_fan_prism_volume(m.coords,Tuple(b.nodes[:,c]))
        elseif b.msh==4
            v+=tet_volume(_pt(m.coords,b.nodes[1,c]),_pt(m.coords,b.nodes[2,c]),
                          _pt(m.coords,b.nodes[3,c]),_pt(m.coords,b.nodes[4,c]))
        end
    end
    return v
end

# Column heights: the last-layer node emitted above input vertex v sits at
# distance sum(offsets); find it among the non-input nodes.
function _fan_column_top(surface::Mesh,m::MixedMesh,v)
    pv=_pt(surface.coords,v)
    dmin=Inf; hit=0
    for i in nnodes(surface)+1:size(m.coords,2)
        d=hypot(m.coords[1,i]-pv[1],m.coords[2,i]-pv[2],m.coords[3,i]-pv[3])
        if abs(d-0.25)<1e-12
            hit+=1
        end
        dmin=min(dmin,d)
    end
    return dmin,hit
end

@testset "boundary layer fan — joined multi-region layers" begin

    cube=box_surface(0.,1.,0.,1.,0.,1.)

    @testset "two-region shell: fans along the shared boundary loop" begin
        regions=([1,2],[3,4,5,6,7,8,9,10,11,12])
        m=mesh_boundary_layer_fan(cube;regions=regions,hwall=0.1,ratio=1.5,
                                  nlayers=2,fan_elements=3)
        @test validate(m).ok
        bnd,man=_fan_boundary(m)
        @test man
        @test bnd==_fan_input_boundary(cube)
        @test _fan_volume(m)≈1.0 atol=1e-9
        m2=mesh_boundary_layer_fan(cube;regions=regions,hwall=0.1,ratio=1.5,
                                   nlayers=2,fan_elements=3)
        @test m2.coords==m.coords
        @test all(b1.nodes==b2.nodes for (b1,b2) in zip(m2.blocks,m.blocks))
    end

    @testset "six-face regions: fans on every edge plus corner blocks" begin
        regions=([1,2],[3,4],[5,6],[7,8],[9,10],[11,12])
        m=mesh_boundary_layer_fan(cube;regions=regions,hwall=0.08,ratio=1.4,
                                  nlayers=2,fan_elements=4)
        @test validate(m).ok
        bnd,man=_fan_boundary(m)
        @test man
        @test bnd==_fan_input_boundary(cube)
        @test _fan_volume(m)≈1.0 atol=1e-9
    end

    @testset "adjacent layered faces with unlayered remainder" begin
        regions=([1,2],[5,6],[3,4,7,8,9,10,11,12])
        m=mesh_boundary_layer_fan(cube;regions=regions,hwall=0.1,ratio=1.5,
                                  nlayers=2,fan_elements=3,unlayered=(3,))
        @test validate(m).ok
        bnd,man=_fan_boundary(m)
        @test man
        @test bnd==_fan_input_boundary(cube)
        @test _fan_volume(m)≈1.0 atol=1e-9
        # vertex 1 belongs to the layered bottom+x0 regions; its offset column
        # top sits at distance hwall+hwall*ratio = 0.25
        dmin,hit=_fan_column_top(cube,m,1)
        @test hit>=1
    end

    @testset "uniform single region: no boundary edges, plain columns" begin
        regions=([i for i in 1:12],)
        m=mesh_boundary_layer_fan(cube;regions=regions,hwall=0.1,ratio=1.5,nlayers=2)
        @test validate(m).ok
        bnd,man=_fan_boundary(m)
        @test man
        @test bnd==_fan_input_boundary(cube)
        @test _fan_volume(m)≈1.0 atol=1e-9
        pr=[b for b in m.blocks if b.msh==6]
        @test length(pr)==1 && size(pr[1].nodes,2)==ntris(cube)*2
    end

    @testset "degenerate and invalid inputs raise precise errors" begin
        @test_throws ArgumentError mesh_boundary_layer_fan(cube;regions=([1,2],),hwall=0.1,ratio=1.5,nlayers=2)
        @test_throws ArgumentError mesh_boundary_layer_fan(cube;regions=([1,2],[2,3,4,5,6,7,8,9,10,11,12]),hwall=0.1,ratio=1.5,nlayers=2)
        @test_throws ArgumentError mesh_boundary_layer_fan(cube;regions=([0,2],[1,3,4,5,6,7,8,9,10,11,12]),hwall=0.1,ratio=1.5,nlayers=2)
        @test_throws ArgumentError mesh_boundary_layer_fan(cube;regions=([1,2],[3,4],[5,6,7,8,9,10,11,12]),hwall=0.1,ratio=1.5,nlayers=2,unlayered=(4,))
        @test_throws ArgumentError mesh_boundary_layer_fan(cube;regions=([1,2],[3,4,5,6,7,8,9,10,11,12]),hwall=0.1,ratio=1.5,nlayers=0)
        @test_throws ArgumentError mesh_boundary_layer_fan(cube;regions=([1,2],[3,4,5,6,7,8,9,10,11,12]),hwall=0.1,ratio=1.5,nlayers=2,fan_elements=0)
        @test_throws ArgumentError mesh_boundary_layer_fan(cube;regions=([1,2],[3,4,5,6,7,8,9,10,11,12]),hwall=-0.1,ratio=1.5,nlayers=2)
        @test_throws ArgumentError mesh_boundary_layer_fan(cube;regions=([1,2],[3,4,5,6,7,8,9,10,11,12]),hwall=0.1,ratio=1.0,nlayers=2)
        @test_throws ArgumentError mesh_boundary_layer_fan(cube;regions=([1,2],[3,4,5,6,7,8,9,10,11,12]),hwall=0.1,ratio=1.5,nlayers=2,unlayered=(1,2))
        # disconnected region
        @test_throws ArgumentError mesh_boundary_layer_fan(cube;regions=([1,12],[2,3,4,5,6,7,8,9,10,11]),hwall=0.1,ratio=1.5,nlayers=2)
        # oversized layer depth: explicit failure, never a defective mesh
        @test_throws Exception mesh_boundary_layer_fan(cube;regions=([1,2],[3,4,5,6,7,8,9,10,11,12]),hwall=10.0,ratio=1.5,nlayers=2)
    end
end
