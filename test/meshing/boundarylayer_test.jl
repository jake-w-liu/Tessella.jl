using Test
using Tessella
using Tessella.MeshTypes: ntris, triangle_area, tet_volume
using Tessella.Elements: validate, mixed_crc

function _bl_quad_area(mesh)
    area=0.0
    for block in mesh.blocks
        if block.msh==3
            for cell in axes(block.nodes,2)
                p1=ntuple(d->mesh.coords[d,block.nodes[1,cell]],3)
                p2=ntuple(d->mesh.coords[d,block.nodes[2,cell]],3)
                p3=ntuple(d->mesh.coords[d,block.nodes[3,cell]],3)
                p4=ntuple(d->mesh.coords[d,block.nodes[4,cell]],3)
                area+=triangle_area(p1,p2,p4)+triangle_area(p4,p2,p3)
            end
        elseif block.msh==2
            for cell in axes(block.nodes,2)
                p1=ntuple(d->mesh.coords[d,block.nodes[1,cell]],3)
                p2=ntuple(d->mesh.coords[d,block.nodes[2,cell]],3)
                p3=ntuple(d->mesh.coords[d,block.nodes[3,cell]],3)
                area+=triangle_area(p1,p2,p3)
            end
        end
    end
    return area
end

_bl_pt(m,i)=ntuple(d->m.coords[d,i],3)

function _bl_mixed_volumes(m)
    vp=vt=0.0
    for block in m.blocks
        if block.msh==6
            for c in axes(block.nodes,2)
                vp+=Tessella.BoundaryLayer._prism_volume6(m.coords,
                    ntuple(i->Int(block.nodes[i,c]),6))
            end
        elseif block.msh==4
            for c in axes(block.nodes,2)
                vt+=tet_volume(_bl_pt(m,block.nodes[1,c]),_bl_pt(m,block.nodes[2,c]),
                               _bl_pt(m,block.nodes[3,c]),_bl_pt(m,block.nodes[4,c]))
            end
        end
    end
    return vp,vt
end

@testset "prismatic boundary layer" begin
    # Unit square in z=0, two triangles.
    coords=Float64[0 1 1 0; 0 0 1 1; 0 0 0 0]
    tris=Int32[1 1; 2 3; 3 4]
    face=Mesh(coords; tris=tris)
    bl=mesh_boundary_layer(face; hwall=0.1, ratio=1.2, nlayers=2)
    @test validate(bl).ok
    @test mixed_crc(bl).sha=="f37a2141b37471b366cdbcaa5b1ede69c9088833b0f54e142aaa945e4ea23651"
    @test size(bl.coords,2)==4*3
    hexes=only(b for b in bl.blocks if b.msh==6)
    @test size(hexes.nodes,2)==ntris(face)*2
    # Independent volume: unit square area 1 times total height
    # H = hwall*((ratio^nl-1)/(ratio-1)) = 0.1*(1.44-1)/0.2 = 0.22
    H=0.1*(1.2^2-1)/(1.2-1)
    @test H≈0.22
    # Prism volume for a straight extrusion of a planar face is area * height.
    # Sample the first layer height hwall=0.1, area=1, two triangles → two prisms
    # of combined volume 0.1; second layer 0.12; total 0.22.
    zmax=maximum(bl.coords[3,:])
    @test zmax≈H atol=1e-12
    @test_throws ArgumentError mesh_boundary_layer(face; hwall=0.1, ratio=1.0, nlayers=1)
    @test_throws ArgumentError mesh_boundary_layer(face; hwall=-0.1, ratio=1.2, nlayers=1)
    @test_throws ArgumentError mesh_boundary_layer(face; hwall=true, ratio=1.2, nlayers=1)
    @test_throws ArgumentError mesh_boundary_layer(face; hwall=0.1, ratio=true, nlayers=1)
    @test_throws ArgumentError mesh_boundary_layer(face; hwall=0.1, ratio=1.2, nlayers=true)
    @test_throws ArgumentError mesh_boundary_layer(face; hwall=0.1, ratio=1.2,
                                                    nlayers=big(typemax(Int))+1)
    @test_throws ArgumentError mesh_boundary_layer(face; hwall=0.1, ratio=1.2,
                                                    nlayers=1, max_prisms=true)
    @test_throws ArgumentError mesh_boundary_layer(face; hwall=0.1, ratio=1.2,
                                                    nlayers=1, max_prisms=0)
    @test_throws ArgumentError mesh_boundary_layer(face; hwall=floatmax(Float64),
                                                    ratio=2.0, nlayers=2)
    @test_throws ArgumentError mesh_boundary_layer(face; hwall=0.1, ratio=1.2,
                                                    nlayers=typemax(Int),
                                                    max_prisms=typemax(Int))
    isolated=Mesh(hcat(coords,[2.0,2.0,0.0]);tris=tris)
    @test_throws ArgumentError mesh_boundary_layer(isolated;hwall=0.1,ratio=1.2,
                                                    nlayers=1)
    invalid=Mesh(coords;tris=reshape(Int32[1,1,2],3,1))
    @test !Tessella.MeshTypes.validate(invalid).ok
    @test_throws ArgumentError mesh_boundary_layer(invalid;hwall=0.1,ratio=1.2,
                                                   nlayers=1)
end

@testset "2-D boundary-layer quads" begin
    # Unit segment along x, left-normal +y. Area = length * H.
    coords=Float64[0 1; 0 0; 0 0]
    segs=Int32[1; 2;;]
    edge=Mesh(coords; segs=segs)
    hw=0.1; ra=1.2; nl=2
    H=hw*(ra^nl-1)/(ra-1)
    bl=mesh_boundary_layer_2d(edge; hwall=hw, ratio=ra, nlayers=nl)
    @test validate(bl).ok
    @test mixed_crc(bl).sha=="bb9e1fb9a0f0e56de42287ff3f85dd93ea3c7115cc9e51d05a6845d97ce8122b"
    quads=only(b for b in bl.blocks if b.msh==3)
    @test size(quads.nodes,2)==nl
    @test _bl_quad_area(bl)≈H atol=1e-12
    @test maximum(bl.coords[2,:])≈H atol=1e-12
    @test_throws ArgumentError mesh_boundary_layer_2d(edge; hwall=hw, ratio=1.0, nlayers=1)
    reversed=Mesh(coords;segs=Int32[2;1;;])
    reverse_bl=mesh_boundary_layer_2d(reversed;hwall=hw,ratio=ra,nlayers=1)
    @test minimum(reverse_bl.coords[2,:])≈-hw atol=1e-12
    @test maximum(reverse_bl.coords[2,:])==0.0
    @test_throws ArgumentError mesh_boundary_layer_2d(edge;hwall=true,ratio=ra,nlayers=1)
    @test_throws ArgumentError mesh_boundary_layer_2d(edge;hwall=hw,ratio=ra,nlayers=true)
    @test_throws ArgumentError mesh_boundary_layer_2d(edge;hwall=hw,ratio=ra,nlayers=1,
                                                      fan_elements=true)
    @test_throws ArgumentError mesh_boundary_layer_2d(edge;hwall=hw,ratio=ra,nlayers=1,
                                                      max_cells=true)
    @test_throws ArgumentError mesh_boundary_layer_2d(edge;hwall=hw,ratio=ra,nlayers=1,
                                                      max_cells=0)
    @test_throws ArgumentError mesh_boundary_layer_2d(edge;hwall=hw,ratio=ra,
                                                      nlayers=typemax(Int),
                                                      max_cells=typemax(Int))

    incoherent=Mesh(Float64[0 1 2;0 0 0;0 0 0];
                    segs=reshape(Int32[1,2,3,2],2,2))
    @test_throws ArgumentError mesh_boundary_layer_2d(incoherent;hwall=hw,ratio=ra,
                                                      nlayers=1)
    isolated=Mesh(Float64[0 1 2;0 0 0;0 0 0];segs=Int32[1;2;;])
    @test_throws ArgumentError mesh_boundary_layer_2d(isolated;hwall=hw,ratio=ra,
                                                      nlayers=1)
end

@testset "2-D boundary-layer closed square" begin
    coords=Float64[0 1 1 0; 0 0 1 1; 0 0 0 0]
    segs=Int32[1 2 3 4; 2 3 4 1]
    square=Mesh(coords; segs=segs)
    hw=0.05; ra=1.4; nl=1
    H=hw*(ra^nl-1)/(ra-1)
    @test H==hw
    bl=mesh_boundary_layer_2d(square; hwall=hw, ratio=ra, nlayers=nl)
    @test validate(bl).ok
    quads=only(b for b in bl.blocks if b.msh==3)
    @test size(quads.nodes,2)==4
    # Vertex-normal extrusion: inner side 1-H√2, strip area 2√2 H - 2 H^2.
    expected=2*sqrt(2)*H - 2*H^2
    @test _bl_quad_area(bl)≈expected atol=1e-12
end

@testset "2-D boundary-layer corner fan" begin
    coords=Float64[0 1 1; 0 0 1; 0 0 0]
    segs=reshape(Int32[1,2,2,3],2,2)
    ell=Mesh(coords; segs=segs)
    hw=0.1; ra=1.5; nl=2; nfan=2
    H=hw*(ra^nl-1)/(ra-1)
    bl=mesh_boundary_layer_2d(ell; hwall=hw, ratio=ra, nlayers=nl, fans=(2,), fan_elements=nfan)
    @test validate(bl).ok
    tris=only(b for b in bl.blocks if b.msh==2)
    quads=only(b for b in bl.blocks if b.msh==3)
    @test size(tris.nodes,2)==nfan
    @test size(quads.nodes,2)==2*nl + nfan*(nl-1)
    # Two unit legs of height H plus a polygonal quarter-disk of radii
    # r_k = hwall*(ratio^k-1)/(ratio-1). Each of nfan equal wedges has
    # area 0.5*sin(π/(2 nfan))*(r_k^2 - r_{k-1}^2).
    α=π/(2*nfan)
    r0=0.0
    fan_area=0.0
    for k in 1:nl
        rk=hw*(ra^k-1)/(ra-1)
        fan_area+=nfan*0.5*sin(α)*(rk^2-r0^2)
        r0=rk
    end
    @test _bl_quad_area(bl)≈2*H + fan_area atol=1e-12
    @test_throws ArgumentError mesh_boundary_layer_2d(ell; hwall=hw, ratio=ra, nlayers=nl,
                                                     fans=(1,), fan_elements=nfan)
    @test_throws ArgumentError mesh_boundary_layer_2d(ell; hwall=hw, ratio=ra, nlayers=nl,
                                                     fans=(2,), fan_elements=1)
    @test_throws ArgumentError mesh_boundary_layer_2d(ell; hwall=hw, ratio=ra, nlayers=nl,
                                                     fans=(true,), fan_elements=nfan)
    @test_throws ArgumentError mesh_boundary_layer_2d(ell; hwall=hw, ratio=ra, nlayers=nl,
                                                     fans=(2.0,), fan_elements=nfan)
    @test_throws ArgumentError mesh_boundary_layer_2d(ell; hwall=hw, ratio=ra, nlayers=nl,
                                                     fans=(big(typemax(Int))+1,),
                                                     fan_elements=nfan)
end

@testset "filled prismatic boundary layer" begin
    # Unit cube wall: two inward layers leave a core that must close the volume.
    cube=Tessella.Geometry.box_surface(0.0,1.0,0.0,1.0,0.0,1.0)
    m=mesh_boundary_layer_filled(cube; hwall=0.1, ratio=1.2, nlayers=2)
    @test validate(m).ok
    nv=8
    prisms=only(b for b in m.blocks if b.msh==6)
    tets=only(b for b in m.blocks if b.msh==4)
    @test size(prisms.nodes,2)==ntris(cube)*2
    @test size(tets.nodes,2)>0
    vp,vt=_bl_mixed_volumes(m)
    @test vp+vt≈1.0 rtol=1e-12
    # Interface conformity: every last-layer node is shared with the tet block.
    capids=Set(2nv+1:3nv)
    tetnodes=Set(Int(tets.nodes[i,c]) for c in axes(tets.nodes,2), i in 1:4)
    @test capids ⊆ tetnodes

    # Annulus: outer solid wall + inner cavity wall.
    outer=Tessella.Geometry.sphere_surface((0.0,0.0,0.0),1.0)
    inner=Tessella.Geometry.sphere_surface((0.0,0.0,0.0),0.4)
    n1=size(outer.coords,2)
    shell=Mesh(hcat(outer.coords,inner.coords);
               tris=hcat(outer.tris, inner.tris .+ Int32(n1)))
    ann=mesh_boundary_layer_filled(shell; hwall=0.04, ratio=1.25, nlayers=2,
                                   cavities=(2,))
    @test validate(ann).ok
    _,at=_bl_mixed_volumes(ann)
    Vout=4/3*pi; Vin=4/3*pi*0.4^3
    # The vertex-normal offset eats more than an ideal offset, so the filled
    # core must sit strictly between zero and the analytic annulus volume.
    @test 0 < at < (Vout-Vin)

    # Blockers: oversized layers fold the shell and violate its identity.
    @test_throws ErrorException mesh_boundary_layer_filled(cube; hwall=0.6,
                                                           ratio=1.05, nlayers=2)
    # Open surface has no enclosed interior to fill.
    open_sq=Mesh(Float64[0 1 1; 0 0 1; 0 0 0]; tris=[1 2 3]')
    @test_throws ArgumentError mesh_boundary_layer_filled(open_sq; hwall=0.1,
                                                          ratio=1.2, nlayers=2)
    # Cavity bookkeeping.
    @test_throws ArgumentError mesh_boundary_layer_filled(shell; hwall=0.04,
        ratio=1.25, nlayers=2, cavities=(1,2))
    @test_throws ArgumentError mesh_boundary_layer_filled(shell; hwall=0.04,
        ratio=1.25, nlayers=2, cavities=(5,))
    @test_throws ArgumentError mesh_boundary_layer_filled(shell; hwall=0.04,
        ratio=1.25, nlayers=2, cavities=(2,2))
    @test_throws ArgumentError mesh_boundary_layer_filled(cube;hwall=true,ratio=1.2,
                                                          nlayers=1)
    @test_throws ArgumentError mesh_boundary_layer_filled(cube;hwall=0.1,ratio=1.2,
                                                          nlayers=true)
    @test_throws ArgumentError mesh_boundary_layer_filled(cube;hwall=0.1,ratio=1.2,
                                                          nlayers=1,max_prisms=true)
    @test_throws ArgumentError mesh_boundary_layer_filled(cube;hwall=0.1,ratio=1.2,
                                                          nlayers=1,max_tets=true)
    @test_throws ArgumentError mesh_boundary_layer_filled(cube;hwall=0.1,ratio=1.2,
                                                          nlayers=1,cavities=(true,))
    @test_throws ArgumentError mesh_boundary_layer_filled(cube;hwall=0.1,ratio=1.2,
                                                          nlayers=1,
                                                          cavities=(big(typemax(Int))+1,))
end

@testset "boundary-layer public documentation" begin
    @test isempty(Docs.undocumented_names(Tessella.BoundaryLayer;private=false))
end
