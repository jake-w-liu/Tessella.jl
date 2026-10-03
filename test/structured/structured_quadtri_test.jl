using Test, Tessella
using Tessella.Elements: ElementBlock, MixedMesh
using Tessella.MeshTypes: tet_signed_volume
using Tessella.StructuredQuadTri: _transfinite_quadtri

# A single logical cube with independent 1--3 boundary diagonals. The
# transfinite kernel supplies only the base hex; the center fan is tested
# against an analytic solid rather than a second call to the fan kernel.
function _structured_quadtri_fixture()
    corners=[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),
             (0.,0.,1.),(1.,0.,1.),(1.,1.,1.),(0.,1.,1.)]
    base=mesh_transfinite_volume(corners,(1,1,1);recombine=true)
    quads=only(filter(b->b.msh==3,base.blocks))
    tris=Matrix{Int32}(undef,3,2size(quads.nodes,2))
    for i in axes(quads.nodes,2)
        a,b,c,d=quads.nodes[:,i]
        tris[1,2i-1]=a;tris[2,2i-1]=b;tris[3,2i-1]=c
        tris[1,2i]=a;tris[2,2i]=c;tris[3,2i]=d
    end
    return base,tris
end

function _structured_quadtri_boundary(base,tris)
    return [(1,MixedMesh(base.coords,
        [ElementBlock(2,tris,fill(Int32(1),size(tris,2)))]))]
end

@testset "Structured QuadTri centroid certification" begin
    # Raising the (1,1,1) corner by 0.2 creates two equal-area roof
    # triangles. Each triangle's mean height increases by 0.2/3, so the
    # solid's exact analytic volume is 1+0.2/3.
    base,tris=_structured_quadtri_fixture()
    for i in axes(base.coords,2)
        Tuple(base.coords[:,i])==(1.,1.,1.) && (base.coords[3,i]=1.2)
    end
    output=_transfinite_quadtri(base,_structured_quadtri_boundary(base,tris))
    tets=only(filter(b->b.msh==4,output.blocks))
    volume=0.
    for cell in axes(tets.nodes,2)
        a,b,c,d=ntuple(k->Tuple(output.coords[:,tets.nodes[k,cell]]),4)
        part=tet_signed_volume(a,b,c,d)
        @test part>0
        volume+=part
    end
    @test volume≈1+0.2/3 atol=1e-14
    @test size(output.coords,2)==9

    # MixedMesh validation checks structure, so folded and flat logical
    # cells reach this kernel with valid indices. Its exact orientation
    # certification must reject them before returning any center fan.
    for mode in (:folded,:flat)
        base,tris=_structured_quadtri_fixture()
        for i in axes(base.coords,2)
            if mode===:flat
                base.coords[3,i]=0.
            elseif Tuple(base.coords[:,i])==(1.,1.,1.)
                base.coords[3,i]=-1.
            end
        end
        @test validate(base).ok
        err=try
            _transfinite_quadtri(base,_structured_quadtri_boundary(base,tris))
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError
        @test occursin("folded",sprint(showerror,err))
    end

    # A repeated or missing surface cell cannot be silently accepted as a
    # complete boundary, even when the center fan itself is well oriented.
    for mode in (:duplicate,:missing)
        base,tris=_structured_quadtri_fixture()
        malformed=mode===:duplicate ? hcat(tris,tris[:,1]) : tris[:,2:end]
        @test_throws ArgumentError _transfinite_quadtri(
            base,_structured_quadtri_boundary(base,malformed))
    end
end
