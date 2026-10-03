using Test, Tessella
using Tessella.Elements: ElementBlock, MixedMesh, msh_dimension
using Tessella.StructuredQuadTri: _transfinite_quadtri
using Tessella.Model: mesh_model_volume, model_to_mixed
using Tessella.MeshTypes: tet_signed_volume, nnodes

const _QT_BOX=read(joinpath(@__DIR__,"..","artifacts","quadtri_box.geo"),String)
const _QT_PRISM=read(joinpath(@__DIR__,"..","artifacts","quadtri_prism.geo"),String)

function _qt_execute(source)
    mktemp() do path,io
        write(io,source);close(io)
        execute_geo(path;mesh_dim=3)
    end
end

function _qt_source(geometry,mask;n=3,arrangement="Left",extra="",compact=false)
    recombine=isempty(mask) ? "" : "Recombine Surface{"*join(mask,",")*"};"
    return geometry*"""
        Mesh.TransfiniteTri=$(compact ? 1 : 0);
        Transfinite Curve{:}=$(n+1);
        Transfinite Surface{:} $arrangement;
        Transfinite Volume{1};
        TransfQuadTri{1};
        $recombine
        $extra
        """
end

_qt_counts(mesh)=Dict(b.msh=>size(b.nodes,2) for b in mesh.blocks
                      if b isa ElementBlock && msh_dimension(b.msh)==3)

# Independent analytic volumes and first-order MSH cell decompositions.
const _QT_DECOMPS=Dict(
    4=>((1,2,3,4),),
    5=>((1,2,3,7),(1,3,4,7),(1,4,8,7),(1,8,5,7),(1,5,6,7),(1,6,2,7)),
    6=>((1,2,3,6),(1,2,6,5),(1,5,6,4)),
    7=>((1,2,3,5),(1,3,4,5)))

function _qt_volume(mesh)
    result=0.0
    for b in mesh.blocks
        b isa ElementBlock && haskey(_QT_DECOMPS,b.msh) || continue
        for i in axes(b.nodes,2),tet in _QT_DECOMPS[b.msh]
            ps=ntuple(k->Tuple(mesh.coords[:,b.nodes[tet[k],i]]),4)
            volume=tet_signed_volume(ps...)
            @test isfinite(volume) && volume>0
            result+=volume
        end
    end
    return result
end

@testset "TransfQuadTri transfinite transitions" begin
    @testset "all six-face masks and arrangements" begin
        for maskbits in 0:63
            mask=[s for s in 1:6 if !iszero(maskbits & (1<<(s-1)))]
            execution=_qt_execute(_qt_source(_QT_BOX,mask))
            volume=geo_entity_mesh(execution,3,1)
            @test volume isa MixedMesh
            @test validate(volume).ok
            @test validate(execution.mesh).ok
            @test _qt_volume(volume)≈1.0 atol=1e-13
            # Classification includes each tetrahedron/pyramid block and
            # its source volume ownership, plus the actual boundary faces.
            mixed=model_to_mixed(execution.model,volume,3,1)
            @test validate(mixed).ok
            @test _qt_counts(mixed)==_qt_counts(volume)
            @test all(all(==(1),owners) for (block,owners) in
                zip(mixed.blocks,mixed.elementary_entities) if
                block isa ElementBlock && msh_dimension(block.msh)==3)
        end
        for arrangement in ("Left","Right","AlternateLeft","AlternateRight"),n in (1,2,4)
            execution=_qt_execute(_qt_source(_QT_BOX,[1,4];n,arrangement))
            @test validate(execution.mesh).ok
            @test _qt_volume(geo_entity_mesh(execution,3,1))≈1.0 atol=1e-13
        end
    end
    @testset "five-face prism oracle fixtures" begin
        # Gmsh 4.15.2 n=3 oracle counts: nodes, tetrahedra, pyramids, hexes, prisms.
        for (mask,counts) in ((Int[],(75,90,81,1,3)),
                (collect(1:5),(52,0,0,18,9)),
                ([1],(74,78,81,2,3)),([3],(71,66,75,2,6)),
                ([4],(74,72,84,2,3)),([3,4,5],(64,24,60,6,9)),
                ([1,3,4,5],(58,12,30,12,9)))
            execution=_qt_execute(_qt_source(_QT_PRISM,mask))
            part=geo_entity_mesh(execution,3,1)
            actual=_qt_counts(part)
            @test (nnodes(part),get(actual,4,0),get(actual,7,0),
                   get(actual,5,0),get(actual,6,0))==counts
            @test validate(execution.mesh).ok
            @test _qt_volume(part)≈0.5 atol=1e-13
            if mask==[3,4,5]
                # Full coordinate/connectivity/ownership pins are recorded
                # in test/artifacts/quadtri_crc.txt after the pinned oracle.
                @test mixed_crc(execution.mesh)==(n_nodes=64,n_blocks=18,n_cells=183,
                    bbox=((0.0,0.0,0.0),(1.0,1.0,1.0)),
                    sha="17d4afbbf28ef179760aff9bce7e1c44feb722f623af05a14b8cab76dbe277c1")
                @test mixed_crc(part).sha==
                    "3188df0104d79397f470294da0fa7cd3898fd6a4b9ec141f4f19ceb021708efc"
            end
        end
        execution=_qt_execute(_qt_source(_QT_PRISM,collect(1:5);compact=true))
        @test nnodes(execution.mesh)==46
        @test _qt_counts(geo_entity_mesh(execution,3,1))==Dict(6=>27)
        @test _qt_volume(geo_entity_mesh(execution,3,1))≈0.5 atol=1e-13
    end
    @testset "volume attributes and deterministic output" begin
        plain=_qt_execute(_qt_source(_QT_BOX,[1,2,3]))
        again=_qt_execute(_qt_source(_QT_BOX,[1,2,3]))
        @test mixed_crc(plain.mesh)==mixed_crc(again.mesh)
        @test mixed_crc(plain.mesh)==(n_nodes=85,n_blocks=21,n_cells=276,
            bbox=((0.0,0.0,0.0),(1.0,1.0,1.0)),
            sha="2938c7188f2ccb1159c9df3db102fa7e1c3ca43d07ce4cb86b9729de93973d69")
        # Every first-order volume family, including new pyramids, reverses.
        a=geo_entity_mesh(plain,3,1)
        @test mixed_crc(a).sha==
            "5ade95ed527ba619ff150c536b27d2c863efc8bba4ed933fd5ea85986c2200e1"
        Tessella.Model.set_reverse!(plain.model,3,1,true)
        b=mesh_model_volume(plain.model,1)
        @test _qt_counts(a)==_qt_counts(b)
        for (ab,bb) in zip(a.blocks,b.blocks)
            ab isa ElementBlock && ab.msh==7 || continue
            @test ab.nodes[1,:]==bb.nodes[3,:]
            @test ab.nodes[3,:]==bb.nodes[1,:]
        end
    end
    @testset "resource limits and linear allocation growth" begin
        bases=MixedMesh[];boundaries=[]
        for n in (3,6)
            ex=_qt_execute(_qt_source(_QT_BOX,Int[];n))
            m=ex.model
            delete!(m.meshing.quad_tri,1)
            Tessella.Model.set_recombine!(m,2,1)
            for s in 2:6;Tessella.Model.set_recombine!(m,2,s);end
            base=mesh_model_volume(m,1)
            for s in 1:6;delete!(m.meshing.recombine,(2,s));end
            parts=[(s,Tessella.Model.mesh_model_surface(m,s)) for s in 1:6]
            push!(bases,base);push!(boundaries,parts)
        end
        _transfinite_quadtri(bases[1],boundaries[1])
        _transfinite_quadtri(bases[2],boundaries[2])
        small=@allocated _transfinite_quadtri(bases[1],boundaries[1])
        large=@allocated _transfinite_quadtri(bases[2],boundaries[2])
        # Doubling each axis grows source cells by 8; memory must stay linear.
        @test large<9small
        @test_throws ArgumentError _transfinite_quadtri(bases[1],boundaries[1];max_nodes=64)
        @test_throws ArgumentError _transfinite_quadtri(bases[1],boundaries[1];max_tets=1)
    end
end
