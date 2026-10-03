using Test
using Tessella
using Tessella.Elements: MixedMesh, ElementBlock, _VOLUME_CELL_FACES
using Tessella.Model: mesh_model_volume, mesh_model_surface, model_to_mixed
using Tessella.MeshTypes: nnodes, ntris, ntets, validate

const _QT_EXTRUDE_ALLOCATION_RESULT=Ref((0,0))

function _qt_extrude_square(n;quad=true)
    return """
    Point(1)={0,0,0,1}; Point(2)={1,0,0,1};
    Point(3)={1,1,0,1}; Point(4)={0,1,0,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
    Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
    Transfinite Curve{:}=$(n+1);Transfinite Surface{1};
    $(quad ? "Recombine Surface{1};" : "")
    """
end

function _qt_extrude_execute(source;dim=3)
    mktempdir() do directory
        path=joinpath(directory,"quadtri.geo")
        write(path,source)
        return execute_geo(path;mesh_dim=dim)
    end
end

function _qt_extrude_counts(mesh)
    mesh isa Mesh && return Dict(2=>ntris(mesh),4=>ntets(mesh))
    return Dict(Int(b.msh)=>size(b.nodes,2) for b in mesh.blocks if b isa ElementBlock)
end

function _qt_extrude_two_squares(n;both_quad=true)
    return _qt_extrude_square(n)*"""
    Point(5)={2,0,0,1};Point(6)={2,1,0,1};
    Line(5)={2,5};Line(6)={5,6};Line(7)={6,3};
    Curve Loop(2)={5,6,7,-2};Plane Surface(2)={2};
    Transfinite Curve{5,6,7}=$(n+1);Transfinite Surface{2};
    $(both_quad ? "Recombine Surface{2};" : "")
    """
end

function _qt_extrude_boundary(mesh)
    faces=Dict{Tuple,Int}()
    blocks=mesh isa Mesh ? [ElementBlock(4,mesh.tets)] : mesh.blocks
    for b in blocks
        b isa ElementBlock && haskey(_VOLUME_CELL_FACES,Int(b.msh)) || continue
        for i in axes(b.nodes,2),face in _VOLUME_CELL_FACES[Int(b.msh)]
            key=Tuple(sort!([Tuple(mesh.coords[:,b.nodes[k,i]]) for k in face]))
            faces[key]=get(faces,key,0)+1
        end
    end
    @test all(c in (1,2) for c in values(faces))
    return Set(k for (k,c) in faces if c==1)
end

function _qt_extrude_surface_faces(execution;tags=nothing)
    faces=Set{Tuple}()
    for (dim,tag,part) in execution.mesh_parts
        dim==2 || continue
        tags===nothing || tag in tags || continue
        blocks=part isa Mesh ? [ElementBlock(2,part.tris)] : part.blocks
        for b in blocks
            b isa ElementBlock && b.msh in (2,3) || continue
            for i in axes(b.nodes,2)
                push!(faces,Tuple(sort!([Tuple(part.coords[:,v]) for v in b.nodes[:,i]])))
            end
        end
    end
    return faces
end

@testset "QuadTriAddVerts extrusion" begin
    # Native counts pinned to Gmsh 4.15.2, including transition pyramids and
    # the untouched core. The boundary must match without hidden diagonals.
    for (n,layers,laterals,nnode,volume_counts) in (
        (1,1,false,9,Dict(4=>10,7=>1)),
        (1,1,true,9,Dict(4=>2,7=>5)),
        (2,2,false,35,Dict(4=>40,7=>28)),
        (2,2,true,31,Dict(4=>8,5=>4,7=>20)),
        (3,3,false,89,Dict(4=>90,5=>2,7=>105)),
        (3,3,true,73,Dict(4=>18,5=>18,7=>45)))
        source=_qt_extrude_square(n)*"""
        Extrude{0,0,1}{Surface{1};Layers{$layers};Recombine;
        QuadTriAddVerts$(laterals ? " RecombLaterals" : "");}
        """
        execution=_qt_extrude_execute(source)
        volume=geo_entity_mesh(execution,3,1)
        @test nnodes(execution.mesh)==nnode
        @test _qt_extrude_counts(volume)==volume_counts
        if n==2 && layers==2 && laterals
            @test mixed_crc(execution.mesh)==(n_nodes=31,n_blocks=21,n_cells=84,
                bbox=((0.0,0.0,0.0),(1.0,1.0,1.0)),
                sha="731dff627c14f4982ac1dcab342c69f49c07ed13e1d78b35b64f72d79eb5a2d4")
            @test mixed_crc(volume).sha==
                "510be5d0090dcf135ac3f26f8dd9b89e86a186acd643ed9289d1226abf3e4377"
        end
        @test validate(execution.mesh).ok
        @test validate(volume).ok
        @test _qt_extrude_boundary(volume)==_qt_extrude_surface_faces(execution)
        direct=mesh_model_volume(execution.model,1)
        @test direct.coords==volume.coords
        @test _qt_extrude_counts(direct)==volume_counts
        projected=model_to_mixed(execution.model,volume,3,1)
        @test validate(projected).ok
        @test _qt_extrude_counts(projected)[7]==volume_counts[7]
    end
    for (laterals,nnode,counts) in (
        (false,10,Dict(4=>12,7=>2)),(true,8,Dict(6=>2)))
        source=_qt_extrude_square(1;quad=false)*"""
        Extrude{0,0,1}{Surface{1};Layers{1};Recombine;
        QuadTriAddVerts$(laterals ? " RecombLaterals" : "");}
        """
        execution=_qt_extrude_execute(source)
        volume=geo_entity_mesh(execution,3,1)
        @test nnodes(execution.mesh)==nnode
        @test _qt_extrude_counts(volume)==counts
        @test _qt_extrude_boundary(volume)==_qt_extrude_surface_faces(execution)
        @test validate(volume).ok
    end
    @testset "inert on curve extrusions and unrecombined volumes" begin
        curve="""
        Point(1)={0,0,0};Point(2)={1,0,0};Line(1)={1,2};
        Transfinite Curve{1}=3;
        Extrude{0,0,1}{Curve{1};Layers{2};Recombine;QuadTriAddVerts;}
        """
        execution=_qt_extrude_execute(curve;dim=2)
        @test sum(size(b.nodes,2) for b in execution.mesh.blocks if b.msh==3)==4
        source=_qt_extrude_square(1;quad=false)*"""
        Extrude{0,0,1}{Surface{1};Layers{1};QuadTriAddVerts RecombLaterals;}
        """
        execution=_qt_extrude_execute(source)
        @test ntets(execution.mesh)==6
        @test nnodes(execution.mesh)==8
        @test validate(execution.mesh).ok
    end
    @testset "invalid sweeps fail explicitly" begin
        source=_qt_extrude_square(1)*"Extrude{0,0,1}{Surface{1};Layers{2};}\n"
        @test_throws r"Cannot extrude quadrangles without Recombine" _qt_extrude_execute(source)
        source=_qt_extrude_square(1;quad=false)*"Extrude{{0,0,1},{0,0,0},Pi/4}{Surface{1};Layers{2};}\n"
        @test_throws r"degenerate or nonfinite volume" _qt_extrude_execute(source)
    end
    @testset "shared laterals respect neighboring volume cells" begin
        for neighbor in (:quadtri,:recombined,:tetrahedral)
            source=_qt_extrude_two_squares(1;both_quad=neighbor!==:tetrahedral)
            source*="""
            Extrude{0,0,1}{Surface{1};Layers{1};Recombine;QuadTriAddVerts;}
            Extrude{0,0,1}{Surface{2};Layers{1};
            $(neighbor!==:tetrahedral ? "Recombine;" : "")
            $(neighbor===:quadtri ? "QuadTriAddVerts;" : "")}
            """
            execution=_qt_extrude_execute(source)
            left=geo_entity_mesh(execution,3,1)
            right=geo_entity_mesh(execution,3,2)
            expected=neighbor===:tetrahedral ? Dict(4=>10,7=>1) : Dict(4=>8,7=>2)
            @test _qt_extrude_counts(left)==expected
            @test validate(execution.mesh).ok
            shared=intersect(_qt_extrude_boundary(left),_qt_extrude_boundary(right))
            @test length(shared)==(neighbor===:tetrahedral ? 2 : 1)
            @test all(all(p[1]==1.0 for p in face) for face in shared)
            for (dim,tag,part) in execution.mesh_parts
                dim==2 || continue
                all(part.coords[1,:].==1.0) || continue
                @test Set(_qt_extrude_surface_faces((mesh_parts=[(dim,tag,part)],)))==shared
            end
        end
    end
    @testset "graded translation, revolution, and twist" begin
        transforms=(
            "Extrude{0,0,2}",
            "Extrude{{0,1,0},{-1,0,0},Pi/3}",
            "Extrude{{0,0,2},{0,0,1},{0,0,0},Pi/3}")
        for transform in transforms,laterals in (false,true)
            source=_qt_extrude_square(2)*transform*"""
            {Surface{1};Layers{{2,3},{0.4,1.0}};Recombine;
            QuadTriAddVerts$(laterals ? " RecombLaterals" : "");}
            """
            execution=_qt_extrude_execute(source)
            volume=geo_entity_mesh(execution,3,1)
            if occursin("Extrude{{0,0,2}",transform)
                # Gmsh keeps each orphan spline control-point entity as a
                # separate mesh node, even at identical connector positions.
                @test nnodes(execution.mesh)==(laterals ? 74 : 90)
            end
            @test validate(execution.mesh).ok
            @test validate(volume).ok
            @test _qt_extrude_boundary(volume)==_qt_extrude_surface_faces(execution)
            again=_qt_extrude_execute(source)
            @test again.mesh.coords==execution.mesh.coords
            @test _qt_extrude_counts(again.mesh)==_qt_extrude_counts(execution.mesh)
        end
    end
    @testset "output-sized allocation growth" begin
        source(n)=_qt_extrude_square(n)*"""
        Extrude{0,0,1}{Surface{1};Layers{3};Recombine;QuadTriAddVerts;}
        """
        small=_qt_extrude_execute(source(4);dim=0).model
        large=_qt_extrude_execute(source(8);dim=0).model
        mesh_model_volume(small,1);mesh_model_volume(large,1)
        a=minimum(@allocated(mesh_model_volume(small,1)) for _ in 1:3)
        b=minimum(@allocated(mesh_model_volume(large,1)) for _ in 1:3)
        _QT_EXTRUDE_ALLOCATION_RESULT[]=(a,b)
        @test b<=5a
    end
    @testset "centroid fans certify faces before orientation repair" begin
        v=((0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,1.0,0.0),(0.0,1.0,0.0),
           (0.0,0.0,1.0),(1.0,0.0,1.0),(1.0,1.0,1.0),(0.0,1.0,1.0))
        empty_edges=Set{NTuple{2,NTuple{3,Float64}}}()
        center=Tessella.Model._extrude_quadtri_centroid(v)
        @test Tessella.Model._extrude_quadtri_certify(v,center,empty_edges)===nothing
        folded=ntuple(k->k==7 ? (1.0,1.0,-1.0) : v[k],8)
        center=Tessella.Model._extrude_quadtri_centroid(folded)
        @test_throws r"folded or degenerate" Tessella.Model._extrude_quadtri_certify(
            folded,center,empty_edges)
        conflicts=Set((Tessella.Model._extrude_ekey(v[5],v[7]),
                       Tessella.Model._extrude_ekey(v[6],v[8])))
        @test_throws r"conflicting diagonals" Tessella.Model._extrude_quadtri_certify(
            v,Tessella.Model._extrude_quadtri_centroid(v),conflicts)
        reverse=ntuple(k->(v[k][1],v[k][2],-v[k][3]),8)
        @test Tessella.Model._extrude_quadtri_certify(
            reverse,Tessella.Model._extrude_quadtri_centroid(reverse),empty_edges)===nothing
        distorted=[(0.0,0.0,0.0),(3.0,0.0,0.0),(3.0,2.0,0.0),(0.0,1.0,0.0)]
        @test !Tessella.Model._extrude_quadtri_top_diagonal(distorted,(1,2,3,4),false)
        @test Tessella.Model._extrude_quadtri_top_diagonal(distorted,(1,2,3,4),true)
    end
    @testset "revolution with fixed columns" begin
        for axis in ("0,1,0","1,-1,0"),laterals in (false,true)
            source=_qt_extrude_square(1)*"""
            Extrude{{$axis},{0,0,0},Pi/3}{Surface{1};Layers{3};Recombine;
            QuadTriAddVerts$(laterals ? " RecombLaterals" : "");}
            """
            execution=_qt_extrude_execute(source)
            volume=geo_entity_mesh(execution,3,1)
            @test validate(execution.mesh).ok
            @test validate(volume).ok
            @test _qt_extrude_boundary(volume)==_qt_extrude_surface_faces(execution)
        end
    end
    @testset "closed revolutions preserve section quads and source seam" begin
        for (n,layers,stages,laterals,nnode) in (
                (1,1,(1,2,3,4),false,22),(1,1,(1,2,3,4),true,22),
                (2,2,(1,2,3,4),false,106),(2,2,(1,2,3,4),true,90),
                (2,2,(1,3),false,90),(2,2,(1,3),true,82))
            source=_qt_extrude_square(n)
            source=replace(source,"{0,0,0,1}"=>"{1,0,0,1}","{1,0,0,1}"=>"{2,0,0,1}",
                "{1,1,0,1}"=>"{2,0,1,1}","{0,1,0,1}"=>"{1,0,1,1}")
            # The square helper includes whitespace, so transform its points
            # by replacing the literal coordinates rather than point tags.
            for stage in 1:4
                from=stage==1 ? "1" : "q$(stage-1)[0]"
                quadtri=stage in stages ? "QuadTriAddVerts$(laterals ? " RecombLaterals" : "");" : ""
                source*="q$stage[]=Extrude{{0,0,1},{0,0,0},Pi/2}{Surface{$from};Layers{$layers};Recombine;$quadtri};\n"
            end
            execution=_qt_extrude_execute(source)
            @test nnodes(execution.mesh)==nnode
            @test validate(execution.mesh).ok
            for tag in (1,26,48,70)
                @test _qt_extrude_counts(geo_entity_mesh(execution,2,tag))==Dict(3=>n*n)
            end
            for region in 1:4
                volume=geo_entity_mesh(execution,3,region)
                counts=region in stages && !laterals ?
                    (n==1 ? Dict(4=>8,7=>2) : Dict(4=>32,7=>32)) : Dict(5=>n*n*layers)
                @test _qt_extrude_counts(volume)==counts
                @test validate(volume).ok
                tags=Set(abs.(Tessella.Model._model_volume_boundary_surfaces(
                    execution.model,region,"QuadTri toroidal regression")))
                @test _qt_extrude_boundary(volume)==_qt_extrude_surface_faces(execution;tags=tags)
                direct=mesh_model_volume(execution.model,region)
                @test direct.coords==volume.coords
            end
            first=geo_entity_mesh(execution,3,1)
            last=geo_entity_mesh(execution,3,4)
            root=_qt_extrude_surface_faces(execution;tags=Set([1]))
            @test intersect(_qt_extrude_boundary(first),_qt_extrude_boundary(last))==root
        end
    end
end
