module MultiPeriodicPairsChecks
using Test,Tessella
using Tessella.MeshTypes: validate

# PER02 multi-direction slice: two independent opposite-face translation
# pairs on a unit box. Gmsh 4.15.2 oracle on the same fixture accepts both
# declarations, generates the volume, and reports 2 surface links x 12 node
# pairs, 7 curve links x 3 pairs, and 6 point links x 1 pair in $Periodic —
# every link carrying its relation affine. Native slave copies reparametrize
# on the slave geometry, so coordinates satisfy the affine only up to the
# ~1e-8·lc `parFromPoint` residual rather than bitwise.

const AFF_X=Float64[1,0,0,1, 0,1,0,0, 0,0,1,0, 0,0,0,1]
const AFF_Y=Float64[1,0,0,0, 0,1,0,1, 0,0,1,0, 0,0,0,1]

function box(api)
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
            (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,point) in enumerate(points)
        api.model.add_point(point...;tag,meshSize=0.5)
    end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
           (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(edges)
        api.model.add_line(ends...;tag)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,lines) in enumerate(loops)
        api.model.add_curve_loop(collect(lines);tag)
        api.model.add_plane_surface([tag];tag)
    end
    api.model.add_surface_loop(collect(1:6);tag=1)
    api.model.add_volume([1];tag=1)
end

function affine_image(api,affine,point)
    coefficients,translation=
        Tessella.Transform._periodic_affine_3x4(Tuple(affine))
    return Tessella.Model._model_affine_point(
        coefficients,translation,Tuple(point),"test",0)
end

function check_pairs(api,dim,slave,count)
    relation=api.mesh.get_periodic_nodes(dim,slave)
    @test relation.master_entity>0
    @test length(relation.slave_nodes)==count
    @test length(relation.master_nodes)==count
    @test allunique(relation.slave_nodes)&&allunique(relation.master_nodes)
    isempty(relation.affine) && return relation
    for (slave_uid,master_uid) in
        zip(relation.slave_nodes,relation.master_nodes)
        slave_point=Tuple(first(api.mesh.get_node(slave_uid)))
        master_point=Tuple(first(api.mesh.get_node(master_uid)))
        expected=affine_image(api,relation.affine,master_point)
        @test isapprox(collect(slave_point),collect(expected);atol=1e-7)
    end
    return relation
end

@testset "multiple independent periodic translation pairs" begin
    api=Tessella.API
    api.initialize()
    try
        api.option("General.Terminal",0)
        api.option("Mesh.Renumber",0)
        box(api)
        api.mesh.set_periodic(2,[4],[6],AFF_X)
        api.mesh.set_periodic(2,[5],[3],AFF_Y)
        api.mesh.generate(3)
        @test validate(api.mesh.get()).ok
        mesh=Tessella.API.LAST_MESH[]

        # Declared surface links: Gmsh reports 12 pairs on each.
        x=check_pairs(api,2,4,12)
        @test x.master_entity==6
        y=check_pairs(api,2,5,12)
        @test y.master_entity==3

        # Induced boundary relations: 7 curve links x 3 pairs and 6 point
        # links x 1 pair, each slave resolving to exactly one master even
        # where the two directions share corners and edges.
        links=mesh.periodic_links
        @test links!==nothing
        bydim=Dict{Int,Vector{eltype(links)}}()
        for link in links
            push!(get!(()->eltype(links)[],bydim,link.dim),link)
        end
        @test length(get(bydim,0,[]))==6
        @test length(get(bydim,1,[]))==7
        @test length(get(bydim,2,[]))==2
        @test all(link->length(link.slave_nodes)==1,bydim[0])
        @test all(link->length(link.slave_nodes)==3,bydim[1])
        @test all(link->length(link.slave_nodes)==12,bydim[2])
        @test all(link->link.affine!==nothing,links)
        slave_points=sort!(Int[Int(link.slave_entity) for link in bydim[0]])
        @test allunique(slave_points)
        slave_curves=sort!(Int[Int(link.slave_entity) for link in bydim[1]])
        @test allunique(slave_curves)
        @test sort!(Int[Int(link.slave_entity) for link in bydim[2]])==[4,5]

        # Public node maps exist for every induced curve and point slave —
        # same counts the Gmsh oracle reports for the equivalent fixture.
        for curve in slave_curves
            check_pairs(api,1,curve,3)
        end
        for point in slave_points
            check_pairs(api,0,point,1)
        end

        # Serialization emits every relation: Gmsh's $Periodic block for
        # this fixture carries 15 links in the same dim split.
        path=joinpath(mktempdir(),"box2per.msh")
        api.write(path)
        text=read(path,String)
        section=match(r"\$Periodic\n(\d+)\n(.*)\$EndPeriodic"s,text)
        @test section!==nothing
        @test parse(Int,section.captures[1])==15
        body=split(strip(section.captures[2]),'\n';keepempty=false)
        dims=Int[]
        index=1
        while index<=length(body)
            fields=split(strip(body[index]))
            push!(dims,parse(Int,fields[1]))
            index+=1
            parse(Int,split(strip(body[index]))[1])==16 && (index+=1)
            index+=1+parse(Int,strip(body[index]))
        end
        @test index==length(body)+1
        @test sort!(dims)==[0,0,0,0,0,0,1,1,1,1,1,1,1,2,2]
    finally
        api.finalize()
    end
end
end
