module GeoPartitionMeshFunctionalChecks
# Functional checks for the `.geo` `PartitionMesh(n)` statement: the native
# OCC box mesh is the smallest input whose lower-BRep transport exercises
# native topology normalization (OCC wire/shell use flags vs actual cell
# directions), the multi-sheet interface junction (one nonmanifold valence-4
# edge where two partition-2 wedges meet the same partition-1 carrier arc),
# junction seam chains and Point15 junction entities, ghost layers and the
# MSH 4.1 `$PartitionedEntities` writer that `.msh` Save upgrades into.
using Test, Tessella
const E = Tessella.Elements
const API = Tessella.API

const BOX = """
SetFactory("OpenCASCADE");
Box(1) = {0,0,0, 1,1,1};
Mesh 3;
"""

function _pm_execute(extra;dir=missing)
    src = joinpath(dir===missing ? tempdir() : dir,"pm_case.geo")
    write(src,BOX*extra)
    return Tessella.GeoExec.execute_geo(src)
end

function _pm_cell_partitions(mesh)
    data = mesh.entity_data
    out = Dict{UInt64,Int32}()
    for (block,tags,owners) in zip(mesh.blocks,data.external_element_tags,
                                   data.block_entities),
        (uid,owner) in zip(tags,owners)
        block.msh==4 || continue
        entity = mesh.partition_data.entities[(3,Int(owner))]
        out[uid] = only(entity.partitions)
    end
    return out
end

_counts(mesh) = Dict(Int(block.msh)=>count(b->b.msh==block.msh,mesh.blocks) for
    block in mesh.blocks)
_ncells(mesh,msh) = sum(size(b.nodes,2) for b in mesh.blocks if b.msh==msh)

@testset "GEO PartitionMesh partitions an OCC box mesh" begin
    mktempdir() do directory
        out = joinpath(directory,"out.msh")
        ex = _pm_execute("PartitionMesh(3);\nSave \"$(replace(out,"\\"=>"/"))\";\n")
        mesh = ex.mesh
        @test mesh isa MixedMesh
        @test mesh.partition_data!==nothing
        @test mesh.partition_data.num_partitions == 3
        dim_counts = Dict(d=>count(k->k[1]==d,keys(mesh.partition_data.entities)) for d in 0:3)
        @test dim_counts[3]==3                       # one Region child per partition
        @test dim_counts[2]>0 && dim_counts[1]>0    # surface + interface/seam children
        # The assignment matches api.mesh.partition(3) semantics exactly —
        # compare against a fresh unpartitioned reference run.
        base = _pm_execute("")
        _,ref,_ = API._higher_frontend_input(base.model,base.mesh,nothing,
            nothing,nothing,"check")
        assignment = API._partition_lower_inherit(API._partition_elements(ref,3),ref)
        geo = _pm_cell_partitions(mesh)
        let offset = 0, refparts = Dict{UInt64,Int32}()
            for (bi,block) in enumerate(ref.blocks)
                tags = ref.entity_data.external_element_tags[bi]
                if block.msh==4
                    for (i,uid) in enumerate(tags)
                        refparts[uid] = assignment[offset+i]
                    end
                end
                offset += size(block.nodes,2)
            end
            @test sort(collect(keys(geo))) == sort(collect(keys(refparts)))
            @test all(uid->geo[uid]==refparts[uid],keys(refparts))
        end
        @test all(k->count(==(k),values(geo))>0,Int32(1):Int32(3))
        # Saved file carries MSH 4.1 partition metadata even through `.msh`.
        txt = read(out,String)
        @test occursin("\$PartitionedEntities",txt)
        @test occursin("4.1",split(txt,'\n')[2])
        @test !occursin("\$GhostElements",txt)
        pe = split(txt[first(findfirst("\$PartitionedEntities",txt)):end],'\n')
        @test parse(Int,pe[2]) == 3                # numPartitions
    end
end

@testset "GEO PartitionMesh option combinations" begin
    seen = Dict{Tuple{Int,Int},NamedTuple}()
    for (topo,ghost) in ((0,0),(0,1),(1,0),(1,1))
        ex = _pm_execute("Mesh.PartitionCreateTopology=$topo;\n" *
                         "Mesh.PartitionCreateGhostCells=$ghost;\n" *
                         "PartitionMesh(3);\n")
        mesh = ex.mesh
        @test mesh.partition_data.num_partitions == 3
        n_interface = count(keys(mesh.partition_data.entities)) do key
            key[1]==2 && length(mesh.partition_data.entities[key].partitions)>1
        end
        nghost = length(mesh.partition_data.ghost_elements)
        seen[(topo,ghost)] = (interfaces=n_interface,ghosts=nghost,
            tets=_ncells(mesh,4),seams=_ncells(mesh,1))
        @test (n_interface>0) == (topo==1)
        @test (nghost>0) == (ghost==1)
        @test _ncells(mesh,4) > 0
    end
    @test seen[(1,0)].seams > seen[(0,0)].seams    # topology emits seam cells
    @test seen[(1,0)].tets == seen[(0,0)].tets     # topology adds no cells
    @test seen[(0,1)].tets == seen[(0,0)].tets     # ghosts stay aliases
end

@testset "GEO PartitionMesh rejects invalid use" begin
    @test_throws ArgumentError _pm_execute("PartitionMesh(0);\n")
    nomesh = joinpath(tempdir(),"pm_nomesh.geo")
    write(nomesh,BOX*"Delete Meshes;\nPartitionMesh(3);\n")
    @test_throws ArgumentError Tessella.GeoExec.execute_geo(nomesh)
end
end
