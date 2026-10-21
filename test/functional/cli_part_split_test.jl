module CLIPartSplitFunctionalChecks
# Functional checks for `-part_split` and its aliases `-part_split_files` /
# `-oneFilePerPart` (upstream `Mesh.PartitionSplitMeshFiles`): one
# `<stem>_<p>.msh` per partition carrying that partition's own cells plus its
# ghost copies as regular cells on its ghost entity. Before the flag was
# implemented these flags hit the _REJECTED_FLAGS table and threw
# `ArgumentError("per-partition file output is not implemented")`.
using Test, Tessella
const E=Tessella.Elements

const SQUARE="""
Point(1)={0,0,0,0.4};Point(2)={1,0,0,0.4};
Point(3)={1,1,0,0.4};Point(4)={0,1,0,0.4};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
"""
# Two Tet4 sharing a face on volume entity 1 — the smallest ghost-producing
# partition: each side owns one tet and ghosts the other.
const TWO_TETS="""
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Entities
0 0 0 1
1 0 0 0 1 1 0 0 0
\$EndEntities
\$Nodes
1 5 1 5
3 1 0 5
1
2
3
4
5
0 0 0
1 0 0
0 1 0
0 0 1
1 1 1
\$EndNodes
\$Elements
1 2 1 2
3 1 4 2
1 1 2 3 4
2 2 3 4 5
\$EndElements
"""

_block_ncells(block::E.ElementBlock)=size(block.nodes,2)
_block_ncells(block::E.SpecialElementBlock)=length(block.tags)
_block_nodes(block::E.ElementBlock,j)=block.nodes[:,j]
_block_nodes(block::E.SpecialElementBlock,j)=
    view(block.connectivity,Int(block.offsets[j]):Int(block.offsets[j+1])-1)
_block_dim(block::E.ElementBlock)=E.msh_dimension(block.msh)
_block_dim(block::E.SpecialElementBlock)=E.MSH_SPECIAL_TYPES[block.msh].dim
total_cells(m)=sum(_block_ncells,m.blocks)

# Cell identity across files: the sorted public node tags — independent of
# each file's node renumbering.
function cell_key(m,bi,j)
    tags=m.entity_data.external_node_tags
    return sort!([tags[Int(node)] for node in _block_nodes(m.blocks[bi],j)])
end

# Split a partition file's cells into own cells (entities whose `partitions`
# list holds `partition`) and ghost copies (cells on `partition`'s ghost
# entity, the `(tag, partition)` records of `$PartitionedEntities`).
function partition_sets(m,partition)
    data=m.entity_data;pdata=m.partition_data
    ghost_tags=Set(Int(tag) for (tag,part) in pdata.ghost_entities
                   if Int(part)==partition)
    own=Set{Vector{UInt64}}();ghost=Set{Vector{UInt64}}()
    for (bi,block) in enumerate(m.blocks)
        entities=data.block_entities[bi];dim=_block_dim(block)
        for j in 1:_block_ncells(block)
            key=cell_key(m,bi,j);tag=Int(entities[j])
            record=get(pdata.entities,(dim,tag),nothing)
            if record!==nothing && Int32(partition) in record.partitions
                push!(own,key)
            elseif tag in ghost_tags
                push!(ghost,key)
            else
                error("split file cell on unexpected entity ($dim,$tag)")
            end
        end
    end
    return own,ghost
end

cell_key_set(m)=Set(cell_key(m,bi,j) for (bi,block) in enumerate(m.blocks)
                    for j in 1:_block_ncells(block))

@testset "CLI -part_split/-oneFilePerPart/-part_split_files" begin
    mktempdir() do directory
        two_tets=joinpath(directory,"two_tets.msh");write(two_tets,TWO_TETS)
        original=E.read_mixed_msh(two_tets)
        original_set=Set(cell_key(original,bi,j)
                         for (bi,block) in enumerate(original.blocks)
                         for j in 1:_block_ncells(block))
        @test length(original_set)==2

        # Ghost-producing split: each file holds its own tet plus a regular
        # cell copy of the neighbour's tet on its ghost entity.
        out=joinpath(directory,"split.msh")
        @test Tessella.CLI.main([two_tets,"-part","2","-part_ghosts",
            "-part_split","-o",out])==out
        @test !isfile(out)
        files=[joinpath(directory,"split_$p.msh") for p in 1:2]
        @test all(isfile,files)
        views=[E.read_mixed_msh(f) for f in files]
        sets=[partition_sets(v,p) for (p,v) in enumerate(views)]
        for (p,v) in enumerate(views)
            @test Tessella.validate(v).ok
            @test v.partition_data!==nothing
            @test Int(v.partition_data.num_partitions)==2
            own,ghost=sets[p]
            @test length(own)==1 && length(ghost)==1
            @test only(ghost)==only(sets[3-p][1])   # ghost = the other's own tet
            # One record per ghost copy: owner partition, single-destination
            # list naming this file's partition, index into this file's cells.
            records=v.partition_data.ghost_elements
            @test length(records)==1
            record=only(records)
            @test Int(record.partition)==3-p
            @test record.ghost_partitions==Int32[p]
            data=v.entity_data;index=0;ghost_entity_hit=false
            for (bi,block) in enumerate(v.blocks),j in 1:_block_ncells(block)
                index+=1
                index==Int(record.element) || continue
                tag=Int(data.block_entities[bi][j])
                ghost_entity_hit=any(
                    t==tag && Int(part)==p for (t,part) in
                    v.partition_data.ghost_entities)
            end
            @test ghost_entity_hit
            # The partition's ghost entity is declared in $Entities — the
            # Tessella convention that keeps split files self-consistent.
            gtag=only(Int(tag) for (tag,part) in v.partition_data.ghost_entities
                      if Int(part)==p)
            @test haskey(data.entities,(3,gtag))
        end
        @test sets[1][1] ∪ sets[2][1] == original_set   # own cells partition the input

        # Both aliases produce identical split files.
        for (alias,stem) in (("-oneFilePerPart","al1"),("-part_split_files","al2"))
            Tessella.CLI.main([two_tets,"-part","2","-part_ghosts",
                alias,"-o",joinpath(directory,"$stem.msh")])
            for p in 1:2
                a=joinpath(directory,"$(stem)_$p.msh")
                @test isfile(a)
                @test cell_key_set(E.read_mixed_msh(a))==
                      cell_key_set(views[p])
            end
        end

        # Label-mode `-part` (no topology/ghosts): each file holds the cells
        # its partition owns on the original entity — the union reproduces
        # the input exactly, with no partition metadata left over.
        lab=joinpath(directory,"lab.msh")
        Tessella.CLI.main([two_tets,"-part","2","-part_split","-o",lab])
        @test !isfile(lab)
        label_views=[E.read_mixed_msh(joinpath(directory,"lab_$p.msh"))
                     for p in 1:2]
        @test all(v->Tessella.validate(v).ok,label_views)
        @test all(v->v.partition_data===nothing,label_views)
        label_union=union((Set(cell_key(v,bi,j)
            for (bi,block) in enumerate(v.blocks)
            for j in 1:_block_ncells(block)) for v in label_views)...)
        @test label_union==original_set
        @test sum(total_cells,label_views)==2

        # .geo input without -o derives the stem from the input file. The
        # mesher re-rolls between runs, so there is no cross-run baseline —
        # the cross-file contract is disjoint coverage of one fresh mesh.
        square_geo=joinpath(directory,"square.geo");write(square_geo,SQUARE)
        Tessella.CLI.main([square_geo,"-2","-part","3","-part_split"])
        @test !isfile(joinpath(directory,"square.msh"))
        square_files=[joinpath(directory,"square_$p.msh") for p in 1:3]
        @test all(isfile,square_files)
        square_views=[E.read_mixed_msh(f) for f in square_files]
        @test all(v->Tessella.validate(v).ok,square_views)
        square_sets=[cell_key_set(v) for v in square_views]
        @test all(!isempty,square_sets)
        @test isempty(square_sets[1] ∩ square_sets[2])
        @test isempty(square_sets[1] ∩ square_sets[3])
        @test isempty(square_sets[2] ∩ square_sets[3])
        @test length(union(square_sets...))==sum(total_cells,square_views)

        # `-part_split` without `-part` is inert: the normal single output.
        plain=joinpath(directory,"plain.msh")
        @test Tessella.CLI.main([two_tets,"-save","-part_split","-o",plain])==plain
        @test isfile(plain)
        @test !isfile(joinpath(directory,"plain_1.msh"))
        @test !isfile(joinpath(directory,"plain_2.msh"))
        @test Tessella.validate(E.read_mixed_msh(plain)).ok

        # Flag contract: split flags conflict with each other like the other
        # partition toggles; -part alone still writes the combined file.
        @test_throws ArgumentError Tessella.CLI.main(
            [two_tets,"-part","2","-part_split","-oneFilePerPart","-o",lab])
        @test_throws ArgumentError Tessella.CLI.main(
            [two_tets,"-part","2","-part_split_files","-part_split","-o",lab])
        combined=joinpath(directory,"combined.msh")
        Tessella.CLI.main([two_tets,"-part","2","-o",combined])
        @test isfile(combined)
        @test !isfile(joinpath(directory,"combined_1.msh"))

        # A single partition degenerates to the normal write (upstream keeps
        # the unsuffixed name when there is nothing to split).
        single=joinpath(directory,"single.msh")
        Tessella.CLI.main([two_tets,"-part","1","-part_ghosts","-part_split",
            "-o",single])
        @test isfile(single)
        @test !isfile(joinpath(directory,"single_1.msh"))
    end
end
end
