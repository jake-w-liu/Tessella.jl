module ComputeRenumberingChecks
# Functional checks for `API.mesh.compute_renumbering` covering the native
# method tiers: "RCMK" (reverse Cuthill–McKee), "Hilbert" (space-filling
# order + deterministic 2-opt polish) and "Metis" (native nested
# dissection — the upstream METIS_NodeND tier without a METIS dependency).
using Test, Tessella, Random
const A=Tessella.API

# Structured N×N triangle grid with node tags assigned in a deterministically
# shuffled order so identity numbering has a wide front.
function _scrambled_tri_grid!(N::Int)
    A.model.add_discrete_entity(2,1)
    n=N+1
    # Deterministic permutation: spatial neighbours carry far-apart tags.
    order=sortperm(rand(MersenneTwister(1234),n*n))
    coords=Dict{Int,Vector{Float64}}()
    for j in 1:n,i in 1:n
        coords[(j-1)*n+i]=[i-1.0,j-1.0,0.0]
    end
    tags=collect(1:length(order))
    flat=Float64[]
    for (new,old) in enumerate(order)
        append!(flat,coords[old])
    end
    # nodes[old_position]=new tag position lookup
    A.mesh.add_nodes(2,1,tags,flat)
    tris=Int[];conn=Int32[]
    # connectivity must reference the scrambled node tags
    newpos=Dict{Int,Int}(old=>new for (new,old) in enumerate(order))
    for j in 1:N,i in 1:N
        a=(j-1)*n+i;b=a+1;c=a+n;d=c+1
        push!(tris,length(tris)+1)
        append!(conn,Int32[newpos[a],newpos[b],newpos[d]])
        push!(tris,length(tris)+1)
        append!(conn,Int32[newpos[a],newpos[d],newpos[c]])
    end
    A.mesh.add_elements_by_type(1,2,tris,conn)
    A.mesh.generate(2)
    return nothing
end

# Bandwidth of the current mesh under a node permutation: `newpos` maps
# old public tag -> new position (as returned by compute_renumbering).
# `get_elements` returns (types, element_tags, flat_node_tags) blocks.
function _bandwidth(newpos)
    mx=0
    families,_,cells=A.mesh.get_elements(2,1)
    for (f,block) in enumerate(cells)
        width=families[f]==Int32(2) ? 3 : 4
        for c in 1:width:length(block)
            vs=[newpos[Int(block[c+r])] for r in 0:width-1]
            mx=max(mx,maximum(vs)-minimum(vs))
        end
    end
    return mx
end

@testset "compute_renumbering methods" begin
    A.finalize();A.initialize();A.model.add("renumber")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    try
        _scrambled_tri_grid!(8)    # 9x9 nodes, 128 tris, scrambled tags
        tags,_,_=A.mesh.get_nodes()
        identity=Dict{Int,Int}(Int(t)=>i for (i,t) in enumerate(sort(Int.(tags))))
        base=_bandwidth(identity)
        for method in ("RCMK","Hilbert","Metis")
            old,new=A.mesh.compute_renumbering(method)
            @test length(old)==length(new)==81
            # `new` is a permutation of 1..n over the involved tags.
            @test sort(collect(new))==collect(1:81)
            @test sort(Int.(old))==sort(Int.(tags))
            pos=Dict{Int,Int}(Int(t)=>Int(nt)
                              for (t,nt) in zip(old,new))
            bw=_bandwidth(pos)
            @test bw<base
            # Deterministic.
            @test A.mesh.compute_renumbering(method)==(old,new)
        end
        # Restricted element subset: only nodes of selected elements move.
        families,uids,cells=A.mesh.get_elements(2,1)
        tri_i=findfirst(==(Int32(2)),families)
        subset=Int.(uids[tri_i][1:16])
        old,new=A.mesh.compute_renumbering("RCMK",subset)
        @test length(old)<81
        @test sort(collect(new))==collect(1:length(old))
        @test_throws ArgumentError A.mesh.compute_renumbering("Bogus")
        @test_throws ArgumentError A.mesh.compute_renumbering("RCMK",[9999])
    finally
        A.finalize()
    end
end

@testset "renumber_nodes applies a computed permutation" begin
    A.finalize();A.initialize();A.model.add("renumber apply")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    try
        _scrambled_tri_grid!(6)
        old,new=A.mesh.compute_renumbering("RCMK")
        A.mesh.renumber_nodes(old,new)
        tags,_,_=A.mesh.get_nodes()
        @test sort(Int.(tags))==collect(Int.(old))
        # Element connectivity still references existing node tags.
        _,_,cells=A.mesh.get_elements(2,1)
        node_set=Set(Int.(tags))
        for block in cells
            @test all(t->Int(t) in node_set,block)
        end
    finally
        A.finalize()
    end
end
end
