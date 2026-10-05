# Accepted physical input witnesses, retained from the independent 30,264
# product authority. Every rank is realized by reordering actual source IDs;
# every processing order is realized by reordering actual source cells.
const RETAINED_WITNESSES=(
    (name="two_centers",m=7,n=1,
     cells=((1,2,10,9),(10,2,3,11),(11,3,4,12),(12,4,5,13),(13,5,6,14),(7,15,14,6),(8,16,15,7)),
     permutation=(6,14,12,13,5,11,1,3,7,4,15,16,9,2,10,8),
     order=(4,3,6,5,1,7,2),directions=((1,2),(1,9),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,16),(10,9),(11,10),(12,11),(13,12),(14,13),(15,14),(16,15)),
     problems=((2,1),(5,1)),counts=(40,0,0,7)),
    (name="nonterminal_center",m=9,n=2,
     cells=((11,1,2,12),(12,2,3,13),(4,14,13,3),(14,4,5,15),(15,5,6,16),(7,17,16,6),(7,8,18,17),(9,19,18,8),(10,20,19,9)),
     permutation=(7,8,6,15,10,20,3,13,4,17,5,16,18,9,19,11,2,12,14,1),
     order=(9,8,2,1,7,5,3,6,4),directions=((2,1),(1,11),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,9),(10,20),(11,12),(12,13),(13,14),(14,15),(15,16),(16,17),(17,18),(18,19),(19,20)),
     problems=((4,1),),counts=(80,0,0,17)))

function retained_product(witness;direction=1)
    M=Tessella.Model;Free=M._ExtrudeNoNewB4Free;m=witness.m;n=witness.n
    directions=Dict(minmax(pair...)=>(pair[1],pair[2]) for pair in witness.directions)
    chains=(collect(1:m+1),[m+1,2m+2],collect(2m+2:-1:m+2),[m+2,1])
    reversals=0
    for curve in 1:4
        same=[directions[minmax(chains[curve][k],chains[curve][k+1])]==(chains[curve][k],chains[curve][k+1]) for k in 1:length(chains[curve])-1]
        require(all(same) || all(!,same),"witness generatrix direction changes within one CAD curve")
        all(same) || (reversals|=1<<(curve-1))
    end
    f=fixture("$(witness.name)_$(direction)";strip_length=m,layers=:one,
        height=Float64(direction),pins=(1,2,3,4),curve_reverse=reversals,tags=:sparse)
    f=merge(f,(source=replace(f.source,"Layers{1}"=>"Layers{$n}"),
        levels=collect(range(0.,1.;length=n+1)),intervals=n))
    geometry=execute(f.source;dim=0);model=geometry.model
    original=M.mesh_model_surface(model,f.surface);independent=source_complex(original,f)
    canonical=vcat(collect(independent.chains[1]),reverse(collect(independent.chains[3])))
    permutation=collect(witness.permutation);inverse=invperm(permutation)
    coordinates=original.coords[:,canonical[permutation]]
    cells=hcat([Int32[inverse[node] for node in witness.cells[cell]] for cell in witness.order]...)
    source=MixedMesh(coordinates,[ElementBlock(3,cells)])
    source_complex(source,f)
    out=Int.(geometry.lists["sweep"]);tag=out[2];caller="actual retained free B4 witness"
    spec=model.meshing.extrude_specs[(3,tag)]
    data=M._extrude_nonew_rect_grid_source(model,f.surface,source,spec,caller;mode=:b4_strip)
    refs=fill((Int32(1),Int32(1)),n)
    catalog=Free.plan(data,f.levels,refs,caller)
    positions=[(Int32(row),Int32(interval)) for (row,cell) in enumerate(witness.order)
        for interval in 1:n if (cell,interval) in witness.problems]
    require(catalog.problem_positions==positions && catalog.cell_counts==witness.counts,"actual source no longer realizes its accepted witness")
    params=M._extrude_gate(model,3,tag)
    cols=M._extrude_volume_columns(model,tag,f.surface,source,params,spec,f.levels,caller)
    plan=Free.finish(model,tag,params,f.surface,source,catalog,cols,out[1],Tuple(out[3:end]),caller)
    volume=M._extrude_volume_part(model,plan.sweep,caller)
    # Build an explicit operation scope from this actual completed source,
    # using the production top/lateral adapters and projection unchanged.
    surfaces=Dict{Int,Union{Tessella.Mesh,MixedMesh}}(f.surface=>source)
    top_params,_,_=M._extrude_entity_params(model,2,out[1],caller)
    surfaces[out[1]]=M._extrude_nonew_rect_grid_top(model,out[1],top_params,cols,catalog,caller)
    for lateral in out[3:end]
        lateral_params,_,link=M._extrude_entity_params(model,2,lateral,caller)
        surfaces[lateral]=M._extrude_nonew_rect_grid_lateral(model,lateral,lateral_params,link[2],cols,catalog,plan.edges,caller)
    end
    M._extrude_nonew_certify_boundary(volume,surfaces,caller;oriented_internal=true,face_capacity=catalog.face_capacity)
    completed=M._ExtrudeNoNewCompletePlan(volume,surfaces,cols,catalog)
    scope=M._ExtrudeNoNewScope(model,Dict(tag=>completed),surfaces)
    projected=M.model_to_mixed(model,volume,3,tag;_extrude_scope=scope)
    certificate=certify(volume,f,source)
    require(Set(values(certificate.center_parent))==Set(Tuple(Int.(p)) for p in positions),"actual center owners differ from retained source witnesses")
    return (;f,source,volume,projected,geometry,scope,catalog,certificate,
        expected_centers=length(positions),nonterminal=any(last(p)<n for p in positions),
        expected_positions=positions)
end

retained_products()=[retained_product(witness;direction) for witness in RETAINED_WITNESSES for direction in (-1,1)]
