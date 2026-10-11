# A periodic declaration records node identities at an operation boundary.
# A retired MVertex tag can later coincide with another live node's tag;
# these records therefore keep liveness separately from their literal UIDs.
struct _ModelPeriodicUIDPair
    slave::UInt64
    master::UInt64
    slave_live::Bool
    master_live::Bool
    function _ModelPeriodicUIDPair(slave::UInt64,master::UInt64,
                                   slave_live::Bool,master_live::Bool)
        return new(slave,master,slave_live,master_live)
    end
end
Base.:(==)(a::_ModelPeriodicUIDPair,b::_ModelPeriodicUIDPair)=
    all(name->getfield(a,name)==getfield(b,name),fieldnames(_ModelPeriodicUIDPair))

struct _ModelPeriodicUIDHistory
    master::Int32
    affine::Union{Nothing,NTuple{16,Float64}}
    primary::Vector{_ModelPeriodicUIDPair}
    higher::Vector{_ModelPeriodicUIDPair}
    function _ModelPeriodicUIDHistory(master::Int32,
            affine::Union{Nothing,NTuple{16,Float64}},
            primary::Vector{_ModelPeriodicUIDPair},higher::Vector{_ModelPeriodicUIDPair})
        master>0 || throw(ArgumentError("periodic UID history requires a positive master"))
        affine===nothing || all(isfinite,affine) || throw(ArgumentError(
            "periodic UID history requires a finite affine transform"))
        return new(master,affine,copy(primary),copy(higher))
    end
end
Base.:(==)(a::_ModelPeriodicUIDHistory,b::_ModelPeriodicUIDHistory)=
    all(name->getfield(a,name)==getfield(b,name),fieldnames(_ModelPeriodicUIDHistory))

function _model_periodic_uid_query(history::_ModelPeriodicUIDHistory,
                                    include_high_order::Bool)
    pairs=include_high_order ? vcat(history.primary,history.higher) : history.primary
    return (master_entity=history.master,
        slave_nodes=UInt64[pair.slave for pair in pairs],
        master_nodes=UInt64[pair.master for pair in pairs],
        affine=history.affine===nothing ? Float64[] : collect(history.affine))
end

function _model_periodic_uid_remap(pair::_ModelPeriodicUIDPair,
                                   mapping::AbstractDict{UInt64,UInt64})
    function side(uid,live)
        live || return uid,false
        # Absence from the certified operation map retires this identity.
        # It does not revive when a later node reuses the same literal UID.
        return haskey(mapping,uid) ? (mapping[uid],true) : (uid,false)
    end
    slave,slave_live=side(pair.slave,pair.slave_live)
    master,master_live=side(pair.master,pair.master_live)
    return _ModelPeriodicUIDPair(slave,master,slave_live,master_live)
end

function _model_periodic_uid_remap(history::_ModelPeriodicUIDHistory,
                                   mapping::AbstractDict{UInt64,UInt64})
    return _ModelPeriodicUIDHistory(history.master,history.affine,
        [_model_periodic_uid_remap(pair,mapping) for pair in history.primary],
        [_model_periodic_uid_remap(pair,mapping) for pair in history.higher])
end

function _model_periodic_uid_set_order(history::_ModelPeriodicUIDHistory,
        mapping::AbstractDict{UInt64,UInt64},
        higher::Vector{_ModelPeriodicUIDPair})
    primary=[_model_periodic_uid_remap(pair,mapping) for pair in history.primary]
    return _ModelPeriodicUIDHistory(history.master,history.affine,primary,higher)
end

function _model_periodic_uid_refine(history::_ModelPeriodicUIDHistory,
        mapping::AbstractDict{UInt64,UInt64},
        children::Vector{_ModelPeriodicUIDPair})
    primary=[_model_periodic_uid_remap(pair,mapping)
             for pair in vcat(history.primary,history.higher)]
    # Only live identities may coincide with current child rows. Retired
    # records retain their separate event even if the absolute UIDs collide.
    for pair in children
        pair.slave_live && pair.master_live || throw(ArgumentError(
            "periodic refinement child receipt must bind live node identities"))
        any(old->old.slave_live && old.master_live &&
            old.slave==pair.slave && old.master==pair.master,primary) && continue
        push!(primary,pair)
    end
    return _ModelPeriodicUIDHistory(history.master,history.affine,
        primary,_ModelPeriodicUIDPair[])
end

function _model_periodic_uid_retag(histories,dimension::Int,old_tag::Int,new_tag::Int)
    replacement=Dict{Tuple{Int,Int},_ModelPeriodicUIDHistory}()
    for ((dim,slave),history) in histories
        moved_slave=dim==dimension && slave==old_tag ? new_tag : slave
        moved_master=dim==dimension && history.master==old_tag ? Int32(new_tag) : history.master
        key=(dim,moved_slave)
        haskey(replacement,key) && throw(ArgumentError(
            "periodic UID history retagging would merge distinct declarations"))
        replacement[key]=_ModelPeriodicUIDHistory(moved_master,history.affine,
            history.primary,history.higher)
    end
    return replacement
end

function _model_periodic_uid_clear!(histories,selection)
    for key in selection
        history=get(histories,key,nothing)
        history===nothing && continue
        histories[key]=_ModelPeriodicUIDHistory(history.master,history.affine,
            _ModelPeriodicUIDPair[],_ModelPeriodicUIDPair[])
    end
    return nothing
end

function _model_periodic_uid_retire(histories,removed_uids::AbstractSet{UInt64},
                                   removed_entities=Set{Tuple{Int,Int}}())
    replacement=Dict{Tuple{Int,Int},_ModelPeriodicUIDHistory}()
    for (key,history) in histories
        key in removed_entities && continue
        retired(pair)=_ModelPeriodicUIDPair(pair.slave,pair.master,
            pair.slave_live && !(pair.slave in removed_uids),
            pair.master_live && !(pair.master in removed_uids))
        replacement[key]=_ModelPeriodicUIDHistory(history.master,history.affine,
            retired.(history.primary),retired.(history.higher))
    end
    return replacement
end

function _model_periodic_uid_remap!(m,mapping::AbstractDict{UInt64,UInt64})
    replacement=Dict(key=>_model_periodic_uid_remap(history,mapping)
        for (key,history) in m.meshing.periodic_uid_history)
    m.meshing.periodic_uid_history=replacement
    return nothing
end

function _model_periodic_surface_reset!(m,slaves::Vector{Int},caller::AbstractString)
    all(slave->_model_entity_known(m,2,slave),slaves) || throw(ArgumentError(
        "$caller: unknown Surface in periodic self relation"))
    periodic=copy(m.periodic)
    histories=copy(m.meshing.periodic_uid_history)
    for slave in slaves
        delete!(periodic,(2,slave))
        for (dimension,tag) in model_boundary(m,[(2,slave)],false,false,false)
            dimension==1 || throw(ArgumentError(
                "$caller: periodic Surface has an invalid literal Curve boundary"))
            delete!(periodic,(1,abs(Int(tag))))
        end
        histories[(2,slave)]=_ModelPeriodicUIDHistory(Int32(slave),nothing,
            _ModelPeriodicUIDPair[],_ModelPeriodicUIDPair[])
    end
    m.periodic=periodic
    m.meshing.periodic_uid_history=histories
    return nothing
end

function _model_periodic_surface_native_declaration(m,slave_entities,master_entities,
        affine,atol,caller::AbstractString)
    slaves=_periodic_entity_tags(slave_entities,2,caller,"slave")
    masters=_periodic_entity_tags(master_entities,2,caller,"master")
    length(slaves)==length(masters) && !isempty(slaves) || throw(ArgumentError(
        "$caller: Surface slave/master counts must agree and be nonempty"))
    length(unique(abs.(slaves)))==length(slaves) || throw(ArgumentError(
        "$caller: repeated Surface slave identity"))
    _periodic_affine_input(affine,caller;name="affine transform")
    _model_periodic_tolerance(atol,caller)
    staged=_model_mesh_transaction_copy(m)
    for (slave,master) in zip(slaves,masters)
        if abs(slave)==abs(master)
            _model_periodic_surface_reset!(staged,[abs(slave)],caller)
        else
            set_periodic!(staged,2,[slave],[master],affine;atol,overwrite=true)
        end
    end
    return staged
end

# This is the explicit declared-affine correspondence operation. It records
# pair events and never merges nodes, changes ownership, or establishes CAD
# identity from coordinates. The caller supplies the actual owned node rows
# of the declared Surface closure after its CAD relation has been certified.
struct _ModelPeriodicUIDVertex
    uid::UInt64
    xyz::NTuple{3,Float64}
    owner::Tuple{Int,Int}
    function _ModelPeriodicUIDVertex(uid::UInt64,xyz::NTuple{3,Float64},
                                     owner::Tuple{Int,Int})
        all(isfinite,xyz) && 0<=owner[1]<=2 && owner[2]>=0 ||
            throw(ArgumentError("periodic correspondence requires a finite owned node row"))
        return new(uid,xyz,owner)
    end
end

function _model_periodic_uid_inverse(affine::NTuple{16,Float64})
    all(isfinite,affine) && affine[13:16]==(0.,0.,0.,1.) || throw(ArgumentError(
        "periodic correspondence requires a finite affine transform"))
    a,b,c,d,e,f,g,h,i=affine[1],affine[2],affine[3],affine[5],affine[6],
        affine[7],affine[9],affine[10],affine[11]
    cof=(e*i-f*h,c*h-b*i,b*f-c*e,f*g-d*i,a*i-c*g,c*d-a*f,
         d*h-e*g,b*g-a*h,a*e-b*d)
    determinant=a*cof[1]+b*cof[4]+c*cof[7]
    isfinite(determinant) && determinant!=0. || throw(ArgumentError(
        "periodic correspondence requires a nonsingular affine transform"))
    inverse=ntuple(j->cof[j]/determinant,9)
    all(isfinite,inverse) || throw(ArgumentError(
        "periodic correspondence inverse is not representable"))
    return inverse
end

function _model_periodic_uid_correspondence(
        slaves::Vector{_ModelPeriodicUIDVertex},masters::Vector{_ModelPeriodicUIDVertex},
        affine::NTuple{16,Float64},caller::AbstractString)
    length(unique(vertex.uid for vertex in slaves))==length(slaves) &&
        length(unique(vertex.uid for vertex in masters))==length(masters) ||
        throw(ArgumentError("$caller: periodic closure has repeated live node identities"))
    inverse=_model_periodic_uid_inverse(affine)
    pairs=_ModelPeriodicUIDPair[]
    isempty(masters) && return pairs
    for slave in slaves
        shifted=(slave.xyz[1]-affine[4],slave.xyz[2]-affine[8],slave.xyz[3]-affine[12])
        point=ntuple(axis->inverse[3axis-2]*shifted[1]+
            inverse[3axis-1]*shifted[2]+inverse[3axis]*shifted[3],3)
        all(isfinite,point) || throw(ArgumentError(
            "$caller: inverse-affine periodic node is not representable"))
        closest=0;minimum_distance=Inf;ambiguous=false
        for (index,master) in enumerate(masters)
            delta=ntuple(axis->point[axis]-master.xyz[axis],3)
            distance=delta[1]*delta[1]+delta[2]*delta[2]+delta[3]*delta[3]
            isfinite(distance) || throw(ArgumentError(
                "$caller: periodic correspondence distance is not representable"))
            if distance<minimum_distance
                closest=index;minimum_distance=distance;ambiguous=false
            elseif distance==minimum_distance
                ambiguous=true
            end
        end
        ambiguous && throw(ArgumentError(
            "$caller: tied master-node correspondence is not certified"))
        closest>0 || throw(ArgumentError("$caller: no periodic master-node correspondence"))
        push!(pairs,_ModelPeriodicUIDPair(slave.uid,masters[closest].uid,true,true))
    end
    return pairs
end

function _model_periodic_uid_owned_vertices(m,surface::Int,caller::AbstractString)
    _model_entity_known(m,2,surface) || throw(ArgumentError(
        "$caller: unknown periodic Surface[$surface]"))
    owners=Set{Tuple{Int,Int}}([(2,surface)])
    pending=Tuple{Int,Int}[(2,surface)]
    while !isempty(pending)
        entity=pop!(pending)
        entity[1]==0 && continue
        for (dimension,tag) in model_boundary(m,[entity],false,false,false)
            child=(Int(dimension),abs(Int(tag)))
            child in owners && continue
            push!(owners,child);push!(pending,child)
        end
    end
    vertices=_ModelPeriodicUIDVertex[]
    identities=Dict{UInt64,_ModelPeriodicUIDVertex}()
    for (owner,record) in sort!(collect(_discrete_mesh_records_model(m));by=first)
        owner in owners || continue
        size(record.node_coords)==(3,length(record.node_tags)) || throw(ArgumentError(
            "$caller: periodic owned record has an invalid coordinate table"))
        for (row,uid) in enumerate(record.node_tags)
            vertex=_ModelPeriodicUIDVertex(UInt64(uid),
                ntuple(axis->Float64(record.node_coords[axis,row]),3),owner)
            previous=get(identities,vertex.uid,nothing)
            previous===nothing || throw(ArgumentError(
                "$caller: periodic closure repeats owned node UID $(vertex.uid)"))
            identities[vertex.uid]=vertex;push!(vertices,vertex)
        end
    end
    return vertices
end

function _model_periodic_uid_declaration_history(m,constraint,caller::AbstractString)
    constraint.dim==2 && constraint.affine!==nothing || throw(ArgumentError(
        "$caller: Surface UID declaration requires its certified affine CAD relation"))
    slave=Int(constraint.slave_entity);master=Int(constraint.master_entity)
    slaves=_model_periodic_uid_owned_vertices(m,slave,caller)
    masters=_model_periodic_uid_owned_vertices(m,master,caller)
    pairs=_model_periodic_uid_correspondence(slaves,masters,constraint.affine,caller)
    return _ModelPeriodicUIDHistory(Int32(master),constraint.affine,pairs,_ModelPeriodicUIDPair[])
end
