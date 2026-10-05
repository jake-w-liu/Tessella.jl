# Free all-boundary source strips: complete physical face propagation,
# deferred literal factory dispatch, and retained actual-center emission.
module _ExtrudeNoNewB4Free
const Model = parentmodule(@__MODULE__)

module SourcePhase
const Model = getfield(parentmodule(@__MODULE__), :Model)
module LocalPhase
const FACES=((1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8),
             (1,2,3,4),(5,6,7,8))
const PAIRS=((1,3),(2,4),(5,6))
@inline bit(face::Int)=UInt8(1)<<(face-1)
@inline has(mask::UInt8,face::Int)=!iszero(mask&bit(face))
@inline function through(face::Int,vertex::Int)
    cycle=FACES[face]
    return vertex==cycle[1] || vertex==cycle[3] ? UInt8(1) : UInt8(2)
end
@inline function contains(face::Int,state::UInt8,vertex::Int)
    state==0 && return false
    cycle=FACES[face]
    return state==1 ? vertex==cycle[1] || vertex==cycle[3] :
                      vertex==cycle[2] || vertex==cycle[4]
end
@inline function lowest(face::Int,ranks::NTuple{8,Int})
    cycle=FACES[face];winner=cycle[1]
    for index in 2:4
        vertex=cycle[index]
        ranks[vertex]<ranks[winner] && (winner=vertex)
    end
    return through(face,winner)
end
@inline function incident(vertex::Int)
    local_vertex=mod1(vertex,4)
    return (local_vertex,mod1(local_vertex-1,4),vertex<=4 ? 5 : 6)
end
@inline accepts(states,wild,face,vertex)=
    has(wild,face) || contains(face,states[face],vertex)
@inline matching(pair_index::Int,state::UInt8)=pair_index==3 ? state : UInt8(3)-state

@inline function resolve_pair(states::NTuple{6,UInt8},wild::UInt8,
        pair_index::Int,ranks::NTuple{8,Int})
    first,second=PAIRS[pair_index]
    (states[first]!=0 || has(wild,first)) &&
        (states[second]!=0 || has(wild,second)) || return states,wild,false
    if states[first]!=0 && states[second]!=0
        aligned=pair_index==3 ? states[first]==states[second] : states[first]!=states[second]
        aligned || return states,wild,false
    end
    if has(wild,first) && has(wild,second)
        states=Base.setindex(states,lowest(first,ranks),first)
        wild&=~bit(first)
    end
    if states[first]!=0 && has(wild,second)
        states=Base.setindex(states,matching(pair_index,states[first]),second)
        wild&=~bit(second)
    elseif states[second]!=0 && has(wild,first)
        states=Base.setindex(states,matching(pair_index,states[second]),first)
        wild&=~bit(first)
    end
    return states,wild,true
end

@inline function initial_state(fixed,adjustable,free::UInt8,ranks,tier,face)
    fixed[face]!=0 && return fixed[face]
    tier==0 && adjustable[face]!=0 && return adjustable[face]
    ((tier in (1,2) && adjustable[face]!=0) ||
        (tier==2 && has(free,face))) && return lowest(face,ranks)
    return UInt8(0)
end
@inline initial_states(fixed,adjustable,free,ranks,tier)=
    ntuple(face->initial_state(fixed,adjustable,free,ranks,tier,face),6)
@inline function finish_states(fixed,adjustable,selected,done)
    return ntuple(face->has(done,face) ? selected[face] :
        fixed[face]!=0 ? fixed[face] : adjustable[face],6)
end

function full_hex(fixed::NTuple{6,UInt8},adjustable::NTuple{6,UInt8},
        free::UInt8,ranks::NTuple{8,Int})
    forbidden=UInt8(0)
    fixed_count=0;adjustable_count=0
    for face in 1:6
        fixed_count+=Int(fixed[face]!=0)
        adjustable_count+=Int(adjustable[face]!=0)
        fixed[face]==0 && adjustable[face]==0 && !has(free,face) &&
            (forbidden|=bit(face))
    end
    if count_ones(free)==6
        winner=1
        for vertex in 2:8
            ranks[vertex]<ranks[winner] && (winner=vertex)
        end
        output=fixed
        for face in incident(winner)
            output=Base.setindex(output,through(face,winner),face)
        end
        return output,forbidden,false
    end
    if 6-count_ones(forbidden)<2
        problem=fixed_count+adjustable_count>0
        if !problem
            for face in 1:6
                if has(free,face)
                    forbidden|=bit(face)
                    break
                end
            end
        end
        return finish_states(fixed,adjustable,fixed,UInt8(0)),forbidden,problem
    end

    required=fixed_count+adjustable_count
    held_pair=0;held_first=UInt8(0);held_second=UInt8(0)
    done=UInt8(0);found=false;selected=fixed
    for tier in 0:3
        states=initial_states(fixed,adjustable,free,ranks,tier)
        wild=UInt8(0)
        if tier==3
            for face in 1:6
                (adjustable[face]!=0 || has(free,face)) && (wild|=bit(face))
            end
        end
        for pair_index in 1:3
            states,wild,available=resolve_pair(states,wild,pair_index,ranks)
            available || continue
            first,second=PAIRS[pair_index]
            pair_mask=bit(first)|bit(second)
            if held_pair==0 && required<=Int(fixed[first]!=0 || adjustable[first]!=0)+
                    Int(fixed[second]!=0 || adjustable[second]!=0)
                held_pair=pair_index;held_first=states[first];held_second=states[second]
            end
            for other_index in 1:3
                other_index==pair_index && continue
                states,wild,available=resolve_pair(states,wild,other_index,ranks)
                if available
                    a,b=PAIRS[other_index]
                    done=pair_mask|bit(a)|bit(b)
                    selected=states;found=true
                    break
                end
            end
            found && break
            for vertex in 1:8
                count=0;corner_mask=UInt8(0)
                for face in incident(vertex)
                    if !has(pair_mask,face) && accepts(states,wild,face,vertex)
                        corner_mask|=bit(face);count+=1
                        count==2 && break
                    end
                end
                count==2 || continue
                done=pair_mask|corner_mask
                for face in incident(vertex)
                    has(wild,face) && (states=Base.setindex(states,through(face,vertex),face))
                end
                selected=states;found=true
                break
            end
            found && break
        end
        if !found
            for vertex in 1:8
                a,b,c=incident(vertex)
                accepts(states,wild,a,vertex) && accepts(states,wild,b,vertex) &&
                    accepts(states,wild,c,vertex) || continue
                done=bit(a)|bit(b)|bit(c)
                for face in (a,b,c)
                    has(wild,face) && (states=Base.setindex(states,through(face,vertex),face))
                end
                selected=states;found=true
                break
            end
        end
        found && break
    end
    if !found && held_pair!=0
        first,second=PAIRS[held_pair]
        done=bit(first)|bit(second)
        selected=Base.setindex(fixed,held_first,first)
        selected=Base.setindex(selected,held_second,second)
        forbidden|=free&~done
    end
    return finish_states(fixed,adjustable,selected,done),forbidden,!found && held_pair==0
end
end
const Local = LocalPhase
const Source = Model._ExtrudeNoNewRectGridSource

struct Phase
    source::Source
    intervals::Int
    top_states::Vector{UInt8}
    lateral::Matrix{UInt8}
    cap::Matrix{UInt8}
    lateral_forbidden::BitMatrix
    cap_forbidden::BitMatrix
    preferred::Vector{UInt8}
    problems::BitMatrix
end

@noinline fail(reason)=throw(ArgumentError("free B4 source phase: "*reason))
@inline reverse_state(state::UInt8)=state==0 ? UInt8(0) : UInt8(3)-state

function checked_sizes(source::Source,intervals::Int)
    intervals>0 || fail("positive intervals required")
    a,b=source.grid_shape
    min(a,b)==1 && max(a,b)>=5 || fail("certified B4 strip required")
    sizes=try
        m=Base.checked_mul(a,b)
        v=Base.checked_mul(Base.checked_add(a,1),Base.checked_add(b,1))
        e=Base.checked_add(Base.checked_mul(3,m),1)
        levels=Base.checked_add(intervals,1)
        macros=Base.checked_mul(m,intervals)
        nodes=Base.checked_add(Base.checked_mul(v,levels),macros)
        cells=Base.checked_mul(12,macros)
        Base.checked_mul(e,intervals)
        Base.checked_mul(m,levels)
        (m,v,e,nodes,cells)
    catch err
        err isa OverflowError || rethrow()
        fail("size arithmetic overflows Int")
    end
    m,v,e,nodes,cells=sizes
    nodes<=Model._EXTRUDE_NONEW_MAX_NODES && nodes<=typemax(Int32) && cells<=typemax(Int32) ||
        fail("checked output node or cell limit exceeded")
    length(source.source_cells)==length(source.edge_indices)==length(source.category)==length(source.boundary_masks)==m ||
        fail("certified cell dimensions differ")
    length(source.source_coordinates)==length(source.boundary_vertices)==v && all(source.boundary_vertices) ||
        fail("certified vertex dimensions differ")
    length(source.edges)==length(source.edge_cells)==length(source.edge_sides)==length(source.edge_curve)==length(source.edge_curve_direction)==e ||
        fail("certified edge dimensions differ")
    return m,v,e
end

function verify_incidence(source::Source,m::Int,v::Int,e::Int)
    for cell in 1:m
        nodes=source.source_cells[cell]
        allunique(nodes) && all(node->1<=node<=v,nodes) || fail("invalid original cell IDs")
        source.category[cell]==4 && source.boundary_masks[cell]==0x0f || fail("non-B4 source cell")
        for side in 1:4
            edge=Int(source.edge_indices[cell][side])
            1<=edge<=e || fail("invalid physical edge index")
            a,b=nodes[side],nodes[mod1(side+1,4)]
            source.edges[edge]==minmax(a,b) || fail("physical edge differs from cell side")
            owners=source.edge_cells[edge];sides=source.edge_sides[edge]
            (owners[1]==cell && sides[1]==side) || (owners[2]==cell && sides[2]==side) ||
                fail("cell side missing from physical incidence")
        end
    end
    for edge in 1:e
        edge_left,edge_right=source.edges[edge]
        1<=edge_left<edge_right<=length(source.source_coordinates) || fail("unnormalized physical edge")
        owners=source.edge_cells[edge];sides=source.edge_sides[edge]
        1<=owners[1]<=m && 1<=sides[1]<=4 || fail("invalid first physical incidence")
        source.edge_indices[owners[1]][sides[1]]==edge || fail("first incidence does not point back")
        if owners[2]==0
            sides[2]==0 && source.edge_curve[edge] in 0x01:0x04 && source.edge_curve_direction[edge] in (-1,1) ||
                fail("invalid exterior carrier or direction")
        else
            1<=owners[2]<=m && owners[2]!=owners[1] && 1<=sides[2]<=4 || fail("invalid second physical incidence")
            source.edge_indices[owners[2]][sides[2]]==edge || fail("second incidence does not point back")
            first=source.source_cells[owners[1]];second=source.source_cells[owners[2]]
            first[sides[1]]==second[mod1(Int(sides[2])+1,4)] &&
                first[mod1(Int(sides[1])+1,4)]==second[sides[2]] || fail("shared physical edge has equal orientation")
            source.edge_curve[edge]==0 && source.edge_curve_direction[edge]==0 || fail("internal edge has exterior carrier")
        end
    end
    return nothing
end

function Phase(source::Source,intervals::Int)
    m,v,e=checked_sizes(source,intervals)
    verify_incidence(source,m,v,e)
    top=Vector{UInt8}(undef,m)
    for cell in 1:m
        top[cell]=Model._extrude_nonew_b4_strip_top_state(source.source_cells[cell])
    end
    caps=zeros(UInt8,m,intervals+1)
    for cell in 1:m;caps[cell,end]=top[cell];end
    preferred=Vector{UInt8}(undef,e)
    for edge in 1:e
        direction=source.edge_curve_direction[edge]
        preferred[edge]=direction==0 ? 0x00 : direction>0 ? 0x01 : 0x02
    end
    return Phase(source,intervals,top,zeros(UInt8,e,intervals),caps,
        falses(e,intervals),falses(m,intervals+1),preferred,falses(m,intervals))
end

@inline function physical_face(phase::Phase,cell::Int,interval::Int,face::Int)
    if face<=4
        edge=Int(phase.source.edge_indices[cell][face])
        reversed=phase.source.source_cells[cell][face]!=phase.source.edges[edge][1]
        return true,edge,interval,reversed
    end
    return false,cell,interval+face-5,false
end

@inline function face_state(phase::Phase,cell::Int,interval::Int,face::Int)
    lateral,row,column,reversed=physical_face(phase,cell,interval,face)
    forbidden=lateral ? phase.lateral_forbidden[row,column] : phase.cap_forbidden[row,column]
    forbidden && return UInt8(0)
    state=lateral ? phase.lateral[row,column] : phase.cap[row,column]
    return reversed ? reverse_state(state) : state
end
@inline mask(phase::Phase,cell::Int,interval::Int)=ntuple(face->face_state(phase,cell,interval,face),Val(6))

function classification(phase::Phase,cell::Int,interval::Int)
    fixed=ntuple(_->UInt8(0),Val(6));adjustable=fixed;free=UInt8(0)
    for face in 1:6
        lateral,row,column,reversed=physical_face(phase,cell,interval,face)
        forbidden=lateral ? phase.lateral_forbidden[row,column] : phase.cap_forbidden[row,column]
        (forbidden || (face==5 && interval==1)) && continue
        state=lateral ? phase.lateral[row,column] : phase.cap[row,column]
        preferred=lateral ? phase.preferred[row] : UInt8(0)
        if state!=0
            fixed=Base.setindex(fixed,reversed ? reverse_state(state) : state,face)
        elseif preferred!=0
            adjustable=Base.setindex(adjustable,reversed ? reverse_state(preferred) : preferred,face)
        else
            free|=Local.bit(face)
        end
    end
    nodes=phase.source.source_cells[cell];levels=phase.intervals+1
    a=(Int(nodes[1])-1)*levels;b=(Int(nodes[2])-1)*levels
    c=(Int(nodes[3])-1)*levels;d=(Int(nodes[4])-1)*levels
    ranks=(a+interval,b+interval,c+interval,d+interval,a+interval+1,b+interval+1,c+interval+1,d+interval+1)
    return fixed,adjustable,free,ranks
end

function commit!(phase::Phase,cell::Int,interval::Int,output,forbidden::UInt8)
    for face in 1:6
        lateral,row,column,reversed=physical_face(phase,cell,interval,face)
        states=lateral ? phase.lateral : phase.cap
        restrictions=lateral ? phase.lateral_forbidden : phase.cap_forbidden
        if output[face]!=0
            !restrictions[row,column] || fail("committing a forbidden physical face")
            state=reversed ? reverse_state(output[face]) : output[face]
            (states[row,column]==0 || states[row,column]==state) || fail("conflicting physical diagonal")
            states[row,column]=state
        end
        if Local.has(forbidden,face)
            states[row,column]==0 || fail("forbidding a fixed physical diagonal")
            restrictions[row,column]=true
        end
    end
    return nothing
end

function brute_force(phase::Phase,cell::Int,interval::Int,classes)
    if phase.problems[cell,interval]
        fixed,adjustable,free,_=classes
        output=fixed;forbidden=UInt8(0)
        for face in 1:6
            fixed[face]==0 && (output=Base.setindex(output,adjustable[face],face))
            fixed[face]==0 && adjustable[face]==0 && !Local.has(free,face) && (forbidden|=Local.bit(face))
        end
        return output,forbidden,true
    end
    return Local.full_hex(classes...)
end

function run!(phase::Phase)
    for cell in eachindex(phase.source.source_cells)
        terminal=phase.top_states[cell]
        for interval in 1:phase.intervals
            interval==1 && phase.cap[cell,interval]==0 && (phase.cap_forbidden[cell,interval]=true)
            seeded=false
            for attempt in 1:2
                if attempt==1 && phase.cap[cell,interval+1]==0 && !phase.cap_forbidden[cell,interval+1]
                    phase.cap[cell,interval+1]=terminal
                    seeded=true
                elseif attempt==2
                    (seeded && phase.problems[cell,interval]) || break
                    phase.problems[cell,interval]=false
                    phase.cap[cell,interval+1]=0
                end
                output,forbidden,problem=brute_force(phase,cell,interval,classification(phase,cell,interval))
                commit!(phase,cell,interval,output,forbidden)
                problem && (phase.problems[cell,interval]=true)
            end
        end
    end
    return phase
end
end

module FinalFactories
# Deferred final-mask dispatch, independently reviewed literal CPP authority.
# Existing immutable315 records supply geometry; retained problems use their actual center fan.
const INDICES=(UInt16(1),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(239),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(214),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(281),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(262),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(229),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(264),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(254),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(237),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(56),UInt16(57),UInt16(58),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(234),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(288),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(59),UInt16(60),UInt16(61),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(231),UInt16(0),UInt16(0),UInt16(2),UInt16(0),UInt16(0),UInt16(3),UInt16(0),UInt16(0),UInt16(4),UInt16(0),UInt16(0),UInt16(5),UInt16(259),UInt16(271),UInt16(6),UInt16(256),UInt16(62),UInt16(7),UInt16(64),UInt16(0),UInt16(8),UInt16(0),UInt16(0),UInt16(9),UInt16(0),UInt16(0),UInt16(10),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(29),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(30),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(221),UInt16(31),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(32),UInt16(0),UInt16(0),UInt16(295),UInt16(0),UInt16(0),UInt16(246),UInt16(0),UInt16(0),UInt16(33),UInt16(0),UInt16(0),UInt16(251),UInt16(83),UInt16(84),UInt16(85),UInt16(86),UInt16(87),UInt16(34),UInt16(89),UInt16(90),UInt16(91),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(35),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(36),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(278),UInt16(226),UInt16(37),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(240),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(162),UInt16(163),UInt16(164),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(110),UInt16(0),UInt16(0),UInt16(111),UInt16(0),UInt16(0),UInt16(112),UInt16(0),UInt16(0),UInt16(113),UInt16(0),UInt16(0),UInt16(114),UInt16(0),UInt16(0),UInt16(115),UInt16(0),UInt16(215),UInt16(116),UInt16(165),UInt16(166),UInt16(117),UInt16(283),UInt16(228),UInt16(118),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(266),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(253),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(168),UInt16(169),UInt16(170),UInt16(0),UInt16(0),UInt16(0),UInt16(298),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(301),UInt16(0),UInt16(304),UInt16(0),UInt16(0),UInt16(0),UInt16(243),UInt16(0),UInt16(0),UInt16(242),UInt16(65),UInt16(66),UInt16(67),UInt16(0),UInt16(218),UInt16(0),UInt16(171),UInt16(172),UInt16(173),UInt16(0),UInt16(220),UInt16(0),UInt16(0),UInt16(0),UInt16(119),UInt16(0),UInt16(0),UInt16(120),UInt16(293),UInt16(291),UInt16(121),UInt16(0),UInt16(0),UInt16(122),UInt16(0),UInt16(0),UInt16(123),UInt16(68),UInt16(69),UInt16(70),UInt16(299),UInt16(217),UInt16(124),UInt16(174),UInt16(175),UInt16(125),UInt16(294),UInt16(235),UInt16(126),UInt16(0),UInt16(11),UInt16(0),UInt16(276),UInt16(12),UInt16(277),UInt16(0),UInt16(13),UInt16(0),UInt16(300),UInt16(14),UInt16(245),UInt16(274),UInt16(15),UInt16(260),UInt16(71),UInt16(16),UInt16(73),UInt16(0),UInt16(17),UInt16(0),UInt16(177),UInt16(18),UInt16(178),UInt16(0),UInt16(19),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(38),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(39),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(179),UInt16(180),UInt16(40),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(127),UInt16(0),UInt16(0),UInt16(41),UInt16(0),UInt16(0),UInt16(129),UInt16(0),UInt16(0),UInt16(130),UInt16(0),UInt16(0),UInt16(42),UInt16(0),UInt16(0),UInt16(132),UInt16(92),UInt16(93),UInt16(94),UInt16(95),UInt16(96),UInt16(43),UInt16(98),UInt16(99),UInt16(100),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(44),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(45),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(185),UInt16(186),UInt16(46),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(136),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(137),UInt16(241),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(138),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(139),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(282),UInt16(140),UInt16(263),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(141),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(265),UInt16(0),UInt16(0),UInt16(142),UInt16(0),UInt16(188),UInt16(189),UInt16(190),UInt16(191),UInt16(192),UInt16(193),UInt16(194),UInt16(143),UInt16(196),UInt16(0),UInt16(216),UInt16(0),UInt16(0),UInt16(238),UInt16(0),UInt16(0),UInt16(144),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(145),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(74),UInt16(75),UInt16(76),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(147),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(148),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(77),UInt16(78),UInt16(79),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(150),UInt16(0),UInt16(0),UInt16(20),UInt16(0),UInt16(0),UInt16(21),UInt16(0),UInt16(0),UInt16(22),UInt16(0),UInt16(197),UInt16(23),UInt16(199),UInt16(200),UInt16(24),UInt16(202),UInt16(80),UInt16(25),UInt16(82),UInt16(0),UInt16(26),UInt16(0),UInt16(0),UInt16(27),UInt16(0),UInt16(0),UInt16(28),UInt16(0),UInt16(307),UInt16(0),UInt16(0),UInt16(0),UInt16(0),UInt16(47),UInt16(0),UInt16(154),UInt16(0),UInt16(0),UInt16(0),UInt16(257),UInt16(0),UInt16(0),UInt16(48),UInt16(310),UInt16(155),UInt16(258),UInt16(0),UInt16(232),UInt16(0),UInt16(313),UInt16(224),UInt16(49),UInt16(0),UInt16(156),UInt16(0),UInt16(0),UInt16(0),UInt16(309),UInt16(0),UInt16(0),UInt16(50),UInt16(285),UInt16(157),UInt16(287),UInt16(0),UInt16(0),UInt16(249),UInt16(0),UInt16(0),UInt16(51),UInt16(284),UInt16(158),UInt16(248),UInt16(101),UInt16(102),UInt16(103),UInt16(104),UInt16(105),UInt16(52),UInt16(107),UInt16(108),UInt16(109),UInt16(0),UInt16(308),UInt16(0),UInt16(268),UInt16(267),UInt16(53),UInt16(0),UInt16(159),UInt16(0),UInt16(206),UInt16(207),UInt16(208),UInt16(209),UInt16(210),UInt16(54),UInt16(211),UInt16(160),UInt16(213),UInt16(0),UInt16(233),UInt16(0),UInt16(270),UInt16(223),UInt16(55),UInt16(0),UInt16(161),UInt16(0))

@inline function factory_index(states::NTuple{6,UInt8})
    for state in states
        state<=0x02 || throw(ArgumentError("invalid final physical face state"))
    end
    code=1+Int(states[1])+3Int(states[2])+9Int(states[3])+27Int(states[4])+81Int(states[5])+243Int(states[6])
    index=INDICES[code]
    index>0 || throw(ArgumentError("no existing-corner factory for final physical face state"))
    return index
end
end

module CenterFan
# Actual retained-center fan. Positive reference orientation points each
# base toward the actual interior mean; a problem retains its center even if
# subsequent propagation makes its final mask eligible for a corner factory.
const FACES=((1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8),
             (1,2,3,4),(5,6,7,8))
struct Cell
    msh::UInt8
    nodes::NTuple{8,UInt8}
end
struct Fan
    faces::NTuple{6,UInt8}
    ncells::UInt8
    cells::NTuple{12,Cell}
end
const EMPTY=Cell(0x00,ntuple(_->UInt8(0),Val(8)))
@inline tetra(a,b,c)=Cell(0x04,(UInt8(a),UInt8(b),UInt8(c),0x09,0x00,0x00,0x00,0x00))
@inline pyramid(a,b,c,d)=Cell(0x07,(UInt8(a),UInt8(b),UInt8(c),UInt8(d),0x09,0x00,0x00,0x00))

function fan(states::NTuple{6,UInt8})
    cells=ntuple(_->EMPTY,Val(12));filled=0
    for face in 1:6
        state=states[face]
        state in (0x00,0x01,0x02) || throw(ArgumentError("invalid physical face state"))
        a,b,c,d=FACES[face]
        inward=face==5
        if state==0
            cell=inward ? pyramid(a,b,c,d) : pyramid(b,a,d,c)
            filled+=1;cells=Base.setindex(cells,cell,filled)
        else
            if state==1
                first=inward ? tetra(a,b,c) : tetra(b,a,c)
                second=inward ? tetra(a,c,d) : tetra(c,a,d)
            else
                first=inward ? tetra(a,b,d) : tetra(b,a,d)
                second=inward ? tetra(b,c,d) : tetra(c,b,d)
            end
            filled+=1;cells=Base.setindex(cells,first,filled)
            filled+=1;cells=Base.setindex(cells,second,filled)
        end
    end
    return Fan(states,UInt8(filled),cells)
end
end

const Phase = SourcePhase
const Factory = FinalFactories
const Fan = CenterFan
const Source = Model._ExtrudeNoNewRectGridSource

# Recorded problems remain centers even if their final mask has a factory.
struct Catalog <: Model._ExtrudeNoNewDynamicGridCatalog
    source::Source
    top_states::Vector{UInt8}
    template_indices::Matrix{UInt16}
    problem_positions::Vector{NTuple{2,Int32}}
    problem_masks::Vector{NTuple{6,UInt8}}
    cell_counts::NTuple{4,Int}
    face_capacity::Int
    diagonal_count::Int
    recombined::Bool
    levels::Vector{Float64}
    layer_refs::Vector{NTuple{2,Int32}}
end

@noinline fail(reason)=throw(ArgumentError("free all-boundary strip: "*reason))

function catalog(phase::Phase.Phase,levels::Vector{Float64},refs::Vector{NTuple{2,Int32}})
    source=phase.source;m=length(source.source_cells);n=phase.intervals
    length(refs)==n && length(levels)==n+1 && levels[1]==0.0 && levels[end]==1.0 ||
        fail("inconsistent normalized layers")
    for level in eachindex(levels)
        isfinite(levels[level]) && (level==1 || levels[level]>levels[level-1]) ||
            fail("collapsed or nonfinite levels")
    end
    for cell in 1:m
        phase.cap[cell,1]==0 && phase.cap[cell,end]==phase.top_states[cell] ||
            fail("source/top final physical face state differs")
    end
    for edge in eachindex(source.edges),interval in 1:n
        source.edge_cells[edge][2]!=0 || phase.lateral[edge,interval] in (0x01,0x02) ||
            fail("unfinished exterior physical face")
    end
    indices=Matrix{UInt16}(undef,m,n)
    positions=NTuple{2,Int32}[]
    masks=NTuple{6,UInt8}[]
    centers=count(phase.problems)
    sizehint!(positions,centers);sizehint!(masks,centers)
    counts=(0,0,0,0)
    templates=Model._extrude_nonew_templates()
    for cell in 1:m,interval in 1:n
        states=Phase.mask(phase,cell,interval)
        if phase.problems[cell,interval]
            indices[cell,interval]=0
            push!(positions,(Int32(cell),Int32(interval)));push!(masks,states)
            diagonal=count(!iszero,states)
            counts=Base.setindex(counts,Base.checked_add(counts[1],2*diagonal),1)
            counts=Base.setindex(counts,Base.checked_add(counts[4],6-diagonal),4)
        else
            index=Factory.factory_index(states);indices[cell,interval]=index
            template=templates[Int(index)]
            template.faces==states || fail("factory mask differs from completed physical faces")
            for position in 1:Int(template.ncells)
                family=Int(template.cells[position].msh)-3
                1<=family<=4 || fail("unsupported existing-corner family")
                counts=Base.setindex(counts,Base.checked_add(counts[family],1),family)
            end
        end
    end
    length(positions)==length(masks)==centers || fail("retained center accounting differs")
    nodes=Base.checked_add(Base.checked_mul(length(source.source_coordinates),n+1),centers)
    cells=foldl(Base.checked_add,counts;init=0)
    nodes<=Model._EXTRUDE_NONEW_MAX_NODES && nodes<=typemax(Int32) && cells<=typemax(Int32) ||
        fail("actual output exceeds checked limits")
    incidence=0
    for family in 1:4
        incidence=Base.checked_add(incidence,Base.checked_mul((4,6,5,5)[family],counts[family]))
    end
    exterior=Base.checked_add(Base.checked_mul(3,m),Base.checked_mul(Base.checked_mul(4,m+1),n))
    doubled=Base.checked_add(incidence,exterior)
    iseven(doubled) || fail("nonintegral actual typed-face capacity")
    diagonals=Base.checked_add(count(!iszero,phase.lateral),count(!iszero,phase.cap))
    return Catalog(source,copy(phase.top_states),indices,positions,masks,counts,
        doubled÷2,diagonals,false,levels,refs)
end

function emit(completed::Catalog,cols,spec,caller)
    Model._extrude_nonew_rect_grid_product_certify(cols,completed,spec,caller)
    source=completed.source
    nlevels,intervals=Model._extrude_nonew_rect_grid_dimensions(cols,completed,caller)
    ncells=length(source.source_cells)
    size(completed.template_indices)==(ncells,intervals) || fail("final index dimensions differ")
    centers=length(completed.problem_masks)
    length(completed.problem_positions)==centers || fail("problem position dimensions differ")
    column_count=Base.checked_mul(length(source.source_coordinates),nlevels)
    node_count=Base.checked_add(column_count,centers)
    coordinates=Matrix{Float64}(undef,3,node_count)
    for node in eachindex(source.source_coordinates),level in 1:nlevels
        index=Model._extrude_nonew_rect_grid_node(Int32(node),level,nlevels)
        point=cols[level,node]
        coordinates[1,index]=point[1];coordinates[2,index]=point[2];coordinates[3,index]=point[3]
    end
    arrays=ntuple(family->Matrix{Int32}(undef,(4,8,6,5)[family],completed.cell_counts[family]),4)
    filled=zeros(Int,4)
    edges=Set{NTuple{2,NTuple{3,Float64}}}()
    sizehint!(edges,completed.diagonal_count)
    orientation=Int(source.source_orientation)*Int(source.normal_direction)
    templates=Model._extrude_nonew_templates()
    center_position=0
    for source_cell in 1:ncells,interval in 1:intervals
        original=source.source_cells[source_cell]
        corners=Model._extrude_nonew_corners(cols,original,interval)
        index=completed.template_indices[source_cell,interval]
        if index==0
            center_position+=1
            center_position<=centers || fail("unreserved actual center")
            completed.problem_positions[center_position]==(Int32(source_cell),Int32(interval)) ||
                fail("actual center order differs from original cell/interval order")
            states=completed.problem_masks[center_position]
            template=Fan.fan(states)
            Model._extrude_nonew_add_diagonals!(edges,corners,states)
            Model._extrude_nonew_certify(corners,edges,nothing,caller)
            center=Model._extrude_quadtri_centroid(corners)
            Model._extrude_nonew_certify_template((corners...,center),template,orientation,caller)
            center_index=Int32(column_count+center_position)
            coordinates[1,center_index]=center[1]
            coordinates[2,center_index]=center[2]
            coordinates[3,center_index]=center[3]
            Model._extrude_nonew_b4_strip_emit!(arrays,filled,template,original,
                interval,nlevels,center_index,caller)
        else
            template=templates[Int(index)]
            Model._extrude_nonew_add_diagonals!(edges,corners,template.faces)
            Model._extrude_nonew_certify(corners,edges,template,caller)
            Model._extrude_nonew_b4_strip_emit!(arrays,filled,template,original,
                interval,nlevels,Int32(0),caller)
        end
    end
    center_position==centers && Tuple(filled)==completed.cell_counts &&
        length(edges)==completed.diagonal_count || fail("actual final capacities differ")
    blocks=Model.ElementBlock[];sizehint!(blocks,4)
    for family in 1:4
        completed.cell_counts[family]>0 || continue
        Model._extrude_make_positive!(coordinates,arrays[family],family+3)
        push!(blocks,Model.ElementBlock(family+3,arrays[family]))
    end
    return (volume=Model.MixedMesh(coordinates,blocks),edges=edges,catalog=completed)
end

@noinline function _error(caller, reason)
    throw(ArgumentError("$caller: QuadTriNoNewVerts free all-boundary strip $reason"))
end

# At most one real retained center and twelve cells per source cell/interval.
# Layer construction performs the interval-dependent checked size admission.
function preflight(shape, caller)
    try
        vertices, cells, _, _ = Model._extrude_nonew_rect_grid_sizes(
            shape[1], shape[2], caller; mode=:b4_strip)
        cells_per_interval = Base.checked_mul(12, cells)
        Base.checked_add(vertices, cells)
        cells_per_interval <= typemax(Int32) ||
            _error(caller, "output cell capacity exceeds Int32")
        return vertices, cells_per_interval, 0, cells
    catch err
        err isa OverflowError && _error(caller, "source dimensions overflow Int")
        rethrow()
    end
end

function plan(source::Source, levels::Vector{Float64},
        refs::Vector{NTuple{2,Int32}}, caller)
    intervals = length(refs)
    intervals > 0 && length(levels) == Base.checked_add(intervals, 1) ||
        _error(caller, "requires one more level than positive intervals")
    levels[1] == 0.0 && levels[end] == 1.0 ||
        _error(caller, "requires normalized layer levels")
    for level in eachindex(levels)
        isfinite(levels[level]) || _error(caller, "requires finite layer levels")
        level == 1 || levels[level] > levels[level-1] ||
            _error(caller, "layer levels collapse in Float64")
    end
    try
        phase = SourcePhase.Phase(source, intervals)
        SourcePhase.run!(phase)
        # Factory dispatch is deferred until every physical column is complete.
        return catalog(phase, levels, refs)
    catch err
        err isa OverflowError && _error(caller, "size arithmetic overflows Int")
        if err isa ArgumentError && !startswith(err.msg, string(caller, ":"))
            _error(caller, err.msg)
        end
        rethrow()
    end
end

function finish(m::Model.GeoModel, t::Int, params, source_tag::Int,
        source_mesh, completed::Catalog, cols, top, laterals, caller)
    !params.recomb_laterals && !completed.recombined ||
        _error(caller, "requires free lateral faces")
    spec = get(m.meshing.extrude_specs, (3, t), nothing)
    spec !== nothing || _error(caller, "volume lacks extrusion provenance")
    actual = try
        emit(completed, cols, spec, caller)
    catch err
        err isa OverflowError && _error(caller, "size arithmetic overflows Int")
        if err isa ArgumentError && !startswith(err.msg, string(caller, ":"))
            _error(caller, err.msg)
        end
        rethrow()
    end
    sweep = Model._ExtrudeNoNewIndexedSweep(t, cols, actual.volume)
    return (sweep=sweep, edges=actual.edges, source_mesh=source_mesh,
        source_tag=source_tag, top_tag=top, lateral_tags=laterals, catalog=completed)
end

end
