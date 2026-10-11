module HilbertIntervals

# The geometry uses binary reflected Gray codes and an entry/exit frame.
# Background: Hamilton and Rau-Chaplin, "Compact Hilbert Indices for
# Multi-Dimensional Data", CISIS 2007, section 2.
# https://web.cs.dal.ca/~arc/publications/2-43/paper.pdf
# Contiguous Gray-rank intervals form a balanced binary partition tree.
# The implementation below generates that tree; it does not encode a fixed
# sequence of seven partition assignments.

struct Frame
    entry::Int
    exit_axis::Int
end

const CHILD_ENTRY=(0,0,0,3,3,6,6,5)
const CHILD_EXIT_AXIS=(0,1,1,2,2,1,1,0)

@inline function rotate3(bits::Int,amount::Int)
    shift=mod(amount,3)
    return ((bits<<shift)|(bits>>(3-shift)))&7
end
@inline function corner(frame::Frame,rank::Int)
    reflected=rank⊻(rank>>1)
    return frame.entry⊻rotate3(reflected,frame.exit_axis+1)
end
@inline function child_frame(frame::Frame,rank::Int)
    return Frame(frame.entry⊻rotate3(CHILD_ENTRY[rank+1],frame.exit_axis+1),
                 mod(CHILD_EXIT_AXIS[rank+1]+frame.exit_axis+1,3))
end

# A deterministic two-finger partition. The half containing the first
# corner uses a strict comparison; points on the midpoint join the other
# half. Swapping those fingers preserves the specified input-row tie order.
function partition!(order,points,start::Int,count::Int,axis::Int,
                    midpoint::Float64,left_is_lower::Bool)
    first=start
    last=start+count-1
    while first<=last
        while first<=last
            coordinate=points[order[first]][axis]
            (left_is_lower ? coordinate<midpoint : coordinate>midpoint) || break
            first+=1
        end
        while last>=first
            coordinate=points[order[last]][axis]
            (left_is_lower ? coordinate>=midpoint : coordinate<=midpoint) || break
            last-=1
        end
        first>last && break
        order[first],order[last]=order[last],order[first]
    end
    return first-start
end

function bisect_ranks!(order,points,start::Int,offset::Int,count::Int,
                       first_rank::Int,last_rank::Int,frame::Frame,lo,hi,
                       boundaries::NTuple{9,Int})
    middle=(first_rank+last_rank)÷2
    left_corner=corner(frame,middle-1)
    right_corner=corner(frame,middle)
    axis=trailing_zeros(left_corner⊻right_corner)+1
    midpoint=0.5*(lo[axis]+hi[axis])
    left_is_lower=(left_corner&(1<<(axis-1)))==0
    left_count=partition!(order,points,start+offset,count,axis,midpoint,left_is_lower)
    middle_offset=offset+left_count
    boundaries=Base.setindex(boundaries,middle_offset,middle+1)
    if middle-first_rank>1
        boundaries=bisect_ranks!(order,points,start,offset,left_count,
            first_rank,middle,frame,lo,hi,boundaries)
    end
    if last_rank-middle>1
        boundaries=bisect_ranks!(order,points,start,middle_offset,count-left_count,
            middle,last_rank,frame,lo,hi,boundaries)
    end
    return boundaries
end

function rank_boundaries!(order,points,start::Int,count::Int,frame::Frame,lo,hi)
    boundaries=(0,0,0,0,0,0,0,0,count)
    return bisect_ranks!(order,points,start,0,count,0,8,frame,lo,hi,boundaries)
end

@inline function child_bounds(frame::Frame,rank::Int,lo,hi)
    code=corner(frame,rank)
    midpoint=ntuple(axis->0.5*(lo[axis]+hi[axis]),3)
    lower=ntuple(axis->(code&(1<<(axis-1)))==0 ? lo[axis] : midpoint[axis],3)
    upper=ntuple(axis->(code&(1<<(axis-1)))==0 ? midpoint[axis] : hi[axis],3)
    return lower,upper
end

function sort_cell!(order,points,start::Int,count::Int,frame::Frame,lo,hi,level::Int)
    boundaries=rank_boundaries!(order,points,start,count,frame,lo,hi)
    level+1==1000 && return nothing
    for rank in 0:7
        offset=boundaries[rank+1]
        child_count=boundaries[rank+2]-offset
        child_count>2 || continue
        lower,upper=child_bounds(frame,rank,lo,hi)
        sort_cell!(order,points,start+offset,child_count,child_frame(frame,rank),
                   lower,upper,level+1)
    end
    return nothing
end

function multiscale!(order,points,start::Int,count::Int,lo,hi)
    # Production row counts are checked Int32 values, for which division by
    # eight equals truncation of the binary64 product count*0.125 exactly.
    prefix=count>=64 ? count÷8 : 0
    prefix>0 && multiscale!(order,points,start,prefix,lo,hi)
    sort_cell!(order,points,start+prefix,count-prefix,Frame(0,0),lo,hi,0)
    return nothing
end

function order!(order,points,lo,hi)
    multiscale!(order,points,1,length(order),lo,hi)
    return order
end

end
