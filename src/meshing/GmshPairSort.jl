# Adapted from the original HP/SGI STL stl_algo.h and stl_heap.h,
# https://github.com/karottc/sgi-stl/blob/master/stl_algo.h
# https://github.com/karottc/sgi-stl/blob/master/stl_heap.h
# Full notices from both permissive source headers are retained below.

#
# Copyright (c) 1994
# Hewlett-Packard Company
#
# Permission to use, copy, modify, distribute and sell this software
# and its documentation for any purpose is hereby granted without fee,
# provided that the above copyright notice appear in all copies and
# that both that copyright notice and this permission notice appear
# in supporting documentation.  Hewlett-Packard Company makes no
# representations about the suitability of this software for any
# purpose.  It is provided "as is" without express or implied warranty.
#
#
# Copyright (c) 1996
# Silicon Graphics Computer Systems, Inc.
#
# Permission to use, copy, modify, distribute and sell this software
# and its documentation for any purpose is hereby granted without fee,
# provided that the above copyright notice appear in all copies and
# that both that copyright notice and this permission notice appear
# in supporting documentation.  Silicon Graphics makes no
# representations about the suitability of this software for any
# purpose.  It is provided "as is" without express or implied warranty.

#
# Copyright (c) 1994
# Hewlett-Packard Company
#
# Permission to use, copy, modify, distribute and sell this software
# and its documentation for any purpose is hereby granted without fee,
# provided that the above copyright notice appear in all copies and
# that both that copyright notice and this permission notice appear
# in supporting documentation.  Hewlett-Packard Company makes no
# representations about the suitability of this software for any
# purpose.  It is provided "as is" without express or implied warranty.
#
# Copyright (c) 1997
# Silicon Graphics Computer Systems, Inc.
#
# Permission to use, copy, modify, distribute and sell this software
# and its documentation for any purpose is hereby granted without fee,
# provided that the above copyright notice appear in all copies and
# that both that copyright notice and this permission notice appear
# in supporting documentation.  Silicon Graphics makes no
# representations about the suitability of this software for any
# purpose.  It is provided "as is" without express or implied warranty.

# The HP/SGI STL introsort/heap/insertion algorithms below are adapted to Julia.
# Native Windows Gmsh comparison/pivot/tie behavior is established by actual
# public pair/element-label oracles. Its observed policy uses median-of-three,
# a depth budget of twice floor(log2(n)), a sixteen-entry insertion prefix,
# and heapsort at exhausted depth. Ties and unordered NaNs retain this policy.
# The implementation is independently checked through public primary pair and
# element-label identities, including lists above the partition threshold.
function _gmsh_pair_heap_adjust!(pairs,start,hole,count,value)
    top=hole;child=2hole+2
    while child<count
        _gmsh_pair_less(pairs[start+child],pairs[start+child-1]) && (child-=1)
        pairs[start+hole]=pairs[start+child];hole=child
        child=2child+2
    end
    if child==count
        pairs[start+hole]=pairs[start+child-1];hole=child-1
    end
    parent=(hole-1)÷2
    while hole>top && _gmsh_pair_less(pairs[start+parent],value)
        pairs[start+hole]=pairs[start+parent];hole=parent;parent=(hole-1)÷2
    end
    pairs[start+hole]=value
    return nothing
end

function _gmsh_pair_heapsort!(pairs,start,stop)
    count=stop-start
    count<2 && return nothing
    for parent in (count-2)÷2:-1:0
        _gmsh_pair_heap_adjust!(pairs,start,parent,count,pairs[start+parent])
    end
    while count>1
        count-=1
        saved=pairs[start+count];pairs[start+count]=pairs[start]
        _gmsh_pair_heap_adjust!(pairs,start,0,count,saved)
    end
    return nothing
end

function _gmsh_pair_partition!(pairs,start,stop)
    first=start+1;middle=start+(stop-start)÷2;last=stop-1
    # Encode the eight median comparison outcomes explicitly. False NaN
    # comparisons select the same pivot index as the pinned public oracle.
    outcomes=4Int(_gmsh_pair_less(pairs[first],pairs[middle]))+
        2Int(_gmsh_pair_less(pairs[middle],pairs[last]))+
        Int(_gmsh_pair_less(pairs[first],pairs[last]))
    median=(middle,first,last,first,first,last,middle,middle)[outcomes+1]
    pairs[start],pairs[median]=pairs[median],pairs[start]
    left=start+1;right=stop
    while true
        while _gmsh_pair_less(pairs[left],pairs[start])
            left+=1
        end
        right-=1
        while _gmsh_pair_less(pairs[start],pairs[right])
            right-=1
        end
        left>=right && return left
        pairs[left],pairs[right]=pairs[right],pairs[left];left+=1
    end
end

function _gmsh_pair_introsort!(pairs,start,stop,depth)
    while stop-start>16
        if depth==0
            _gmsh_pair_heapsort!(pairs,start,stop)
            return nothing
        end
        depth-=1
        split=_gmsh_pair_partition!(pairs,start,stop)
        _gmsh_pair_introsort!(pairs,split,stop,depth)
        stop=split
    end
    return nothing
end

function _gmsh_pair_insert!(pairs,index)
    saved=pairs[index];previous=index-1
    while _gmsh_pair_less(saved,pairs[previous])
        pairs[index]=pairs[previous];index=previous;previous-=1
    end
    pairs[index]=saved
    return nothing
end

function _gmsh_pair_sort!(pairs)
    count=length(pairs)
    count<2 && return pairs
    depth=2*(8sizeof(Int)-1-leading_zeros(count))
    _gmsh_pair_introsort!(pairs,1,count+1,depth)
    prefix=min(count,16)
    for index in 2:prefix
        if _gmsh_pair_less(pairs[index],pairs[1])
            saved=pairs[index]
            for previous in index-1:-1:1
                pairs[previous+1]=pairs[previous]
            end
            pairs[1]=saved
        else
            _gmsh_pair_insert!(pairs,index)
        end
    end
    for index in prefix+1:count
        _gmsh_pair_insert!(pairs,index)
    end
    return pairs
end
