# GNU11.3 libstdc++ _Prime_rehash_policy and unique-key hashtable list order.
# Copyright (C) 2007-2021 Free Software Foundation, Inc.
# Derived policy constants follow the GCC Runtime Library Exception.
# This is an explicit candidate policy; public native topology proof is separate.
const _GMSH_INITIAL_BUCKET_PRIMES=(
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 103, 109, 113, 127, 137,
    139, 149, 157, 167, 179, 193, 199, 211, 227, 241,
    257, 277, 293, 313, 337, 359, 383, 409, 439, 467,
    503, 541, 577, 619, 661, 709, 761, 823, 887, 953,
    1031, 1109, 1193, 1289, 1381, 1493, 1613, 1741, 1879, 2029,
    2179, 2357, 2549, 2753, 2971, 3209, 3469, 3739, 4027, 4349,
    4703, 5087, 5503, 5953, 6427, 6949, 7517, 8123, 8783, 9497,
    10273, 11113, 12011, 12983, 14033, 15173, 16411, 17749, 19183, 20753,
    22447, 24281, 26267, 28411, 30727, 33223, 35933, 38873, 42043, 45481,
    49201, 53201, 57557, 62233, 67307, 72817, 78779, 85229, 92203, 99733,
    107897, 116731, 126271, 136607, 147793, 159871, 172933, 187091, 202409, 218971,
    236897, 256279, 277261, 299951, 324503, 351061, 379787, 410857, 444487, 480881,
    520241, 562841, 608903, 658753, 712697, 771049, 834181, 902483, 976369, 1056323,
    1142821, 1236397, 1337629, 1447153, 1565659, 1693859, 1832561, 1982627, 2144977, 2320627,
    2510653, 2716249, 2938679, 3179303, 3439651, 3721303, 4026031, 4355707, 4712381, 5098259,
    5515729, 5967347, 6456007, 6984629, 7556579, 8175383, 8844859, 9569143, 10352717, 11200489,
    12117689, 13109983, 14183539, 15345007, 16601593, 17961079, 19431899, 21023161, 22744717, 24607243,
    26622317, 28802401, 31160981, 33712729, 36473443, 39460231, 42691603, 46187573, 49969847, 54061849,
    58488943, 63278561, 68460391, 74066549, 80131819, 86693767, 93793069, 101473717, 109783337, 118773397,
    128499677, 139022417, 150406843, 162723577, 176048909, 190465427, 206062531, 222936881, 241193053, 260944219,
    282312799, 305431229, 330442829, 357502601, 386778277, 418451333, 452718089, 489790921, 529899637, 573292817,
    620239453, 671030513, 725980837, 785430967, 849749479, 919334987, 994618837, 1076067617, 1164186217, 1259520799,
    1362662261, 1474249943, 1594975441, 1725587117, 1866894511, 2019773507, 2185171673, 2364114217, 2557710269, 2767159799,
    2993761039, 3238918481, 3504151727, 3791104843, 4101556399, 4294967291,
)

function _gmsh_initial_bucket_count(n)
    lo=1;hi=length(_GMSH_INITIAL_BUCKET_PRIMES)
    while lo<hi
        mid=lo+div(hi-lo,2)
        if _GMSH_INITIAL_BUCKET_PRIMES[mid]<n
            lo=mid+1
        else
            hi=mid
        end
    end
    value=_GMSH_INITIAL_BUCKET_PRIMES[lo]
    value<=typemax(Int32) || throw(ArgumentError("native initial mesh: bucket capacity exceeds Int32"))
    return value
end

function _gmsh_initial_rehash!(next,tags,buckets,count)
    replacement=fill(Int32(-1),count)
    node=next[1];next[1]=0;begin_bucket=0
    while node!=0
        following=next[node+1]
        bucket=Int(mod(tags[node],UInt64(count)))+1
        if replacement[bucket]==-1
            next[node+1]=next[1];next[1]=node;replacement[bucket]=0
            following_head=next[node+1]
            following_head!=0 && (replacement[begin_bucket]=node)
            begin_bucket=bucket
        else
            predecessor=replacement[bucket]
            next[node+1]=next[predecessor+1]
            next[predecessor+1]=node
        end
        node=following
    end
    return replacement
end

function _gmsh_initial_uid_order(input)
    tags=UInt64[];source_rows=Int[];indices=Dict{UInt64,Int32}()
    next=Int32[0];buckets=Int32[-1];resize_at=0
    for (row,raw_tag) in enumerate(input)
        raw_tag isa Integer && !(raw_tag isa Bool) && raw_tag>0 ||
            throw(ArgumentError("native initial mesh: source UID must be a positive integer"))
        tag=UInt64(raw_tag)
        haskey(indices,tag) && continue
        count=length(tags)
        count<typemax(Int32)-1 || throw(ArgumentError("native initial mesh: UID storage exceeds Int32"))
        if count+1>resize_at
            minimum_buckets=max(count+1,resize_at==0 ? 11 : 0)
            if minimum_buckets>=length(buckets)
                requested=max(minimum_buckets+1,Base.checked_mul(length(buckets),2))
                new_count=_gmsh_initial_bucket_count(requested)
                buckets=_gmsh_initial_rehash!(next,tags,buckets,new_count)
                resize_at=new_count
            else
                resize_at=length(buckets)
            end
        end
        push!(tags,tag);push!(source_rows,row);push!(next,0)
        node=Int32(count+1);indices[tag]=node
        bucket=Int(mod(tag,UInt64(length(buckets))))+1
        predecessor=buckets[bucket]
        if predecessor!=-1
            next[node+1]=next[predecessor+1]
            next[predecessor+1]=node
        else
            next[node+1]=next[1];next[1]=node
            if next[node+1]!=0
                old_bucket=Int(mod(tags[next[node+1]],UInt64(length(buckets))))+1
                buckets[old_bucket]=node
            end
            buckets[bucket]=0
        end
    end
    result=Int[];sizehint!(result,length(tags));node=next[1]
    while node!=0
        push!(result,source_rows[node]);node=next[node+1]
    end
    return result
end
