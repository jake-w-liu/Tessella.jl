# Native interior quality policy, evaluated with an independently derived
# fixed-size Gaussian/barycentric-gradient kernel. Pool topology and birth
# decisions remain in the bounded arena pass below.
include("ModelNativeTetrahedronQualityFactors.jl")
const QUALITY_EDGE_VERSIONS=(0,1,2,3,8,5)
const QUALITY_NORMAL_PAIRS=((1,2),(2,3),(3,4),(1,4),(3,1),(2,4))
const QUALITY_MAX_COS=cos(179.0/180.0*pi)

struct QualityWorkspace end
@inline quality_dihedrals(points,::QualityWorkspace) = QualityFactors.tetrahedron_cosines(points)

@inline function _quality_improves(next,previous)
    difference=next-previous
    abs(difference/previous)<1e-8 && (difference=0.0)
    return difference>0.0
end

# The kernel scan avoids constructing a second incidence arena when the native
# pre-constructor quality stage has no bad real tetrahedron.
function has_bad_quality(kernel,work::QualityWorkspace)
    for cell in kernel.cells
        vertices=cell.vertices
        vertices[1]!=0 && all(node->node<=kernel.count,vertices) || continue
        values=quality_dihedrals(ntuple(index->kernel.points[vertices[index]],4),work)
        minimum(values)<QUALITY_MAX_COS && return true
    end
    return false
end

function improve_quality32!(arena::Arena,work::QualityWorkspace=QualityWorkspace())
    queue=Tuple{Int32,NTuple{4,Int32}}[]
    for (slot,cell) in enumerate(arena.cells)
        cell.alive && all(node->node<=arena.nreal,cell.vertices) || continue
        quality=quality_dihedrals(ntuple(index->arena.points[cell.vertices[index]],4),work)
        minimum(quality)<QUALITY_MAX_COS && push!(queue,(Int32(slot),cell.vertices))
    end
    flips=0
    while !isempty(queue)
        nextqueue=Tuple{Int32,NTuple{4,Int32}}[]
        for (slot,vertices) in queue
            cell=arena.cells[slot]
            cell.alive && cell.vertices==vertices || continue
            quality=quality_dihedrals(ntuple(index->arena.points[vertices[index]],4),work)
            for edge in 1:6
                quality[edge]<QUALITY_MAX_COS || continue
                handle=Handle(slot,QUALITY_EDGE_VERSIONS[edge])
                ring=try
                    spin(arena,handle)
                catch error
                    error isa Unsupported || rethrow()
                    continue
                end
                length(ring)==3 || continue
                all(h->all(node->node<=arena.nreal,arena.cells[h.cell].vertices),ring) || continue
                first,second,third=ring
                e=org(arena,first);d=dest(arena,first)
                p=apex(arena,first);q=apex(arena,second);r=apex(arena,third)
                candidates=((p,q,r,d),(q,p,r,e))
                newquality=ntuple(2) do index
                    minimum(quality_dihedrals(
                        ntuple(vertex->arena.points[candidates[index][vertex]],4),work))
                end
                all(value->_quality_improves(value,quality[edge]),newquality) || continue
                flip32!(arena,ring) || continue
                flips+=1
                # Every admitted operation deletes one real pool cell. Thus
                # both successful work and queue growth are bounded by input.
                for h in (first,second)
                    current=arena.cells[h.cell].vertices
                    values=quality_dihedrals(ntuple(index->arena.points[current[index]],4),work)
                    minimum(values)<QUALITY_MAX_COS && push!(nextqueue,(h.cell,current))
                end
                break
            end
        end
        queue=nextqueue
    end
    return flips
end
