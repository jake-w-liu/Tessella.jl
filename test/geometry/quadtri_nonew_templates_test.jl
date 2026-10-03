module QuadTriNoNewTemplateTests

using Tessella, Test
using Tessella.Model
using Tessella.MeshTypes: tet_signed_volume
using Tessella.Predicates: orient3

const CUBE = ((0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,1.0,0.0),(0.0,1.0,0.0),
              (0.0,0.0,1.0),(1.0,0.0,1.0),(1.0,1.0,1.0),(0.0,1.0,1.0))
const OUTER = ((1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8),
               (1,2,3,4),(5,6,7,8))
const FACES = Dict(4=>((1,3,2),(1,2,4),(2,3,4),(3,1,4)),
    5=>((1,4,3,2),(5,6,7,8),(1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8)),
    6=>((1,3,2),(4,5,6),(1,2,5,4),(2,3,6,5),(3,1,4,6)),
    7=>((1,4,3,2),(1,2,5),(2,3,5),(3,4,5),(4,1,5)))
const WIDTH = Dict(4=>4,5=>8,6=>6,7=>5)

function face_direction(nodes)
    rotations=[Tuple(circshift(collect(nodes),i)) for i in 0:length(nodes)-1]
    reverse_rotations=[Tuple(circshift(reverse(collect(nodes)),i)) for i in 0:length(nodes)-1]
    return minimum(rotations)<minimum(reverse_rotations) ? 1 : -1
end

function expected_boundary(states)
    boundary=Set{Tuple}()
    for (face,state) in zip(OUTER,states)
        a,b,c,d=face
        if state==0
            push!(boundary,Tuple(sort(collect(face))))
        elseif state==1
            push!(boundary,Tuple(sort([a,b,c])))
            push!(boundary,Tuple(sort([a,c,d])))
        else
            push!(boundary,Tuple(sort([a,b,d])))
            push!(boundary,Tuple(sort([b,c,d])))
        end
    end
    return boundary
end

function certificate(template)
    counts=Dict{Tuple,Int}()
    directions=Dict{Tuple,Int}()
    supports=Set{Tuple}()
    total=0.0
    for ci in 1:template.ncells
        cell=template.cells[ci]
        msh=Int(cell.msh)
        width=WIDTH[msh]
        nodes=Int.(cell.nodes[1:width])
        @test all(1<=n<=8 for n in nodes)
        @test length(unique(nodes))==width
        @test all(==(0),cell.nodes[width+1:8])
        @test !((msh,Tuple(sort(collect(nodes)))) in supports)
        push!(supports,(msh,Tuple(sort(collect(nodes)))))
        center=ntuple(k->sum(CUBE[n][k] for n in nodes)/width,3)
        volume=0.0
        for face in FACES[msh]
            global_face=Tuple(nodes[k] for k in face)
            key=Tuple(sort(collect(global_face)))
            counts[key]=get(counts,key,0)+1
            directions[key]=get(directions,key,0)+face_direction(global_face)
            for fi in 2:length(face)-1
                a,b,c=CUBE[nodes[face[1]]],CUBE[nodes[face[fi]]],CUBE[nodes[face[fi+1]]]
                @test orient3(a,b,c,center)>0
                volume+=tet_signed_volume(a,c,b,center)
            end
            if length(face)==4
                a,b,c,d=(CUBE[nodes[k]] for k in face)
                @test orient3(a,b,d,center)>0
                @test orient3(b,c,d,center)>0
            end
        end
        @test volume>0
        total+=volume
    end
    @test isapprox(total,1;atol=8eps(Float64),rtol=0)
    @test all(n in (1,2) for n in values(counts))
    @test all(counts[face]==1 || directions[face]==0 for face in keys(counts))
    @test Set(face for (face,count) in counts if count==1)==expected_boundary(template.faces)
    return supports
end

function empty_sweep()
    return Model._ExtrudeVolumeSweep(1,true,Matrix{NTuple{3,Float64}}(undef,0,0),
        NTuple{6,NTuple{3,Float64}}[],NTuple{4,NTuple{3,Float64}}[],
        NTuple{8,NTuple{3,Float64}}[],NTuple{6,NTuple{3,Float64}}[],
        NTuple{5,NTuple{3,Float64}}[],NTuple{3,Float64}[],Bool[false])
end

function certificate_emit(template,supports)
    sweep=empty_sweep()
    Model._extrude_nonew_emit_template!(sweep,CUBE,template)
    emitted=Set{Tuple}()
    for (msh,cells) in ((4,sweep.tets),(5,sweep.hexes),(6,sweep.prisms),(7,sweep.pyramids))
        for cell in cells
            ids=[findfirst(==(point),CUBE) for point in cell]
            push!(emitted,(msh,Tuple(sort(ids))))
        end
    end
    @test emitted==supports
    @test isempty(sweep.interior)
    @test !sweep.degenerate[]
end

function warmed_emit_bytes(template)
    sweep=empty_sweep()
    for cells in (sweep.tets,sweep.hexes,sweep.prisms,sweep.pyramids)
        sizehint!(cells,6)
    end
    Model._extrude_nonew_emit_template!(sweep,CUBE,template)
    for cells in (sweep.tets,sweep.hexes,sweep.prisms,sweep.pyramids)
        empty!(cells)
    end
    return @allocated Model._extrude_nonew_emit_template!(sweep,CUBE,template)
end

@testset "NoNew full-hex independent template certificate" begin
    templates=Model._extrude_nonew_templates()
    @test isbitstype(Model._ExtrudeNoNewCell)
    @test isbitstype(Model._ExtrudeNoNewTemplate)
    @test !isempty(templates)
    keys=Set{Tuple}()
    for template in templates
        @test 1<=template.ncells<=6
        @test all(state in 0:2 for state in template.faces)
        supports=certificate(template)
        key=(template.faces,Tuple(sort!(collect(supports))))
        @test !(key in keys)
        push!(keys,key)
        certificate_emit(template,supports)
    end
    @test count(t->t.kind==0,templates)==1
    @test Set(Int(t.choice) for t in templates if t.kind==1)==Set(1:8)
    @test Set(Int(t.choice) for t in templates if t.kind==2)==Set(1:6)
    for template in templates
        warmed_emit_bytes(template)
        @test warmed_emit_bytes(template)==0
    end
    boxed=String[]
    for symbol in names(Model;all=true)
        startswith(String(symbol),"_extrude_nonew_") || continue
        value=getfield(Model,symbol)
        value isa Function || continue
        for method in methods(value)
            occursin("Core.Box",sprint(show,Base.uncompressed_ast(method))) &&
                push!(boxed,String(symbol))
        end
    end
    @test isempty(boxed)
    println("TEMPLATE_COUNT ",length(templates))
    println("FACE_MASK_COUNT ",length(Set(t.faces for t in templates)))
    println("KIND_COUNTS ",Dict(kind=>count(t->t.kind==kind,templates) for kind in 0:2))
    println("TEMPLATE_BYTES ",sizeof(Model._ExtrudeNoNewTemplate))
    println("BOXED ",boxed)
end


# Independent one-based copies of prism_v from the pinned full-hex factory.
const SLICES = (((1,2,6,4,3,7),(1,5,6,4,8,7)),
    ((3,4,8,2,1,5),(3,7,8,2,6,5)),
    ((2,3,7,1,4,8),(2,6,7,1,5,8)),
    ((4,1,5,3,2,6),(4,8,5,3,7,6)),
    ((1,2,3,5,6,7),(1,4,3,5,8,7)),
    ((2,3,4,6,7,8),(2,1,4,6,5,8)))
const QUADS = ((1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8),
               (1,2,3,4),(5,6,7,8))
const PRISM_QUADS = ((1,2,5,4),(2,3,6,5),(3,1,4,6))
const PRISM_TRIS = ((1,2,3),(4,5,6))

edge(face,state) = state==0 ? (0,0) :
    state==1 ? minmax(face[1],face[3]) : minmax(face[2],face[4])
support(face) = Tuple(sort(collect(face)))

function corner_possible(faces,states,nvertices)
    for apex in 1:nvertices
        valid=true
        for (face,state) in zip(faces,states)
            if apex in face && !(apex in edge(face,state))
                valid=false;break
            end
        end
        valid && return true
    end
    return false
end

function prism_possible(states)
    all(==(0),states) && return true
    return corner_possible(PRISM_QUADS,states,6)
end

function global_boundary(states)
    boundary=Set{Tuple}()
    for (face,state) in zip(QUADS,states)
        a,b,c,d=face
        if state==0
            push!(boundary,support(face))
        elseif state==1
            push!(boundary,support((a,b,c)));push!(boundary,support((a,c,d)))
        else
            push!(boundary,support((a,b,d)));push!(boundary,support((b,c,d)))
        end
    end
    return boundary
end

function subprism_states(mapping,states,internal_edge)
    result=Int[]
    for face in PRISM_QUADS
        global_face=Tuple(mapping[i] for i in face)
        index=findfirst(q->support(q)==support(global_face),QUADS)
        diagonal=index===nothing ? internal_edge : edge(QUADS[index],states[index])
        if diagonal==(0,0)
            push!(result,0)
        elseif diagonal==edge(global_face,1)
            push!(result,1)
        elseif diagonal==edge(global_face,2)
            push!(result,2)
        else
            return nothing
        end
    end
    return Tuple(result)
end

function slice_possible(mapping_pair,states)
    boundary=global_boundary(states)
    for mapping in mapping_pair,face in PRISM_TRIS
        support(Tuple(mapping[i] for i in face)) in boundary || return false
    end
    first_mapping,second_mapping=mapping_pair
    first_quads=[Tuple(first_mapping[i] for i in q) for q in PRISM_QUADS]
    second_quads=[Tuple(second_mapping[i] for i in q) for q in PRISM_QUADS]
    shared=findfirst(q->any(other->support(q)==support(other),second_quads),first_quads)
    shared===nothing && error("reference slice has no shared quadrangle")
    for state in 0:2
        diagonal=edge(first_quads[shared],state)
        first_states=subprism_states(first_mapping,states,diagonal)
        second_states=subprism_states(second_mapping,states,diagonal)
        first_states===nothing && continue
        second_states===nothing && continue
        prism_possible(first_states) && prism_possible(second_states) && return true
    end
    return false
end

function factory_possible(states)
    all(==(0),states) && return true
    corner_possible(QUADS,states,8) && return true
    return any(pair->slice_possible(pair,states),SLICES)
end

@testset "Complete primary full-hex factory face relation" begin
    @test count(prism_possible,Iterators.product(0:2,0:2,0:2))==13
    masks=Set(Int.(template.faces) for template in Model._extrude_nonew_templates())
    for states in Iterators.product(ntuple(_->0:2,6)...)
        @test (states in masks)==factory_possible(states)
    end
    @test !((0,0,0,0,0,1) in masks)
    @test !((0,0,0,0,0,2) in masks)
    @test ((0,0,0,0,0,0) in masks)
    println("FACTORY_FACE_MASKS ",length(masks))
end


end
