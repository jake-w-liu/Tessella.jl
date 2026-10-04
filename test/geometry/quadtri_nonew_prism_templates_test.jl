module QuadTriNoNewPrismTemplateTests

using Tessella, Test, SHA
using Tessella.Model

const UNIT = ((0.0,0.0,0.0),(1.0,0.0,0.0),(0.0,1.0,0.0),
              (0.0,0.0,1.0),(1.0,0.0,1.0),(0.0,1.0,1.0))
const CAPS = ((1,3,2),(4,5,6))
const LATERALS = ((1,2,5,4),(2,3,6,5),(3,1,4,6))
const FACES = Dict(4=>((1,3,2),(1,2,4),(2,3,4),(3,1,4)),
    6=>((1,3,2),(4,5,6),(1,2,5,4),(2,3,6,5),(3,1,4,6)),
    7=>((1,4,3,2),(1,2,5),(2,3,5),(3,4,5),(4,1,5)))
const WIDTH = Dict(4=>4,6=>6,7=>5)
const RANKED_MASKS = ((0,0,0),(1,0,2),(1,1,2),(1,2,2),(2,1,0),
    (2,1,1),(2,1,2),(0,2,1),(1,2,1),(2,2,1),(2,0,1),(1,2,0),(0,1,2))
const Q = Rational{BigInt}

point(p) = ntuple(axis->Q(p[axis]),3)
subtract(a,b) = ntuple(axis->a[axis]-b[axis],3)
add(a,b) = ntuple(axis->a[axis]+b[axis],3)
scale(a,t) = ntuple(axis->a[axis]*t,3)
determinant(a,b,c) = a[1]*(b[2]*c[3]-b[3]*c[2])-
    a[2]*(b[1]*c[3]-b[3]*c[1])+a[3]*(b[1]*c[2]-b[2]*c[1])
tet_determinant(v) = determinant(subtract(v[2],v[1]),
    subtract(v[3],v[1]),subtract(v[4],v[1]))

# Independent analytic reference gradients; no production evaluator or
# Jacobian certificate is used to establish the template's full-map sign.
function prism_determinant(v,r,s,t)
    lower_r=subtract(v[2],v[1]);lower_s=subtract(v[3],v[1])
    upper_r=subtract(v[5],v[4]);upper_s=subtract(v[6],v[4])
    dr=add(scale(lower_r,1-t),scale(upper_r,t))
    ds=add(scale(lower_s,1-t),scale(upper_s,t))
    dt=add(subtract(v[4],v[1]),add(scale(subtract(upper_r,lower_r),r),
        scale(subtract(upper_s,lower_s),s)))
    return determinant(dr,ds,dt)
end

function pyramid_factor(v,u,w)
    base=add(add(scale(v[1],(1-u)*(1-w)),scale(v[2],u*(1-w))),
        add(scale(v[3],u*w),scale(v[4],(1-u)*w)))
    du=add(scale(subtract(v[2],v[1]),1-w),scale(subtract(v[3],v[4]),w))
    dw=add(scale(subtract(v[4],v[1]),1-u),scale(subtract(v[3],v[2]),u))
    return determinant(du,dw,subtract(v[5],base))
end

function full_map_certificate(msh,v,orientation)
    if msh==4
        @test orientation*tet_determinant(v)>0
    elseif msh==6
        for (r,s) in ((Q(0),Q(0)),(Q(1),Q(0)),(Q(0),Q(1)))
            b0=orientation*prism_determinant(v,r,s,Q(0))
            b2=orientation*prism_determinant(v,r,s,Q(1))
            h=4orientation*prism_determinant(v,r,s,Q(1,2))-b0-b2
            @test b0>0 && b2>0
            @test h>=0 || 4b0*b2>h*h
        end
    else
        factors=ntuple(index->begin
            u=(index-1)%2;w=(index-1)÷2
            orientation*pyramid_factor(v,Q(u),Q(w))
        end,4)
        @test all(>(0),factors)
        # The collapsed determinant factor is bilinear on the whole base.
        # Midpoint identities independently verify the asserted degree on
        # these actual cells; positive corner factors prove the open map.
        for u in (Q(0),Q(1,2),Q(1)),w in (Q(0),Q(1,2),Q(1))
            expected=(1-u)*(1-w)*factors[1]+u*(1-w)*factors[2]+
                (1-u)*w*factors[3]+u*w*factors[4]
            @test orientation*pyramid_factor(v,u,w)==expected
        end
    end
end

support(face)=Tuple(sort!(collect(face)))
function direction(face)
    forward=minimum(Tuple(circshift(collect(face),shift)) for shift in 0:length(face)-1)
    backward=minimum(Tuple(circshift(reverse(collect(face)),shift)) for shift in 0:length(face)-1)
    return forward<backward ? 1 : -1
end

function external_boundary(states)
    result=Set{Tuple}(support(cap) for cap in CAPS)
    for (face,state) in zip(LATERALS,states)
        a,b,c,d=face
        if state==0
            push!(result,support(face))
        elseif state==1
            push!(result,support((a,b,c)));push!(result,support((a,c,d)))
        else
            push!(result,support((a,b,d)));push!(result,support((b,c,d)))
        end
    end
    return result
end

function complex_certificate(template)
    counts=Dict{Tuple,Int}();directions=Dict{Tuple,Int}();typed=Set{Tuple}()
    total=Q(0)
    for index in 1:Int(template.ncells)
        cell=template.cells[index];msh=Int(cell.msh);width=WIDTH[msh]
        ids=Int.(cell.nodes[1:width])
        @test all(1<=id<=6 for id in ids)
        @test length(unique(ids))==width
        @test all(==(0),cell.nodes[width+1:8])
        key=(msh,support(ids))
        @test key ∉ typed
        push!(typed,key)
        v=Tuple(point(UNIT[id]) for id in ids)
        full_map_certificate(msh,v,1)
        center=ntuple(axis->sum(p[axis] for p in v)/width,3)
        volume=Q(0)
        for face in FACES[msh]
            global_face=Tuple(ids[local_id] for local_id in face)
            key=support(global_face)
            counts[key]=get(counts,key,0)+1
            directions[key]=get(directions,key,0)+direction(global_face)
            for local_id in 2:length(face)-1
                a,b,c=v[face[1]],v[face[local_id]],v[face[local_id+1]]
                # Outward oriented face has negative determinant to center.
                value=determinant(subtract(b,a),subtract(c,a),subtract(center,a))
                @test value<0
                volume-=value/6
            end
            if length(face)==4
                a,b,c,d=(v[local_id] for local_id in face)
                @test determinant(subtract(b,a),subtract(d,a),subtract(center,a))<0
                @test determinant(subtract(c,b),subtract(d,b),subtract(center,b))<0
            end
        end
        @test volume>0
        total+=volume
    end
    @test total==Q(1,2)
    @test all(count in (1,2) for count in values(counts))
    @test all(counts[key]==1 || directions[key]==0 for key in keys(counts))
    @test Set(key for (key,count) in counts if count==1)==external_boundary(template.faces)
    for index in Int(template.ncells)+1:3
        @test template.cells[index].msh==0
        @test all(iszero,template.cells[index].nodes)
    end
    return typed
end

function empty_sweep()
    return Model._ExtrudeVolumeSweep(1,true,Matrix{NTuple{3,Float64}}(undef,0,0),
        NTuple{6,NTuple{3,Float64}}[],NTuple{4,NTuple{3,Float64}}[],
        NTuple{8,NTuple{3,Float64}}[],NTuple{6,NTuple{3,Float64}}[],
        NTuple{5,NTuple{3,Float64}}[],NTuple{3,Float64}[],Bool[false])
end

function emitter_certificate(template,expected)
    sweep=empty_sweep()
    Model._extrude_nonew_emit_template!(sweep,UNIT,template)
    actual=Set{Tuple}()
    for (msh,cells) in ((4,sweep.tets),(6,sweep.prisms),(7,sweep.pyramids))
        for cell in cells
            ids=Tuple(findfirst(==(p),UNIT) for p in cell)
            push!(actual,(msh,support(ids)))
        end
    end
    @test actual==expected
    @test isempty(sweep.hexes)
    @test isempty(sweep.interior)
    @test !sweep.degenerate[]
end

function warmed_emit_bytes(template)
    sweep=empty_sweep()
    for cells in (sweep.tets,sweep.prisms,sweep.pyramids)
        sizehint!(cells,3)
    end
    Model._extrude_nonew_emit_template!(sweep,UNIT,template)
    for cells in (sweep.tets,sweep.prisms,sweep.pyramids)
        empty!(cells)
    end
    return @allocated Model._extrude_nonew_emit_template!(sweep,UNIT,template)
end

# Independent logical predicate: an apex must lie on both diagonals of its
# incident lateral faces. The opposite face remains free; only000 permits a
# whole prism. This predicts all27 masks without calling either factory.
function mask_possible(states)
    all(iszero,states) && return true
    for apex in 1:6
        possible=true
        for (face,state) in zip(LATERALS,states)
            apex in face || continue
            diagonal=state==1 ? (face[1],face[3]) :
                state==2 ? (face[2],face[4]) : (0,0)
            if apex ∉ diagonal
                possible=false;break
            end
        end
        possible && return true
    end
    return false
end

function contains_box(value)
    value===Core.Box && return true
    value isa GlobalRef && return value.mod===Core && value.name===:Box
    value isa QuoteNode && return contains_box(value.value)
    value isa Expr && return any(contains_box,value.args)
    value isa Core.CodeInfo && return any(contains_box,value.code)
    return false
end

function deliberate_box(flag)
    value=0
    captured=()->value
    value=flag ? 1 : 2
    return captured()
end

@testset "NoNew true-prism immutable complete template relation" begin
    templates=Model._extrude_nonew_prism_templates()
    hexes=Model._extrude_nonew_templates()
    hex_before=sha256(reinterpret(UInt8,hexes))
    @test length(hexes)==315
    @test templates isa NTuple{13,Model._ExtrudeNoNewPrismTemplate}
    @test isbitstype(Model._ExtrudeNoNewPrismTemplate)
    @test isbitstype(typeof(templates))
    @test sizeof(Model._ExtrudeNoNewPrismTemplate)==33
    @test sizeof(templates)==429
    @test Tuple(Int.(template.faces) for template in templates)==RANKED_MASKS
    @test templates==Model._extrude_nonew_prism_build_templates()
    @test count(template->template.kind==0,templates)==1
    @test Set(Int(template.choice) for template in templates if template.kind==1)==Set(1:6)
    unique_keys=Set{Tuple}()
    for template in templates
        @test 1<=template.ncells<=3
        @test all(state in 0:2 for state in template.faces)
        typed=complex_certificate(template)
        key=(template.faces,Tuple(sort!(collect(typed))))
        @test key ∉ unique_keys
        push!(unique_keys,key)
        emitter_certificate(template,typed)
        warmed_emit_bytes(template)
        @test warmed_emit_bytes(template)==0
    end
    masks=Set(Tuple(Int.(template.faces)) for template in templates)
    for a in 0:2,b in 0:2,c in 0:2
        states=(a,b,c)
        @test (states in masks)==mask_possible(states)
    end
    @test length(masks)==13
    @test (1,1,1) ∉ masks && (2,2,2) ∉ masks
    @test count(template->template.ncells==1,templates)==1
    @test count(template->template.ncells==2,templates)==6
    @test count(template->template.ncells==3,templates)==6
    @test contains_box(Base.uncompressed_ast(only(methods(deliberate_box))))
    boxed=Symbol[]
    for name in names(Model;all=true)
        startswith(String(name),"_extrude_nonew_prism_") || continue
        value=getfield(Model,name)
        value isa Function || continue
        for method in methods(value)
            contains_box(Base.uncompressed_ast(method)) && push!(boxed,name)
        end
    end
    for method in methods(Model._extrude_nonew_emit_template!)
        contains_box(Base.uncompressed_ast(method)) && push!(boxed,:_extrude_nonew_emit_template!)
    end
    @test isempty(boxed)
    @test hex_before==sha256(reinterpret(UInt8,Model._extrude_nonew_templates()))
    println("PRISM_TEMPLATE_RELATION_OK templates=13 masks=13/27 bytes=429 boxed=",boxed)
end

@testset "Prism factory actual maps retain affine orientation and graded heights" begin
    fixtures=(
        (ntuple(i->(UNIT[i][1],UNIT[i][2],UNIT[i][3]/8),6),1),
        (ntuple(i->(-UNIT[i][1],UNIT[i][2],UNIT[i][3]),6),-1),
        ((UNIT[1],UNIT[3],UNIT[2],UNIT[4],UNIT[6],UNIT[5]),-1),
        (ntuple(i->(UNIT[i][1]+UNIT[i][3]/4,UNIT[i][2]+UNIT[i][3]/8,
            UNIT[i][3]/2),6),1))
    for (corners,orientation) in fixtures,template in Model._extrude_nonew_prism_templates()
        for index in 1:Int(template.ncells)
            cell=template.cells[index];width=WIDTH[Int(cell.msh)]
            v=ntuple(k->point(corners[cell.nodes[k]]),width)
            full_map_certificate(Int(cell.msh),v,orientation)
        end
    end
end

end
