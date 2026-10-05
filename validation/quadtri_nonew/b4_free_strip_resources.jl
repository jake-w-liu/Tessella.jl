# Normal compilation and --check-bounds=yes; no Gmsh dependency.
# Actual free M5 strip, N500/1000/2000, separately charged entry stages.
# Three actual M-growth rows are measured separately; three small P2 audits are
# unmeasured query/proof checks.
# Allocations/times below charge construction only; geometry/query proofs and
# public-volume extraction are outside @timed.
using Tessella, Test, SHA, TOML
using Tessella.Model: mesh_model_volume, mesh_model_surface, model_to_mixed
using Tessella.Elements: ElementBlock, MixedMesh, msh_dimension
include(joinpath(@__DIR__, "public_label_resources.jl"))
include(joinpath(@__DIR__, "../../test/geometry/quadtri_nonew_b4_rec_strip_certificates.jl"))
const RectGridBase = QuadTriNoNewB4RecStripCertificates
const B4Cert = QuadTriNoNewB4RecStripCertificates
const RectGridQuad = QuadTriNoNewCertificates
const RectGridModel = Tessella.Model
const FreeB4 = RectGridModel._ExtrudeNoNewB4Free
const B4_FREE_STRIP_ROOT = normpath(joinpath(@__DIR__, "../.."))
b4_free_strip_inputs() = sort!(vcat([joinpath(B4_FREE_STRIP_ROOT,"Project.toml")],
    [joinpath(directory,file) for (directory,_,files) in walkdir(joinpath(B4_FREE_STRIP_ROOT,"src"))
     for file in files if endswith(file,".jl")], [@__FILE__,joinpath(@__DIR__,"public_label_resources.jl")],
    [joinpath(B4_FREE_STRIP_ROOT,"test/geometry",file) for file in
        ("quadtri_nonew_b4_rec_strip_certificates.jl","quadtri_nonew_rect_grid_certificates.jl",
         "quadtri_nonew_quad_strip_certificates.jl","quadtri_nonew_quad_patch_certificates.jl",
         "quadtri_nonew_certificates.jl")]))
b4_free_strip_snapshot() = Dict(relpath(path,B4_FREE_STRIP_ROOT)=>bytes2hex(sha256(read(path))) for path in b4_free_strip_inputs())
const B4_FREE_STRIP_START = b4_free_strip_snapshot()
println("START B4_FREE_STRIP_RESOURCE Julia=",VERSION," inputs=",length(B4_FREE_STRIP_START))
flush(stdout)

function b4_free_strip_boxes(value)
    value isa Core.Box && return 1
    value isa GlobalRef && return value==GlobalRef(Core,:Box) ? 1 : 0
    value isa Core.CodeInfo && return sum(b4_free_strip_boxes,value.code;init=0)
    value isa Expr && return sum(b4_free_strip_boxes,value.args;init=0)
    value isa Core.NewvarNode && return b4_free_strip_boxes(value.slot)
    value isa QuoteNode && return b4_free_strip_boxes(value.value)
    return value===Core.Box ? 1 : 0
end

function b4_free_strip_box_control()
    captured=0
    read=()->captured
    captured=1
    return read
end

function b4_free_strip_keyword_box_control(value=1;initial=0)
    captured=initial
    read=()->captured
    captured=value
    return read
end

function b4_free_strip_lowered_methods(method::Method)
    result=Method[method]
    # Optional positional wrappers can lack Base.bodyfunction metadata.
    # Generated keyword bindings retain the complete implementation bodies.
    prefix="#"*String(method.name)*"#"
    for name in names(method.module;all=true)
        startswith(String(name),prefix) || continue
        value=getfield(method.module,name)
        if value isa Function
            append!(result,methods(value))
        elseif value isa Type && value<:Function
            matches=Base._methods_by_ftype(Tuple{value,Vararg{Any}},-1,Base.get_world_counter())
            matches===nothing || foreach(match->push!(result,match.method),matches)
        end
    end
    return result
end

function b4_free_strip_module_methods(mod)
    checked=Set{Method}()
    for name in names(mod;all=true)
        value=getfield(mod,name)
        if value isa Function
            for method in methods(value)
                method.module===mod || continue
                union!(checked,b4_free_strip_lowered_methods(method))
            end
        elseif value isa Type && value<:Function
            # Closure-call methods are actual implementation bodies too.
            matches=Base._methods_by_ftype(Tuple{value,Vararg{Any}},-1,Base.get_world_counter())
            matches===nothing || foreach(match->push!(checked,match.method),matches)
        end
    end
    return checked
end

function b4_free_strip_source_copy(source;kwargs...)
    fields=fieldnames(typeof(source))
    data=NamedTuple{fields}(Tuple(getfield(source,name) for name in fields))
    changed=merge(data,(;kwargs...))
    return typeof(source)(Tuple(getfield(changed,name) for name in fields)...)
end

# Saved accepted rank-order witnesses. These are real Source inputs, with
# native certification and physical propagation; no problem flag is injected.
const B4_FREE_RETAINED_WITNESSES=(
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

function b4_free_strip_retained_witness(witness)
    m=witness.m;n=witness.n
    directions=Dict(minmax(pair...)=>(pair[1],pair[2]) for pair in witness.directions)
    chains=(collect(1:m+1),[m+1,2m+2],collect(2m+2:-1:m+2),[m+2,1])
    reversals=0
    for curve in 1:4
        same=[directions[minmax(chains[curve][k],chains[curve][k+1])]==
            (chains[curve][k],chains[curve][k+1]) for k in 1:length(chains[curve])-1]
        @test all(same) || all(!,same)
        all(same) || (reversals|=1<<(curve-1))
    end
    f=B4Cert.fixture(witness.name;strip_length=m,laterals=false,layers=:one,
        pins=(1,2,3,4),curve_reverse=reversals,tags=:sparse)
    f=merge(f,(source=replace(f.source,"Layers{1}"=>"Layers{$n}"),
        levels=collect(range(0.,1.;length=n+1)),intervals=n))
    execution=B4Cert.execute(f.source;dim=0);model=execution.model
    original=mesh_model_surface(model,f.surface)
    independent=B4Cert.source_complex(original,f)
    canonical=vcat(collect(independent.chains[1]),reverse(collect(independent.chains[3])))
    permutation=collect(witness.permutation);inverse=invperm(permutation)
    coordinates=original.coords[:,canonical[permutation]]
    cells=hcat([Int32[inverse[node] for node in witness.cells[cell]] for cell in witness.order]...)
    source_mesh=MixedMesh(coordinates,[ElementBlock(3,cells)])
    # Both independent geometry and production Source admission see the
    # reordered actual source; they do not share the planner's decisions.
    independent=B4Cert.source_complex(source_mesh,f)
    @test length(independent.cells)==m
    volume=Int(execution.lists["sweep"][2]);spec=model.meshing.extrude_specs[(3,volume)]
    source=RectGridModel._extrude_nonew_rect_grid_source(model,f.surface,source_mesh,spec,
        "free resource witness";mode=:b4_strip)
    refs=fill((Int32(1),Int32(1)),n)
    completed=FreeB4.plan(source,f.levels,refs,"free resource witness")
    expected=[(Int32(row),Int32(interval)) for (row,cell) in enumerate(witness.order)
        for interval in 1:n if (cell,interval) in witness.problems]
    @test completed.problem_positions==expected
    @test completed.cell_counts==witness.counts
    params=RectGridModel._extrude_gate(model,3,volume)
    cols=RectGridModel._extrude_volume_columns(model,volume,f.surface,source_mesh,params,spec,
        f.levels,"free resource witness")
    result=FreeB4.emit(completed,cols,spec,"free resource witness")
    audit=b4_free_strip_audit(result.volume,source_mesh,n,false;grid=(m,1))
    @test audit.centers==length(witness.problems) && audit.cell_counts==collect(witness.counts)
    @test audit.typed_faces==completed.face_capacity
    @test witness.name=="nonterminal_center" ? audit.nonterminal_centers==1 : audit.nonterminal_centers==0
    return Dict{String,Any}("name"=>witness.name,"source_cells"=>m,"layers"=>n,
        "scope"=>"actual_certified_source_with_saved_rank_order",
        (String(key)=>value for (key,value) in pairs(audit))...)
end

function b4_free_strip_rejection_controls()
    caller="free resource guards"
    @test FreeB4.preflight((5,1),caller)==(12,60,0,5)
    for shape in ((4,1),(2,2),(typemax(Int),1))
        @test_throws ArgumentError FreeB4.preflight(shape,caller)
    end
    f=B4Cert.fixture("resource_guards";strip_length=5,laterals=false,layers=:one,pins=(1,2,3,4))
    execution=B4Cert.execute(f.source;dim=0);model=execution.model
    mesh=mesh_model_surface(model,f.surface);volume=Int(execution.lists["sweep"][2])
    spec=model.meshing.extrude_specs[(3,volume)]
    source=RectGridModel._extrude_nonew_rect_grid_source(model,f.surface,mesh,spec,caller;mode=:b4_strip)
    refs=[(Int32(1),Int32(1))]
    for levels in ([0.,.5],[0.,NaN,1.],[0.,0.,1.],[1.,2.])
        actual_refs=fill((Int32(1),Int32(1)),length(levels)-1)
        @test_throws ArgumentError FreeB4.plan(source,levels,actual_refs,caller)
    end
    @test_throws ArgumentError FreeB4.plan(source,[0.,1.],NTuple{2,Int32}[],caller)
    for intervals in (0,-1,typemax(Int))
        @test_throws ArgumentError FreeB4.SourcePhase.Phase(source,intervals)
    end
    bad_edges=copy(source.edges);bad_edges[1]=reverse(bad_edges[1])
    @test_throws ArgumentError FreeB4.SourcePhase.Phase(b4_free_strip_source_copy(source;edges=bad_edges),1)
    bad_category=copy(source.category);bad_category[1]=0x03
    @test_throws ArgumentError FreeB4.SourcePhase.Phase(b4_free_strip_source_copy(source;category=bad_category),1)
    @test_throws ArgumentError FreeB4.SourcePhase.Phase(b4_free_strip_source_copy(source;source_cells=source.source_cells[1:end-1]),1)
    bad_indices=copy(source.edge_indices);bad_indices[1]=(Int32(0),bad_indices[1][2:4]...)
    @test_throws ArgumentError FreeB4.SourcePhase.Phase(b4_free_strip_source_copy(source;edge_indices=bad_indices),1)
    coordinates=copy(mesh.coords);coordinates[:,2]=coordinates[:,1]
    malformed=MixedMesh(coordinates,[ElementBlock(3,copy(only(mesh.blocks).nodes))])
    @test_throws ArgumentError B4Cert.source_complex(malformed,f)
    @test_throws ArgumentError RectGridModel._extrude_nonew_rect_grid_source(model,f.surface,malformed,spec,caller;mode=:b4_strip)
    repeated=copy(only(mesh.blocks).nodes);repeated[2,1]=repeated[1,1]
    malformed=MixedMesh(copy(mesh.coords),[ElementBlock(3,repeated)])
    @test_throws ArgumentError B4Cert.source_complex(malformed,f)
    @test_throws ArgumentError RectGridModel._extrude_nonew_rect_grid_source(model,f.surface,malformed,spec,caller;mode=:b4_strip)
    # Empty/foreign masks must give a diagnostic, never an empty mesh.
    @test_throws ArgumentError FreeB4.FinalFactories.factory_index((0x03,0x00,0x00,0x00,0x00,0x00))
    @test_throws ArgumentError FreeB4.CenterFan.fan((0x00,0x00,0x00,0x00,0x00,0x03))
end

# Outward interval arithmetic encloses the exact arithmetic of the represented
# Float64 corners. The resource fixtures stay at ordinary finite unit scales.
# A positive lower endpoint is a whole-map proof, not a point-sample criterion.
b4_free_strip_interval(x::Float64)=(x,x)
b4_free_strip_add(a,b)=(prevfloat(a[1]+b[1]),nextfloat(a[2]+b[2]))
b4_free_strip_sub(a,b)=(prevfloat(a[1]-b[2]),nextfloat(a[2]-b[1]))
function b4_free_strip_mul(a,b)
    products=(a[1]*b[1],a[1]*b[2],a[2]*b[1],a[2]*b[2])
    return (prevfloat(minimum(products)),nextfloat(maximum(products)))
end
b4_free_strip_scale(a,x::Float64)=b4_free_strip_mul(a,b4_free_strip_interval(x))
b4_free_strip_difference(p,q)=ntuple(k->b4_free_strip_sub(b4_free_strip_interval(p[k]),b4_free_strip_interval(q[k])),3)
b4_free_strip_average(a,b)=ntuple(k->b4_free_strip_scale(b4_free_strip_add(a[k],b[k]),.5),3)
function b4_free_strip_determinant(a,b,c)
    first=b4_free_strip_mul(a[1],b4_free_strip_sub(b4_free_strip_mul(b[2],c[3]),b4_free_strip_mul(b[3],c[2])))
    second=b4_free_strip_mul(a[2],b4_free_strip_sub(b4_free_strip_mul(b[1],c[3]),b4_free_strip_mul(b[3],c[1])))
    third=b4_free_strip_mul(a[3],b4_free_strip_sub(b4_free_strip_mul(b[1],c[2]),b4_free_strip_mul(b[2],c[1])))
    return b4_free_strip_add(b4_free_strip_sub(first,second),third)
end

# Differentiate the reference interpolation directly. A translated Hex8 has
# a bilinear base determinant; a collapsed Pyramid5 has four base coefficients.
# A Prism6 is affine on its triangle and quadratic along its interval: three
# positive Bernstein coefficients at each triangle vertex certify that domain.
function b4_free_strip_map(points,msh)
    if msh==4
        bound=b4_free_strip_determinant(b4_free_strip_difference(points[2],points[1]),
            b4_free_strip_difference(points[3],points[1]),b4_free_strip_difference(points[4],points[1]))
        @test bound[1]>0
    elseif msh==5
        delta=points[5].-points[1]
        @test all(index->points[index+4].-points[index]==delta,1:4)
        for index in 1:4
            bound=b4_free_strip_determinant(b4_free_strip_difference(points[mod1(index+1,4)],points[index]),
                b4_free_strip_difference(points[mod1(index-1,4)],points[index]),b4_free_strip_difference(points[5],points[1]))
            @test bound[1]>0
        end
    elseif msh==6
        low_u=b4_free_strip_difference(points[2],points[1]);high_u=b4_free_strip_difference(points[5],points[4])
        low_v=b4_free_strip_difference(points[3],points[1]);high_v=b4_free_strip_difference(points[6],points[4])
        for vertex in 1:3
            direction=b4_free_strip_difference(points[vertex+3],points[vertex])
            first=b4_free_strip_determinant(low_u,low_v,direction)
            last=b4_free_strip_determinant(high_u,high_v,direction)
            middle=b4_free_strip_determinant(b4_free_strip_average(low_u,high_u),b4_free_strip_average(low_v,high_v),direction)
            coefficient=b4_free_strip_sub(b4_free_strip_scale(middle,2.),b4_free_strip_scale(b4_free_strip_add(first,last),.5))
            @test first[1]>0 && last[1]>0 && coefficient[1]>0
        end
    elseif msh==7
        for index in 1:4
            bound=b4_free_strip_determinant(b4_free_strip_difference(points[mod1(index+1,4)],points[index]),
                b4_free_strip_difference(points[mod1(index-1,4)],points[index]),b4_free_strip_difference(points[5],points[index]))
            @test bound[1]>0
        end
    else
        @test false # A foreign family violates this declared product.
    end
end

function b4_free_strip_audit(volume,source,layers,recombined;grid=(5,1))
    a,b=grid;vertices=(a+1)*(b+1);macros=a*b;boundary=2(a+b)
    @test b==1 && a>=5 && !recombined
    primary=vertices*(layers+1)
    @test primary<=size(volume.coords,2)<=primary+macros*layers
    @test all(isfinite,volume.coords)
    source_cells=Tuple(Tuple(Int.(cell)) for block in source.blocks
        if block isa ElementBlock && block.msh==3 for cell in eachcol(block.nodes))
    @test length(source_cells)==macros && size(source.coords)==(3,vertices)
    points=[ntuple(d->source.coords[d,node],3) for node in 1:vertices]
    columns=Dict((p[1],p[2])=>node for (node,p) in enumerate(points))
    @test length(columns)==vertices
    planes=sort!(unique(volume.coords[3,node] for node in axes(volume.coords,2)
        if haskey(columns,(volume.coords[1,node],volume.coords[2,node]))))
    @test length(planes)==layers+1 && first(planes)==0. && last(planes)==1.
    heights=Dict(height=>level for (level,height) in enumerate(planes))
    logical=fill((0,0),size(volume.coords,2))
    grid=Set{Tuple{Int,Int}}();extra=Int[]
    for node in axes(volume.coords,2)
        xy=(volume.coords[1,node],volume.coords[2,node])
        if haskey(columns,xy)
            @test haskey(heights,volume.coords[3,node])
            logical[node]=(columns[xy],heights[volume.coords[3,node]])
            @test !(logical[node] in grid)
            push!(grid,logical[node])
        else
            push!(extra,node)
        end
    end
    @test length(grid)==primary
    # At most four parents touch a source vertex. Intersect actual column
    # incidence, rather than scan F parents for every output cell.
    parents=[Int[] for _ in 1:vertices]
    for (parent,cell) in enumerate(source_cells),node in cell
        push!(parents[node],parent)
    end
    source_edges=Dict{Tuple{Int,Int},Int}()
    for cell in source_cells,k in 1:4
        edge=minmax(cell[k],cell[mod1(k+1,4)])
        source_edges[edge]=get(source_edges,edge,0)+1
    end
    @test all(value->value in (1,2),values(source_edges))
    exterior_edges=Set(edge for (edge,uses) in source_edges if uses==1)
    @test length(exterior_edges)==boundary
    areas=Tuple(begin
        origin=points[first(cell)]
        abs(sum((points[cell[index]][1]-origin[1])*(points[cell[mod1(index+1,4)]][2]-origin[2])-
            (points[cell[index]][2]-origin[2])*(points[cell[mod1(index+1,4)]][1]-origin[1]) for index in 1:4))/2
    end for cell in source_cells)
    # The represented source is a complete unit-area disk. Summing thousands
    # of separately rounded Float64 cell areas can exceed a fixed few-ulp
    # budget, so check this global cancellation in exact dyadic arithmetic.
    exact_area=sum(begin
        origin=Tuple(Rational{BigInt}(value) for value in points[first(cell)])
        p=Tuple(Tuple(Rational{BigInt}(value) for value in points[node]) for node in cell)
        abs(sum((p[k][1]-origin[1])*(p[mod1(k+1,4)][2]-origin[2])-
            (p[k][2]-origin[2])*(p[mod1(k+1,4)][1]-origin[1]) for k in 1:4))/2
    end for cell in source_cells)
    @test all(>(0),areas) && exact_area==1
    # The represented one-row source has two strictly ordered boundary chains.
    # Discover every actual center and its actual interval. No terminal-center
    # count or position from the recombined policy participates in this audit.
    parent_order=sortperm([minimum(points[node][1] for node in cell) for cell in source_cells])
    lower=[minimum(points[node][1] for node in source_cells[parent]) for parent in parent_order]
    upper=[maximum(points[node][1] for node in source_cells[parent]) for parent in parent_order]
    @test first(lower)==0. && last(upper)==1.
    # Native opposing chain samples can differ in their represented X values.
    # Their tiny X bounding-box overlap is not geometric overlap. At y=0 and
    # y=1 the exact shared IDs and strict X order establish disjoint strips
    # between affine joining edges over the entire represented unit domain.
    bottom=Vector{Int}[];top=Vector{Int}[]
    for parent in parent_order
        cell=source_cells[parent]
        low=sort!([node for node in cell if points[node][2]==0.];by=node->points[node][1])
        high=sort!([node for node in cell if points[node][2]==1.];by=node->points[node][1])
        @test length(low)==length(high)==2
        @test points[low[1]][1]<points[low[2]][1] && points[high[1]][1]<points[high[2]][1]
        push!(bottom,low);push!(top,high)
    end
    @test all(bottom[index][2]==bottom[index+1][1] && top[index][2]==top[index+1][1]
        for index in 1:macros-1)
    center_parent=zeros(Int,size(volume.coords,2));center_level=zeros(Int,size(volume.coords,2))
    locations=Set{Tuple{Int,Int}}()
    for node in extra
        p=ntuple(d->volume.coords[d,node],3)
        slot=searchsortedlast(lower,p[1]);@test 1<=slot<=macros
        parent=parent_order[slot];cell=source_cells[parent]
        level=searchsortedlast(planes,p[3]);@test 1<=level<=layers
        @test lower[slot]<p[1]<upper[slot] && planes[level]<p[3]<planes[level+1]
        for side in 1:4
            x,y=points[cell[side]],points[cell[mod1(side+1,4)]]
            @test (y[1]-x[1])*(p[2]-x[2])-(y[2]-x[2])*(p[1]-x[1])>0
        end
        for d in 1:3
            total=(0.,0.)
            for height in planes[level:level+1],corner in cell
                value=d==3 ? height : points[corner][d]
                total=b4_free_strip_add(total,b4_free_strip_interval(value))
            end
            mean=b4_free_strip_scale(total,.125)
            @test mean[1]<=p[d]<=mean[2]
        end
        @test !((parent,level) in locations)
        push!(locations,(parent,level));center_parent[node]=parent;center_level[node]=level
    end
    domain_counts=zeros(Int,macros,layers);domain_volumes=zeros(Float64,macros,layers)
    families=Dict{Int,Int}();faces=Dict{Tuple,Tuple{Int,Tuple}}()
    used=falses(size(volume.coords,2));signatures=Tuple[];actual_edges=Set{Tuple}()
    for block in volume.blocks
        block isa ElementBlock && msh_dimension(block.msh)==3 || continue
        @test block.msh in (4,5,6,7)
        families[Int(block.msh)]=get(families,Int(block.msh),0)+size(block.nodes,2)
        for cell in eachcol(block.nodes)
            @test length(unique(cell))==length(cell)
            used[cell].=true
            body=[node for node in cell if center_parent[node]!=0]
            ids=unique(logical[node][1] for node in cell if logical[node][1]!=0)
            if isempty(body)
                candidates=copy(parents[first(ids)])
                for node in ids;intersect!(candidates,parents[node]);end
                @test length(candidates)==1
                parent=only(candidates)
                levels=[logical[node][2] for node in cell];level=minimum(levels)
                @test maximum(levels)==level+1
            else
                @test length(body)==1 && block.msh in (4,7)
                parent=center_parent[only(body)];level=center_level[only(body)]
                @test all(node->node in source_cells[parent],ids)
                @test all(node->logical[node][2] in (level,level+1),cell[cell.!=only(body)])
            end
            domain_counts[parent,level]+=1
            p=Tuple(ntuple(d->volume.coords[d,node],3) for node in cell)
            b4_free_strip_map(p,Int(block.msh))
            # Bilinear quadrangle flux integrates the represented curved face.
            # A fixed hidden diagonal would mismeasure warped Prism6 faces.
            signed=sum(RectGridQuad.face_volume(volume,Tuple(cell[index] for index in face),p[1])
                for face in RectGridQuad.FACES[Int(block.msh)])
            @test isfinite(signed) && signed>0
            domain_volumes[parent,level]+=signed
            for pattern in RectGridQuad.FACES[Int(block.msh)]
                face=Tuple(Int(cell[index]) for index in pattern)
                key=Tuple(sort!(collect(face)));cycle=RectGridBase.cycle(face)
                previous=get(faces,key,nothing)
                if previous===nothing
                    faces[key]=(1,cycle)
                else
                    @test previous[1]==1 && previous[2]==RectGridBase.cycle(reverse(face))
                    faces[key]=(2,previous[2])
                end
            end
            for pattern in RectGridQuad.FACES[Int(block.msh)]
                face=Tuple(Int(cell[index]) for index in pattern)
                for k in eachindex(face)
                    push!(actual_edges,minmax(face[k],face[mod1(k+1,length(face))]))
                end
            end
            push!(signatures,(Int(block.msh),Tuple(sort!(collect(p)))))
        end
    end
    @test all(used)
    @test allunique(signatures)
    @test sum(values(families))<=12macros*layers
    @test length(used)-length(actual_edges)+length(faces)-sum(values(families))==1
    exterior=(recombined ? boundary : 2boundary)*layers+3macros
    @test 2length(faces)==sum(length(RectGridQuad.FACES[msh])*count for (msh,count) in families)+exterior
    groups=Dict{Tuple,Vector{Tuple}}()
    for (uses,cycle) in values(faces)
        uses==1 || continue
        @test all(node->center_parent[node]==0,cycle)
        ids=unique(logical[node][1] for node in cycle)
        levels=unique(logical[node][2] for node in cycle)
        if length(levels)==1
            level=only(levels)
            @test level in (1,layers+1)
            candidates=copy(parents[first(ids)])
            for node in ids;intersect!(candidates,parents[node]);end
            @test length(candidates)==1
            carrier=(:cap,level,only(candidates))
        else
            @test length(ids)==2 && length(levels)==2 && maximum(levels)==minimum(levels)+1
            edge=minmax(ids...)
            @test edge in exterior_edges
            carrier=(:side,edge,minimum(levels))
        end
        push!(get!(groups,carrier,Tuple[]),cycle)
    end
    @test length(groups)==2macros+boundary*layers
    for (carrier,cycles) in groups
        if carrier[1]===:cap
            _,level,parent=carrier
            @test length(cycles)==(level==1 ? 1 : 2)
            @test all(face->length(face)==(level==1 ? 4 : 3),cycles)
            @test Set(logical[node] for face in cycles for node in face)==
                Set((node,level) for node in source_cells[parent])
        else
            _,edge,level=carrier
            @test length(cycles)==(recombined ? 1 : 2)
            @test all(face->length(face)==(recombined ? 4 : 3),cycles)
            @test Set(logical[node] for face in cycles for node in face)==
                Set((node,height) for node in edge for height in (level,level+1))
        end
    end
    for parent in 1:macros,level in 1:layers
        @test 1<=domain_counts[parent,level]<=12
        @test isapprox(domain_volumes[parent,level],areas[parent]*(planes[level+1]-planes[level]);rtol=64eps(),atol=0)
    end
    @test count(value->value[1]==1,values(faces))==exterior
    @test isapprox(sum(domain_volumes),1.;rtol=128eps())
    io=IOBuffer();show(io,sort!(signatures))
    typed=IOBuffer()
    for block in volume.blocks
        msh_dimension(block.msh)==3 || continue
        write(typed,Int32(block.msh));write(typed,Int64(size(block.nodes,2)));write(typed,vec(block.nodes))
    end
    return (;digest=bytes2hex(sha256(take!(io))),
        coordinate_sha256=bytes2hex(sha256(reinterpret(UInt8,vec(volume.coords)))),
        typed_connectivity_sha256=bytes2hex(sha256(take!(typed))),
        nodes=size(volume.coords,2),centers=length(extra),
        nonterminal_centers=count(location->location[2]<layers,locations),
        cell_counts=[get(families,msh,0) for msh in (4,5,6,7)],
        cells=sum(values(families)),typed_faces=length(faces),edges=length(actual_edges),boundary=exterior)
end

function b4_free_strip_profile(make,extract,source,layers,recombined;grid=(5,1))
    result=make()
    samples=[@timed(make()) for _ in 1:3]
    selected=argmin(sample->sample.bytes,samples)
    allocated=selected.bytes
    allocation_count=Base.gc_alloc_count(selected.gcstats)
    seconds=minimum(sample.time for sample in samples)
    volume=extract(result)
    audit=b4_free_strip_audit(volume,source,layers,recombined;grid)
    return merge((;layers,allocated,allocation_count,seconds,
        allocation_samples=[sample.bytes for sample in samples],
        allocation_count_samples=[Base.gc_alloc_count(sample.gcstats) for sample in samples]),audit)
end

# The measured entry products are P1. Actual P2 support identity, carrier
# ownership and public map/parameter queries are audited separately below.
# Lower-cell carriers come from the actual classified P1 projection, rather
# than an endpoint-owner heuristic or a guessed subdivision of CAD faces.
b4_free_strip_support_key(nodes)=Tuple(sort!(Int.(collect(nodes))))

function b4_free_strip_actual_carriers(api,projected)
    tags,xyz,_=api.mesh.get_nodes(3,1,true)
    points=reshape(xyz,3,:)
    positions=Dict(Tuple(p)=>tag for (tag,p) in zip(tags,eachcol(points)))
    @test length(positions)==length(tags)
    remap=[positions[Tuple(projected.coords[:,node])] for node in axes(projected.coords,2)]
    edges=Dict{Tuple,Tuple{Int,Int}}();quads=Dict{Tuple,Tuple{Int,Int}}()
    for (index,block) in enumerate(projected.blocks)
        dimension=msh_dimension(block.msh)
        dimension in (1,2) || continue
        owners=projected.elementary_entities[index]
        for (column,cell) in enumerate(eachcol(block.nodes))
            carrier=(dimension,Int(owners[column]))
            if dimension==1
                edges[b4_free_strip_support_key(remap[node] for node in cell)]=carrier
            else
                for i in eachindex(cell)
                    support=b4_free_strip_support_key((remap[cell[i]],remap[cell[mod1(i+1,length(cell))]]))
                    haskey(edges,support) && edges[support][1]<dimension && continue
                    @test !haskey(edges,support) || edges[support]==carrier
                    edges[support]=carrier
                end
                block.msh==3 && (quads[b4_free_strip_support_key(remap[node] for node in cell)]=carrier)
            end
        end
    end
    primary_owners=Dict(tag=>begin
        _,_,dimension,owner=api.mesh.get_node(tag);(Int(dimension),Int(owner))
    end for tag in tags)
    faces=Dict{Int,Vector{Tuple}}()
    for (index,block) in enumerate(projected.blocks)
        msh_dimension(block.msh)==2 || continue
        for (column,cell) in enumerate(eachcol(block.nodes))
            owner=Int(projected.elementary_entities[index][column])
            push!(get!(faces,owner,Tuple[]),Tuple(remap[node] for node in cell))
        end
    end
    return (;edges,quads,primary_owners,faces,remap)
end

function b4_free_strip_surface_query(api,surface,boundary,expected)
    # Native cache queries compute UV for all selected nodes. This deliberately
    # does not assert parity with Gmsh's partially empty stored P2 parameters.
    nodes,coordinates,parameters=api.mesh.get_nodes(2,surface,boundary,true)
    @test length(nodes)==expected && allunique(nodes)
    @test length(coordinates)==3expected && length(parameters)==2expected
    @test maximum(abs.(api.model.get_value(2,surface,parameters).-coordinates))<=2e-11
    return Set(nodes)
end

function b4_free_strip_quadratic_audit(api,projected,layers,recombined;grid=(5,1))
    @test !recombined
    expected=b4_free_strip_actual_carriers(api,projected)
    model=api.CURRENT[]
    surfaces=RectGridModel._model_volume_boundary_surfaces(model,1,"B4 free-strip resource P2 audit")
    p1=RectGridBase.public_volume(api,1)
    required=Set{Tuple}()
    primary_tags,primary_xyz,_=api.mesh.get_nodes(3,1,true)
    positions=Dict(Tuple(p)=>tag for (tag,p) in zip(primary_tags,eachcol(reshape(primary_xyz,3,:))))
    for block in p1.blocks,cell in eachcol(block.nodes)
        msh_dimension(block.msh)==3 || continue
        # Public-volume extraction uses local IDs, whereas carrier keys use the
        # original API tags. Match by coordinates before raising the order.
        nodes=Tuple(positions[Tuple(p1.coords[:,node])] for node in cell)
        for pattern in RectGridQuad.FACES[Int(block.msh)]
            face=Tuple(nodes[k] for k in pattern)
            for k in eachindex(face)
                push!(required,b4_free_strip_support_key((face[k],face[mod1(k+1,length(face))])))
            end
            length(face)==4 && push!(required,b4_free_strip_support_key(face))
        end
        block.msh==5 && push!(required,b4_free_strip_support_key(nodes))
    end
    api.mesh.set_order(2)
    nodes,xyz,_=api.mesh.get_nodes(3,1,true)
    @test length(nodes)==size(p1.coords,2)+length(required)
    @test length(xyz)==3length(nodes) && allunique(nodes)
    points=Dict(node=>Tuple(p) for (node,p) in zip(nodes,eachcol(reshape(xyz,3,:))))
    owners=Dict{UInt64,Tuple{Int,Int}}();counts=zeros(Int,4)
    for node in nodes
        p,parameters,dimension,owner=api.mesh.get_node(node)
        @test dimension in 0:3 && Tuple(p)==points[node]
        owners[node]=(dimension,owner);counts[dimension+1]+=1
        @test dimension in (1,2) ? length(parameters)==dimension : isempty(parameters)
        dimension in (1,2) && (@test maximum(abs.(api.model.get_value(dimension,owner,parameters).-p))<=2e-11)
    end
    types,element_tags,connectivity=api.mesh.get_elements(3,1)
    p1_families=Dict(Int(block.msh)=>size(block.nodes,2) for block in p1.blocks if msh_dimension(block.msh)==3)
    quadratic_families=Dict(4=>11,5=>12,6=>13,7=>14)
    @test Dict(Int(msh)=>length(tags) for (msh,tags) in zip(types,element_tags))==
        Dict(quadratic_families[msh]=>count for (msh,count) in p1_families)
    supports=Dict{UInt64,Tuple}()
    for (msh,flat) in zip(types,connectivity)
        for edge in eachcol(reshape(api.mesh.get_element_edge_nodes(msh,1,false),3,:))
            support=b4_free_strip_support_key(edge[1:2]);node=edge[3]
            @test !haskey(supports,node) || supports[node]==support
            supports[node]=support
        end
        if msh in (12,13,14)
            for face in eachcol(reshape(api.mesh.get_element_face_nodes(msh,4,1,false),9,:))
                support=b4_free_strip_support_key(face[1:4]);node=face[9]
                @test !haskey(supports,node) || supports[node]==support
                supports[node]=support
            end
        end
        if msh==12
            for cell in eachcol(reshape(flat,27,:))
                support=b4_free_strip_support_key(cell[1:8]);node=cell[27]
                @test !haskey(supports,node) || supports[node]==support
                supports[node]=support
            end
        end
        reference,_=api.mesh.get_integration_points(msh,"Gauss3")
        jacobians,determinants,coordinates=api.mesh.get_jacobians(msh,reference,1)
        @test all(isfinite,jacobians) && all(>(0),determinants) && all(isfinite,coordinates)
        @test length(jacobians)==9length(determinants) && length(coordinates)==3length(determinants)
    end
    @test length(supports)==length(required) && Set(values(supports))==required
    @test allunique(values(supports))
    expected_owners=copy(expected.primary_owners)
    for (node,support) in supports
        carrier=length(support)==2 ? get(expected.edges,support,(3,1)) :
            length(support)==4 ? get(expected.quads,support,(3,1)) : (3,1)
        expected_owners[node]=carrier
        @test owners[node]==carrier
        mean=ntuple(d->sum(points[primary][d] for primary in support)/length(support),3)
        @test maximum(abs.(points[node].-mean))<=2e-11
    end
    @test owners==expected_owners
    support_nodes=Dict(support=>node for (node,support) in supports)
    surface_records=Dict{String,Any}[]
    for signed_surface in surfaces
        surface=abs(signed_surface)
        owned=Set(node for (node,carrier) in expected_owners if carrier==(2,surface))
        closure=Set{UInt64}()
        for face in expected.faces[surface]
            union!(closure,face)
            for k in eachindex(face)
                edge=b4_free_strip_support_key((face[k],face[mod1(k+1,length(face))]))
                push!(closure,support_nodes[edge])
            end
            length(face)==4 && push!(closure,support_nodes[b4_free_strip_support_key(face)])
        end
        own=b4_free_strip_surface_query(api,surface,false,length(owned))
        boundary=b4_free_strip_surface_query(api,surface,true,length(closure))
        @test own==owned && boundary==closure && own⊆boundary
        push!(surface_records,Dict("surface"=>surface,"owned"=>length(own),"closure"=>length(boundary)))
    end
    # Whole-domain P2 positivity and exact integrals are stronger than the
    # sampled public quadrature protocol above. This helper differentiates
    # the represented quadratic maps independently of production Jacobians.
    volume=RectGridBase.public_volume(api,1);total=zero(B4Cert.Q);cells=0
    for block in volume.blocks,cell in eachcol(block.nodes)
        msh_dimension(block.msh)==3 || continue
        proof=B4Cert.quadratic_map_certificate(Int(block.msh),Tuple(B4Cert.point(volume,node) for node in cell))
        @test proof.volume>0
        total+=proof.volume;cells+=1
    end
    @test abs(total-1)<=B4Cert.Q(2e-11)
    @test cells==sum(values(p1_families))
    println("B4_FREE_STRIP_P2_AUDIT grid=",grid," layers=",layers,
            " nodes=",length(nodes)," supports=",length(supports)," owners=",counts)
    signatures=[(length(support),Tuple(sort!([points[primary] for primary in support])),
                 points[node],owners[node]) for (node,support) in supports]
    buffer=IOBuffer();show(buffer,sort!(signatures))
    return Dict{String,Any}("grid"=>collect(grid),"layers"=>layers,"policy"=>"free",
        "nodes"=>length(nodes),"primary_nodes"=>size(p1.coords,2),"cells"=>cells,
        "supports"=>length(supports),"owners"=>counts,"surfaces"=>surface_records,
        "digest"=>bytes2hex(sha256(take!(buffer))))
end

b4_free_strip_records=Dict{String,Any}[]
b4_free_strip_p2_records=Dict{String,Any}[]
b4_free_strip_source_records=Dict{String,Any}[]
b4_free_strip_witness_records=Dict{String,Any}[]
b4_free_strip_ast_records=Dict{String,Any}[]
b4_free_strip_growth_records=Dict{String,Any}[]
public_label_records=Dict{String,Any}[]

function b4_free_strip_fixture(grid,layers,recombined)
    a,b=grid;@assert a>=5 && b==1 && !recombined
    f=B4Cert.fixture("resource";strip_length=a,layers=:one,pins=(1,2,3,4),laterals=false)
    return replace(f.source,"Layers{1}"=>"Layers{$layers}")
end
const B4_TINY=get(ENV,"TESSELLA_B4_RESOURCE_TINY","0")=="1"
const B4_LAYERS=B4_TINY ? (1,3,4) : (500,1000,2000)
const B4_SOURCES=B4_TINY ? (5,7,9) : (1000,2000,4000)

@testset "B4 free-strip allocation and actual complete product" begin
    @test b4_free_strip_boxes(Core.Box)==1
    @test b4_free_strip_boxes(Core.Box(1))==1
    @test b4_free_strip_boxes(GlobalRef(Core,:Box))==1
    @test b4_free_strip_boxes(Base.uncompressed_ast(first(methods(b4_free_strip_box_control))))>0
    for method in methods(b4_free_strip_keyword_box_control)
        @test sum(body->b4_free_strip_boxes(Base.uncompressed_ast(body)),
            b4_free_strip_lowered_methods(method);init=0)>0
    end
    b4_free_strip_rejection_controls()
    for witness in B4_FREE_RETAINED_WITNESSES
        push!(b4_free_strip_witness_records,b4_free_strip_retained_witness(witness))
    end
    # An independent non-column prism control reaches the resource helper's
    # warped-face arm independently of the selected unit-grid factories.
    # Its reference determinant is 1-s^2/64 on triangle times s in [0,1],
    # giving volume (1/2)*(1-1/192)=191/384. A quad diagonal cannot supply it.
    warped=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),
            (0.,0.,1.),(1.,.125,1.),(.125,1.,1.))
    b4_free_strip_map(warped,6)
    bound=b4_free_strip_determinant(b4_free_strip_difference(warped[5],warped[4]),
        b4_free_strip_difference(warped[6],warped[4]),b4_free_strip_difference(warped[4],warped[1]))
    exact=big(63)//64
    @test Rational{BigInt}(bound[1])<=exact<=Rational{BigInt}(bound[2])
    control=MixedMesh(hcat((collect(p) for p in warped)...),
        [ElementBlock(6,reshape(Int32.(1:6),6,1))])
    flux=sum(RectGridQuad.face_volume(control,face,warped[1]) for face in RectGridQuad.FACES[6])
    @test isapprox(flux,191/384;rtol=16eps(),atol=0)
    @test any(length(face)==4 && RectGridBase.det(warped[face[2]].-warped[face[1]],
        warped[face[3]].-warped[face[1]],warped[face[4]].-warped[face[1]])!=0
        for face in RectGridQuad.FACES[6])
    for recombined in (false,)
        previous=Dict{Symbol,Any}()
        for layers in B4_LAYERS
            text=b4_free_strip_fixture((5,1),layers,recombined)
            mktempdir() do directory
                path=joinpath(directory,"grid.geo");write(path,text)
                geometry=execute_geo(path;mesh_dim=0)
                source=mesh_model_surface(geometry.model,1)
                base=mesh_model_volume(geometry.model,1)
                plan=RectGridModel._extrude_nonew_plan(geometry.model,1,"free resource stage")
                @test plan.catalog isa FreeB4.Catalog
                @test !plan.catalog.recombined
                @test Tuple(size(block.nodes,2) for block in base.blocks)==Tuple(
                    count for count in plan.catalog.cell_counts if count>0)
                caller="free resource stages"
                catalog=plan.catalog;cols=plan.sweep.cols
                spec=geometry.model.meshing.extrude_specs[(3,1)]
                projected=model_to_mixed(geometry.model,base,3,1)
                generated=execute_geo(path;mesh_dim=3)
                merge_parts=generated.mesh_parts
                merge_owners=Tessella.GeoExec._geo_mesh_part_node_entities(geometry.model,merge_parts,caller)
                paths=((:plan,()->FreeB4.plan(catalog.source,catalog.levels,catalog.layer_refs,caller),
                        completed->FreeB4.emit(completed,cols,spec,caller).volume),
                    (:emitter,()->FreeB4.emit(catalog,cols,spec,caller),result->result.volume),
                    (:standalone,()->mesh_model_volume(geometry.model,1),identity),
                    (:geo,()->execute_geo(path;mesh_dim=3),result->geo_entity_mesh(result,3,1)),
                    (:projection,()->model_to_mixed(geometry.model,base,3,1),identity),
                    (:classification,()->Tessella.API._classify_cached_mesh(geometry.model,base,3,1,base),
                        classification->begin
                            @test classification.node_entities==projected.entity_data.node_entities
                            @test sum(length,values(classification.cell_entities))==sum(size(block.nodes,2) for block in base.blocks)
                            classification.mesh
                        end),
                    (:merge,()->Tessella.GeoExec._geo_merge_entity_meshes_mixed(merge_parts,caller;node_entities=merge_owners),
                        result->first(result)))
                expected=nothing
                for (name,make,extract) in paths
                    row=b4_free_strip_profile(make,extract,source,layers,recombined)
                    expected===nothing ? (expected=row.digest) : (@test row.digest==expected)
                    if haskey(previous,name) && !B4_TINY
                        bound=2.15previous[name].allocated+65536
                        @test row.allocated<=bound
                        push!(b4_free_strip_growth_records,Dict("axis"=>"layers","path"=>String(name),
                            "from"=>previous[name].layers,"to"=>layers,"from_bytes"=>previous[name].allocated,
                            "to_bytes"=>row.allocated,"bound"=>bound,"passed"=>row.allocated<=bound))
                    end
                    previous[name]=row
                    push!(b4_free_strip_records,Dict{String,Any}("policy"=>recombined ? "recombined" : "free",
                        "path"=>String(name),(String(key)=>value for (key,value) in pairs(row))...))
                    println("B4_FREE_STRIP_RESOURCE policy=",recombined," path=",name," ",row)
                    flush(stdout)
                end
                api=Tessella.API;api.initialize()
                try
                    api.open_geo!(path;mesh_dim=0)
                    row=b4_free_strip_profile(()->api.mesh.generate(3),_->RectGridBase.public_volume(api,1),source,layers,recombined)
                    @test row.digest==expected
                    if haskey(previous,:api) && !B4_TINY
                        bound=2.15previous[:api].allocated+65536
                        @test row.allocated<=bound
                        push!(b4_free_strip_growth_records,Dict("axis"=>"layers","path"=>"api",
                            "from"=>previous[:api].layers,"to"=>layers,"from_bytes"=>previous[:api].allocated,
                            "to_bytes"=>row.allocated,"bound"=>bound,"passed"=>row.allocated<=bound))
                    end
                    previous[:api]=row
                    push!(b4_free_strip_records,Dict{String,Any}("policy"=>recombined ? "recombined" : "free",
                        "path"=>"api",(String(key)=>value for (key,value) in pairs(row))...))
                    println("B4_FREE_STRIP_RESOURCE policy=",recombined," path=api ",row)
                    flush(stdout)
                finally
                    api.finalize()
                end
            end
        end
    end
    # Vary actual source size separately from layer count. This includes source
    # sampling, B4 topology/Jordan admission, physical propagation and emission.
    # No pure planner measurement is substituted for this actual entry path.
    for recombined in (false,)
        previous=nothing
        for a in B4_SOURCES
            grid=(a,1);text=b4_free_strip_fixture(grid,3,recombined)
            mktempdir() do directory
                path=joinpath(directory,"source-growth.geo");write(path,text)
                geometry=execute_geo(path;mesh_dim=0)
                source=mesh_model_surface(geometry.model,1)
                row=b4_free_strip_profile(()->mesh_model_volume(geometry.model,1),identity,source,3,recombined;grid)
                if previous!==nothing && !B4_TINY
                    bound=2.15previous.allocated+65536
                    @test row.allocated<=bound
                    push!(b4_free_strip_growth_records,Dict("axis"=>"source","path"=>"standalone",
                        "from"=>a÷2,"to"=>a,"from_bytes"=>previous.allocated,"to_bytes"=>row.allocated,
                        "bound"=>bound,"passed"=>row.allocated<=bound))
                end
                previous=row
                push!(b4_free_strip_source_records,Dict{String,Any}("policy"=>recombined ? "recombined" : "free",
                    "source_cells"=>a,"source_nodes"=>2(a+1),
                    (String(key)=>value for (key,value) in pairs(row))...))
                println("B4_FREE_STRIP_SOURCE_RESOURCE policy=",recombined," grid=",grid," ",row)
                flush(stdout)
            end
        end
    end
    # Three small P2 audits exercise actual support identity/owners/closure and
    # computed native UV, without charging query/proof allocations to creation.
    for recombined in (false,),(grid,layers) in (((5,1),1),((5,1),3),((7,1),3))
        mktempdir() do directory
            path=joinpath(directory,"quadratic.geo");write(path,b4_free_strip_fixture(grid,layers,recombined))
            api=Tessella.API;api.initialize()
            try
                api.open_geo!(path;mesh_dim=0);api.mesh.generate(3)
                projected=model_to_mixed(api.CURRENT[],RectGridBase.public_volume(api,1),3,1)
                push!(b4_free_strip_p2_records,b4_free_strip_quadratic_audit(api,projected,layers,recombined;grid))
            finally
                api.finalize()
            end
        end
    end
    free_methods=Set{Method}()
    for mod in (FreeB4,FreeB4.SourcePhase,FreeB4.SourcePhase.LocalPhase,
                FreeB4.FinalFactories,FreeB4.CenterFan)
        union!(free_methods,b4_free_strip_module_methods(mod))
    end
    union!(free_methods,methods(FreeB4.Catalog),methods(FreeB4.SourcePhase.Phase))
    private_methods=copy(free_methods)
    # Source admission is shared with the rectangular disk implementation.
    # Inventory those actual bodies alongside the private production module;
    # prototypes' method counts are not a substitute for this final inventory.
    for name in (:_extrude_nonew_rect_grid_sizes,:_extrude_nonew_rect_grid_shape,
                 :_extrude_nonew_rect_grid_axis,:_extrude_nonew_rect_grid_assign!,
                 :_extrude_nonew_rect_grid_topology,:_extrude_nonew_rect_grid_boundary_tree!,
                 :_extrude_nonew_rect_grid_boundary_pair,:_extrude_nonew_rect_grid_boundary_certify,
                 :_extrude_nonew_rect_grid_source)
        for method in methods(getfield(RectGridModel,name))
            union!(free_methods,b4_free_strip_lowered_methods(method))
        end
    end
    @test length(free_methods)>=53
    checked_methods=copy(free_methods)
    for name in names(RectGridModel;all=true)
        ((startswith(String(name),"_extrude_nonew_rect_grid") || startswith(String(name),"_extrude_nonew_b4_strip")) ||
         name in (:_extrude_nonew_plan,:_extrude_nonew_levels,
                  :_extrude_nonew_complete_plan,:_extrude_nonew_dynamic_grid_face_capacity,:_extrude_nonew_face_cycle,:_extrude_nonew_face_orientation,
                  :_extrude_nonew_count_face!,:_extrude_nonew_count_faces!,
                  :_extrude_nonew_projection_context,
                  :_extrude_nonew_cell_counts,:_extrude_nonew_projection_lookup,
                  :_extrude_quadtri_centroid,:_extrude_quadtri_centroid_exact,
                  :_model_curve_transfinite_native_params,:_model_curve_transfinite_density_params,
                  :_model_native_line_sampling_speed,:_model_native_line_sampling_length,
                  :_transfinite_native_line_speed,:_model_curve_transfinite_mass,:_transfinite_val,
                  :_model_surface_curve_writeback!,
                  :_model_surface_curve_parameter_match)) || continue
        value=getfield(RectGridModel,name);value isa Function || continue
        for method in methods(value)
            union!(checked_methods,b4_free_strip_lowered_methods(method))
        end
    end
    for method in methods(Tessella.GeoExec._geo_merge_entity_meshes_mixed)
        union!(checked_methods,b4_free_strip_lowered_methods(method))
    end
    for function_value in (Tessella.API._classify_cached_mesh,Tessella.API._mixed_classification,
                           Tessella.API._merge_classified_parts,
                           Tessella.API._reconcile_generated_records!,
                           Tessella.API._validate_record_cache_labels,
                           Tessella.API._apply_mesh_order,
                           Tessella.API._dim01_add_node!,
                           Tessella.API._dim01_quadratic!,
                           Tessella.API._dim01_reconcile_records!,
                           Tessella.API._append_native_boundary_nodes!,
                           Tessella.API._generate,Tessella.API._get_nodes)
        for method in methods(function_value)
            union!(checked_methods,b4_free_strip_lowered_methods(method))
        end
    end
    for name in PUBLIC_LABEL_RESOURCE_API_TARGETS
        @test isdefined(Tessella.API,name)
        for method in methods(getfield(Tessella.API,name))
            union!(checked_methods,b4_free_strip_lowered_methods(method))
        end
    end
    for function_value in (public_label_resource_fixture,public_label_resource_rows,
                           Tessella.Model.add_discrete_elements!,
                           Tessella.Recombine._gmsh_recombine_pair_measure,
                           Tessella.Recombine._recombine_angle,Tessella.Recombine.recombine_triangles)
        for method in methods(function_value)
            union!(checked_methods,b4_free_strip_lowered_methods(method))
        end
    end
    for method in sort!(collect(checked_methods);by=string)
        boxes=b4_free_strip_boxes(Base.uncompressed_ast(method))
        boxes==0 || println("B4_FREE_STRIP_RESOURCE_BOX_METHOD=",method," boxes=",boxes)
        @test boxes==0
        push!(b4_free_strip_ast_records,Dict("method"=>string(method),"boxes"=>boxes,
            "module"=>string(method.module),"native_free_body"=>method in free_methods,
            "private_free_body"=>method in private_methods))
    end
    println("B4_FREE_STRIP_RESOURCE_BOX_METHODS=",length(checked_methods),
        " native_free_and_source_bodies=",length(free_methods)," private_free_bodies=",length(private_methods))
    append!(public_label_records,public_label_resource_rows())
    final=b4_free_strip_snapshot()
    changed=sort!(collect(path for path in union(keys(B4_FREE_STRIP_START),keys(final))
        if get(B4_FREE_STRIP_START,path,nothing)!=get(final,path,nothing)))
    isempty(changed) || println("B4_FREE_STRIP_RESOURCE_CHANGED_INPUTS=",changed)
    @test isempty(changed)
end
if !isempty(ARGS)
    ispath(ARGS[1]) && error("refusing to overwrite resource evidence")
    open(ARGS[1],"w") do io
        TOML.print(io,Dict("julia"=>string(VERSION),"records"=>b4_free_strip_records,
            "p2_records"=>b4_free_strip_p2_records,"source_records"=>b4_free_strip_source_records,
            "witness_records"=>b4_free_strip_witness_records,"ast_records"=>b4_free_strip_ast_records,
            "public_label_records"=>public_label_records,
            "growth_records"=>b4_free_strip_growth_records,"growth_factor"=>2.15,"growth_slack_bytes"=>65536,
            "tiny"=>B4_TINY,"inputs"=>B4_FREE_STRIP_START,"inputs_after"=>b4_free_strip_snapshot()))
    end
end
println("B4_FREE_STRIP_RESOURCE_OK records=",length(b4_free_strip_records)," source_records=",length(b4_free_strip_source_records),
        " p2_records=",length(b4_free_strip_p2_records)," inputs_stable=",b4_free_strip_snapshot()==B4_FREE_STRIP_START)
