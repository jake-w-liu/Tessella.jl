# Normal compilation and --check-bounds=yes; no Gmsh dependency.
# 24 measured rows: 2x3 grid, N500/1000/2000, both policies, four entry paths.
# Six actual F-growth rows are measured separately; six small P2 audits are
# unmeasured query/proof checks.
# Allocations/times below charge construction only; geometry/query proofs and
# public-volume extraction are outside @allocated/@elapsed.
using Tessella, Test, SHA, TOML
using Tessella.Model: mesh_model_volume, mesh_model_surface, model_to_mixed
using Tessella.Elements: ElementBlock, MixedMesh, msh_dimension
include(joinpath(@__DIR__, "../../test/geometry/quadtri_nonew_two_tri_certificates.jl"))
const RectGridBase = QuadTriNoNewTwoTriCertificates
const RectGridQuad = QuadTriNoNewCertificates
const RectGridModel = Tessella.Model
const RECT_GRID_ROOT = normpath(joinpath(@__DIR__, "../.."))
rect_grid_inputs() = sort!(vcat([joinpath(RECT_GRID_ROOT,"Project.toml")],
    [joinpath(directory,file) for (directory,_,files) in walkdir(joinpath(RECT_GRID_ROOT,"src"))
     for file in files if endswith(file,".jl")],
    [@__FILE__,joinpath(RECT_GRID_ROOT,"test/geometry/quadtri_nonew_two_tri_certificates.jl"),
     joinpath(RECT_GRID_ROOT,"test/geometry/quadtri_nonew_certificates.jl")]))
rect_grid_snapshot() = Dict(relpath(path,RECT_GRID_ROOT)=>bytes2hex(sha256(read(path))) for path in rect_grid_inputs())
const RECT_GRID_START = rect_grid_snapshot()
println("START RECT_GRID_RESOURCE Julia=",VERSION," inputs=",length(RECT_GRID_START))

function rect_grid_boxes(value)
    value isa Core.Box && return 1
    value isa GlobalRef && return value==GlobalRef(Core,:Box) ? 1 : 0
    value isa Core.CodeInfo && return sum(rect_grid_boxes,value.code;init=0)
    value isa Expr && return sum(rect_grid_boxes,value.args;init=0)
    value isa Core.NewvarNode && return rect_grid_boxes(value.slot)
    value isa QuoteNode && return rect_grid_boxes(value.value)
    return value===Core.Box ? 1 : 0
end

function rect_grid_box_control()
    captured=0
    read=()->captured
    captured=1
    return read
end

function rect_grid_keyword_box_control(value=1;initial=0)
    captured=initial
    read=()->captured
    captured=value
    return read
end

function rect_grid_lowered_methods(method::Method)
    result=Method[method]
    # Optional positional wrappers can lack Base.bodyfunction metadata.
    # Generated keyword bindings retain the complete implementation bodies.
    prefix="#"*String(method.name)*"#"
    for name in names(method.module;all=true)
        startswith(String(name),prefix) || continue
        value=getfield(method.module,name)
        value isa Function || continue
        append!(result,methods(value))
    end
    return result
end

# Outward interval arithmetic encloses the exact arithmetic of the represented
# Float64 corners. The resource fixtures stay at ordinary finite unit scales.
# A positive lower endpoint is a whole-map proof, not a point-sample criterion.
rect_grid_interval(x::Float64)=(x,x)
rect_grid_add(a,b)=(prevfloat(a[1]+b[1]),nextfloat(a[2]+b[2]))
rect_grid_sub(a,b)=(prevfloat(a[1]-b[2]),nextfloat(a[2]-b[1]))
function rect_grid_mul(a,b)
    products=(a[1]*b[1],a[1]*b[2],a[2]*b[1],a[2]*b[2])
    return (prevfloat(minimum(products)),nextfloat(maximum(products)))
end
rect_grid_scale(a,x::Float64)=rect_grid_mul(a,rect_grid_interval(x))
rect_grid_difference(p,q)=ntuple(k->rect_grid_sub(rect_grid_interval(p[k]),rect_grid_interval(q[k])),3)
rect_grid_average(a,b)=ntuple(k->rect_grid_scale(rect_grid_add(a[k],b[k]),.5),3)
function rect_grid_determinant(a,b,c)
    first=rect_grid_mul(a[1],rect_grid_sub(rect_grid_mul(b[2],c[3]),rect_grid_mul(b[3],c[2])))
    second=rect_grid_mul(a[2],rect_grid_sub(rect_grid_mul(b[1],c[3]),rect_grid_mul(b[3],c[1])))
    third=rect_grid_mul(a[3],rect_grid_sub(rect_grid_mul(b[1],c[2]),rect_grid_mul(b[2],c[1])))
    return rect_grid_add(rect_grid_sub(first,second),third)
end

# Differentiate the reference interpolation directly. A translated Hex8 has
# a bilinear base determinant; a collapsed Pyramid5 has four base coefficients.
# A Prism6 is affine on its triangle and quadratic along its interval: three
# positive Bernstein coefficients at each triangle vertex certify that domain.
function rect_grid_map(points,msh)
    if msh==4
        bound=rect_grid_determinant(rect_grid_difference(points[2],points[1]),
            rect_grid_difference(points[3],points[1]),rect_grid_difference(points[4],points[1]))
        @test bound[1]>0
    elseif msh==5
        delta=points[5].-points[1]
        @test all(index->points[index+4].-points[index]==delta,1:4)
        for index in 1:4
            bound=rect_grid_determinant(rect_grid_difference(points[mod1(index+1,4)],points[index]),
                rect_grid_difference(points[mod1(index-1,4)],points[index]),rect_grid_difference(points[5],points[1]))
            @test bound[1]>0
        end
    elseif msh==6
        low_u=rect_grid_difference(points[2],points[1]);high_u=rect_grid_difference(points[5],points[4])
        low_v=rect_grid_difference(points[3],points[1]);high_v=rect_grid_difference(points[6],points[4])
        for vertex in 1:3
            direction=rect_grid_difference(points[vertex+3],points[vertex])
            first=rect_grid_determinant(low_u,low_v,direction)
            last=rect_grid_determinant(high_u,high_v,direction)
            middle=rect_grid_determinant(rect_grid_average(low_u,high_u),rect_grid_average(low_v,high_v),direction)
            coefficient=rect_grid_sub(rect_grid_scale(middle,2.),rect_grid_scale(rect_grid_add(first,last),.5))
            @test first[1]>0 && last[1]>0 && coefficient[1]>0
        end
    elseif msh==7
        for index in 1:4
            bound=rect_grid_determinant(rect_grid_difference(points[mod1(index+1,4)],points[index]),
                rect_grid_difference(points[mod1(index-1,4)],points[index]),rect_grid_difference(points[5],points[index]))
            @test bound[1]>0
        end
    else
        @test false # A foreign family violates this declared product.
    end
end

function rect_grid_audit(volume,source,layers,recombined;grid=(2,3))
    a,b=grid;vertices=(a+1)*(b+1);macros=a*b;boundary=2(a+b)
    primary=vertices*(layers+1)
    @test size(volume.coords)==(3,primary)
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
    logical=Vector{Tuple{Int,Int}}(undef,size(volume.coords,2))
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
    @test length(grid)==primary && isempty(extra)
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
    domain_counts=zeros(Int,macros,layers);domain_volumes=zeros(Float64,macros,layers)
    families=Dict{Int,Int}();faces=Dict{Tuple,Tuple{Int,Tuple}}()
    used=falses(size(volume.coords,2));signatures=Tuple[]
    for block in volume.blocks
        block isa ElementBlock && msh_dimension(block.msh)==3 || continue
        @test block.msh in (4,5,6,7)
        families[Int(block.msh)]=get(families,Int(block.msh),0)+size(block.nodes,2)
        for cell in eachcol(block.nodes)
            @test length(unique(cell))==length(cell)
            used[cell].=true
            ids=unique(logical[node][1] for node in cell)
            @test length(ids)>=3
            candidates=copy(parents[first(ids)])
            for node in ids;intersect!(candidates,parents[node]);end
            @test length(candidates)==1
            parent=only(candidates)
            levels=[logical[node][2] for node in cell]
            level=minimum(levels)
            @test maximum(levels)==level+1 && all(v->v in (level,level+1),levels)
            domain_counts[parent,level]+=1
            p=Tuple(ntuple(d->volume.coords[d,node],3) for node in cell)
            rect_grid_map(p,Int(block.msh))
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
            push!(signatures,(Int(block.msh),Tuple(sort!(collect(p)))))
        end
    end
    @test all(used)
    if recombined
        expected=Dict(4=>4macros-2boundary,7=>macros+boundary)
        layers>1 && (expected[5]=macros*(layers-1))
        @test families==expected
        @test length(faces)==(3macros+boundary÷2)*layers+9macros-3boundary÷2
    else
        @test sum(values(families))<=6macros*layers
    end
    exterior=(recombined ? boundary : 2boundary)*layers+3macros
    @test 2length(faces)==sum(length(RectGridQuad.FACES[msh])*count for (msh,count) in families)+exterior
    groups=Dict{Tuple,Vector{Tuple}}()
    for (uses,cycle) in values(faces)
        uses==1 || continue
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
        @test 1<=domain_counts[parent,level]<=6
        recombined && level<layers && (@test domain_counts[parent,level]==1)
        @test isapprox(domain_volumes[parent,level],areas[parent]*(planes[level+1]-planes[level]);rtol=64eps(),atol=0)
    end
    @test count(value->value[1]==1,values(faces))==exterior
    @test isapprox(sum(domain_volumes),1.;rtol=128eps())
    io=IOBuffer();show(io,sort!(signatures))
    return bytes2hex(sha256(take!(io)))
end

function rect_grid_profile(make,extract,source,layers,recombined;grid=(2,3))
    result=make()
    allocated=minimum(@allocated(make()) for _ in 1:3)
    seconds=minimum(@elapsed(make()) for _ in 1:3)
    volume=extract(result)
    digest=rect_grid_audit(volume,source,layers,recombined;grid)
    return (;layers,allocated,seconds,digest,nodes=size(volume.coords,2))
end

# The measured entry products are P1. Actual P2 support identity, carrier
# ownership and public map/parameter queries are audited separately below.
# Lower-cell carriers come from the actual classified P1 projection, rather
# than an endpoint-owner heuristic or a guessed subdivision of CAD faces.
rect_grid_support_key(nodes)=Tuple(sort!(Int.(collect(nodes))))

function rect_grid_actual_carriers(api,projected)
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
                edges[rect_grid_support_key(remap[node] for node in cell)]=carrier
            else
                for i in eachindex(cell)
                    support=rect_grid_support_key((remap[cell[i]],remap[cell[mod1(i+1,length(cell))]]))
                    haskey(edges,support) && edges[support][1]<dimension && continue
                    @test !haskey(edges,support) || edges[support]==carrier
                    edges[support]=carrier
                end
                block.msh==3 && (quads[rect_grid_support_key(remap[node] for node in cell)]=carrier)
            end
        end
    end
    return (;edges,quads)
end

function rect_grid_surface_query(api,surface,boundary,expected)
    # Native cache queries compute UV for all selected nodes. This deliberately
    # does not assert parity with Gmsh's partially empty stored P2 parameters.
    nodes,coordinates,parameters=api.mesh.get_nodes(2,surface,boundary,true)
    @test length(nodes)==expected && allunique(nodes)
    @test length(coordinates)==3expected && length(parameters)==2expected
    @test maximum(abs.(api.model.get_value(2,surface,parameters).-coordinates))<=2e-11
    return Set(nodes)
end

function rect_grid_quadratic_audit(api,projected,layers,recombined;grid=(2,3))
    expected=rect_grid_actual_carriers(api,projected)
    model=api.CURRENT[]
    surfaces=RectGridModel._model_volume_boundary_surfaces(model,1,"rectangular-grid resource P2 audit")
    api.mesh.set_order(2)
    nodes,xyz,_=api.mesh.get_nodes(3,1,true)
    a,b=grid;macros=a*b;boundary=2(a+b)
    @test length(nodes)==(2a+1)*(2b+1)*(2layers+1)
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
    @test counts==[8,8layers+4boundary-12,
        4(boundary-2)*layers+8macros-4boundary+6,
        (4macros-boundary+1)*(2layers-1)]
    types,element_tags,connectivity=api.mesh.get_elements(3,1)
    if recombined
        families=Dict(Int(msh)=>length(tags) for (msh,tags) in zip(types,element_tags))
        exact=Dict(11=>4macros-2boundary,14=>macros+boundary)
        layers>1 && (exact[12]=macros*(layers-1))
        @test families==exact
    end
    supports=Dict{UInt64,Tuple}()
    for (msh,flat) in zip(types,connectivity)
        for edge in eachcol(reshape(api.mesh.get_element_edge_nodes(msh,1,false),3,:))
            support=rect_grid_support_key(edge[1:2]);node=edge[3]
            @test !haskey(supports,node) || supports[node]==support
            supports[node]=support
        end
        if msh in (12,13,14)
            for face in eachcol(reshape(api.mesh.get_element_face_nodes(msh,4,1,false),9,:))
                support=rect_grid_support_key(face[1:4]);node=face[9]
                @test !haskey(supports,node) || supports[node]==support
                supports[node]=support
            end
        end
        if msh==12
            for cell in eachcol(reshape(flat,27,:))
                supports[cell[27]]=rect_grid_support_key(cell[1:8])
            end
        end
        # This is a public quadrature protocol check. It is distinct from
        # the P1 whole-reference-map proof in rect_grid_map above.
        reference,_=api.mesh.get_integration_points(msh,"Gauss3")
        jacobians,determinants,coordinates=api.mesh.get_jacobians(msh,reference,1)
        @test all(isfinite,jacobians) && all(>(0),determinants) && all(isfinite,coordinates)
        @test length(jacobians)==9length(determinants) && length(coordinates)==3length(determinants)
    end
    @test length(supports)==(7macros+3boundary÷2+1)*layers+3macros+boundary÷2
    @test allunique(values(supports))
    for (node,support) in supports
        carrier=length(support)==2 ? get(expected.edges,support,(3,1)) :
            length(support)==4 ? get(expected.quads,support,(3,1)) : (3,1)
        @test owners[node]==carrier
        mean=ntuple(d->sum(points[primary][d] for primary in support)/length(support),3)
        @test maximum(abs.(points[node].-mean))<=2e-11
    end
    for signed_surface in surfaces
        surface=abs(signed_surface)
        link=get(model.meshing.extrude_sources,(2,surface),nothing)
        lateral=link!==nothing && link[1]==1
        width=lateral ? (abs(link[2]) in (1,3) ? a : b) : 0
        owned=lateral ? (2width-1)*(2layers-1) : (2a-1)*(2b-1)
        closure=lateral ? (2width+1)*(2layers+1) : (2a+1)*(2b+1)
        own=rect_grid_surface_query(api,surface,false,owned)
        boundary=rect_grid_surface_query(api,surface,true,closure)
        @test own⊆boundary
    end
    println("RECT_GRID_P2_AUDIT grid=",grid," layers=",layers," policy=",recombined,
            " nodes=",length(nodes)," supports=",length(supports)," owners=",counts)
    signatures=[(length(support),Tuple(sort!([points[primary] for primary in support])),
                 points[node],owners[node]) for (node,support) in supports]
    buffer=IOBuffer();show(buffer,sort!(signatures))
    return Dict{String,Any}("grid"=>collect(grid),"layers"=>layers,"policy"=>recombined ? "recombined" : "free",
        "nodes"=>length(nodes),"supports"=>length(supports),"owners"=>counts,
        "digest"=>bytes2hex(sha256(take!(buffer))))
end

rect_grid_records=Dict{String,Any}[]
rect_grid_p2_records=Dict{String,Any}[]
rect_grid_source_records=Dict{String,Any}[]

function rect_grid_fixture(grid,layers,recombined)
    a,b=grid
    fixture=RectGridBase.fixture("resource";laterals=recombined,layers=:one,pins=(1,2,3,4))
    return replace(fixture.source,"Transfinite Curve{1,2,3,4}=2;"=>
        "Transfinite Curve{1,3}=$(a+1);Transfinite Curve{2,4}=$(b+1);Recombine Surface{1};",
        "Layers{1}"=>"Layers{$layers}")
end

@testset "Rectangular-grid allocation and actual complete product" begin
    @test rect_grid_boxes(Core.Box)==1
    @test rect_grid_boxes(Core.Box(1))==1
    @test rect_grid_boxes(GlobalRef(Core,:Box))==1
    @test rect_grid_boxes(Base.uncompressed_ast(first(methods(rect_grid_box_control))))>0
    for method in methods(rect_grid_keyword_box_control)
        @test sum(body->rect_grid_boxes(Base.uncompressed_ast(body)),
            rect_grid_lowered_methods(method);init=0)>0
    end
    # An independent non-column prism control reaches the resource helper's
    # warped-face arm independently of the selected unit-grid factories.
    # Its reference determinant is 1-s^2/64 on triangle times s in [0,1],
    # giving volume (1/2)*(1-1/192)=191/384. A quad diagonal cannot supply it.
    warped=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),
            (0.,0.,1.),(1.,.125,1.),(.125,1.,1.))
    rect_grid_map(warped,6)
    bound=rect_grid_determinant(rect_grid_difference(warped[5],warped[4]),
        rect_grid_difference(warped[6],warped[4]),rect_grid_difference(warped[4],warped[1]))
    exact=big(63)//64
    @test Rational{BigInt}(bound[1])<=exact<=Rational{BigInt}(bound[2])
    control=MixedMesh(hcat((collect(p) for p in warped)...),
        [ElementBlock(6,reshape(Int32.(1:6),6,1))])
    flux=sum(RectGridQuad.face_volume(control,face,warped[1]) for face in RectGridQuad.FACES[6])
    @test isapprox(flux,191/384;rtol=16eps(),atol=0)
    @test any(length(face)==4 && RectGridBase.det(warped[face[2]].-warped[face[1]],
        warped[face[3]].-warped[face[1]],warped[face[4]].-warped[face[1]])!=0
        for face in RectGridQuad.FACES[6])
    for recombined in (false,true)
        previous=Dict{Symbol,Any}()
        for layers in (500,1000,2000)
            text=rect_grid_fixture((2,3),layers,recombined)
            mktempdir() do directory
                path=joinpath(directory,"grid.geo");write(path,text)
                geometry=execute_geo(path;mesh_dim=0)
                source=mesh_model_surface(geometry.model,1)
                base=mesh_model_volume(geometry.model,1)
                paths=((:standalone,()->mesh_model_volume(geometry.model,1),identity),
                    (:geo,()->execute_geo(path;mesh_dim=3),result->geo_entity_mesh(result,3,1)),
                    (:projection,()->model_to_mixed(geometry.model,base,3,1),identity))
                expected=nothing
                for (name,make,extract) in paths
                    row=rect_grid_profile(make,extract,source,layers,recombined)
                    expected===nothing ? (expected=row.digest) : (@test row.digest==expected)
                    haskey(previous,name) && (@test row.allocated<=2.15previous[name].allocated+65536)
                    previous[name]=row
                    push!(rect_grid_records,Dict{String,Any}("policy"=>recombined ? "recombined" : "free",
                        "path"=>String(name),(String(key)=>value for (key,value) in pairs(row))...))
                    println("RECT_GRID_RESOURCE policy=",recombined," path=",name," ",row)
                end
                api=Tessella.API;api.initialize()
                try
                    api.open_geo!(path;mesh_dim=0)
                    row=rect_grid_profile(()->api.mesh.generate(3),_->RectGridBase.public_volume(api),source,layers,recombined)
                    @test row.digest==expected
                    haskey(previous,:api) && (@test row.allocated<=2.15previous[:api].allocated+65536)
                    previous[:api]=row
                    push!(rect_grid_records,Dict{String,Any}("policy"=>recombined ? "recombined" : "free",
                        "path"=>"api",(String(key)=>value for (key,value) in pairs(row))...))
                    println("RECT_GRID_RESOURCE policy=",recombined," path=api ",row)
                finally
                    api.finalize()
                end
            end
        end
    end
    # Vary actual source size separately from layer count. This includes source
    # sampling, rectangular topology/Jordan admission, and full emission.
    # The previously measured pure planner is not substituted for this path.
    for recombined in (false,true)
        previous=nothing
        for b in (500,1000,2000)
            grid=(2,b);text=rect_grid_fixture(grid,1,recombined)
            mktempdir() do directory
                path=joinpath(directory,"source-growth.geo");write(path,text)
                geometry=execute_geo(path;mesh_dim=0)
                source=mesh_model_surface(geometry.model,1)
                row=rect_grid_profile(()->mesh_model_volume(geometry.model,1),identity,source,1,recombined;grid)
                previous===nothing || (@test row.allocated<=2.15previous.allocated+65536)
                previous=row
                push!(rect_grid_source_records,Dict{String,Any}("policy"=>recombined ? "recombined" : "free",
                    "source_cells"=>2b,"source_nodes"=>3(b+1),
                    (String(key)=>value for (key,value) in pairs(row))...))
                println("RECT_GRID_SOURCE_RESOURCE policy=",recombined," grid=",grid," ",row)
            end
        end
    end
    # Six small P2 audits exercise actual support identity/owners/closure and
    # computed native UV, without charging query/proof allocations to creation.
    for recombined in (false,true),(grid,layers) in (((2,3),1),((2,3),3),((3,3),3))
        mktempdir() do directory
            path=joinpath(directory,"quadratic.geo");write(path,rect_grid_fixture(grid,layers,recombined))
            api=Tessella.API;api.initialize()
            try
                api.open_geo!(path;mesh_dim=0);api.mesh.generate(3)
                projected=model_to_mixed(api.CURRENT[],RectGridBase.public_volume(api),3,1)
                push!(rect_grid_p2_records,rect_grid_quadratic_audit(api,projected,layers,recombined;grid))
            finally
                api.finalize()
            end
        end
    end
    checked_methods=Set{Method}()
    for name in names(RectGridModel;all=true)
        (startswith(String(name),"_extrude_nonew_rect_grid") ||
         name in (:_extrude_nonew_complete_plan,:_extrude_nonew_face_cycle,:_extrude_nonew_face_orientation,
                  :_extrude_nonew_count_face!,:_extrude_nonew_count_faces!,
                  :_extrude_nonew_projection_context,
                  :_extrude_nonew_cell_counts,:_extrude_nonew_projection_lookup,
                  :_extrude_quadtri_centroid,:_extrude_quadtri_centroid_exact,
                  :_model_curve_transfinite_native_params,:_transfinite_val,
                  :_model_surface_curve_writeback!,
                  :_model_surface_curve_parameter_match)) || continue
        value=getfield(RectGridModel,name);value isa Function || continue
        for method in methods(value)
            union!(checked_methods,rect_grid_lowered_methods(method))
        end
    end
    for method in methods(Tessella.GeoExec._geo_merge_entity_meshes_mixed)
        union!(checked_methods,rect_grid_lowered_methods(method))
    end
    for method in checked_methods
        boxes=rect_grid_boxes(Base.uncompressed_ast(method))
        boxes==0 || println("RECT_GRID_RESOURCE_BOX_METHOD=",method," boxes=",boxes)
        @test boxes==0
    end
    println("RECT_GRID_RESOURCE_BOX_METHODS=",length(checked_methods))
    final=rect_grid_snapshot()
    changed=sort!(collect(path for path in union(keys(RECT_GRID_START),keys(final))
        if get(RECT_GRID_START,path,nothing)!=get(final,path,nothing)))
    isempty(changed) || println("RECT_GRID_RESOURCE_CHANGED_INPUTS=",changed)
    @test isempty(changed)
end
if !isempty(ARGS)
    open(ARGS[1],"w") do io
        TOML.print(io,Dict("julia"=>string(VERSION),"records"=>rect_grid_records,
            "p2_records"=>rect_grid_p2_records,"source_records"=>rect_grid_source_records,"inputs"=>RECT_GRID_START))
    end
end
println("RECT_GRID_RESOURCE_OK records=",length(rect_grid_records)," source_records=",length(rect_grid_source_records),
        " p2_records=",length(rect_grid_p2_records)," inputs_stable=",rect_grid_snapshot()==RECT_GRID_START)
