# Normal compilation and --check-bounds=yes; no Gmsh dependency.
using Tessella, Test, SHA, TOML
using Tessella.Model: mesh_model_volume, mesh_model_surface, model_to_mixed
using Tessella.Elements: ElementBlock, MixedMesh, msh_dimension
include(joinpath(@__DIR__, "../../test/geometry/quadtri_nonew_two_tri_certificates.jl"))
const ThreeStripBase = QuadTriNoNewTwoTriCertificates
const ThreeStripQuad = QuadTriNoNewCertificates
const ThreeStripModel = Tessella.Model
const THREE_STRIP_ROOT = normpath(joinpath(@__DIR__, "../.."))
three_strip_inputs() = sort!(vcat([joinpath(THREE_STRIP_ROOT,"Project.toml")],
    [joinpath(directory,file) for (directory,_,files) in walkdir(joinpath(THREE_STRIP_ROOT,"src"))
     for file in files if endswith(file,".jl")],
    [@__FILE__,joinpath(THREE_STRIP_ROOT,"test/geometry/quadtri_nonew_two_tri_certificates.jl"),
     joinpath(THREE_STRIP_ROOT,"test/geometry/quadtri_nonew_certificates.jl")]))
three_strip_snapshot() = Dict(relpath(path,THREE_STRIP_ROOT)=>bytes2hex(sha256(read(path))) for path in three_strip_inputs())
const THREE_STRIP_START = three_strip_snapshot()
println("START THREE_QUAD_STRIP_RESOURCE Julia=",VERSION," inputs=",length(THREE_STRIP_START))

function three_strip_boxes(value)
    value isa Core.Box && return 1
    value isa GlobalRef && return value==GlobalRef(Core,:Box) ? 1 : 0
    value isa Core.CodeInfo && return sum(three_strip_boxes,value.code;init=0)
    value isa Expr && return sum(three_strip_boxes,value.args;init=0)
    value isa Core.NewvarNode && return three_strip_boxes(value.slot)
    value isa QuoteNode && return three_strip_boxes(value.value)
    return value===Core.Box ? 1 : 0
end

function three_strip_box_control()
    captured=0
    read=()->captured
    captured=1
    return read
end

function three_strip_keyword_box_control(value=1;initial=0)
    captured=initial
    read=()->captured
    captured=value
    return read
end

function three_strip_lowered_methods(method::Method)
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

# These checks differentiate the reference interpolation directly. In a prism,
# each triangular vertex has a quadratic determinant along the interval; its
# three Bernstein coefficients certify the whole triangle-times-interval map.
function three_strip_map(points,msh)
    if msh==4
        @test ThreeStripBase.det(points[2].-points[1],points[3].-points[1],points[4].-points[1])>0
    elseif msh==5
        delta=points[5].-points[1]
        @test all(index->points[index+4].-points[index]==delta,1:4)
        for index in 1:4
            @test ThreeStripBase.det(points[mod1(index+1,4)].-points[index],
                points[mod1(index-1,4)].-points[index],delta)>0
        end
    elseif msh==6
        low_u=points[2].-points[1];high_u=points[5].-points[4]
        low_v=points[3].-points[1];high_v=points[6].-points[4]
        for vertex in 1:3
            direction=points[vertex+3].-points[vertex]
            first=ThreeStripBase.det(low_u,low_v,direction)
            last=ThreeStripBase.det(high_u,high_v,direction)
            middle=ThreeStripBase.det((low_u.+high_u)./2,(low_v.+high_v)./2,direction)
            @test first>0 && last>0 && 2middle-(first+last)/2>0
        end
    elseif msh==7
        @test ThreeStripBase.det(points[2].-points[1],points[3].-points[1],points[4].-points[1])==0
        for index in 1:4
            @test ThreeStripBase.det(points[mod1(index+1,4)].-points[index],
                points[mod1(index-1,4)].-points[index],points[5].-points[index])>0
        end
    else
        @test false # A foreign family violates this declared product.
    end
end

function three_strip_audit(volume,source,layers,recombined)
    primary=8(layers+1)
    @test size(volume.coords)==(3,primary+(recombined ? 3 : 0))
    source_cells=Tuple(Tuple(Int.(cell)) for block in source.blocks
        if block isa ElementBlock && block.msh==3 for cell in eachcol(block.nodes))
    @test length(source_cells)==3 && size(source.coords)==(3,8)
    points=[ntuple(d->source.coords[d,node],3) for node in 1:8]
    columns=Dict((p[1],p[2])=>node for (node,p) in enumerate(points))
    @test length(columns)==8
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
    @test length(grid)==primary && length(extra)==(recombined ? 3 : 0)
    for parent in 1:(recombined ? 3 : 0)
        corners=Tuple(ntuple(d->d==3 ? planes[level] : points[column][d],3)
            for level in (layers,layers+1) for column in source_cells[parent])
        center=ntuple(d->sum(p[d] for p in corners)/8,3)
        matches=[node for node in extra if ntuple(d->volume.coords[d,node],3)==center]
        @test length(matches)==1
        logical[only(matches)]=(-parent,layers)
        @test planes[layers]<center[3]<planes[layers+1]
    end
    parents=Tuple(Set(cell) for cell in source_cells)
    areas=Tuple(abs(sum(points[cell[index]][1]*points[cell[mod1(index+1,4)]][2]-
        points[cell[index]][2]*points[cell[mod1(index+1,4)]][1] for index in 1:4))/2
        for cell in source_cells)
    domain_counts=zeros(Int,3,layers);domain_volumes=zeros(Float64,3,layers)
    families=Dict{Int,Int}();faces=Dict{Tuple,Tuple{Int,Tuple}}()
    used=falses(size(volume.coords,2));signatures=Tuple[]
    for block in volume.blocks
        block isa ElementBlock && msh_dimension(block.msh)==3 || continue
        @test block.msh in (4,5,6,7)
        families[Int(block.msh)]=get(families,Int(block.msh),0)+size(block.nodes,2)
        for cell in eachcol(block.nodes)
            @test length(unique(cell))==length(cell)
            used[cell].=true
            ids=Set(logical[node][1] for node in cell if logical[node][1]>0)
            centers=Set(-logical[node][1] for node in cell if logical[node][1]<0)
            @test length(centers)<=1
            # A fan pyramid on the common face has only shared column corners.
            # Its independently located interior center resolves its macro.
            candidates=[index for index in 1:3 if ids⊆parents[index] &&
                (isempty(centers) || index in centers)]
            @test length(candidates)==1
            parent=only(candidates)
            @test all(logical[node][1]>0 || logical[node][1]==-parent for node in cell)
            levels=[logical[node][2] for node in cell]
            level=isempty(centers) ? minimum(levels) : layers
            @test all(logical[node][1]<0 || logical[node][2] in (level,level+1) for node in cell)
            if isempty(centers)
                @test maximum(levels)==level+1
            else
                # A bottom-face pyramid ends at the actual interior centroid,
                # so it covers part of the terminal interval, not both planes.
                @test minimum(volume.coords[3,node] for node in cell)<
                    maximum(volume.coords[3,node] for node in cell)
            end
            domain_counts[parent,level]+=1
            p=Tuple(ntuple(d->volume.coords[d,node],3) for node in cell)
            three_strip_map(p,Int(block.msh))
            # Bilinear quadrangle flux integrates the represented curved face.
            # A fixed hidden diagonal would mismeasure warped Prism6 faces.
            signed=sum(ThreeStripQuad.face_volume(volume,Tuple(cell[index] for index in face),p[1])
                for face in ThreeStripQuad.FACES[Int(block.msh)])
            @test isfinite(signed) && signed>0
            domain_volumes[parent,level]+=signed
            for pattern in ThreeStripQuad.FACES[Int(block.msh)]
                face=Tuple(Int(cell[index]) for index in pattern)
                key=Tuple(sort!(collect(face)));cycle=ThreeStripBase.cycle(face)
                previous=get(faces,key,nothing)
                if previous===nothing
                    faces[key]=(1,cycle)
                else
                    @test previous[1]==1 && previous[2]==ThreeStripBase.cycle(reverse(face))
                    faces[key]=(2,previous[2])
                end
            end
            push!(signatures,(Int(block.msh),Tuple(sort!(collect(p)))))
        end
    end
    @test all(used)
    if recombined
        expected=Dict(4=>6,7=>15)
        layers>1 && (expected[5]=3(layers-1))
        @test families==expected
        @test length(faces)==13layers+45
    else
        @test !(5 in keys(families)) && sum(values(families))<=18layers-3
    end
    for parent in 1:3,level in 1:layers
        @test recombined ? domain_counts[parent,level]==(level==layers ? 7 : 1) :
            1<=domain_counts[parent,level]<=6
        @test isapprox(domain_volumes[parent,level],areas[parent]*(planes[level+1]-planes[level]);rtol=64eps(),atol=0)
    end
    @test count(value->value[1]==1,values(faces))==(recombined ? 8layers+9 : 16layers+9)
    @test isapprox(sum(domain_volumes),1.;rtol=128eps())
    io=IOBuffer();show(io,sort!(signatures))
    return bytes2hex(sha256(take!(io)))
end

function three_strip_profile(make,extract,source,layers,recombined)
    result=make()
    allocated=minimum(@allocated(make()) for _ in 1:3)
    seconds=minimum(@elapsed(make()) for _ in 1:3)
    volume=extract(result)
    digest=three_strip_audit(volume,source,layers,recombined)
    return (;layers,allocated,seconds,digest,nodes=size(volume.coords,2))
end

# The measured entry products are P1. Actual P2 support identity, carrier
# ownership and public map/parameter queries are audited separately below.
# Lower-cell carriers come from the actual classified P1 projection, rather
# than an endpoint-owner heuristic or a guessed subdivision of CAD faces.
three_strip_support_key(nodes)=Tuple(sort!(Int.(collect(nodes))))

function three_strip_actual_carriers(api,projected)
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
                edges[three_strip_support_key(remap[node] for node in cell)]=carrier
            else
                for i in eachindex(cell)
                    support=three_strip_support_key((remap[cell[i]],remap[cell[mod1(i+1,length(cell))]]))
                    haskey(edges,support) && edges[support][1]<dimension && continue
                    @test !haskey(edges,support) || edges[support]==carrier
                    edges[support]=carrier
                end
                block.msh==3 && (quads[three_strip_support_key(remap[node] for node in cell)]=carrier)
            end
        end
    end
    return (;edges,quads)
end

function three_strip_surface_query(api,surface,boundary,expected)
    nodes,coordinates,parameters=api.mesh.get_nodes(2,surface,boundary,true)
    @test length(nodes)==expected && allunique(nodes)
    @test length(coordinates)==3expected && length(parameters)==2expected
    @test maximum(abs.(api.model.get_value(2,surface,parameters).-coordinates))<=2e-11
    return Set(nodes)
end

function three_strip_quadratic_audit(api,projected,layers,recombined)
    expected=three_strip_actual_carriers(api,projected)
    model=api.CURRENT[]
    surfaces=ThreeStripModel._model_volume_boundary_surfaces(model,1,"three-strip resource P2 audit")
    api.mesh.set_order(2)
    nodes,xyz,_=api.mesh.get_nodes(3,1,true)
    center_count=recombined ? 3 : 0
    @test length(nodes)==42layers+21+8center_count
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
    @test counts==[8,8layers+20,24layers-2,10layers-5+8center_count]
    types,element_tags,connectivity=api.mesh.get_elements(3,1)
    if recombined
        families=Dict(Int(msh)=>length(tags) for (msh,tags) in zip(types,element_tags))
        exact=Dict(11=>6,14=>15)
        layers>1 && (exact[12]=3(layers-1))
        @test families==exact
    end
    supports=Dict{UInt64,Tuple}()
    for (msh,flat) in zip(types,connectivity)
        for edge in eachcol(reshape(api.mesh.get_element_edge_nodes(msh,1,false),3,:))
            support=three_strip_support_key(edge[1:2]);node=edge[3]
            @test !haskey(supports,node) || supports[node]==support
            supports[node]=support
        end
        if msh in (12,13,14)
            for face in eachcol(reshape(api.mesh.get_element_face_nodes(msh,4,1,false),9,:))
                support=three_strip_support_key(face[1:4]);node=face[9]
                @test !haskey(supports,node) || supports[node]==support
                supports[node]=support
            end
        end
        if msh==12
            for cell in eachcol(reshape(flat,27,:))
                supports[cell[27]]=three_strip_support_key(cell[1:8])
            end
        end
        reference,_=api.mesh.get_integration_points(msh,"Gauss3")
        jacobians,determinants,coordinates=api.mesh.get_jacobians(msh,reference,1)
        @test all(isfinite,jacobians) && all(>(0),determinants) && all(isfinite,coordinates)
        @test length(jacobians)==9length(determinants) && length(coordinates)==3length(determinants)
    end
    @test length(supports)==34layers+13+7center_count
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
        long=link!==nothing && link[1]==1 && abs(link[2]) in (1,3)
        lateral=link!==nothing && link[1]==1
        owned=lateral ? (long ? 10layers-5 : 2layers-1) : 5
        closure=lateral ? (long ? 14layers+7 : 6layers+3) : 21
        own=three_strip_surface_query(api,surface,false,owned)
        boundary=three_strip_surface_query(api,surface,true,closure)
        @test own⊆boundary
    end
    println("THREE_QUAD_STRIP_P2_AUDIT layers=",layers," policy=",recombined,
            " nodes=",length(nodes)," supports=",length(supports)," owners=",counts)
    signatures=[(length(support),Tuple(sort!([points[primary] for primary in support])),
                 points[node],owners[node]) for (node,support) in supports]
    buffer=IOBuffer();show(buffer,sort!(signatures))
    return Dict{String,Any}("layers"=>layers,"policy"=>recombined ? "recombined" : "free",
        "nodes"=>length(nodes),"supports"=>length(supports),"owners"=>counts,
        "digest"=>bytes2hex(sha256(take!(buffer))))
end

three_strip_records=Dict{String,Any}[]
three_strip_p2_records=Dict{String,Any}[]
@testset "Three-quadrangle strip allocation and actual complete product" begin
    @test three_strip_boxes(Base.uncompressed_ast(first(methods(three_strip_box_control))))>0
    for method in methods(three_strip_keyword_box_control)
        @test sum(body->three_strip_boxes(Base.uncompressed_ast(body)),
            three_strip_lowered_methods(method);init=0)>0
    end
    # An independent non-column prism control reaches the resource helper's
    # warped-face arm even when the deterministic unit-strip policy uses none.
    # Its reference determinant is 1-s^2/64 on triangle times s in [0,1],
    # giving volume (1/2)*(1-1/192)=191/384. A quad diagonal cannot supply it.
    warped=((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),
            (0.,0.,1.),(1.,.125,1.),(.125,1.,1.))
    three_strip_map(warped,6)
    control=MixedMesh(hcat((collect(p) for p in warped)...),
        [ElementBlock(6,reshape(Int32.(1:6),6,1))])
    flux=sum(ThreeStripQuad.face_volume(control,face,warped[1]) for face in ThreeStripQuad.FACES[6])
    @test isapprox(flux,191/384;rtol=16eps(),atol=0)
    @test any(length(face)==4 && ThreeStripBase.det(warped[face[2]].-warped[face[1]],
        warped[face[3]].-warped[face[1]],warped[face[4]].-warped[face[1]])!=0
        for face in ThreeStripQuad.FACES[6])
    for recombined in (false,true)
        previous=Dict{Symbol,Any}()
        for layers in (1000,2000,4000)
            fixture=ThreeStripBase.fixture("resource";laterals=recombined,layers=:one,pins=(1,2,3,4))
            text=replace(fixture.source,"Transfinite Curve{1,2,3,4}=2;"=>
                "Transfinite Curve{1,3}=4;Transfinite Curve{2,4}=2;Recombine Surface{1};",
                "Layers{1}"=>"Layers{$layers}")
            mktempdir() do directory
                path=joinpath(directory,"strip.geo");write(path,text)
                geometry=execute_geo(path;mesh_dim=0)
                source=mesh_model_surface(geometry.model,1)
                base=mesh_model_volume(geometry.model,1)
                paths=((:standalone,()->mesh_model_volume(geometry.model,1),identity),
                    (:geo,()->execute_geo(path;mesh_dim=3),result->geo_entity_mesh(result,3,1)),
                    (:projection,()->model_to_mixed(geometry.model,base,3,1),identity))
                expected=nothing
                for (name,make,extract) in paths
                    row=three_strip_profile(make,extract,source,layers,recombined)
                    expected===nothing ? (expected=row.digest) : (@test row.digest==expected)
                    haskey(previous,name) && (@test row.allocated<=2.15previous[name].allocated+65536)
                    previous[name]=row
                    push!(three_strip_records,Dict{String,Any}("policy"=>recombined ? "recombined" : "free",
                        "path"=>String(name),(String(key)=>value for (key,value) in pairs(row))...))
                    println("THREE_QUAD_STRIP_RESOURCE policy=",recombined," path=",name," ",row)
                end
                api=Tessella.API;api.initialize()
                try
                    api.open_geo!(path;mesh_dim=0)
                    row=three_strip_profile(()->api.mesh.generate(3),_->ThreeStripBase.public_volume(api),source,layers,recombined)
                    @test row.digest==expected
                    haskey(previous,:api) && (@test row.allocated<=2.15previous[:api].allocated+65536)
                    previous[:api]=row
                    push!(three_strip_records,Dict{String,Any}("policy"=>recombined ? "recombined" : "free",
                        "path"=>"api",(String(key)=>value for (key,value) in pairs(row))...))
                    println("THREE_QUAD_STRIP_RESOURCE policy=",recombined," path=api ",row)
                    projected=model_to_mixed(api.CURRENT[],ThreeStripBase.public_volume(api),3,1)
                    push!(three_strip_p2_records,three_strip_quadratic_audit(api,projected,layers,recombined))
                finally
                    api.finalize()
                end
            end
        end
    end
    checked_methods=Set{Method}()
    for name in names(ThreeStripModel;all=true)
        (startswith(String(name),"_extrude_nonew_three_quad_strip") ||
         name in (:_extrude_nonew_complete_plan,:_extrude_nonew_face_cycle,:_extrude_nonew_face_orientation,
                  :_extrude_nonew_count_faces!,:_extrude_nonew_projection_context,
                  :_extrude_nonew_cell_counts,:_extrude_nonew_projection_lookup,
                  :_extrude_quadtri_centroid,:_extrude_quadtri_centroid_exact,
                  :_model_curve_transfinite_native_params,:_transfinite_val)) || continue
        value=getfield(ThreeStripModel,name);value isa Function || continue
        for method in methods(value)
            union!(checked_methods,three_strip_lowered_methods(method))
        end
    end
    for method in methods(Tessella.GeoExec._geo_merge_entity_meshes_mixed)
        union!(checked_methods,three_strip_lowered_methods(method))
    end
    for method in checked_methods
        boxes=three_strip_boxes(Base.uncompressed_ast(method))
        boxes==0 || println("THREE_QUAD_STRIP_RESOURCE_BOX_METHOD=",method," boxes=",boxes)
        @test boxes==0
    end
    println("THREE_QUAD_STRIP_RESOURCE_BOX_METHODS=",length(checked_methods))
    final=three_strip_snapshot()
    changed=sort!(collect(path for path in union(keys(THREE_STRIP_START),keys(final))
        if get(THREE_STRIP_START,path,nothing)!=get(final,path,nothing)))
    isempty(changed) || println("THREE_QUAD_STRIP_RESOURCE_CHANGED_INPUTS=",changed)
    @test isempty(changed)
end
if !isempty(ARGS)
    open(ARGS[1],"w") do io
        TOML.print(io,Dict("julia"=>string(VERSION),"records"=>three_strip_records,
            "p2_records"=>three_strip_p2_records,"inputs"=>THREE_STRIP_START))
    end
end
println("THREE_QUAD_STRIP_RESOURCE_OK records=",length(three_strip_records)," inputs_stable=",three_strip_snapshot()==THREE_STRIP_START)
