# Normal compilation and --check-bounds=yes; no Gmsh dependency.
using Tessella, Test, SHA, TOML
using Tessella.Model: mesh_model_volume, mesh_model_surface, model_to_mixed
using Tessella.Elements: ElementBlock, MixedMesh, msh_dimension
include(joinpath(@__DIR__, "../../test/geometry/quadtri_nonew_two_tri_certificates.jl"))
const PatchBase = QuadTriNoNewTwoTriCertificates
const PatchFaces = QuadTriNoNewCertificates.FACES
const PatchModel = Tessella.Model
const PATCH_ROOT = normpath(joinpath(@__DIR__, "../.."))
const PATCH_INPUTS = sort!(vcat([joinpath(PATCH_ROOT,"Project.toml")],
    [joinpath(directory,file) for (directory,_,files) in walkdir(joinpath(PATCH_ROOT,"src"))
     for file in files if endswith(file,".jl")],
    [@__FILE__,joinpath(PATCH_ROOT,"test/geometry/quadtri_nonew_two_tri_certificates.jl"),
     joinpath(PATCH_ROOT,"test/geometry/quadtri_nonew_certificates.jl")]))
patch_snapshot() = Dict(relpath(path,PATCH_ROOT)=>bytes2hex(sha256(read(path))) for path in PATCH_INPUTS)
const PATCH_START = patch_snapshot()
println("START QUAD_PATCH_RESOURCE Julia=",VERSION," inputs=",length(PATCH_START))

function patch_boxes(value)
    value isa Core.Box && return 1
    value isa GlobalRef && return value==GlobalRef(Core,:Box) ? 1 : 0
    value isa Core.CodeInfo && return sum(patch_boxes,value.code;init=0)
    value isa Expr && return sum(patch_boxes,value.args;init=0)
    value isa Core.NewvarNode && return patch_boxes(value.slot)
    value isa QuoteNode && return patch_boxes(value.value)
    return value===Core.Box ? 1 : 0
end

function patch_box_control()
    captured=0
    read=()->captured
    captured=1
    return read
end

# Independent outward-face integration, anchored at an actual cell corner.
function patch_signed_volume(points,msh)
    anchor=points[1]
    result=0.0
    for face in PatchFaces[msh],index in 2:length(face)-1
        result+=PatchBase.det(points[face[1]].-anchor,
            points[face[index]].-anchor,points[face[index+1]].-anchor)/6
    end
    return result
end

# The resource geometry consists of four rectangular source macros. These
# reference-map checks use the actual stored corners, independently of factories.
function patch_map(points,msh)
    if msh==4
        @test PatchBase.det(points[2].-points[1],points[3].-points[1],points[4].-points[1])>0
    elseif msh==5
        delta=points[5].-points[1]
        @test all(index->points[index+4].-points[index]==delta,1:4)
        for index in 1:4
            @test PatchBase.det(points[mod1(index+1,4)].-points[index],
                points[mod1(index-1,4)].-points[index],delta)>0
        end
    elseif msh==7
        @test PatchBase.det(points[2].-points[1],points[3].-points[1],points[4].-points[1])==0
        for index in 1:4
            @test PatchBase.det(points[mod1(index+1,4)].-points[index],
                points[mod1(index-1,4)].-points[index],points[5].-points[index])>0
        end
    else
        @test false # A foreign family violates this declared product.
    end
end

function patch_audit(volume,source,layers,recombined)
    @test size(volume.coords)==(3,9(layers+1))
    source_cells=Tuple(Tuple(Int.(cell)) for block in source.blocks
        if block isa ElementBlock && block.msh==3 for cell in eachcol(block.nodes))
    @test length(source_cells)==4 && size(source.coords)==(3,9)
    points=[ntuple(d->source.coords[d,node],3) for node in 1:9]
    columns=Dict((p[1],p[2])=>node for (node,p) in enumerate(points))
    @test length(columns)==9
    planes=sort!(unique(volume.coords[3,:]))
    @test length(planes)==layers+1 && first(planes)==0. && last(planes)==1.
    heights=Dict(height=>level for (level,height) in enumerate(planes))
    logical=Vector{Tuple{Int,Int}}(undef,size(volume.coords,2))
    grid=Set{Tuple{Int,Int}}()
    for node in axes(volume.coords,2)
        xy=(volume.coords[1,node],volume.coords[2,node])
        @test haskey(columns,xy) && haskey(heights,volume.coords[3,node])
        logical[node]=(columns[xy],heights[volume.coords[3,node]])
        @test !(logical[node] in grid)
        push!(grid,logical[node])
    end
    parents=Tuple(Set(cell) for cell in source_cells)
    areas=Tuple(abs(sum(points[cell[index]][1]*points[cell[mod1(index+1,4)]][2]-
        points[cell[index]][2]*points[cell[mod1(index+1,4)]][1] for index in 1:4))/2
        for cell in source_cells)
    domain_counts=zeros(Int,4,layers)
    domain_volumes=zeros(Float64,4,layers)
    families=Dict{Int,Int}()
    faces=Dict{Tuple,Tuple{Int,Tuple}}()
    used=falses(size(volume.coords,2))
    signatures=Tuple[]
    for block in volume.blocks
        block isa ElementBlock && msh_dimension(block.msh)==3 || continue
        @test block.msh in (4,5,7)
        families[Int(block.msh)]=get(families,Int(block.msh),0)+size(block.nodes,2)
        for cell in eachcol(block.nodes)
            @test length(unique(cell))==length(cell)
            used[cell].=true
            ids=Set(logical[node][1] for node in cell)
            parent=only(findall(candidate->ids⊆candidate,parents))
            levels=[logical[node][2] for node in cell]
            level=minimum(levels)
            @test maximum(levels)==level+1
            domain_counts[parent,level]+=1
            p=Tuple(ntuple(d->volume.coords[d,node],3) for node in cell)
            patch_map(p,Int(block.msh))
            signed=patch_signed_volume(p,Int(block.msh))
            @test isfinite(signed) && signed>0
            domain_volumes[parent,level]+=signed
            for pattern in PatchFaces[Int(block.msh)]
                face=Tuple(Int(cell[index]) for index in pattern)
                key=Tuple(sort!(collect(face)))
                cycle=PatchBase.cycle(face)
                previous=get(faces,key,nothing)
                if previous===nothing
                    faces[key]=(1,cycle)
                else
                    @test previous[1]==1 && previous[2]==PatchBase.cycle(reverse(face))
                    faces[key]=(2,previous[2])
                end
            end
            push!(signatures,(Int(block.msh),Tuple(sort!(collect(p)))))
        end
    end
    @test all(used)
    expected=recombined ? Dict(7=>12) : Dict(4=>24layers-8,7=>4)
    recombined && layers>1 && (expected[5]=4(layers-1))
    @test families==expected
    for parent in 1:4,level in 1:layers
        expected_count=recombined ? (level==layers ? 3 : 1) : (level==1 ? 5 : 6)
        @test domain_counts[parent,level]==expected_count
        @test isapprox(domain_volumes[parent,level],areas[parent]*(planes[level+1]-planes[level]);rtol=32eps(),atol=0)
    end
    @test length(faces)==(recombined ? 16layers+24 : 56layers)
    @test count(value->value[1]==1,values(faces))==(recombined ? 8layers+12 : 16layers+12)
    @test count(value->value[1]==2,values(faces))==(recombined ? 8layers+12 : 40layers-12)
    @test isapprox(sum(domain_volumes),1.;rtol=64eps())
    io=IOBuffer();show(io,sort!(signatures))
    return bytes2hex(sha256(take!(io)))
end

function patch_profile(make,extract,source,layers,recombined)
    result=make()
    allocated=minimum(@allocated(make()) for _ in 1:3)
    seconds=minimum(@elapsed(make()) for _ in 1:3)
    volume=extract(result)
    digest=patch_audit(volume,source,layers,recombined)
    return (;layers,allocated,seconds,digest,nodes=size(volume.coords,2))
end

patch_records=Dict{String,Any}[]
@testset "2-by-2 quadrangle allocation and actual complete product" begin
    @test patch_boxes(Base.uncompressed_ast(first(methods(patch_box_control))))>0
    for recombined in (false,true)
        previous=Dict{Symbol,Any}()
        for layers in (1000,2000,4000)
            fixture=PatchBase.fixture("resource";laterals=recombined,layers=:one,pins=(1,2,3,4))
            text=replace(fixture.source,"Transfinite Curve{1,2,3,4}=2;"=>
                "Transfinite Curve{1,2,3,4}=3;Recombine Surface{1};","Layers{1}"=>"Layers{$layers}")
            mktempdir() do directory
                path=joinpath(directory,"patch.geo");write(path,text)
                geometry=execute_geo(path;mesh_dim=0)
                source=mesh_model_surface(geometry.model,1)
                base=mesh_model_volume(geometry.model,1)
                paths=((:standalone,()->mesh_model_volume(geometry.model,1),identity),
                    (:geo,()->execute_geo(path;mesh_dim=3),result->geo_entity_mesh(result,3,1)),
                    (:projection,()->model_to_mixed(geometry.model,base,3,1),identity))
                expected=nothing
                for (name,make,extract) in paths
                    row=patch_profile(make,extract,source,layers,recombined)
                    expected===nothing ? (expected=row.digest) : (@test row.digest==expected)
                    haskey(previous,name) && (@test row.allocated<=2.15previous[name].allocated+65536)
                    previous[name]=row
                    push!(patch_records,Dict{String,Any}("policy"=>recombined ? "recombined" : "free",
                        "path"=>String(name),(String(key)=>value for (key,value) in pairs(row))...))
                    println("QUAD_PATCH_RESOURCE policy=",recombined," path=",name," ",row)
                end
                api=Tessella.API;api.initialize()
                try
                    api.open_geo!(path;mesh_dim=0)
                    row=patch_profile(()->api.mesh.generate(3),_->PatchBase.public_volume(api),source,layers,recombined)
                    @test row.digest==expected
                    haskey(previous,:api) && (@test row.allocated<=2.15previous[:api].allocated+65536)
                    previous[:api]=row
                    push!(patch_records,Dict{String,Any}("policy"=>recombined ? "recombined" : "free",
                        "path"=>"api",(String(key)=>value for (key,value) in pairs(row))...))
                    println("QUAD_PATCH_RESOURCE policy=",recombined," path=api ",row)
                finally
                    api.finalize()
                end
            end
        end
    end
    checked=0
    for name in names(PatchModel;all=true)
        (startswith(String(name),"_extrude_nonew_quad_patch") ||
         name in (:_extrude_nonew_face_cycle,:_extrude_nonew_face_orientation,
                  :_extrude_nonew_count_faces!,:_extrude_nonew_projection_context,
                  :_extrude_nonew_cell_counts,:_extrude_nonew_projection_lookup)) || continue
        value=getfield(PatchModel,name);value isa Function || continue
        for method in methods(value)
            @test patch_boxes(Base.uncompressed_ast(method))==0
            checked+=1
        end
    end
    println("QUAD_PATCH_RESOURCE_BOX_METHODS=",checked)
    final=patch_snapshot()
    changed=sort!(collect(path for path in union(keys(PATCH_START),keys(final))
        if get(PATCH_START,path,nothing)!=get(final,path,nothing)))
    isempty(changed) || println("QUAD_PATCH_RESOURCE_CHANGED_INPUTS=",changed)
    @test isempty(changed)
end
if !isempty(ARGS)
    open(ARGS[1],"w") do io
        TOML.print(io,Dict("julia"=>string(VERSION),"records"=>patch_records,"inputs"=>PATCH_START))
    end
end
println("QUAD_PATCH_RESOURCE_OK records=",length(patch_records)," inputs_stable=",patch_snapshot()==PATCH_START)
