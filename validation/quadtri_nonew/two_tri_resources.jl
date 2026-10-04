# Run with normal compilation and --check-bounds=yes. Gmsh is not required.
using Tessella, Test, SHA, TOML
using Tessella.Model: mesh_model_volume, mesh_model_surface, model_to_mixed
using Tessella.MeshTypes: Mesh
using Tessella.Elements: MixedMesh, ElementBlock, msh_dimension
include(joinpath(@__DIR__, "../../test/geometry/quadtri_nonew_two_tri_certificates.jl"))
const TwoTri = QuadTriNoNewTwoTriCertificates
const Model = Tessella.Model
const ROOT = normpath(joinpath(@__DIR__, "../.."))
const INPUTS = sort!(vcat([joinpath(ROOT,"Project.toml")],
    [joinpath(directory,file) for (directory,_,files) in walkdir(joinpath(ROOT,"src"))
     for file in files if endswith(file,".jl")],
    [@__FILE__, joinpath(ROOT,"test/geometry/quadtri_nonew_two_tri_certificates.jl"),
     joinpath(ROOT,"test/geometry/quadtri_nonew_certificates.jl")]))
snapshot() = Dict(relpath(path,ROOT)=>bytes2hex(sha256(read(path))) for path in INPUTS)
const START = snapshot()
println("START TWO_TRI_RESOURCE Julia=",VERSION," inputs=",length(START))

function count_boxes(value)
    value isa Core.Box && return 1
    value isa Core.CodeInfo && return sum(count_boxes,value.code;init=0)
    value isa Expr && return sum(count_boxes,value.args;init=0)
    value isa Core.NewvarNode && return count_boxes(value.slot)
    value isa QuoteNode && return count_boxes(value.value)
    return value===Core.Box ? 1 : 0
end

# This linear audit uses actual source identities and actual cell coordinates.
# The small geometry suite separately checks exact rational partition volumes.
function audit(volume,source,layers,recombined)
    @test size(volume.coords)==(3,4(layers+1))
    points = [ntuple(d->source.coords[d,node],3) for node in 1:4]
    columns = Dict((p[1],p[2])=>node for (node,p) in enumerate(points))
    planes = sort!(unique(volume.coords[3,:]))
    @test length(planes)==layers+1 && first(planes)==0. && last(planes)==1.
    heights = Dict(height=>level for (level,height) in enumerate(planes))
    logical = Vector{Tuple{Int,Int}}(undef,size(volume.coords,2))
    grid = Set{Tuple{Int,Int}}()
    for node in axes(volume.coords,2)
        xy = (volume.coords[1,node],volume.coords[2,node])
        @test haskey(columns,xy) && haskey(heights,volume.coords[3,node])
        logical[node] = (columns[xy],heights[volume.coords[3,node]])
        @test !(logical[node] in grid)
        push!(grid,logical[node])
    end
    parents = [Set(Int.(cell)) for cell in eachcol(source.tris)]
    domain_counts = zeros(Int,2,layers)
    domain_volumes = zeros(Float64,2,layers)
    areas = [abs((points[cell[2]][1]-points[cell[1]][1])*(points[cell[3]][2]-points[cell[1]][2])-
        (points[cell[2]][2]-points[cell[1]][2])*(points[cell[3]][1]-points[cell[1]][1]))/2
        for cell in eachcol(source.tris)]
    faces = Dict{Tuple,Tuple{Int,Tuple}}()
    used = falses(size(volume.coords,2))
    signatures = Tuple[]
    ncell = 0
    for block in TwoTri.blocks(volume)
        msh_dimension(block.msh)==3 || continue
        @test block.msh==(recombined ? 6 : 4)
        for cell in eachcol(block.nodes)
            ncell+=1
            @test length(unique(cell))==length(cell)
            used[cell].=true
            ids = Set(logical[node][1] for node in cell)
            parent = only(findall(candidate->ids⊆candidate,parents))
            levels = [logical[node][2] for node in cell]
            level = minimum(levels)
            @test maximum(levels)==level+1
            domain_counts[parent,level]+=1
            p = Tuple(ntuple(d->volume.coords[d,node],3) for node in cell)
            if block.msh==4
                determinant = TwoTri.det(p[2].-p[1],p[3].-p[1],p[4].-p[1])
                @test isfinite(determinant) && determinant>0
                domain_volumes[parent,level]+=determinant/6
            else
                delta = p[4].-p[1]
                @test all(i->p[i+3].-p[i]==delta,1:3)
                determinant = TwoTri.det(p[2].-p[1],p[3].-p[1],delta)
                @test isfinite(determinant) && determinant>0
                domain_volumes[parent,level]+=determinant/2
            end
            for pattern in TwoTri.FACES[Int(block.msh)]
                face = Tuple(Int(cell[index]) for index in pattern)
                key = Tuple(sort!(collect(face)))
                cycle = TwoTri.cycle(face)
                previous = get(faces,key,nothing)
                if previous===nothing
                    faces[key]=(1,cycle)
                else
                    @test previous[1]==1 && previous[2]==TwoTri.cycle(reverse(face))
                    faces[key]=(2,previous[2])
                end
            end
            push!(signatures,(Int(block.msh),Tuple(sort!(collect(p)))))
        end
    end
    @test all(used) && ncell==(recombined ? 2layers : 6layers)
    @test all(==(recombined ? 1 : 3),domain_counts)
    for parent in 1:2,level in 1:layers
        @test isapprox(domain_volumes[parent,level],areas[parent]*(planes[level+1]-planes[level]);rtol=16eps(),atol=0)
    end
    exterior=count(value->value[1]==1,values(faces))
    interior=count(value->value[1]==2,values(faces))
    @test exterior==(recombined ? 4layers+4 : 8layers+4)
    @test interior==(recombined ? 3layers-2 : 8layers-2)
    @test isapprox(sum(domain_volumes),1.;rtol=64eps())
    io=IOBuffer();show(io,sort!(signatures))
    return bytes2hex(sha256(take!(io)))
end

function profile(make,extract,source,layers,recombined)
    result=make()
    allocated=minimum(@allocated(make()) for _ in 1:3)
    seconds=minimum(@elapsed(make()) for _ in 1:3)
    volume=extract(result)
    digest=audit(volume,source,layers,recombined)
    return (;layers,allocated,seconds,digest,nodes=size(volume.coords,2))
end

records=Dict{String,Any}[]
@testset "Two-triangle grid allocation and complete payload" begin
    for recombined in (false,true)
        previous=Dict{Symbol,Any}()
        for layers in (1000,2000,4000)
            fixture=TwoTri.fixture("resource";laterals=recombined,layers=:one,pins=(1,2,3,4))
            text=replace(fixture.source,"Layers{1}"=>"Layers{$layers}")
            mktempdir() do directory
                path=joinpath(directory,"grid.geo");write(path,text)
                geometry=execute_geo(path;mesh_dim=0)
                source=mesh_model_surface(geometry.model,1)
                base=mesh_model_volume(geometry.model,1)
                paths=(
                    (:standalone,()->mesh_model_volume(geometry.model,1),identity),
                    (:geo,()->execute_geo(path;mesh_dim=3),result->geo_entity_mesh(result,3,1)),
                    (:projection,()->model_to_mixed(geometry.model,base,3,1),identity))
                expected=nothing
                for (name,make,extract) in paths
                    row=profile(make,extract,source,layers,recombined)
                    expected===nothing ? (expected=row.digest) : (@test row.digest==expected)
                    haskey(previous,name) && (@test row.allocated<=2.15previous[name].allocated+65536)
                    previous[name]=row
                    push!(records,Dict{String,Any}("policy"=>recombined ? "recombined" : "free",
                        "path"=>String(name),(String(key)=>value for (key,value) in pairs(row))...))
                    println("TWO_TRI_RESOURCE policy=",recombined," path=",name," ",row)
                end
                api=Tessella.API;api.initialize()
                try
                    api.open_geo!(path;mesh_dim=0)
                    row=profile(()->api.mesh.generate(3),_->TwoTri.public_volume(api),source,layers,recombined)
                    @test row.digest==expected
                    haskey(previous,:api) && (@test row.allocated<=2.15previous[:api].allocated+65536)
                    previous[:api]=row
                    push!(records,Dict{String,Any}("policy"=>recombined ? "recombined" : "free",
                        "path"=>"api",(String(key)=>value for (key,value) in pairs(row))...))
                    println("TWO_TRI_RESOURCE policy=",recombined," path=api ",row)
                finally
                    api.finalize()
                end
            end
        end
    end
    checked=0
    for name in names(Model;all=true)
        (startswith(String(name),"_extrude_nonew_two_tri") ||
         name in (:_extrude_nonew_face_cycle,:_extrude_nonew_face_orientation,
                  :_extrude_nonew_count_faces!)) || continue
        value=getfield(Model,name);value isa Function || continue
        for method in methods(value)
            @test count_boxes(Base.uncompressed_ast(method))==0
            checked+=1
        end
    end
    println("TWO_TRI_RESOURCE_BOX_METHODS=",checked)
    final=snapshot()
    changed=sort!(collect(path for path in union(keys(START),keys(final))
        if get(START,path,nothing)!=get(final,path,nothing)))
    isempty(changed) || println("TWO_TRI_RESOURCE_CHANGED_INPUTS=",changed)
    @test isempty(changed)
end
if !isempty(ARGS)
    open(ARGS[1],"w") do io
        TOML.print(io,Dict("julia"=>string(VERSION),"records"=>records,"inputs"=>START))
    end
end
println("TWO_TRI_RESOURCE_OK records=",length(records)," inputs_stable=",snapshot()==START)
