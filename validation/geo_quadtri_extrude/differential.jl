#!/usr/bin/env julia
# Gmsh is a pinned differential oracle only; production uses native kernels.
using Tessella
using Tessella.Elements: ElementBlock, MixedMesh
using Tessella.MeshTypes: nnodes, ntris, ntets, validate

function find_gmsh_api()
    configured=get(ENV,"GMSH_JULIA_API","")
    isempty(configured) || return configured
    candidates=("/opt/homebrew/lib/gmsh.jl","/usr/local/lib/gmsh.jl",
        joinpath(get(ENV,"LOCALAPPDATA",""),"Programs","Python","Python312","Lib","gmsh.jl"))
    for path in candidates
        isfile(path) && return path
    end
    return ""
end
const GMSH_API=find_gmsh_api()
isfile(GMSH_API) || error("set GMSH_JULIA_API to the pinned Gmsh 4.15.2 Julia API")
include(GMSH_API)
gmsh.GMSH_API_VERSION=="4.15.2" || error("expected Gmsh 4.15.2 API")

function square(n;quad=true)
    """
    Point(1)={0,0,0,1};Point(2)={1,0,0,1};
    Point(3)={1,1,0,1};Point(4)={0,1,0,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
    Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
    Transfinite Curve{:}=$(n+1);Transfinite Surface{1};
    $(quad ? "Recombine Surface{1};" : "")
    """
end

counts(mesh::MixedMesh)=Dict(Int(b.msh)=>size(b.nodes,2) for b in mesh.blocks if b isa ElementBlock)
counts(mesh::Mesh)=Dict(4=>ntets(mesh))
function coordinate_distance(a,b)
    maximum(minimum(sum((a[k,i]-b[k,j])^2 for k in 1:3) for j in axes(b,2)) for i in axes(a,2))
end
function coordinate_ids(points,oracle)
    ids=Int[]
    for i in axes(points,2)
        distances=[sum((points[k,i]-oracle[k,j])^2 for k in 1:3) for j in axes(oracle,2)]
        index=argmin(distances)
        distances[index]<=4e-22 || error("cell vertex does not match an oracle node")
        push!(ids,index)
    end
    return ids
end
function native_cell_vertices(part,dim,oracle)
    ids=coordinate_ids(part.coords,oracle)
    cells=Dict{Int,Vector{Tuple}}()
    blocks=part isa Mesh ? [ElementBlock(dim==2 ? 2 : 4,dim==2 ? part.tris : part.tets)] : part.blocks
    for b in blocks
        b isa ElementBlock || continue
        target=get!(cells,Int(b.msh),Tuple[])
        for i in axes(b.nodes,2)
            push!(target,Tuple(sort!([ids[node] for node in b.nodes[:,i]])))
        end
        sort!(target)
    end
    filter!(p->!isempty(last(p)),cells)
    return cells
end

cases=Tuple{String,String}[]
for n in (1,2,3),layers in (1,3),laterals in (false,true)
    push!(cases,("quad_n$(n)_l$(layers)_lat$(laterals)",square(n)*"""
    Extrude{0,0,1}{Surface{1};Layers{$layers};Recombine;
    QuadTriAddVerts$(laterals ? " RecombLaterals" : "");}
    """))
end
for laterals in (false,true)
    push!(cases,("tri_lat$(laterals)",square(2;quad=false)*"""
    Extrude{0,0,1}{Surface{1};Layers{3};Recombine;
    QuadTriAddVerts$(laterals ? " RecombLaterals" : "");}
    """))
end
for (name,transform) in (("graded","Extrude{0,0,2}"),
        ("rotation","Extrude{{0,1,0},{-1,0,0},Pi/3}"),
        ("twist","Extrude{{0,0,2},{0,0,1},{0,0,0},Pi/3}")),laterals in (false,true)
    push!(cases,(name*"_lat$(laterals)",square(2)*transform*"""
    {Surface{1};Layers{{2,3},{0.4,1.0}};Recombine;
    QuadTriAddVerts$(laterals ? " RecombLaterals" : "");}
    """))
end
for (name,axis) in (("two_fixed_columns","0,1,0"),("one_fixed_column","1,-1,0")),laterals in (false,true)
    push!(cases,(name*"_lat$(laterals)",square(1)*"""
    Extrude{{$axis},{0,0,0},Pi/3}{Surface{1};Layers{3};Recombine;
    QuadTriAddVerts$(laterals ? " RecombLaterals" : "");}
    """))
end
for (n,layers,stages) in ((1,1,(1,2,3,4)),(2,2,(1,2,3,4)),(2,2,(1,3)),(2,2,(2,))),laterals in (false,true)
    source=square(n)
    source=replace(source,"{0,0,0,1}"=>"{1,0,0,1}","{1,0,0,1}"=>"{2,0,0,1}",
        "{1,1,0,1}"=>"{2,0,1,1}","{0,1,0,1}"=>"{1,0,1,1}")
    for stage in 1:4
        from=stage==1 ? "1" : "q$(stage-1)[0]"
        quadtri=stage in stages ? "QuadTriAddVerts$(laterals ? " RecombLaterals" : "");" : ""
        source*="q$stage[]=Extrude{{0,0,1},{0,0,0},Pi/2}{Surface{$from};Layers{$layers};Recombine;$quadtri};\n"
    end
    push!(cases,("toroidal_n$(n)_l$(layers)_qt$(join(stages))_lat$(laterals)",source))
end

gmsh.initialize(["-nopopup"],false,false)
gmsh.option.setNumber("General.Terminal",0)
try
    runtime=gmsh.option.getString("General.Version")
    startswith(runtime,"4.15.2") || error("expected Gmsh 4.15.2 runtime, got $runtime")
    selected=get(ENV,"QUADTRI_CASE","")
    complete_count=Ref(0)
    for (name,source) in cases
        isempty(selected) || occursin(selected,name) || continue
        mktempdir() do directory
            path=joinpath(directory,name*".geo")
            write(path,source*"Mesh 3;\n")
            native=execute_geo(path)
            gmsh.clear();gmsh.open(path)
            gn,xyz,_=gmsh.model.mesh.getNodes()
            gc=reshape(xyz,3,:)
            if nnodes(native.mesh)!=length(gn)
                groups=Dict{NTuple{3,Float64},Vector{Int}}()
                for i in eachindex(gn)
                    push!(get!(groups,Tuple(round.(gc[:,i];digits=10)),Int[]),i)
                end
                for (p,indices) in groups
                    length(indices)>1 || continue
                    println("MISMATCH coincident $p refs=$([(gn[i],gmsh.model.mesh.getNode(gn[i])[3:4],Tuple(gc[:,i])) for i in indices])")
                end
                for (dim,tag) in gmsh.model.getEntities()
                    dim>=2 || continue
                    tags,_xyz,_=gmsh.model.mesh.getNodes(dim,tag)
                    ty,el,_=gmsh.model.mesh.getElements(dim,tag)
                    part=geo_entity_mesh(native,dim,tag)
                    println("MISMATCH entity=($dim,$tag) gmsh_nodes=$(length(tags)) gmsh_counts=$([(t,length(e)) for (t,e) in zip(ty,el)]) native_counts=$(part isa Mesh && dim==2 ? Dict(2=>ntris(part)) : counts(part))")
                end
                for (label,a,b) in (("gmsh",gc,native.mesh.coords),("native",native.mesh.coords,gc))
                    for i in axes(a,2)
                        distance=minimum(sum((a[k,i]-b[k,j])^2 for k in 1:3) for j in axes(b,2))
                        distance>4e-22 && println("MISMATCH $label coordinate=$(Tuple(a[:,i])) distance=$(sqrt(distance))")
                    end
                end
                error("$name node count mismatch: native=$(nnodes(native.mesh)) gmsh=$(length(gn))")
            end
            nd=coordinate_distance(native.mesh.coords,gc)
            gd=coordinate_distance(gc,native.mesh.coords)
            max(nd,gd)<=4e-22 || error("$name node coordinates differ: $(sqrt(max(nd,gd)))")
            canonical=coordinate_ids(gc,gc)
            oracle_ids=Dict(tag=>canonical[i] for (i,tag) in enumerate(gn))
            for (dim,tag) in gmsh.model.getEntities()
                dim in (2,3) || continue
                ty,el,connectivity=gmsh.model.mesh.getElements(dim,tag)
                expected=Dict(Int(t)=>length(e) for (t,e) in zip(ty,el))
                part=geo_entity_mesh(native,dim,tag)
                actual=part isa Mesh && dim==2 ? Dict(2=>ntris(part)) : counts(part)
                filter!(p->last(p)>0,actual)
                actual==expected || error("$name entity ($dim,$tag): $actual != $expected")
                oracle_cells=Dict{Int,Vector{Tuple}}()
                for (type,nodes) in zip(ty,connectivity)
                    width=gmsh.model.mesh.getElementProperties(type)[4]
                    records=Tuple[]
                    for cell in eachcol(reshape(nodes,width,:))
                        push!(records,Tuple(sort!([oracle_ids[node] for node in cell])))
                    end
                    oracle_cells[Int(type)]=sort!(records)
                end
                native_cell_vertices(part,dim,gc)==oracle_cells ||
                    error("$name entity ($dim,$tag): typed cell vertex sets differ")
            end
            validate(native.mesh).ok || error("$name invalid native output")
            complete_count[]+=1
            println("QUADTRI_CASE_OK name=$name nodes=$(length(gn)) max_distance=$(sqrt(max(nd,gd)))")
        end
    end
    println("GEO_QUADTRI_ADDVERTS_DIFFERENTIAL_OK gmsh=4.15.2 cases=$(complete_count[])")
finally
    gmsh.finalize()
end
