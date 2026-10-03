#!/usr/bin/env julia
# Native API mixed-family refinement against the pinned Gmsh 4.15.2 templates.
using Tessella
using Tessella.Elements: ElementBlock, MixedMesh, lagrange_nodes, msh_dimension
using Tessella.MeshTypes: nnodes, validate

api_path=get(ENV,"GMSH_JULIA_API",joinpath(get(ENV,"LOCALAPPDATA",""),
    "Programs","Python","Python312","Lib","gmsh.jl"))
isfile(api_path) || error("set GMSH_JULIA_API to the pinned Gmsh 4.15.2 API")
include(api_path)
gmsh.GMSH_API_VERSION=="4.15.2" || error("expected Gmsh 4.15.2 API")

function install_native(mesh)
    api=Tessella.API
    api.initialize()
    dim=maximum(msh_dimension(b.msh) for b in mesh.blocks)
    api.model.add_discrete_entity(dim,7)
    owners=Dict(Int(b.msh)=>fill(Int32(7),size(b.nodes,2)) for b in mesh.blocks)
    class=api._mixed_classification(mesh,(dim,Int32(7)),[(dim,Int32(7))],
        fill((dim,Int32(7)),nnodes(mesh)),Dict((dim,Int32(7))=>Int32[]),owners)
    lock(api.STATE_LOCK) do
        api._replace_mesh_cache_locked!(mesh,class)
    end
end

function closest_ids(points,oracle)
    ids=Int[]
    for i in axes(points,2)
        distances=[sum((points[k,i]-oracle[k,j])^2 for k in 1:3) for j in axes(oracle,2)]
        index=argmin(distances)
        distances[index]<=9e-24 || error("node coordinates differ by $(sqrt(distances[index]))")
        push!(ids,index)
    end
    return ids
end

function native_vertices(mesh,oracle)
    ids=closest_ids(mesh.coords,oracle)
    result=Dict{Int,Vector{Tuple}}()
    for b in mesh.blocks
        result[Int(b.msh)]=sort!([Tuple(sort!([ids[n] for n in cell])) for cell in eachcol(b.nodes)])
    end
    return result
end

function single(msh)
    coords=lagrange_nodes(msh)
    return MixedMesh(coords,[ElementBlock(msh,reshape(Int32.(1:size(coords,2)),:,1))])
end

cases=Tuple{String,MixedMesh,Int,Bool}[]
for msh in (1,2,3,4,5,6,7,15),quadratic in (false,true)
    push!(cases,("family$(msh)_quadratic$(quadratic)",single(msh),1,quadratic))
end
for msh in (5,6,7)
    push!(cases,("family$(msh)_twice",single(msh),2,false))
    mesh=single(msh);coords=copy(mesh.coords)
    coords[:,end].+=[0.1,-0.1,0.15]
    push!(cases,("family$(msh)_warped",MixedMesh(coords,mesh.blocks),1,false))
end
for msh in (4,5,6,7)
    source=single(msh)
    nodes=copy(source.blocks[1].nodes)
    column=copy(nodes[:,1])
    Tessella.API._record_reverse_connectivity!(Int32(msh),column)
    nodes[:,1]=column
    push!(cases,("family$(msh)_reversed",MixedMesh(source.coords,[ElementBlock(msh,nodes)]),1,false))
end
cube=[0.0 1 1 0 0 1 1 0;0 0 1 1 0 0 1 1;0 0 0 0 1 1 1 1]
push!(cases,("hex_prism_shared_quad",MixedMesh(hcat(cube,[2.0,0,0],[2.0,0,1]),[
    ElementBlock(5,reshape(Int32.(1:8),8,1)),
    ElementBlock(6,reshape(Int32[2,9,3,6,10,7],6,1))]),1,false))
push!(cases,("pyramid_tet_shared_triangle",MixedMesh(hcat(lagrange_nodes(7),[0.0,-2,0]),[
    ElementBlock(7,reshape(Int32.(1:5),5,1)),
    ElementBlock(4,reshape(Int32[1,2,5,6],4,1))]),1,false))
push!(cases,("quad_coincident_orphan",MixedMesh(hcat(lagrange_nodes(3),zeros(3)),[
    ElementBlock(3,reshape(Int32.(1:4),4,1))]),1,true))
for (msh,primary) in ((10,4),(11,4))
    coords=copy(lagrange_nodes(msh))
    for node in primary+1:size(coords,2)
        coords[3,node]+=0.02*node
    end
    push!(cases,("curved_discrete_P2_$msh",MixedMesh(coords,[
        ElementBlock(msh,reshape(Int32.(1:size(coords,2)),:,1))]),1,false))
end

gmsh.initialize(["-nopopup"],false,false)
gmsh.option.setNumber("General.Terminal",0)
try
    startswith(gmsh.option.getString("General.Version"),"4.15.2") || error("expected Gmsh 4.15.2 runtime")
    for (name,source,levels,quadratic) in cases
        try
            install_native(source)
            api=Tessella.API
            gmsh.clear();gmsh.model.add(name)
            dim=maximum(msh_dimension(b.msh) for b in source.blocks)
            gmsh.model.addDiscreteEntity(dim,7)
            gmsh.model.mesh.addNodes(dim,7,collect(1:nnodes(source)),vec(source.coords))
            offset=0
            for b in source.blocks
                count=size(b.nodes,2)
                gmsh.model.mesh.addElementsByType(7,b.msh,collect(offset+1:offset+count),vec(b.nodes))
                offset+=count
            end
            if quadratic
                api.mesh.set_order(2);gmsh.model.mesh.setOrder(2)
            end
            result=nothing
            for _ in 1:levels
                result=api.mesh.refine();gmsh.model.mesh.refine()
            end
            gn,xyz,_=gmsh.model.mesh.getNodes()
            gc=reshape(xyz,3,:)
            nnodes(result)==length(gn) || error("$name node count mismatch: $(nnodes(result)) != $(length(gn))")
            closest_ids(gc,result.coords)
            canonical=closest_ids(gc,gc)
            ids=Dict(tag=>canonical[i] for (i,tag) in enumerate(gn))
            expected=Dict{Int,Vector{Tuple}}()
            types,_,connectivity=gmsh.model.mesh.getElements()
            for (type,nodes) in zip(types,connectivity)
                width=gmsh.model.mesh.getElementProperties(type)[4]
                expected[Int(type)]=sort!([Tuple(sort!([ids[n] for n in cell])) for cell in eachcol(reshape(nodes,width,:))])
            end
            native_vertices(result,gc)==expected || error("$name typed cell vertex sets differ")
            validate(result).ok || error("$name invalid native output")
            for b in result.blocks
                b.msh in (4,5,6,7) || continue
                all(column->Tessella.Model._extrude_cell_signed_volume(
                    result.coords,b.nodes,column,b.msh)>0,axes(b.nodes,2)) ||
                    error("$name refined volume winding is negative")
            end
            println("API_MIXED_REFINE_CASE_OK name=$name nodes=$(length(gn)) cells=$(sum(length,values(expected)))")
        finally
            Tessella.API.finalize()
        end
    end
    println("API_MIXED_REFINE_DIFFERENTIAL_OK gmsh=4.15.2 cases=$(length(cases))")
finally
    gmsh.finalize()
end
