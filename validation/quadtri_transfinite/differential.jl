#!/usr/bin/env julia
# Pinned Gmsh 4.15.2 oracle: all boundary masks, actual face arrangements,
# node multisets, and per-type/per-entity cell vertex sets. Cell corner
# rotations and global node numbering are intentionally canonicalized.
using Tessella
using Tessella.Elements: ElementBlock
using Tessella.MeshTypes: nnodes
function find_gmsh_api()
    configured=get(ENV,"GMSH_JULIA_API","")
    if !isempty(configured)
        isfile(configured) || error("GMSH_JULIA_API does not name a file")
        return configured
    end
    candidates=["/opt/homebrew/lib/gmsh.jl","/usr/local/lib/gmsh.jl"]
    if Sys.iswindows()
        push!(candidates,joinpath(homedir(),"AppData","Local","Programs",
                                 "Python","Python312","Lib","gmsh.jl"))
    else
        executable=Sys.which("gmsh")
        if executable!==nothing
            prefix=dirname(dirname(realpath(executable)))
            append!(candidates,(joinpath(prefix,"lib","gmsh.jl"),
                                joinpath(prefix,"lib64","gmsh.jl")))
        end
    end
    found=findfirst(isfile,candidates)
    found===nothing && error("set GMSH_JULIA_API to the pinned Gmsh 4.15.2 API")
    return candidates[found]
end
include(find_gmsh_api())
gmsh.GMSH_API_VERSION=="4.15.2" || error("pinned Gmsh 4.15.2 required")

const BOX=read(joinpath(@__DIR__,"..","..","test","artifacts","quadtri_box.geo"),String)
const PRISM=read(joinpath(@__DIR__,"..","..","test","artifacts","quadtri_prism.geo"),String)

function compare(source,name)
    mktemp() do path,io
        write(io,source*"Mesh 3;\n");close(io)
        ex=execute_geo(path)
        gmsh.clear();gmsh.parser.clear()
        gmsh.merge(path)
        gtags,gcoords,_=gmsh.model.mesh.getNodes()
        nnodes(ex.mesh)==length(gtags) || error("$name: node count differs")
        native_map=Dict(Tuple(ex.mesh.coords[:,j])=>j for j in 1:nnodes(ex.mesh))
        owners=fill((4,0),nnodes(ex.mesh));used=falses(nnodes(ex.mesh))
        for (dim,tag,part) in ex.mesh_parts
            remap=[native_map[Tuple(part.coords[:,j])] for j in 1:nnodes(part)]
            for j in remap
                (dim,tag)<owners[j] && (owners[j]=(dim,tag))
            end
            dim==3 || continue
            if part isa Tessella.Elements.MixedMesh
                for b in part.blocks
                    b isa ElementBlock || continue
                    for j in b.nodes;used[remap[j]]=true;end
                end
            else
                for j in part.tets;used[remap[j]]=true;end
            end
        end
        gowners=Dict{UInt64,Tuple{Int,Int}}()
        for (dim,tag) in gmsh.model.getEntities()
            tags,_,_=gmsh.model.mesh.getNodes(dim,tag,false,false)
            for j in tags;gowners[j]=(Int(dim),Int(tag));end
        end
        _,_,volumes=gmsh.model.mesh.getElements(3)
        gused=Set(UInt64[j for list in volumes for j in list])
        # Match coordinates one-to-one with entity ownership and volume-use
        # identity. Compact prisms intentionally store unreferenced volume
        # nodes within floating-point noise of referenced surface/tab nodes;
        # arbitrary nearest-coordinate matches mix those identities.
        free=trues(nnodes(ex.mesh));mapping=Dict{UInt64,Int}()
        for (i,tag) in enumerate(gtags)
            p=(gcoords[3i-2],gcoords[3i-1],gcoords[3i])
            matches=Int[j for j in 1:nnodes(ex.mesh) if
                owners[j]==gowners[tag] && used[j]==(tag in gused) &&
                all(d->abs(ex.mesh.coords[d,j]-p[d])<1e-7,1:3)]
            isempty(matches) && error("$name: unmatched oracle node $p")
            slot=findfirst(j->free[j],matches)
            slot===nothing && error("$name: coordinate multiplicity differs")
            free[matches[slot]]=false
            mapping[tag]=matches[slot]
        end
        native=Dict{Tuple{Int,Int,Int},Vector{Vector{Int}}}()
        for (dim,tag) in gmsh.model.getEntities()
            dim==0 && continue
            part=geo_entity_mesh(ex,Int(dim),Int(tag))
            remap=[native_map[Tuple(part.coords[:,j])] for j in 1:nnodes(part)]
            if part isa Tessella.Elements.MixedMesh
                for b in part.blocks
                    b isa ElementBlock || continue
                    native[(Int(dim),Int(tag),b.msh)]=sort!(
                        [sort!(Int[remap[n] for n in b.nodes[:,i]]) for i in axes(b.nodes,2)])
                end
            else
                for (msh,cells) in ((1,part.segs),(2,part.tris),(4,part.tets))
                    size(cells,2)>0 || continue
                    native[(Int(dim),Int(tag),msh)]=sort!(
                        [sort!(Int[remap[n] for n in cells[:,i]]) for i in axes(cells,2)])
                end
            end
            types,_,cells=gmsh.model.mesh.getElements(dim,tag)
            for (msh,list) in zip(types,cells)
                _,_,_,width,_,_=gmsh.model.mesh.getElementProperties(Int(msh))
                key=(Int(dim),Int(tag),Int(msh))
                expected=sort!([sort!(Int[mapping[list[i+k]] for k in 1:Int(width)])
                    for i in 0:Int(width):length(list)-1])
                get(native,key,nothing)==expected || error("$name: cells differ for $key")
                delete!(native,key)
            end
        end
        isempty(native) || error("$name: native has unmatched cell blocks")
        validate(ex.mesh).ok || error("$name: merged native mesh invalid")
    end
end

function source(geo,mask;n=3,arrangement="Left",extra="")
    rc=isempty(mask) ? "" : "Recombine Surface{"*join(mask,",")*"};"
    return geo*"Transfinite Curve{:}=$(n+1);Transfinite Surface{:} $arrangement;"*
        "Transfinite Volume{1};TransfQuadTri{1};$rc$extra"
end

gmsh.initialize(String[],false,false)
cases=Ref(0)
try
    gmsh.option.setNumber("General.Terminal",0)
    gmsh.option.setNumber("General.NumThreads",1)
    for (name,geo,nface) in (("box",BOX,6),("prism",PRISM,5))
        for bits in 0:(1<<nface)-1
            mask=[s for s in 1:nface if !iszero(bits & (1<<(s-1)))]
            compare(source(geo,mask),"$name-mask-$bits");cases[]+=1
        end
        for arrangement in ("Right","AlternateLeft","AlternateRight"),n in (1,2,4)
            compare(source(geo,[1,3];n,arrangement),"$name-$arrangement-$n");cases[]+=1
        end
    end
    # Nonuniform edge laws and a warped ruled box exercise actual face
    # diagonals after surface permutations and Coons interpolation.
    graded=source(BOX,[1,3,5];extra="Transfinite Curve{1,3,5,7}=4 Using Progression 1.5;")
    compare(graded,"graded-box");cases[]+=1
    warped=replace(BOX,"{1,1,1,0.4}"=>"{1,1,1.2,0.4}","Plane Surface"=>"Surface")
    compare(source(warped,[1,3,5]),"warped-box");cases[]+=1
    compare(source(PRISM,collect(1:5);extra="Mesh.TransfiniteTri=1;"),"compact-prism");cases[]+=1
    println("TRANSFINITE_QUADTRI_DIFFERENTIAL_OK gmsh=4.15.2 cases=$(cases[])")
finally
    gmsh.isInitialized()!=0 && gmsh.finalize()
end
