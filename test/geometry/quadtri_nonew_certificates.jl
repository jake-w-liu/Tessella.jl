module QuadTriNoNewCertificates

using Tessella
using Tessella.Elements: ElementBlock, MixedMesh, _VOLUME_CELL_FACES, msh_dimension,
                        msh_spec, lagrange_nodes
using Tessella.MeshTypes: nnodes, tet_signed_volume
using Tessella.Predicates: orient3
using LinearAlgebra: norm

const CORNERS=((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.))
const SQUARE="""
Point(1)={0,0,0,1};Point(2)={1,0,0,1};
Point(3)={1,1,0,1};Point(4)={0,1,0,1};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
Transfinite Curve{:}=2;Transfinite Surface{1};Recombine Surface{1};
"""

# Outward MSH corner cycles, independently of the extrusion factories.
const FACES=Dict(
    4=>((1,3,2),(1,2,4),(2,3,4),(3,1,4)),
    5=>((1,4,3,2),(5,6,7,8),(1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8)),
    6=>((1,3,2),(4,5,6),(1,2,5,4),(2,3,6,5),(3,1,4,6)),
    7=>((1,4,3,2),(1,2,5),(2,3,5),(3,4,5),(4,1,5)))

require(ok,message)=ok || throw(ArgumentError("NoNew certificate: "*message))
key(nodes)=Tuple(sort!(collect(nodes)))
point(mesh,index)=Tuple(mesh.coords[:,index])
blocks(mesh::MixedMesh)=mesh.blocks
blocks(mesh::Mesh)=[ElementBlock(4,mesh.tets)]
counts(mesh)=Dict(b.msh=>size(b.nodes,2) for b in blocks(mesh)
                 if b isa ElementBlock && msh_dimension(b.msh)==3)

# Reference P1 interpolation is independent of the cache elevation routine.
# It also gives topological support keys for shared P2 nodes.
function primary_weights(msh,u,v,w)
    msh==4 && return (1-u-v-w,u,v,w)
    if msh==5
        signs=((-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
               (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1))
        return Tuple((1+s[1]*u)*(1+s[2]*v)*(1+s[3]*w)/8 for s in signs)
    elseif msh==6
        b=(1-u-v,u,v)
        return (b[1]*(1-w)/2,b[2]*(1-w)/2,b[3]*(1-w)/2,
                b[1]*(1+w)/2,b[2]*(1+w)/2,b[3]*(1+w)/2)
    elseif msh==7
        w==1 && return (0.,0.,0.,0.,1.)
        a,b=u/(1-w),v/(1-w)
        return ((1-w)*(1-a)*(1-b)/4,(1-w)*(1+a)*(1-b)/4,
                (1-w)*(1+a)*(1+b)/4,(1-w)*(1-a)*(1+b)/4,w)
    end
    error("unsupported primary interpolation family $msh")
end

function certify_quadratic(mesh,linear)
    original=Dict{Tuple,Tuple}()
    for block in blocks(linear),cell in eachcol(block.nodes)
        original[(block.msh,Tuple(point(linear,n) for n in cell))]=Tuple(cell)
    end
    back=Dict(11=>4,12=>5,13=>6,14=>7)
    primary=Dict(4=>4,5=>8,6=>6,7=>5)
    support_nodes=Dict{Tuple,Int}();node_supports=Dict{Int,Tuple}()
    seen=Set{Tuple}();used=Set{Int}();ncell=0
    for block in blocks(mesh)
        require(haskey(back,block.msh),"unexpected quadratic volume family $(block.msh)")
        msh=back[block.msh];nprimary=primary[msh]
        reference=lagrange_nodes(block.msh)
        require(size(block.nodes,1)==size(reference,2),"quadratic reference width")
        require(msh_spec(block.msh).order==2,"quadratic order")
        for cell in eachcol(block.nodes)
            ncell+=1;union!(used,cell)
            identity=(msh,Tuple(point(mesh,cell[k]) for k in 1:nprimary))
            require(haskey(original,identity),"quadratic primary cell changed")
            require(!(identity in seen),"quadratic primary cell duplicated")
            push!(seen,identity);source=original[identity]
            for j in eachindex(cell)
                u,v,w=reference[:,j]
                weights=primary_weights(msh,u,v,w)
                expected=ntuple(k->sum(weights[i]*linear.coords[k,source[i]]
                                      for i in 1:nprimary),3)
                require(maximum(abs.(expected.-point(mesh,cell[j])))<=2e-11,
                        "quadratic support placement")
                support=Tuple(sort!([(Int(source[i]),weights[i]) for i in 1:nprimary
                                     if weights[i]!=0]))
                require(get(support_nodes,support,Int(cell[j]))==cell[j],
                        "shared quadratic support split")
                require(get(node_supports,Int(cell[j]),support)==support,
                        "unrelated quadratic support welded")
                support_nodes[support]=Int(cell[j]);node_supports[Int(cell[j])]=support
            end
        end
    end
    require(length(seen)==length(original),"quadratic cell set changed")
    require(used==Set(1:nnodes(mesh)),"unreferenced quadratic nodes")
    require(length(support_nodes)==nnodes(mesh),"quadratic support count")
    return (;ncell,nsupport=length(support_nodes))
end

function execute(source;dim=3)
    mktempdir() do directory
        path=joinpath(directory,"nonew.geo")
        write(path,source)
        execute_geo(path;mesh_dim=dim)
    end
end

function fixture(name,layers,levels,laterals;motion=:translation,
                 delta=(0.,0.,1.),angle=pi/6)
    transform=motion===:translation ? "Extrude{$(join(delta,','))}" :
        motion===:rotation ? "Extrude{{0,1,0},{-1,0,0},Pi/6}" :
                            "Extrude{{0,0,2},{0,0,1},{0,0,0},Pi/6}"
    source=SQUARE*"sweep[]=$transform{Surface{1};$layers;Recombine;"*
        "QuadTriNoNewVerts$(laterals ? " RecombLaterals" : "");};\n"
    return (;name,source,levels=Float64.(levels),laterals,motion,delta,angle)
end

function fixtures()
    result=NamedTuple[]
    groups=(
        ("uniform1","Layers{1}",[0.,1.]),
        ("uniform3","Layers{3}",[0.,1/3,2/3,1.]),
        ("first_single","Layers{{1,2},{0.5,1.0}}",[0.,.5,.75,1.]),
        ("last_single","Layers{{2,1},{0.5,1.0}}",[0.,.25,.5,1.]),
        ("five_layers","Layers{{2,3},{0.25,1.0}}",[0.,.125,.25,.5,.75,1.]),
        ("regrouped","Layers{{1,1,1},{0.5,0.75,1.0}}",[0.,.5,.75,1.]))
    for (name,layers,levels) in groups,laterals in (false,true)
        push!(result,fixture("$(name)_$(laterals)",layers,levels,laterals))
    end
    for (name,motion,delta) in (("tilted",:translation,(.25,-.125,2.)),
                               ("negative",:translation,(0.,0.,-1.)),
                               ("rotation",:rotation,(0.,0.,0.)),
                               ("twist",:twist,(0.,0.,2.))),laterals in (false,true)
        push!(result,fixture("$(name)_$(laterals)","Layers{{1,2},{0.5,1.0}}",
            [0.,.5,.75,1.],laterals;motion,delta))
    end
    return result
end

function paired_source()
    return SQUARE*"""
    Point(101)={3,0,0,1};Point(102)={4,0,0,1};
    Point(103)={4,1,0,1};Point(104)={3,1,0,1};
    Line(101)={101,102};Line(102)={102,103};
    Line(103)={103,104};Line(104)={104,101};
    Curve Loop(101)={101,102,103,104};Plane Surface(101)={101};
    Transfinite Curve{101,102,103,104}=2;
    Transfinite Surface{101};Recombine Surface{101};
    left[]=Extrude{0,0,1}{Surface{1};Layers{3};Recombine;QuadTriNoNewVerts;};
    right[]=Extrude{0,0,1}{Surface{101};Layers{3};Recombine;QuadTriNoNewVerts RecombLaterals;};
    """
end

function helical_fixture(name,pitch,laterals=true)
    levels=collect(range(0.,1.;length=53))
    source="Geometry.AutoCoherence=0;\n"*SQUARE*
        "sweep[]=Extrude{{0,$pitch,0},{0,1,0},{-1,0,0},13*Pi/6}{Surface{1};"*
        "Layers{52};Recombine;QuadTriNoNewVerts$(laterals ? " RecombLaterals" : "");};\n"
    return (;name,source,levels,laterals,motion=:helical,delta=(0.,Float64(pitch),0.),
            angle=13pi/6)
end

const FOLDED_HEX_VERTICES=(
    (1.,1.,0.),(1.,0.,0.),(0.,0.,0.),(0.,1.,0.),
    (5.5224396333287125,-4.486401190603047,-.30524078967639356),
    (5.192952385269027,-5.156913942543361,.3594777995412404),
    (4.522439633328713,-5.486401190603047,-.30524078967639356),
    (4.8519268813883984,-4.815888438662733,-.9699593788940275))

function folded_hex_source()
    return "Geometry.AutoCoherence=0;\n"*SQUARE*"""
    Extrude{{5.608013059266764,-7.535936173815433,-4.3120799162487025},
            {1,1,0},{1.940745333694334,2.08924428794778,2.658811780113151},
            -2.445577491942423}{
        Surface{1};Layers{2};Recombine;QuadTriNoNewVerts RecombLaterals;}
    """
end

function exact_hex_jacobian(vertices,r,s,t)
    bits=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
          (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    q=(r,s,t)
    exact=Tuple(Tuple(Rational{BigInt}(p[k]) for k in 1:3) for p in vertices)
    derivative=[sum((bits[i][axis]==0 ? -1 : 1)*
        prod(bits[i][j]==0 ? 1-q[j] : q[j] for j in 1:3 if j!=axis)*exact[i][k]
        for i in 1:8) for k in 1:3,axis in 1:3]
    a=derivative
    return a[1,1]*(a[2,2]*a[3,3]-a[2,3]*a[3,2])-
           a[1,2]*(a[2,1]*a[3,3]-a[2,3]*a[3,1])+
           a[1,3]*(a[2,1]*a[3,2]-a[2,2]*a[3,1])
end

# Invert the literal Gmsh trilinear Hex8 map, independently of the planner's
# intersection algorithm. A strict reference-interior witness proves overlap.
function trilinear_parameters(vertices,p)
    signs=((-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
           (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1))
    q=zeros(3)
    for _ in 1:20
        weights=primary_weights(5,q...)
        residual=[sum(weights[i]*vertices[i][k] for i in 1:8)-p[k] for k in 1:3]
        norm(residual,Inf)<=2e-13 && return q
        jacobian=[sum(signs[i][axis]/8 *
            prod(1+signs[i][j]*q[j] for j in 1:3 if j!=axis)*vertices[i][k]
            for i in 1:8) for k in 1:3,axis in 1:3]
        q-=jacobian\residual
    end
    error("trilinear witness inversion did not converge")
end

function crc_record(name,mesh)
    digest=mixed_crc(mesh)
    return (;name,n_nodes=digest.n_nodes,n_blocks=digest.n_blocks,
            n_cells=digest.n_cells,sha=digest.sha)
end

# Reproducible native artifacts; Gmsh's pointer-selected cells are certified
# independently and never used as a golden connectivity digest.
function artifact_records()
    records=NamedTuple[]
    for f in fixtures()
        push!(records,crc_record(f.name,execute(f.source).mesh))
    end
    for f in filter(f->startswith(f.name,"uniform"),fixtures())
        api=Tessella.API
        api.initialize()
        try
            mktempdir() do directory
                path=joinpath(directory,"quadratic.geo");write(path,f.source)
                api.open_geo!(path;mesh_dim=0)
                api.mesh.generate(3);api.mesh.set_order(2)
                push!(records,crc_record(f.name*"_api_p2",api.mesh.get()))
            end
        finally
            api.finalize()
        end
    end
    return records
end

function artifact_pins()
    path=joinpath(@__DIR__,"..","artifacts","quadtri_nonew_crc.txt")
    records=Dict{String,NamedTuple}()
    for line in eachline(path)
        isempty(line) || startswith(line,"#") || begin
            fields=split(line,'\t')
            require(length(fields)==5,"artifact row width")
            name=String(fields[1])
            records[name]=(;name,n_nodes=parse(Int,fields[2]),n_blocks=parse(Int,fields[3]),
                            n_cells=parse(Int,fields[4]),sha=String(fields[5]))
        end
    end
    return records
end

function transformed(f,p,u)
    f.motion===:translation && return ntuple(k->p[k]+u*f.delta[k],3)
    c,s=cos(f.angle*u),sin(f.angle*u)
    f.motion in (:rotation,:helical) && return (c*(p[1]+1)+s*p[3]-1,
        p[2]+(f.motion===:helical ? u*f.delta[2] : 0.),-s*(p[1]+1)+c*p[3])
    return (c*p[1]-s*p[2],s*p[1]+c*p[2],p[3]+2u)
end

function canonical_cycle(nodes)
    return minimum(Tuple(circshift(collect(nodes),i)) for i in 0:length(nodes)-1)
end

# Surface flux uses the actual bilinear quad map. Choosing a hidden diagonal
# on a warped quadrangle would change the represented P1 volume.
function face_volume(mesh,face,origin)
    p=[ntuple(k->mesh.coords[k,node]-origin[k],3) for node in face]
    if length(face)==3
        return tet_signed_volume((0.,0.,0.),p...)
    end
    volume=0.
    for x in (-inv(sqrt(3.)),inv(sqrt(3.))),y in (-inv(sqrt(3.)),inv(sqrt(3.)))
        basis=((1-x)*(1-y),(1+x)*(1-y),(1+x)*(1+y),(1-x)*(1+y))
        dx=(-(1-y),1-y,1+y,-(1+y))
        dy=(-(1-x),-(1+x),1+x,1-x)
        q=ntuple(k->sum(basis[i]*p[i][k] for i in 1:4)/4,3)
        a=ntuple(k->sum(dx[i]*p[i][k] for i in 1:4)/4,3)
        b=ntuple(k->sum(dy[i]*p[i][k] for i in 1:4)/4,3)
        volume+=(q[1]*(a[2]*b[3]-a[3]*b[2])+
                 q[2]*(a[3]*b[1]-a[1]*b[3])+
                 q[3]*(a[1]*b[2]-a[2]*b[1]))/3
    end
    return volume
end

function certify_complex(mesh)
    require(all(isfinite,mesh.coords),"nonfinite coordinates")
    faces=Dict{Tuple,Vector{Tuple}}()
    edges=Set{Tuple}();used=Set{Int}();ncell=0;total=0.;minimum_volume=Inf
    for block in blocks(mesh)
        block isa ElementBlock && msh_dimension(block.msh)==3 || continue
        require(haskey(FACES,block.msh),"unexpected volume family $(block.msh)")
        require(Set(key(f) for f in FACES[block.msh])==
                Set(key(f) for f in _VOLUME_CELL_FACES[block.msh]),"MSH face catalog differs")
        for cell in eachcol(block.nodes)
            ncell+=1;union!(used,cell)
            require(length(unique(cell))==length(cell),"repeated cell corner")
            center=ntuple(k->sum(mesh.coords[k,n] for n in cell)/length(cell),3)
            volume=0.
            for local_face in FACES[block.msh]
                face=Tuple(Int(cell[i]) for i in local_face)
                push!(get!(faces,key(face),Tuple[]),face)
                for i in eachindex(face)
                    push!(edges,key((face[i],face[mod1(i+1,length(face))])))
                end
                triangles=length(face)==3 ? (face,) :
                    ((face[1],face[2],face[3]),(face[1],face[3],face[4]),
                     (face[2],face[3],face[4]),(face[2],face[4],face[1]))
                for triangle in triangles
                    ps=map(n->point(mesh,n),triangle)
                    require(orient3(center,ps...)==-1,"nonpositive exact face fan")
                end
                volume+=face_volume(mesh,face,center)
            end
            require(isfinite(volume) && volume>0,"nonpositive cell volume")
            minimum_volume=min(minimum_volume,volume);total+=volume
        end
    end
    require(ncell>0,"empty region")
    boundary=Dict{Tuple,Tuple}()
    for (face,cycles) in faces
        require(length(cycles) in (1,2),"face incidence exceeds two")
        if length(cycles)==2
            require(canonical_cycle(cycles[1])==canonical_cycle(reverse(cycles[2])),
                    "interior face orientations disagree")
        else
            boundary[face]=only(cycles)
        end
    end
    boundary_edges=Dict{Tuple,Int}()
    for face in values(boundary),i in eachindex(face)
        edge=key((face[i],face[mod1(i+1,length(face))]))
        boundary_edges[edge]=get(boundary_edges,edge,0)+1
    end
    require(all(==(2),values(boundary_edges)),"open or nonmanifold boundary")
    require(length(used)-length(edges)+length(faces)-ncell==1,"ball Euler characteristic")
    origin=ntuple(k->sum(mesh.coords[k,n] for n in used)/length(used),3)
    boundary_volume=sum(face_volume(mesh,face,origin) for face in values(boundary))
    require(isapprox(total,boundary_volume;atol=2e-12,rtol=2e-12),"boundary volume differs")
    return (;boundary,total,minimum_volume,ncell)
end

function certify_fixture(mesh,f)
    cert=certify_complex(mesh)
    grid=[transformed(f,p,u) for u in f.levels,p in CORNERS]
    grid_ids=Matrix{Int}(undef,size(grid))
    for i in eachindex(grid)
        matches=[j for j in 1:nnodes(mesh) if maximum(abs.(point(mesh,j).-grid[i]))<2e-11]
        require(length(matches)==1,"swept corner missing or duplicated")
        grid_ids[i]=only(matches)
    end
    require(length(unique(grid_ids))==length(grid_ids),"collapsed sweep")
    extras=setdiff(collect(1:nnodes(mesh)),vec(grid_ids))
    nlayer=length(f.levels)-1
    require(length(extras)==(f.laterals ? 1 : 0),"unexpected extra vertices")
    remaining=Set(keys(cert.boundary))
    bottom=key(grid_ids[1,:]);require(bottom in remaining,"source quad changed")
    delete!(remaining,bottom)
    top=Set(grid_ids[end,:])
    top_faces=[face for face in remaining if issubset(Set(face),top)]
    require(length(top_faces)==2 && all(face->length(face)==3,top_faces),"invalid top split")
    diagonal=key(intersect(collect(top_faces[1]),collect(top_faces[2])))
    require(diagonal in (key(grid_ids[end,[1,3]]),key(grid_ids[end,[2,4]])),"top diagonal is not opposite")
    setdiff!(remaining,top_faces)
    for layer in 1:nlayer,side in 1:4
        next=mod1(side+1,4)
        corners=Set((grid_ids[layer,side],grid_ids[layer,next],
                     grid_ids[layer+1,side],grid_ids[layer+1,next]))
        local_faces=[face for face in remaining if issubset(Set(face),corners)]
        if f.laterals
            require(length(local_faces)==1 && length(only(local_faces))==4,"recombined lateral changed")
        else
            require(length(local_faces)==2 && all(face->length(face)==3,local_faces),"free lateral missing split")
            common=key(intersect(collect(local_faces[1]),collect(local_faces[2])))
            require(common in (key((grid_ids[layer,side],grid_ids[layer+1,next])),
                               key((grid_ids[layer,next],grid_ids[layer+1,side]))),"invalid lateral diagonal")
        end
        setdiff!(remaining,local_faces)
    end
    require(isempty(remaining),"unexpected exterior face")
    family_counts=counts(mesh)
    if f.laterals
        center=ntuple(k->sum(grid[i,j][k] for i in size(grid,1)-1:size(grid,1),j in 1:4)/8,3)
        require(maximum(abs.(point(mesh,only(extras)).-center))<2e-11,"wrong final centroid")
        require(family_counts==filter(p->last(p)>0,Dict(4=>2,5=>nlayer-1,7=>5)),"wrong RecombLaterals families")
        last_corners=Set(grid_ids[end-1:end,:])
        push!(last_corners,only(extras))
        for block in blocks(mesh),cell in eachcol(block.nodes)
            block.msh==5 && continue
            require(only(extras) in cell && issubset(Set(cell),last_corners),"centroid fan outside final interval")
        end
    elseif f.name=="uniform1_false"
        require(family_counts==Dict(4=>4,7=>1),"wrong one-layer free families")
    elseif f.name=="uniform3_false"
        require(family_counts in (Dict(4=>16,7=>1),Dict(4=>12,7=>3),Dict(4=>8,7=>5)),
                "uniform3 free families outside repeated oracle outcomes")
    end
    if f.motion===:translation
        require(isapprox(cert.total,abs(f.delta[3])*f.levels[end];atol=2e-12,rtol=2e-12),"analytic translation volume differs")
    end
    return merge(cert,(;grid_ids,family_counts,extras))
end

function surface_faces(parts,volume;tags=nothing)
    lookup=Dict(point(volume,i)=>i for i in 1:nnodes(volume))
    result=Set{Tuple}()
    for (dim,tag,part) in parts
        dim==2 && (tags===nothing || tag in tags) || continue
        part_blocks=part isa Mesh ? [ElementBlock(2,part.tris)] : part.blocks
        for b in part_blocks
            b isa ElementBlock && b.msh in (2,3) || continue
            for cell in eachcol(b.nodes)
                push!(result,key(lookup[point(part,n)] for n in cell))
            end
        end
    end
    return result
end

function typed_signature(mesh)
    return sort!([(b.msh,key(point(mesh,n) for n in cell)) for b in blocks(mesh)
                  if b isa ElementBlock && msh_dimension(b.msh)==3 for cell in eachcol(b.nodes)])
end

end
