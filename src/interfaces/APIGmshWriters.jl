# Gmsh 4.15.2's mesh writers for the foreign formats, ported line by line
# (GModel::writeSU2, writeBDF, writeUNV, writeMESH, writeINP, writeVTK,
# writeSTL and the MVertex/MElement write* methods they call), so that
# `API.write` produces the bytes Gmsh writes for the same model. They read the
# export projection `_model_export_mesh` builds, which already holds only the
# saved cells (physical selection or Mesh.SaveAll) and carries the entity
# records, node owners and the actual node and element tags.
#
# Entities are visited in GModel::getEntities order (points, curves, surfaces,
# volumes; by tag), each entity's nodes in their stored order (mesh_vertices)
# and its elements in GEntity::getMeshElement order (per element type: lines;
# triangles, quadrangles; tetrahedra, hexahedra, prisms, pyramids). The
# writers cover first-order cells; `_gmsh_writer_supported` sends anything
# else to the older writers.
using Printf: Format, format

const _GW_RANK=Dict(:pnt=>1,:lin=>1,:tri=>1,:qua=>2,:tet=>1,:hex=>2,:pri=>3,:pyr=>4)
const _GW_FIRST_ORDER=(15,1,2,3,4,5,6,7)

struct _GmshWriterView
    mesh::MixedMesh
    entities::Vector{Tuple{Int,Int}}
    nodes::Vector{Vector{Int32}}                 # rows per entity
    cells::Vector{Vector{Tuple{Int32,Int32}}}    # (block, column) per entity
    physicals::Vector{Vector{Int32}}             # entity physicals in order
    referenced::BitVector                        # node rows used by a saved cell
    node_tags::Vector{UInt64}
end

function _gmsh_writer_supported(m)
    m isa MixedMesh && m.entity_data!==nothing || return false
    return all(b->b isa ElementBlock && b.msh in _GW_FIRST_ORDER,m.blocks)
end

function _gmsh_writer_view(m::MixedMesh)
    data=m.entity_data
    declared=Set{Tuple{Int,Int}}(keys(data.entities))
    m.partition_data===nothing || union!(declared,keys(m.partition_data.entities))
    for owner in data.node_entities;push!(declared,(owner[1],Int(owner[2])));end
    for (bi,block) in enumerate(m.blocks),owner in data.block_entities[bi]
        push!(declared,(msh_dimension(block.msh),Int(owner)))
    end
    entities=sort!(collect(declared))
    index=Dict(e=>i for (i,e) in enumerate(entities))
    nodes=[Int32[] for _ in entities]
    for (row,owner) in enumerate(data.node_entities)
        push!(nodes[index[(owner[1],Int(owner[2]))]],Int32(row))
    end
    cells=[Tuple{Int32,Int32}[] for _ in entities]
    referenced=falses(size(m.coords,2))
    for (bi,block) in enumerate(m.blocks)
        dim=msh_dimension(block.msh);owners=data.block_entities[bi]
        for j in axes(block.nodes,2)
            push!(cells[index[(dim,Int(owners[j]))]],(Int32(bi),Int32(j)))
            for k in axes(block.nodes,1);referenced[block.nodes[k,j]]=true;end
        end
    end
    for list in cells
        sort!(list;by=c->_GW_RANK[msh_family(m.blocks[c[1]].msh)],alg=Base.Sort.MergeSort)
    end
    physicals=Vector{Vector{Int32}}(undef,length(entities))
    for (i,e) in enumerate(entities)
        record=get(data.entities,e,nothing)
        record===nothing && m.partition_data!==nothing &&
            (record=get(m.partition_data.entities,e,nothing))
        physicals[i]=record===nothing ? Int32[] : record.physical_tags
    end
    return _GmshWriterView(m,entities,nodes,cells,physicals,referenced,data.external_node_tags)
end

_gw_type(m::MixedMesh,c)=msh_family(m.blocks[c[1]].msh)
_gw_vertex(m::MixedMesh,c,k)=m.blocks[c[1]].nodes[k,c[2]]
_gw_nvertices(m::MixedMesh,c)=size(m.blocks[c[1]].nodes,1)
_gw_elem_tag(m::MixedMesh,c)=m.entity_data.external_element_tags[c[1]][c[2]]

# indexMeshVertices(all, 0, renumber): rows of saved cells get an index, in
# entity order; the node tag itself unless `renumber`.
function _gw_index(v::_GmshWriterView;renumber::Bool)
    index=fill(-1,length(v.node_tags));count=0
    for rows in v.nodes,row in rows
        v.referenced[row] || continue
        count+=1;index[row]=renumber ? count : Int(v.node_tags[row])
    end
    return index,count
end

function _gw_open(f,path::AbstractString,caller)
    target=abspath(path);parent=dirname(target)
    isdir(parent) || throw(ArgumentError("$caller: parent directory does not exist: $parent"))
    io=IOBuffer(;sizehint=1<<20)
    f(io)
    Base.write(target,take!(io))
    return nothing
end

const _GW_G16=Format("%.16g")

# ── allocation-free number output ───────────────────────────────────────────
const _GW_DIGITS=Vector{UInt8}(undef,32)

# Decimal digits of a non-negative integer into _GW_DIGITS, least significant
# first; returns the count.
function _gw_digits(n::UInt64)
    k=0
    while true
        k+=1;@inbounds _GW_DIGITS[k]=UInt8('0')+UInt8(n%10);n÷=10
        n==0 && break
    end
    return k
end

function _gw_int(io::IOBuffer,n::Integer)
    if n<0;Base.write(io,UInt8('-'));n=-n;end
    k=_gw_digits(UInt64(n))
    @inbounds for i in k:-1:1;Base.write(io,_GW_DIGITS[i]);end
    return nothing
end
_gw_ndigits(n::Integer)=n<0 ? 1+_gw_digits(UInt64(-n)) : _gw_digits(UInt64(n))
function _gw_lpad_int(io::IOBuffer,n::Integer,width::Int)
    for _ in 1:width-_gw_ndigits(n);Base.write(io,UInt8(' '));end
    _gw_int(io,n)
end
function _gw_rpad_int(io::IOBuffer,n::Integer,width::Int)
    _gw_int(io,n)
    for _ in 1:width-_gw_ndigits(n);Base.write(io,UInt8(' '));end
end

const _GW_EXP=Vector{UInt8}(undef,64)
const _GW_SIG=Vector{UInt8}(undef,16)
# C's printf("%.16g", x): 16 significant digits, trailing zeros removed,
# exponent style below 1e-4 and from 1e16.
function _gw_g16(io::IOBuffer,x::Float64)
    isfinite(x) || return (format(io,_GW_G16,x);nothing)
    if x==0
        signbit(x) && Base.write(io,UInt8('-'))
        Base.write(io,UInt8('0'));return nothing
    end
    n=Base.Ryu.writeexp(_GW_EXP,1,x,15)-1
    b=_GW_EXP;pos=1
    neg=b[1]==UInt8('-');neg && (pos=2)
    epos=pos+17                       # 'e'
    esign=b[epos+1]==UInt8('-') ? -1 : 1
    X=0
    for i in epos+2:n;X=10X+Int(b[i]-UInt8('0'));end
    X*=esign
    d=_GW_SIG
    d[1]=b[pos]
    for i in 1:15;d[i+1]=b[pos+1+i];end
    nd=16
    while nd>1 && d[nd]==UInt8('0');nd-=1;end
    neg && Base.write(io,UInt8('-'))
    if X< -4 || X>=16
        Base.write(io,d[1])
        if nd>1
            Base.write(io,UInt8('.'))
            for i in 2:nd;Base.write(io,d[i]);end
        end
        Base.write(io,UInt8('e'));Base.write(io,X<0 ? UInt8('-') : UInt8('+'))
        a=abs(X);a<10 && Base.write(io,UInt8('0'));_gw_int(io,a)
    elseif X>=0
        for i in 1:X+1;Base.write(io,i<=nd ? d[i] : UInt8('0'));end
        if nd>X+1
            Base.write(io,UInt8('.'))
            for i in X+2:nd;Base.write(io,d[i]);end
        end
    else
        Base.write(io,UInt8('0'));Base.write(io,UInt8('.'))
        for _ in 1:-X-1;Base.write(io,UInt8('0'));end
        for i in 1:nd;Base.write(io,d[i]);end
    end
    return nothing
end

# ── SU2 ──────────────────────────────────────────────────────────────────────
const _GW_VTK_TYPE=Dict(:pnt=>1,:lin=>3,:tri=>5,:qua=>9,:tet=>10,:hex=>12,:pri=>13,:pyr=>14)
_gw_vtk_map(t::Symbol,k::Int)=t===:pri ? (1,3,2,4,6,5)[k] : k

function _gmsh_write_su2(path,m::MixedMesh,physical_names)
    v=_gmsh_writer_view(m)
    ndime=maximum((msh_dimension(b.msh) for b in m.blocks);init=0)
    ndime in (2,3) || throw(ArgumentError("API.write: SU2 mesh output valid only for 2D or 3D models (not $(ndime)D)"))
    # Elements of the top-dimension entities, ordered by element number;
    # nodes they use, ordered by node number, indexed from 0.
    elems=Tuple{Int32,Int32}[]
    for (i,e) in enumerate(v.entities)
        e[1]==ndime || continue
        append!(elems,v.cells[i])
    end
    sort!(elems;by=c->_gw_elem_tag(m,c))
    used=Int32[]
    for c in elems,k in 1:_gw_nvertices(m,c);push!(used,_gw_vertex(m,c,k));end
    sort!(unique!(used);by=row->v.node_tags[row])
    index=zeros(Int,size(m.coords,2))
    for (i,row) in enumerate(used);index[row]=i;end
    _gw_open(path,"API.write") do io
        println(io,"NDIME= ",ndime)
        println(io,"NELEM= ",length(elems))
        for (num,c) in enumerate(elems)
            t=_gw_type(m,c);_gw_int(io,_GW_VTK_TYPE[t]);Base.write(io,UInt8(' '))
            for k in 1:_gw_nvertices(m,c)
                _gw_int(io,index[_gw_vertex(m,c,_gw_vtk_map(t,k))]-1);Base.write(io,UInt8(' '))
            end
            _gw_int(io,num-1);Base.write(io,UInt8('\n'))
        end
        println(io,"NPOIN= ",length(used))
        for (i,row) in enumerate(used)
            _gw_g16(io,m.coords[1,row]);Base.write(io,UInt8(' '));_gw_g16(io,m.coords[2,row])
            ndime==3 && (Base.write(io,UInt8(' '));_gw_g16(io,m.coords[3,row]))
            Base.write(io,UInt8(' '));_gw_int(io,i-1);Base.write(io,UInt8('\n'))
        end
        groups=_gw_physical_groups(v,ndime-1)
        if !isempty(groups)
            println(io,"NMARK= ",length(groups))
            for (tag,members) in groups
                n=sum(length(v.cells[i]) for i in members;init=0)
                n==0 && continue
                name=get(physical_names,(ndime-1,Int(tag)),"")
                isempty(name) && (name=(ndime==3 ? "PhysicalSurface" : "PhysicalLine")*string(tag))
                println(io,"MARKER_TAG= ",replace(name,' '=>'_'))
                println(io,"MARKER_ELEMS= ",n)
                for i in members,c in v.cells[i]
                    t=_gw_type(m,c);print(io,_GW_VTK_TYPE[t]," ")
                    for k in 1:_gw_nvertices(m,c)
                        print(io,index[_gw_vertex(m,c,_gw_vtk_map(t,k))]-1," ")
                    end
                    println(io)
                end
            end
        end
    end
end

# GModel::getPhysicalGroups(dim): physical tag => member entities by tag.
function _gw_physical_groups(v::_GmshWriterView,dim::Int)
    groups=Dict{Int32,Vector{Int}}()
    for (i,e) in enumerate(v.entities)
        e[1]==dim || continue
        for p in v.physicals[i]
            members=get!(()->Int[],groups,p)
            i in members || push!(members,i)
        end
    end
    return sort!(collect(groups);by=first)
end

# ── Nastran BDF (small field format, elementary tags) ───────────────────────
const _GW_BDF=Dict(:lin=>"CBAR",:tri=>"CTRIA3",:qua=>"CQUAD4",:tet=>"CTETRA",:hex=>"CHEXA",
    :pri=>"CPENTA",:pyr=>"CPYRAM")
const _GW_E2=Format("%.2E");const _GW_E1=Format("%.1E");const _GW_F=Format("%f")
function _gw_char8(x::Float64)
    s=x>=1e6 ? sprint(format,_GW_E2,x) : x>=1e-3 ? sprint(format,_GW_F,x) :
      x>=0 ? sprint(format,_GW_E2,x) : x>=-1e-3 ? sprint(format,_GW_E1,x) :
      x>=-1e6 ? sprint(format,_GW_F,x) : sprint(format,_GW_E1,x)
    return length(s)>8 ? s[1:8] : s
end
_gw_left(s,n)=rpad(s,n)

function _gmsh_write_bdf(path,m::MixedMesh)
    v=_gmsh_writer_view(m)
    index,_=_gw_index(v;renumber=false)
    _gw_open(path,"API.write") do io
        println(io,"\$ Created by Gmsh")
        for rows in v.nodes,row in rows
            index[row]<0 && continue
            print(io,"GRID    ",_gw_left(string(index[row]),8),_gw_left("0",8),
                _gw_left(_gw_char8(m.coords[1,row]),8),_gw_left(_gw_char8(m.coords[2,row]),8),
                _gw_left(_gw_char8(m.coords[3,row]),8),'\n')
        end
        cont=("E","F","G","H")
        for (i,e) in enumerate(v.entities)
            physical=isempty(v.physicals[i]) ? 0 : Int(v.physicals[i][1])
            for c in v.cells[i]
                t=_gw_type(m,c);str=get(_GW_BDF,t,nothing);str===nothing && continue
                num=_gw_elem_tag(m,c);n=_gw_nvertices(m,c)
                order=physical<0 ? reverse(1:n) : (1:n)
                print(io,_gw_left(str,8));_gw_rpad_int(io,num,8);_gw_rpad_int(io,e[2],8)
                ncont=0
                for (i0,k) in enumerate(order)
                    _gw_rpad_int(io,index[_gw_vertex(m,c,k)],8)
                    i=i0-1
                    if i!=n-1 && (i+3)%8==0
                        print(io,"+",cont[ncont+1],rpad(string(num),6),"\n+",cont[ncont+1],rpad(string(num),6))
                        ncont+=1
                    end
                end
                n==2 && print(io,"0.      0.      0.      ")
                print(io,'\n')
            end
        end
        println(io,"ENDDATA")
    end
end

# ── I-DEAS UNV (Mesh.UnvStrictFormat 1, Mesh.SaveGroupsOfElements/Nodes) ─────
const _GW_UNV=Dict(:lin=>21,:tri=>91,:qua=>94,:tet=>111,:hex=>115,:pri=>112)
const _GW_E25=Format("%25.16E")

function _gmsh_write_unv(path,m::MixedMesh,physical_names;groups_of_elements::Int=0,
                         groups_of_nodes::Int=0)
    v=_gmsh_writer_view(m)
    index,_=_gw_index(v;renumber=false)
    _gw_open(path,"API.write") do io
        print(io,lpad("-1",6),'\n',lpad("2411",6),'\n')
        for rows in v.nodes,row in rows
            index[row]<0 && continue
            _gw_lpad_int(io,index[row],10);_gw_lpad_int(io,1,10);_gw_lpad_int(io,1,10)
            _gw_lpad_int(io,11,10);Base.write(io,UInt8('\n'))
            line=sprint() do s
                format(s,_GW_E25,m.coords[1,row]);format(s,_GW_E25,m.coords[2,row])
                format(s,_GW_E25,m.coords[3,row])
            end
            print(io,replace(line,'E'=>'D'),'\n')
        end
        print(io,lpad("-1",6),'\n',lpad("-1",6),'\n',lpad("2412",6),'\n')
        for (i,e) in enumerate(v.entities),c in v.cells[i]
            t=_gw_type(m,c);type=get(_GW_UNV,t,0);type==0 && continue
            n=_gw_nvertices(m,c)
            _gw_lpad_int(io,_gw_elem_tag(m,c),10);_gw_lpad_int(io,type,10)
            _gw_lpad_int(io,e[2],10);_gw_lpad_int(io,0,10);_gw_lpad_int(io,7,10)
            _gw_lpad_int(io,n,10);Base.write(io,UInt8('\n'))
            (type==21 || type==24) && print(io,lpad("0",10),lpad("0",10),lpad("0",10),'\n')
            for k in 1:n
                _gw_lpad_int(io,index[_gw_vertex(m,c,k)],10)
                k%8==0 && Base.write(io,UInt8('\n'))
            end
            # `if(n - 1 % 8 != 7)`: always true, so always a newline.
            print(io,'\n')
        end
        print(io,lpad("-1",6),'\n')
        if groups_of_nodes!=0 || groups_of_elements!=0
            print(io,lpad("-1",6),'\n',lpad("2477",6),'\n')
            for dim in 0:3
                save_nodes=groups_of_nodes>0 || (groups_of_nodes<0 && ((-groups_of_nodes)÷10^dim)%10==1)
                save_elems=groups_of_elements>0 || (groups_of_elements<0 && ((-groups_of_elements)÷10^dim)%10==1)
                for (tag,members) in _gw_physical_groups(v,dim)
                    nodes=Int32[]
                    if save_nodes
                        for i in members,c in v.cells[i],k in 1:_gw_nvertices(m,c)
                            push!(nodes,_gw_vertex(m,c,k))
                        end
                        sort!(unique!(nodes);by=row->v.node_tags[row])
                    end
                    nele=save_elems ? sum(length(v.cells[i]) for i in members;init=0) : 0
                    print(io,lpad(string(tag),10))
                    for _ in 1:6;print(io,lpad("0",10));end
                    print(io,lpad(string(length(nodes)+nele),10),'\n')
                    name=get(physical_names,(dim,Int(tag)),"")
                    isempty(name) && (name=("PhysicalPoint","PhysicalLine","PhysicalSurface","PhysicalVolume")[dim+1]*string(tag))
                    print(io,name,'\n')
                    row=0
                    for node in nodes
                        row==2 && (print(io,'\n');row=0)
                        print(io,lpad("7",10),lpad(string(index[node]),10),lpad("0",10),lpad("0",10))
                        row+=1
                    end
                    row==2 && (print(io,'\n');row=0)
                    if save_elems
                        for i in members,c in v.cells[i]
                            row==2 && (print(io,'\n');row=0)
                            _gw_lpad_int(io,8,10);_gw_lpad_int(io,_gw_elem_tag(m,c),10)
                            _gw_lpad_int(io,0,10);_gw_lpad_int(io,0,10)
                            row+=1
                        end
                        print(io,'\n')
                    else
                        row==1 && print(io,'\n')
                    end
                end
            end
            print(io,lpad("-1",6),'\n')
        end
    end
end

# ── INRIA MEDIT .mesh (Mesh.SaveElementTagType 1: elementary) ───────────────
const _GW_G14=Format("%20.14G")
function _gmsh_write_mesh(path,m::MixedMesh;element_tag_type::Int=1)
    v=_gmsh_writer_view(m)
    index,nvertices=_gw_index(v;renumber=true)
    owner=fill(0,length(index))
    for (i,e) in enumerate(v.entities),row in v.nodes[i];owner[row]=e[2];end
    _gw_open(path,"API.write") do io
        print(io," MeshVersionFormatted 2\n Dimension\n 3\n Vertices\n ",nvertices,'\n')
        for rows in v.nodes,row in rows
            index[row]<0 && continue
            print(io,' ')
            format(io,_GW_G14,m.coords[1,row]);print(io,"      ")
            format(io,_GW_G14,m.coords[2,row]);print(io,"      ")
            format(io,_GW_G14,m.coords[3,row]);print(io,"      ",owner[row],'\n')
        end
        for (family,title) in ((:lin," Edges"),(:tri," Triangles"),(:qua," Quadrilaterals"),
                               (:tet," Tetrahedra"),(:hex," Hexahedra"))
            dim=family===:lin ? 1 : family in (:tri,:qua) ? 2 : 3
            n=0
            for (i,e) in enumerate(v.entities)
                e[1]==dim || continue
                n+=count(c->_gw_type(m,c)===family,v.cells[i])
            end
            n==0 && continue
            print(io,title,"\n ",n,'\n')
            for (i,e) in enumerate(v.entities)
                e[1]==dim || continue
                physical=isempty(v.physicals[i]) ? 0 : Int(v.physicals[i][1])
                tag=element_tag_type==2 ? abs(physical) : e[2]
                for c in v.cells[i]
                    _gw_type(m,c)===family || continue
                    nv=_gw_nvertices(m,c)
                    for k in (physical<0 ? reverse(1:nv) : (1:nv))
                        Base.write(io,UInt8(' '));_gw_int(io,index[_gw_vertex(m,c,k)])
                    end
                    Base.write(io,UInt8(' '));_gw_int(io,tag);Base.write(io,UInt8('\n'))
                end
            end
        end
        print(io," End\n")
    end
end

# ── legacy VTK (ASCII) ──────────────────────────────────────────────────────
function _gmsh_write_vtk(path,m::MixedMesh,model_name::AbstractString)
    v=_gmsh_writer_view(m)
    index,nvertices=_gw_index(v;renumber=true)
    _gw_open(path,"API.write") do io
        print(io,"# vtk DataFile Version 2.0\n",model_name,", Created by Gmsh 4.15.2 \nASCII\n",
            "DATASET UNSTRUCTURED_GRID\nPOINTS ",nvertices," double\n")
        for rows in v.nodes,row in rows
            index[row]<0 && continue
            _gw_g16(io,m.coords[1,row]);Base.write(io,UInt8(' '))
            _gw_g16(io,m.coords[2,row]);Base.write(io,UInt8(' '))
            _gw_g16(io,m.coords[3,row]);Base.write(io,UInt8('\n'))
        end
        print(io,'\n')
        ncells=0;nint=0
        for list in v.cells,c in list
            ncells+=1;nint+=_gw_nvertices(m,c)+1
        end
        print(io,"CELLS ",ncells," ",nint,'\n')
        for list in v.cells,c in list
            t=_gw_type(m,c);n=_gw_nvertices(m,c);_gw_int(io,n)
            for k in 1:n
                Base.write(io,UInt8(' '));_gw_int(io,index[_gw_vertex(m,c,_gw_vtk_map(t,k))]-1)
            end
            Base.write(io,UInt8('\n'))
        end
        print(io,'\n',"CELL_TYPES ",ncells,'\n')
        physicals=Int[];have=false
        for (i,list) in enumerate(v.cells),c in list
            _gw_int(io,_GW_VTK_TYPE[_gw_type(m,c)]);Base.write(io,UInt8('\n'))
            if isempty(v.physicals[i]);push!(physicals,-1)
            else push!(physicals,Int(v.physicals[i][1]));have=true;end
        end
        if have
            print(io,"\nCELL_DATA ",ncells,"\nSCALARS CellEntityIds int 1\nLOOKUP_TABLE default\n")
            for p in physicals;_gw_int(io,p);Base.write(io,UInt8('\n'));end
        end
    end
end

# ── STL (one solid, Mesh.StlOneSolidPerSurface 0) ───────────────────────────
const _GW_GN=Format("%g")
function _gw_normal(m::MixedMesh,c)
    a=_gw_vertex(m,c,1);b=_gw_vertex(m,c,2);d=_gw_vertex(m,c,3);x=m.coords
    t1=(x[1,b]-x[1,a],x[2,b]-x[2,a],x[3,b]-x[3,a]);t2=(x[1,d]-x[1,a],x[2,d]-x[2,a],x[3,d]-x[3,a])
    n2=t1[1]*t2[2]-t1[2]*t2[1];n1=-t1[1]*t2[3]+t1[3]*t2[1];n0=t1[2]*t2[3]-t1[3]*t2[2]
    mod=sqrt(n0*n0+n1*n1+n2*n2)
    if mod!=0.0
        inv=1.0/mod;n0*=inv;n1*=inv;n2*=inv
    end
    return (n0,n1,n2)
end

function _gmsh_write_stl(path,m::MixedMesh;binary::Bool=false)
    v=_gmsh_writer_view(m)
    faces=[i for (i,e) in enumerate(v.entities) if e[1]==2]
    facets=Tuple{Int32,Int32,Int}[]       # (block, column, sub-triangle)
    for i in faces,c in v.cells[i]
        t=_gw_type(m,c)
        t===:tri && push!(facets,(c[1],c[2],1))
        t===:qua && (push!(facets,(c[1],c[2],1));push!(facets,(c[1],c[2],2)))
    end
    _gw_open(path,"API.write") do io
        name="Created by Gmsh"
        if binary
            header=zeros(UInt8,80);copyto!(header,codeunits(name));Base.write(io,header)
            Base.write(io,UInt32(length(facets)))
        else
            print(io,"solid ",name,'\n')
        end
        for (bi,j,sub) in facets
            c=(bi,j);n=_gw_normal(m,c)
            corners=sub==1 ? (1,2,3) : (1,3,4)
            if binary
                Base.write(io,Float32(n[1]),Float32(n[2]),Float32(n[3]))
                for k in corners,a in 1:3;Base.write(io,Float32(m.coords[a,_gw_vertex(m,c,k)]));end
                Base.write(io,UInt16(0))
            else
                print(io,"facet normal ");format(io,_GW_GN,n[1]);print(io,' ')
                format(io,_GW_GN,n[2]);print(io,' ');format(io,_GW_GN,n[3])
                print(io,"\n  outer loop\n")
                for k in corners
                    row=_gw_vertex(m,c,k);print(io,"    vertex ")
                    _gw_g16(io,m.coords[1,row]);Base.write(io,UInt8(' '))
                    _gw_g16(io,m.coords[2,row]);Base.write(io,UInt8(' '))
                    _gw_g16(io,m.coords[3,row]);Base.write(io,UInt8('\n'))
                end
                print(io,"  endloop\nendfacet\n")
            end
        end
        binary || print(io,"endsolid ",name,'\n')
    end
end

# ── Abaqus INP ──────────────────────────────────────────────────────────────
const _GW_INP=Dict(:lin=>"T3D2",:tri=>"CPS3",:qua=>"CPS4",:tet=>"C3D4",:hex=>"C3D8",
    :pri=>"C3D6",:pyr=>"C3D5")
const _GW_G14S=Format("%.14g")
function _gmsh_write_inp(path,m::MixedMesh,physical_names;groups_of_elements::Int=1,
                         groups_of_nodes::Int=0)
    v=_gmsh_writer_view(m)
    index,_=_gw_index(v;renumber=false)
    _gw_open(path,"API.write") do io
        print(io,"*Heading\n ",path,"\n*NODE\n")
        for rows in v.nodes,row in rows
            index[row]<0 && continue
            print(io,index[row],", ")
            format(io,_GW_G14S,m.coords[1,row]);print(io,", ")
            format(io,_GW_G14S,m.coords[2,row]);print(io,", ")
            format(io,_GW_G14S,m.coords[3,row]);print(io,'\n')
        end
        print(io,"******* E L E M E N T S *************\n")
        for dim in 0:3
            (groups_of_elements>=0 || ((-groups_of_elements)÷10^dim)%10!=0) || continue
            for (i,e) in enumerate(v.entities)
                e[1]==dim || continue
                for family in (:pnt,:lin,:tri,:qua,:tet,:hex,:pri,:pyr)
                    cells=[c for c in v.cells[i] if _gw_type(m,c)===family]
                    isempty(cells) && continue
                    typ=get(_GW_INP,family,nothing);typ===nothing && continue
                    str=("Point","Line","Surface","Volume")[dim+1]
                    print(io,"*ELEMENT, type=",typ,", ELSET=",str,e[2],'\n')
                    for c in cells
                        _gw_int(io,_gw_elem_tag(m,c));print(io,", ");n=_gw_nvertices(m,c)
                        for k in 1:n
                            _gw_int(io,index[_gw_vertex(m,c,k)])
                            if k!=n
                                print(io,", ")
                                (k-1)!=0 && (k+1)%16==0 && print(io,'\n')
                            end
                        end
                        print(io,'\n')
                    end
                end
            end
        end
        name(dim,tag)=begin
            s=get(physical_names,(dim,Int(tag)),"")
            isempty(s) && (s=("PhysicalPoint","PhysicalLine","PhysicalSurface","PhysicalVolume")[dim+1]*string(tag))
            replace(s,' '=>'_')
        end
        if groups_of_elements!=0
            for dim in 1:3
                (groups_of_elements>0 || (groups_of_elements<0 && ((-groups_of_elements)÷10^dim)%10==1)) || continue
                for (tag,members) in _gw_physical_groups(v,dim)
                    print(io,"*ELSET,ELSET=",name(dim,tag),'\n');n=0
                    for i in members,c in v.cells[i]
                        n>0 && n%10==0 && print(io,'\n')
                        print(io,_gw_elem_tag(m,c),", ");n+=1
                    end
                    print(io,'\n')
                end
            end
        end
        if groups_of_nodes!=0
            for dim in 0:3
                (groups_of_nodes>0 || (groups_of_nodes<0 && ((-groups_of_nodes)÷10^dim)%10==1)) || continue
                for (tag,members) in _gw_physical_groups(v,dim)
                    nodes=Int32[]
                    for i in members,c in v.cells[i],k in 1:_gw_nvertices(m,c);push!(nodes,_gw_vertex(m,c,k));end
                    sort!(unique!(nodes);by=row->v.node_tags[row])
                    print(io,"*NSET,NSET=",name(dim,tag),'\n');n=0
                    for row in nodes
                        n>0 && n%10==0 && print(io,'\n')
                        print(io,index[row],", ");n+=1
                    end
                    print(io,'\n')
                end
            end
        end
    end
end
