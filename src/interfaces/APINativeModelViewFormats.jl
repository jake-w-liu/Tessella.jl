# Model CSV is a rectangular table of owned node/cell/data records. Model VTK
# retains full supported cells and explicit sparse values; POS/TXT keep their
# existing visualization contract. All payloads are assembled before opening.
function _native_model_view_export_plan(view,plan,caller;vtk=false)
    mesh=plan.mesh;data=mesh.entity_data
    cells=_native_view_text_cells(mesh;primary_only=false)
    nodes=Dict(uid=>row for (row,uid) in enumerate(data.external_node_tags))
    cell_rows=Dict(cell.uid=>cell for cell in cells)
    tables=[_native_view_text_table(view,key) for key in plan.group.keys]
    for table in tables
        bound=table.kind=="NodeData" ? nodes : cell_rows
        for (uid,values) in table.rows
            haskey(bound,uid) || throw(ArgumentError(
                "$caller: model view contains unknown $(table.kind=="NodeData" ? "node" : "element") UID $uid"))
            if vtk && table.kind=="ElementNodeData"
                length(values)==table.components*length(cell_rows[uid].nodes) || throw(ArgumentError(
                    "$caller: VTK ElementNodeData must provide every node of its complete element"))
            end
        end
    end
    ranges=[range for range in view.model_step_ranges if range.model_slot===plan.group.slot]
    return (;mesh,data,cells,nodes,cell_rows,tables,ranges,slot=plan.group.slot)
end

function _native_model_csv_value(io,value)
    if value isa AbstractFloat
        @printf(io,"%.17g",value)
    else
        string_value=string(value)
        if any(character->character in (',','"','\r','\n'),string_value)
            print(io,'"',replace(string_value,"\""=>"\"\""),'"')
        else
            print(io,string_value)
        end
    end
end

function _native_model_view_csv_payload(exports,kind)
    components=maximum((table.components for dataset in exports for table in dataset.tables);init=0)
    header=vcat(["record","dataset","model","step","step_last","time","kind","components","uid",
        "value_node","node_uid","entity_dim","entity_tag","x","y","z","msh","physical_tag","connectivity"],
        ["value_$component" for component in 1:components])
    columns=Dict(name=>index for (index,name) in enumerate(header))
    io=IOBuffer();println(io,join(header,','))
    function row(record,fields;values=Float64[])
        entries=Any["" for _ in header];entries[1]=record
        for (key,value) in pairs(fields);entries[columns[string(key)]]=value;end
        for (component,value) in enumerate(values);entries[length(header)-components+component]=value;end
        for (column,value) in enumerate(entries)
            column==1 || print(io,',');_native_model_csv_value(io,value)
        end
        println(io)
    end
    for (dataset_index,dataset) in enumerate(exports)
        mesh=dataset.mesh;data=dataset.data;name=dataset.slot.name
        ordinal=dataset_index-1
        for range in dataset.ranges
            row("range",(dataset=ordinal,model=name,step=range.first,step_last=range.last,kind,
                components=range.components))
        end
        for (position,uid) in enumerate(data.external_node_tags)
            owner=data.node_entities[position]
            row("node",(dataset=ordinal,model=name,uid,node_uid=uid,entity_dim=owner[1],entity_tag=owner[2],
                x=mesh.coords[1,position],y=mesh.coords[2,position],z=mesh.coords[3,position]))
        end
        for cell in dataset.cells
            row("cell",(dataset=ordinal,model=name,uid=cell.uid,entity_dim=cell.owner[1],entity_tag=cell.owner[2],
                msh=cell.msh,physical_tag=cell.physical,connectivity=join(cell.node_uids,';')))
        end
        for table in dataset.tables,uid in sort!(collect(keys(table.rows)))
            values=table.rows[uid];count=length(values)÷table.components
            for node in 1:count
                common=(dataset=ordinal,model=name,step=table.step,time=table.time,kind=table.kind,
                    components=table.components,uid,value_node=node-1)
                if table.kind=="NodeData"
                    position=dataset.nodes[uid];owner=data.node_entities[position]
                    geometry=(node_uid=uid,entity_dim=owner[1],entity_tag=owner[2],
                        x=mesh.coords[1,position],y=mesh.coords[2,position],z=mesh.coords[3,position])
                else
                    cell=dataset.cell_rows[uid]
                    geometry=(entity_dim=cell.owner[1],entity_tag=cell.owner[2],msh=cell.msh)
                    if table.kind=="ElementNodeData" && node<=length(cell.nodes)
                        position=cell.nodes[node]
                        geometry=merge(geometry,(node_uid=cell.node_uids[node],
                            x=mesh.coords[1,position],y=mesh.coords[2,position],z=mesh.coords[3,position]))
                    end
                end
                first=table.components*(node-1)+1
                row("data",merge(common,geometry);values=@view(values[first:first+table.components-1]))
            end
        end
    end
    return take!(io)
end

# Quadratic VTK cells use these canonical node identities. The
# reference lattice is independent of model coordinates, external UIDs and
# discontinuous nodal values. See the corresponding Kitware VTK cell classes.
function _native_model_vtk_quadratic_layout(msh,caller)
    type,reference = if msh in (10,16)
        (msh==10 ? 28 : 23,((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),
             (.5,0.,0.),(1.,.5,0.),(.5,1.,0.),(0.,.5,0.),(.5,.5,0.)))
    elseif msh in (12,17)
        (msh==12 ? 29 : 25,((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),
             (0.,0.,1.),(1.,0.,1.),(1.,1.,1.),(0.,1.,1.),
             (.5,0.,0.),(1.,.5,0.),(.5,1.,0.),(0.,.5,0.),
             (.5,0.,1.),(1.,.5,1.),(.5,1.,1.),(0.,.5,1.),
             (0.,0.,.5),(1.,0.,.5),(1.,1.,.5),(0.,1.,.5),
             (0.,.5,.5),(1.,.5,.5),(.5,0.,.5),(.5,1.,.5),
             (.5,.5,0.),(.5,.5,1.),(.5,.5,.5)))
    elseif msh in (13,18)
        (msh==13 ? 32 : 26,((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),
             (0.,0.,1.),(1.,0.,1.),(0.,1.,1.),
             (.5,0.,0.),(.5,.5,0.),(0.,.5,0.),
             (.5,0.,1.),(.5,.5,1.),(0.,.5,1.),
             (0.,0.,.5),(1.,0.,.5),(0.,1.,.5),
             (.5,0.,.5),(.5,.5,.5),(0.,.5,.5)))
    elseif msh==19
        # Compare corner/edge identities in the Euclidean pyramid lattice.
        # VTK's collapsed parametric apex coordinates are not this lattice.
        (27,((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),(.5,.5,1.),
             (.5,0.,0.),(1.,.5,0.),(.5,1.,0.),(0.,.5,0.),
             (.25,.25,.5),(.75,.25,.5),(.75,.75,.5),(.25,.75,.5)))
    else
        throw(ArgumentError("$caller: unsupported complete quadratic VTK cell type $msh"))
    end
    msh==16 && (reference=reference[1:8])
    msh==17 && (reference=reference[1:20])
    msh==18 && (reference=reference[1:15])
    family=Elements.msh_family(msh)
    local_nodes=Elements.lagrange_nodes(msh)
    identities=Dict{NTuple{3,Float64},Int}()
    for node in axes(local_nodes,2)
        identity=ntuple(3) do axis
            value=local_nodes[axis,node]
            family in (:qua,:pyr) ? (axis==3 ? value : (value+1)/2) :
                (family==:hex || axis==3 ? (value+1)/2 : value)
        end
        haskey(identities,identity) && throw(ArgumentError(
            "$caller: duplicate canonical node identity for Gmsh element type $msh"))
        identities[identity]=node
    end
    length(identities)==length(reference) && all(haskey(identities,node) for node in reference) ||
        throw(ArgumentError("$caller: incomplete canonical node mapping for Gmsh element type $msh"))
    return (type=type,order=[identities[node] for node in reference])
end

function _native_model_view_vtk_cell(msh,caller)
    linear=msh==15 ? 1 : msh==1 ? 3 : msh==2 ? 5 : msh==3 ? 9 :
        msh==4 ? 10 : msh==5 ? 12 : msh==6 ? 13 : msh==7 ? 14 : nothing
    if linear!==nothing
        return (type=linear,order=collect(1:Elements.msh_num_nodes(msh)))
    end
    msh==8 && return (type=21,order=collect(1:3))
    msh==9 && return (type=22,order=collect(1:6))
    # Gmsh's last two tetrahedral mid-edge nodes are (3,4),(2,4);
    # VTK's are (2,4),(3,4). Geometry and nodal data stay bound to the same UID.
    msh==11 && return (type=24,order=[1,2,3,4,5,6,7,8,10,9])
    msh in (10,12,13,16,17,18,19) && return _native_model_vtk_quadratic_layout(msh,caller)
    # VTK Lagrange curves infer their degree from the complete node count.
    # Both formats place the two endpoints first, then increasing edge nodes.
    spec=Elements.msh_spec(msh)
    if spec.family===:lin && 3<=spec.order<=10
        return (type=68,order=collect(1:spec.nnodes))
    end
    throw(ArgumentError("$caller: VTK model export has no complete node ordering for Gmsh element type $msh"))
end

function _native_model_vtk_array(io,name,components,tuples,type,values)
    println(io,name,' ',components,' ',tuples,' ',type)
    for (index,value) in enumerate(values)
        if value isa AbstractFloat
            @printf(io,"%.17g",value)
        else
            print(io,value)
        end
        index%components==0 ? println(io) : print(io,' ')
    end
    isempty(values) && println(io)
end

function _native_model_view_vtk_payload(dataset,kind,caller)
    mesh=dataset.mesh;data=dataset.data;cells=dataset.cells
    layouts=[_native_model_view_vtk_cell(cell.msh,caller) for cell in cells]
    nodal=kind=="ElementNodeData"
    if nodal
        # ElementNodeData can be discontinuous across a shared mesh node.
        # Duplicate cell points and keep unreferenced model nodes explicitly.
        points=Int32[];connectivity=Vector{Int}[];used=Set{Int32}()
        for (cell,layout) in zip(cells,layouts)
            offset=length(points);append!(points,cell.nodes);union!(used,cell.nodes)
            push!(connectivity,[offset+slot-1 for slot in layout.order])
        end
        append!(points,(Int32(row) for row in axes(mesh.coords,2) if !(Int32(row) in used)))
    else
        points=Int32.(axes(mesh.coords,2))
        connectivity=[[Int(cell.nodes[slot])-1 for slot in layout.order] for (cell,layout) in zip(cells,layouts)]
    end
    io=IOBuffer();println(io,"# vtk DataFile Version 2.0\nTessella model view\nASCII\nDATASET UNSTRUCTURED_GRID")
    println(io,"POINTS ",length(points)," double")
    for row in points;@printf(io,"%.17g %.17g %.17g\n",mesh.coords[1,row],mesh.coords[2,row],mesh.coords[3,row]);end
    println(io,"CELLS ",length(cells),' ',sum(length(nodes)+1 for nodes in connectivity;init=0))
    for nodes in connectivity;println(io,length(nodes),' ',join(nodes,' '));end
    println(io,"CELL_TYPES ",length(cells));foreach(layout->println(io,layout.type),layouts)
    println(io,"FIELD FieldData 3")
    _native_model_vtk_array(io,"STEP",1,length(dataset.tables),"vtktypeint64",[table.step for table in dataset.tables])
    _native_model_vtk_array(io,"TIME",1,length(dataset.tables),"double",[table.time for table in dataset.tables])
    _native_model_vtk_array(io,"ALLOCATED_RANGES",2,length(dataset.ranges),"vtktypeint64",
        [value for range in dataset.ranges for value in (range.first,range.last)])
    println(io,"POINT_DATA ",length(points),"\nFIELD PointData ",3+(kind=="ElementData" ? 0 : 2length(dataset.tables)))
    _native_model_vtk_array(io,"node_uid",1,length(points),"vtktypeuint64",data.external_node_tags[points])
    _native_model_vtk_array(io,"node_entity_dim",1,length(points),"int",[data.node_entities[row][1] for row in points])
    _native_model_vtk_array(io,"node_entity_tag",1,length(points),"int",[data.node_entities[row][2] for row in points])
    if kind!="ElementData"
        for table in dataset.tables
            values=fill(NaN,table.components,length(points));present=zeros(Int,length(points))
            if !nodal
                for (position,row) in enumerate(points)
                    stored=get(table.rows,data.external_node_tags[row],nothing)
                    stored===nothing && continue
                    values[:,position].=stored;present[position]=1
                end
            else
                offset=0
                for cell in cells
                    width=length(cell.nodes);stored=get(table.rows,cell.uid,nothing)
                    if stored!==nothing
                        values[:,offset+1:offset+width].=reshape(stored,table.components,width)
                        present[offset+1:offset+width].=1
                    end
                    offset+=width
                end
            end
            _native_model_vtk_array(io,"step_$(table.step)",table.components,length(points),"double",vec(values))
            _native_model_vtk_array(io,"present_$(table.step)",1,length(points),"int",present)
        end
    end
    println(io,"CELL_DATA ",length(cells),"\nFIELD CellData ",4+(kind=="ElementData" ? 2length(dataset.tables) : 0))
    _native_model_vtk_array(io,"element_uid",1,length(cells),"vtktypeuint64",[cell.uid for cell in cells])
    _native_model_vtk_array(io,"element_entity_dim",1,length(cells),"int",[cell.owner[1] for cell in cells])
    _native_model_vtk_array(io,"element_entity_tag",1,length(cells),"int",[cell.owner[2] for cell in cells])
    _native_model_vtk_array(io,"physical_tag",1,length(cells),"int",[cell.physical for cell in cells])
    if kind=="ElementData"
        for table in dataset.tables
            values=fill(NaN,table.components,length(cells));present=zeros(Int,length(cells))
            for (position,cell) in enumerate(cells)
                stored=get(table.rows,cell.uid,nothing);stored===nothing && continue
                values[:,position].=stored;present[position]=1
            end
            _native_model_vtk_array(io,"step_$(table.step)",table.components,length(cells),"double",vec(values))
            _native_model_vtk_array(io,"present_$(table.step)",1,length(cells),"int",present)
        end
    end
    return take!(io)
end

function _native_model_view_publications(view,plans,path,extension,caller)
    exports=[_native_model_view_export_plan(view,plan,caller;vtk=extension==".vtk") for plan in plans]
    extension==".csv" && return [(path=String(path),bytes=_native_model_view_csv_payload(exports,view.model_kind))]
    # Validate/build every held model before opening any split destination.
    return [(path=plan.group.path,bytes=_native_model_view_vtk_payload(dataset,view.model_kind,caller))
        for (plan,dataset) in zip(plans,exports)]
end
