# Native HighOrder.cpp Curve factories. Geometry placement and original mesh
# interpolation are independent, including conversion from an existing LineN.
function _dim01_curve_order_mesh_point(assembly,cell,coordinate,coefficient)
    weights=_dim01_native_line_weights(msh_spec(cell.msh).order,coordinate,coefficient)
    point=(0.0,0.0,0.0)
    for (position,weight) in zip(cell.nodes,weights)
        source=assembly.nodes[position].point
        point=(point[1]+weight*source[1],point[2]+weight*source[2],
               point[3]+weight*source[3])
    end
    return point
end

function _dim01_curve_order!(assembly,m,entities,order;only_visible=false,linear=false)
    first_orders=Dict{Tuple{Int,Int32},Int}()
    point_positions=Dict{Tuple{Int,Int32},Int32}()
    for (position,node) in enumerate(assembly.nodes)
        node.owner[1]==0 && get!(point_positions,node.owner,Int32(position))
    end
    for cell in assembly.cells
        cell.owner[1]==1 || continue
        get!(first_orders,cell.owner,msh_spec(cell.msh).order)
    end
    geometries=Dict{Int,Any}()
    coefficients=Dict{Int,Matrix{Float64}}()
    changed=false
    for cell in assembly.cells
        cell.owner[1]==1 || continue
        only_visible && model_entity_visibility(m,1,Int(cell.owner[2]))==0 && continue
        first_orders[cell.owner]==order && continue
        target=msh_type(:lin,order)
        curve=Int(cell.owner[2])
        geometry=linear ? nothing : get!(geometries,curve) do
            Model._model_high_order_curve_geometry(m,curve,assembly.caller)
        end
        source=copy(cell.nodes)
        coefficient=if geometry===nothing
            source_order=msh_spec(cell.msh).order
            get!(coefficients,source_order) do
                _dim01_native_line_coefficients(source_order)
            end
        else
            nothing
        end
        supports=Int32[]
        sizehint!(supports,order-1)
        parameters=if geometry===nothing
            Float64[]
        else
            first=assembly.nodes[source[1]]
            last=assembly.nodes[source[2]]
            lower=_api_p2_curve_parameter(m,curve,geometry,source[1],
                first.point,first.parameter,point_positions)
            upper=_api_p2_curve_parameter(m,curve,geometry,source[2],
                last.point,last.parameter,point_positions;periodic_end=true)
            Model._model_high_order_curve_parameters(geometry,lower,upper,order,assembly.caller)
        end
        for index in 1:order-1
            coordinate=Float64(index)*(2.0/order)-1.0
            parameter=geometry===nothing ? nothing : Float64[parameters[index]]
            point=geometry===nothing ?
                _dim01_curve_order_mesh_point(assembly,cell,coordinate,coefficient) :
                Model._model_high_order_curve_point(geometry,parameters[index],assembly.caller)
            push!(supports,_dim01_add_node!(assembly,_dim01_next_tag!(assembly,true),
                point,cell.owner,parameter))
        end
        cell.nodes=vcat(source[1:2],supports)
        cell.msh=target
        cell.tag=_dim01_next_tag!(assembly,false)
        changed=true
    end
    return changed
end
