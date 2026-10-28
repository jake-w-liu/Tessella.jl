# Public higher-order inversion follows MElement::xyz2uvw. Its stopping rule
# and original-vertex sums are independent of the stable standalone FEM maps.
function _native_cell_query_shape(msh,uvw,data,caller)
    data===nothing || return _native_cell_shape(data,uvw)
    coordinates=Float64[uvw...]
    _,values,_=mesh_basis_functions(msh,coordinates,"Lagrange";caller)
    _,gradients,_=mesh_basis_functions(msh,coordinates,"GradLagrange";caller)
    return values,[(gradients[3i-2],gradients[3i-1],gradients[3i]) for i in eachindex(values)]
end

@inline function _native_cell_inverse_step(frame,coordinates,rhs)
    a,b,c,d,e,f,g,h,i=frame
    determinant=a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
    determinant==0.0 && return coordinates
    inverse=1.0/determinant
    u=coordinates[1]+((e*i-f*h)*inverse)*rhs[1]+
        (-(d*i-f*g)*inverse)*rhs[2]+((d*h-e*g)*inverse)*rhs[3]
    v=coordinates[2]+(-(b*i-c*h)*inverse)*rhs[1]+
        ((a*i-c*g)*inverse)*rhs[2]+(-(a*h-b*g)*inverse)*rhs[3]
    w=coordinates[3]+((b*f-c*e)*inverse)*rhs[1]+
        (-(a*f-c*d)*inverse)*rhs[2]+((a*e-b*d)*inverse)*rhs[3]
    return u,v,w
end

function _native_cell_local_coordinates(mesh,cells,cell,msh,p,caller,tag;location=false)
    spec=msh_spec(msh)
    _native_cell_numeric_contract(spec,caller)
    references=_native_cell_reference_nodes(spec)
    nearest=1
    first=_mixed_node(mesh,cells,cell,1)
    dx=first[1]-p[1];dy=first[2]-p[2];dz=first[3]-p[3]
    distance=dx*dx+dy*dy+dz*dz
    for node in 2:spec.nnodes
        vertex=_mixed_node(mesh,cells,cell,node)
        dx=vertex[1]-p[1];dy=vertex[2]-p[2];dz=vertex[3]-p[3]
        candidate=dx*dx+dy*dy+dz*dz
        if candidate<distance
            nearest=node;distance=candidate
        end
    end
    coordinates=(references[1,nearest],references[2,nearest],references[3,nearest])
    data=spec.family in (:tri,:qua,:tet,:hex,:pri) ? _native_cell_polynomial(msh,caller) : nothing
    iteration=1;error=1.0
    while error>1.e-6 && iteration<20
        shape=_native_cell_query_shape(msh,coordinates,data,caller)
        frame,determinant,mapped=_native_cell_public_frame(mesh,cells,cell,spec.dim,shape,caller,tag)
        determinant==0.0 && break
        updated=_native_cell_inverse_step(frame,coordinates,
            (p[1]-mapped[1],p[2]-mapped[2],p[3]-mapped[3]))
        all(isfinite,updated) || throw(ArgumentError(
            "$caller: local coordinates for element $tag are not Float64-representable"))
        du=updated[1]-coordinates[1];dv=updated[2]-coordinates[2];dw=updated[3]-coordinates[3]
        error=sqrt(du*du+dv*dv+dw*dw)
        coordinates=updated;iteration+=1
    end
    shape=_native_cell_query_shape(msh,coordinates,data,caller)
    _,_,mapped=_native_cell_public_frame(mesh,cells,cell,spec.dim,shape,caller,tag)
    scale=_mixed_scale(mesh,cells,cell)
    residual=hypot(p[1]-mapped[1],p[2]-mapped[2],p[3]-mapped[3])
    return coordinates,scale>0.0 ? residual/scale : residual,spec.dim
end
