# Private sampled native quality providers for public higher-order queries.
# Public minSICN samples ideal-gradient condition numbers. Public minSIGE
# samples physical-gradient edge-length ratios. Neither subdivides Bernstein
# coefficients and neither is an alias for adaptive minIsotropy.
struct _NativeCellSampleRule
    gradients::NTuple{3,Matrix{Float64}}
    ideal_gradients::NTuple{3,Matrix{Float64}}
    primary_gradients::Matrix{Float64}
    dimension::Int
    family::Symbol
end

const _NATIVE_CELL_SAMPLE_RULES=Dict{Int,_NativeCellSampleRule}()

function _native_cell_quality_points(family,p)
    points=NTuple{3,Float64}[]
    if family===:tri
        scale=1.0/p
        for j in 0:p,i in 0:p-j
            push!(points,(Float64(i)*scale,Float64(j)*scale,0.0))
        end
    elseif family===:tet
        scale=1.0/p
        for k in 0:p,j in 0:p-k,i in 0:p-j-k
            push!(points,(Float64(i)*scale,Float64(j)*scale,Float64(k)*scale))
        end
    elseif family===:qua
        scale=2.0/p
        for j in 0:p,i in 0:p
            push!(points,(Float64(i)*scale-1.0,Float64(j)*scale-1.0,0.0))
        end
    elseif family===:hex
        scale=2.0/p
        for k in 0:p,j in 0:p,i in 0:p
            push!(points,(Float64(i)*scale-1.0,Float64(j)*scale-1.0,Float64(k)*scale-1.0))
        end
    elseif family===:pri
        # FuncSpaceData(TYPE_PRI) samples the tensor product of the order-p
        # triangle lattice in (u,v) with the order-p [-1,1] line lattice in w.
        scale=1.0/p
        for k in 0:p,j in 0:p,i in 0:p-j
            push!(points,(Float64(i)*scale,Float64(j)*scale,Float64(2k)*scale-1.0))
        end
    elseif family===:pyr
        # FuncSpaceData(TYPE_PYR,true,1,p-1,false), unlike the tensorial
        # determinant/bounding space.
        scale=2.0/p
        for k in p-1:-1:0,j in 0:k+1,i in 0:k+1
            z=Float64(p-1-k)/p;width=1.0-z
            push!(points,(-width+Float64(i)*scale,-width+Float64(j)*scale,z))
        end
    else
        throw(ArgumentError("native sampled quality is unavailable for $family"))
    end
    return points
end

function _native_cell_ideal_gradients!(gradients,family)
    x,y,z=gradients
    if family in (:qua,:hex,:pyr)
        for position in eachindex(x)
            x[position]*=2.0;y[position]*=2.0
        end
    elseif family in (:tri,:tet,:pri)
        factor=2.0/sqrt(3.0);axpy=-1.0/sqrt(3.0)
        for position in eachindex(y)
            y[position]*=factor;y[position]+=x[position]*axpy
        end
    end
    if family in (:hex,:pri)
        for position in eachindex(z)
            z[position]*=2.0
        end
    elseif family===:pyr
        factor=sqrt(2.0)
        for position in eachindex(z)
            z[position]*=factor
        end
    elseif family===:tet
        factor=sqrt(1.5);ax=-3.0/2.0/sqrt(6.0);ay=-1.0/2.0/sqrt(2.0)
        for position in eachindex(z)
            z[position]*=factor;z[position]+=x[position]*ax;z[position]+=y[position]*ay
        end
    end
    return gradients
end

function _native_cell_sample_rule_build(msh,caller)
    spec=msh_spec(msh);family=spec.family
    _native_cell_numeric_contract(spec,caller)
    points=_native_cell_quality_points(family,spec.order)
    data=family===:pyr ? _native_bezier_pyramid_reference(msh,caller) : _native_cell_polynomial(msh,caller)
    gradients=ntuple(_->Matrix{Float64}(undef,length(points),spec.nnodes),3)
    for (sample,point) in enumerate(points)
        _,shape=family===:pyr ? _native_bezier_pyramid_shape(data,point) :
            _native_cell_query_shape(msh,point,data,caller)
        for node in 1:spec.nnodes,axis in 1:3
            gradients[axis][sample,node]=shape[node][axis]
        end
    end
    ideal=map(copy,gradients);_native_cell_ideal_gradients!(ideal,family)
    primary=msh_type(family,1)
    primary_data=family===:pyr ? _native_bezier_pyramid_reference(primary,caller) : _native_cell_polynomial(primary,caller)
    _,primary_shape=family===:pyr ? _native_bezier_pyramid_shape(primary_data,(0.0,0.0,0.0)) :
        _native_cell_query_shape(primary,(0.0,0.0,0.0),primary_data,caller)
    primary_gradients=Matrix{Float64}(undef,3,length(primary_shape))
    for node in eachindex(primary_shape),axis in 1:3
        primary_gradients[axis,node]=primary_shape[node][axis]
    end
    return _NativeCellSampleRule(gradients,ideal,primary_gradients,spec.dim,family)
end

function _native_cell_sample_rule(msh,caller)
    lock(_NATIVE_BEZIER_LOCK)
    try
        haskey(_NATIVE_CELL_SAMPLE_RULES,msh) && return _NATIVE_CELL_SAMPLE_RULES[msh]
        rule=_native_cell_sample_rule_build(msh,caller)
        _NATIVE_CELL_SAMPLE_RULES[msh]=rule
        return rule
    finally
        unlock(_NATIVE_BEZIER_LOCK)
    end
end

@inline _native_cell_sample_norm(a)=sqrt(a[1]*a[1]+a[2]*a[2]+a[3]*a[3])
@inline _native_cell_sample_sub(a,b)=(a[1]-b[1],a[2]-b[2],a[3]-b[3])
@inline _native_cell_sample_sqrt(x)=x<0.0 ? NaN : sqrt(x)

function _native_cell_sample_condition2(a,b,normal)
    xx=a[1]*a[1];yx=a[2]*a[2];zx=a[3]*a[3]
    xy=b[1]*b[1];yy=b[2]*b[2];zy=b[3]*b[3]
    dx=xx-xy;dy=yx-yy;cx=a[1]*b[1];cy=a[2]*b[2]
    s1=zy*zy;s2=(zx-dy-dx)*zy;s3=(cy+cx)*a[3]*b[3]
    s4=zx*zx;s5=(dy+dx)*zx;s6=cx*cy;s7=yx*yx;s8=dy*yx
    s9=xx*xx;s10=dx*xx;s11=dy*dy;s12=dx*dy;s13=dx*dx
    s=2.0*(s2+s5+s12)+4.0*(s7-s8+s9-s10)+8.0*(s3+s6)+s1+s4+s11+s13
    n=xx+xy+yx+yy+zx+zy
    root=s>0.0 ? sqrt(s) : 0.0
    first=0.5*(n+root);second=0.5*(n-root)
    condition=2.0*_native_cell_sample_sqrt(first*second)/(first+second)
    nx=a[2]*b[3]-a[3]*b[2];ny=a[3]*b[1]-a[1]*b[3];nz=a[1]*b[2]-a[2]*b[1]
    dot=nx*normal[1]+ny*normal[2]+nz*normal[3]
    return dot>=0.0 ? condition : -condition
end

function _native_cell_sample_condition3(a,b,c)
    x1,y1,z1=a;x2,y2,z2=b;x3,y3,z3=c
    determinant=x1*(y2*z3-z2*y3)-x2*(y1*z3-z1*y3)+x3*(y1*z2-z1*y2)
    determinant==0.0 && return 0.0
    i11=y2*z3-y3*z2;i12=x3*z2-x2*z3;i13=x2*y3-x3*y2
    i21=y3*z1-y1*z3;i22=x1*z3-x3*z1;i23=x3*y1-x1*y3
    i31=y1*z2-y2*z1;i32=x2*z1-x1*z2;i33=x1*y2-x2*y1
    normj=x1*x1+x2*x2+x3*x3+y1*y1+y2*y2+y3*y3+z1*z1+z2*z2+z3*z3
    norminv=i11*i11+i12*i12+i13*i13+i21*i21+i22*i22+i23*i23+i31*i31+i32*i32+i33*i33
    return 3.0*determinant/sqrt(normj*norminv)
end

function _native_cell_sample_ige(a,b,c,determinant,family)
    if family===:pyr
        l0=_native_cell_sample_norm((2a[1],2a[2],2a[3]))
        l1=_native_cell_sample_norm((2b[1],2b[2],2b[3]))
        l2=_native_cell_sample_norm((c[1]+a[1]+b[1],c[2]+a[2]+b[2],c[3]+a[3]+b[3]))
        l3=_native_cell_sample_norm((c[1]-a[1]+b[1],c[2]-a[2]+b[2],c[3]-a[3]+b[3]))
        l4=_native_cell_sample_norm((c[1]-a[1]-b[1],c[2]-a[2]-b[2],c[3]-a[3]-b[3]))
        l5=_native_cell_sample_norm((c[1]+a[1]-b[1],c[2]+a[2]-b[2],c[3]+a[3]-b[3]))
        return 4.0*sqrt(2.0)*determinant*(1.0/l0/l1/l2+1.0/l0/l1/l3+
            1.0/l0/l1/l4+1.0/l0/l1/l5+1.0/l2/l3/l4+1.0/l2/l3/l5+
            1.0/l4/l5/l2+1.0/l4/l5/l3)/8.0
    end
    l0=_native_cell_sample_norm(a);l1=_native_cell_sample_norm(b)
    family===:qua && return determinant/l0/l1
    family===:tri && return (2.0/sqrt(3.0))*determinant*(1.0/l0/l1+
        1.0/l0/_native_cell_sample_norm(_native_cell_sample_sub(b,a))+
        1.0/l1/_native_cell_sample_norm(_native_cell_sample_sub(b,a)))/3.0
    l2=_native_cell_sample_norm(c)
    family===:hex && return determinant/l0/l1/l2
    l3=_native_cell_sample_norm(_native_cell_sample_sub(b,a))
    l4=_native_cell_sample_norm(_native_cell_sample_sub(c,a))
    l5=_native_cell_sample_norm(_native_cell_sample_sub(c,b))
    return sqrt(2.0)*determinant*(1.0/l0/l5/l1+1.0/l0/l5/l2+1.0/l0/l5/l3+
        1.0/l0/l5/l4+1.0/l1/l4/l0+1.0/l1/l4/l2+1.0/l1/l4/l3+
        1.0/l1/l4/l5+1.0/l2/l3/l0+1.0/l2/l3/l1+1.0/l2/l3/l4+1.0/l2/l3/l5)/12.0
end

function _native_cell_quality_sampled(mesh,cells,cell,msh,name,caller)
    rule=_native_cell_sample_rule(msh,caller);nodes=size(cells,1)
    coordinates=Matrix{Float64}(undef,nodes,3)
    for node in 1:nodes
        point=_mixed_node(mesh,cells,cell,node)
        for axis in 1:3
            coordinates[node,axis]=point[axis]
        end
    end
    gradients=name=="minSICN" ? rule.ideal_gradients : rule.gradients
    derivatives=ntuple(_->zeros(Float64,size(gradients[1],1),3),3)
    for axis in 1:rule.dimension
        _native_bezier_product!(derivatives[axis],gradients[axis],coordinates)
    end
    normal=(0.0,0.0,0.0)
    if rule.dimension==2
        if name=="minSICN"
            first=_mixed_node(mesh,cells,cell,1);second=_mixed_node(mesh,cells,cell,2)
            third=_mixed_node(mesh,cells,cell,3)
            a=_native_cell_sample_sub(second,first);b=_native_cell_sample_sub(third,first)
        else
            primary=zeros(Float64,2,3)
            for node in axes(rule.primary_gradients,2),axis in 1:2,coordinate in 1:3
                primary[axis,coordinate]+=rule.primary_gradients[axis,node]*coordinates[node,coordinate]
            end
            a=(primary[1,1],primary[1,2],primary[1,3]);b=(primary[2,1],primary[2,2],primary[2,3])
        end
        raw=(a[2]*b[3]-a[3]*b[2],-a[1]*b[3]+a[3]*b[1],a[1]*b[2]-a[2]*b[1])
        scale=1.0/_native_cell_sample_norm(raw)
        normal=(raw[1]*scale,raw[2]*scale,raw[3]*scale)
    end
    minimum_value=floatmax(Float64)
    for sample in axes(derivatives[1],1)
        a=(derivatives[1][sample,1],derivatives[1][sample,2],derivatives[1][sample,3])
        b=(derivatives[2][sample,1],derivatives[2][sample,2],derivatives[2][sample,3])
        c=(derivatives[3][sample,1],derivatives[3][sample,2],derivatives[3][sample,3])
        value=if name=="minSICN"
            rule.dimension==2 ? _native_cell_sample_condition2(a,b,normal) : _native_cell_sample_condition3(a,b,c)
        else
            determinant=_native_bezier_det(a,b,rule.dimension==2 ? normal : c)
            _native_cell_sample_ige(a,b,c,determinant,rule.family)
        end
        # SICN uses std::min_element; sampled SIGE scans from finite DBL_MAX.
        # Both ignore later NaNs, while SICN retains a NaN in its first row.
        name=="minSICN" && sample==1 && (minimum_value=value)
        value<minimum_value && (minimum_value=value)
    end
    return minimum_value
end
