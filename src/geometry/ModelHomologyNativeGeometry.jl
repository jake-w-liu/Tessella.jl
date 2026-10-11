import ..MeshQuadrature: mesh_integration_points
import ..MeshFunctionSpaces: mesh_basis_functions
import ..Elements: msh_spec

function _native_homology_distance(first,last)
    x=first[1]-last[1];y=first[2]-last[2];z=first[3]-last[3]
    return sqrt(x*x+y*y+z*z)
end

function _native_homology_measure(msh_type,nodes,coordinates)
    spec=msh_spec(msh_type);dim=spec.dim
    dim==0 && return 0.0
    points=NTuple{3,Float64}[]
    for node in nodes
        point=get(coordinates,node,nothing)
        point===nothing && throw(ArgumentError(
            "compute_homology!: source vertex $node has no coordinates"))
        push!(points,point)
    end
    # MLine and MTetrahedron inherit these corner-only volume methods at
    # higher orders as well. Surface high-order cells use MElement quadrature.
    spec.family===:lin && return _native_homology_distance(points[1],points[2])
    a,b,c=points[1],points[2],points[3]
    if spec.family===:tet
        d=points[4]
        ax,ay,az=b[1]-a[1],b[2]-a[2],b[3]-a[3]
        bx,by,bz=c[1]-a[1],c[2]-a[2],c[3]-a[3]
        cx,cy,cz=d[1]-a[1],d[2]-a[2],d[3]-a[3]
        determinant=ax*(by*cz-bz*cy)-bx*(ay*cz-az*cy)+cx*(ay*bz-az*by)
        return abs(determinant*0.166666666666666666)
    elseif spec.family===:tri && length(nodes)==3
        ux,uy,uz=b[1]-a[1],b[2]-a[2],b[3]-a[3]
        vx,vy,vz=c[1]-a[1],c[2]-a[2],c[3]-a[3]
        x,y,z=uy*vz-uz*vy,uz*vx-ux*vz,ux*vy-uy*vx
        return sqrt(x*x+y*y+z*z)/2.0
    elseif spec.family===:qua && length(nodes)==4
        d=points[4]
        aa=_native_homology_distance(a,b);bb=_native_homology_distance(b,c)
        cc=_native_homology_distance(c,d);dd=_native_homology_distance(d,a)
        m=_native_homology_distance(a,c);n=_native_homology_distance(b,d)
        mn=2.0*m*n;abcd=aa*aa-bb*bb+cc*cc-dd*dd
        radicand=mn*mn-abcd*abcd
        return radicand<0 ? NaN : sqrt(radicand)/4.0
    end
    local_points,weights=mesh_integration_points(msh_type,"Gauss"*string(spec.order+3))
    _,gradients,_=mesh_basis_functions(msh_type,local_points,"GradLagrange")
    volume=0.0;count=length(nodes)
    for point in eachindex(weights)
        jac=zeros(Float64,3,3)
        for slot in 1:count
            offset=3*(count*(point-1)+slot)-3
            for direction in 1:3, component in 1:3
                jac[direction,component]+=points[slot][component]*gradients[offset+direction]
            end
        end
        determinant=if dim==2
            x=jac[1,1]*jac[2,2]-jac[1,2]*jac[2,1]
            y=jac[1,3]*jac[2,1]-jac[1,1]*jac[2,3]
            z=jac[1,2]*jac[2,3]-jac[1,3]*jac[2,2]
            sqrt(x*x+y*y+z*z)
        else
            jac[1,1]*jac[2,2]*jac[3,3]+jac[1,3]*jac[2,1]*jac[3,2]+
            jac[1,2]*jac[2,3]*jac[3,1]-jac[1,3]*jac[2,2]*jac[3,1]-
            jac[1,1]*jac[2,3]*jac[3,2]-jac[1,2]*jac[2,1]*jac[3,3]
        end
        volume+=determinant*weights[point]
    end
    return abs(volume)
end
