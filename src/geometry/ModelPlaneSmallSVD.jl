# Independent small-matrix Plane fit. Column-pivoted Householder QR reduces
# the centered sample matrix to a 3x3 triangular matrix. Two-sided Jacobi
# rotations then supply its right singular vectors. Explicit reductions
# retain the arithmetic policy qualified against actual native Plane frames;
# neither a native library nor a restricted-source implementation is used.

function _gmsh_plane_squared_norm(matrix,column,first_row)
    count=size(matrix,1)-first_row+1
    packets=count÷2
    first_lane=0.0;second_lane=0.0
    if packets>0
        first_lane=matrix[first_row,column]^2
        second_lane=matrix[first_row+1,column]^2
        if packets>=2
            other_first=matrix[first_row+2,column]^2
            other_second=matrix[first_row+3,column]^2
            for packet in 3:2:packets-packets%2
                row=first_row+2*(packet-1)
                first_lane+=matrix[row,column]^2
                second_lane+=matrix[row+1,column]^2
                other_first+=matrix[row+2,column]^2
                other_second+=matrix[row+3,column]^2
            end
            first_lane+=other_first;second_lane+=other_second
            if isodd(packets)
                row=first_row+2*(packets-1)
                first_lane+=matrix[row,column]^2
                second_lane+=matrix[row+1,column]^2
            end
        end
    end
    result=first_lane+second_lane
    isodd(count) && (result+=matrix[first_row+count-1,column]^2)
    return result
end

function _gmsh_plane_reflector_dot(matrix,reflector,column,first_row;
                                  scalar_result::Bool=false)
    # A single remaining right-hand side is a scalar inner product. Its
    # two packet accumulators reduce before the last scalar, while the
    # matrix/vector product below streams into one pair of lanes.
    if scalar_result
        count=size(matrix,1)-first_row+1
        packets=count÷2
        first_lane=0.0;second_lane=0.0
        if packets>0
            first_lane=reflector[1]*matrix[first_row,column]
            second_lane=reflector[2]*matrix[first_row+1,column]
            if packets>=2
                other_first=reflector[3]*matrix[first_row+2,column]
                other_second=reflector[4]*matrix[first_row+3,column]
                for packet in 3:2:packets-packets%2
                    offset=2*(packet-1);row=first_row+offset
                    first_lane+=reflector[offset+1]*matrix[row,column]
                    second_lane+=reflector[offset+2]*matrix[row+1,column]
                    other_first+=reflector[offset+3]*matrix[row+2,column]
                    other_second+=reflector[offset+4]*matrix[row+3,column]
                end
                first_lane+=other_first;second_lane+=other_second
                if isodd(packets)
                    offset=2*(packets-1);row=first_row+offset
                    first_lane+=reflector[offset+1]*matrix[row,column]
                    second_lane+=reflector[offset+2]*matrix[row+1,column]
                end
            end
        end
        result=first_lane+second_lane
        isodd(count) && (result+=reflector[count]*matrix[first_row+count-1,column])
        return result
    end
    first_lane=0.0;second_lane=0.0
    count=size(matrix,1)-first_row+1
    paired=count-count%2
    for position in 0:2:paired-1
        row=first_row+position
        first_lane+=reflector[position+1]*matrix[row,column]
        second_lane+=reflector[position+2]*matrix[row+1,column]
    end
    result=first_lane+second_lane
    if isodd(count)
        result+=reflector[count]*matrix[first_row+count-1,column]
    end
    return result
end

function _gmsh_plane_right_singular_vectors(matrix)
    rows=size(matrix,1)
    rows>=3 && size(matrix,2)==3 || throw(ArgumentError(
        "native mean Plane fit requires at least three 3D samples"))
    scale=maximum(abs,matrix)
    isfinite(scale) || throw(ArgumentError("native mean Plane samples are not finite"))
    factor=scale==0.0 ? 1.0 : scale
    work=matrix./factor
    norms=[sqrt(_gmsh_plane_squared_norm(work,column,1)) for column in 1:3]
    direct=copy(norms);permutation=[1,2,3]
    for pivot_column in 1:3
        selected=pivot_column
        for column in pivot_column+1:3
            norms[column]>norms[selected] && (selected=column)
        end
        if selected!=pivot_column
            for row in 1:rows
                work[row,pivot_column],work[row,selected]=work[row,selected],work[row,pivot_column]
            end
            norms[pivot_column],norms[selected]=norms[selected],norms[pivot_column]
            direct[pivot_column],direct[selected]=direct[selected],direct[pivot_column]
            permutation[pivot_column],permutation[selected]=permutation[selected],permutation[pivot_column]
        end
        initial=work[pivot_column,pivot_column]
        tail=_gmsh_plane_squared_norm(work,pivot_column,pivot_column+1)
        reflector=zeros(rows-pivot_column)
        if tail==0.0
            diagonal=initial;coefficient=0.0
        else
            diagonal=sqrt(initial*initial+tail)
            initial>=0.0 && (diagonal=-diagonal)
            for row in pivot_column+1:rows
                reflector[row-pivot_column]=work[row,pivot_column]/(initial-diagonal)
            end
            coefficient=(diagonal-initial)/diagonal
        end
        work[pivot_column,pivot_column]=diagonal
        for row in pivot_column+1:rows
            work[row,pivot_column]=reflector[row-pivot_column]
        end
        for column in pivot_column+1:3
            projection=_gmsh_plane_reflector_dot(work,reflector,column,pivot_column+1;
                scalar_result=size(work,2)-pivot_column==1)+work[pivot_column,column]
            work[pivot_column,column]-=coefficient*projection
            for row in pivot_column+1:rows
                work[row,column]-=(coefficient*reflector[row-pivot_column])*projection
            end
            if norms[column]!=0.0
                ratio=abs(work[pivot_column,column])/norms[column]
                remaining=max(0.0,(1.0+ratio)*(1.0-ratio))
                relative=norms[column]/direct[column]
                if remaining*(relative*relative)<=sqrt(eps(Float64))
                    direct[column]=sqrt(_gmsh_plane_squared_norm(work,column,pivot_column+1))
                    norms[column]=direct[column]
                else
                    norms[column]*=sqrt(remaining)
                end
            end
        end
    end
    triangular=[row<=column ? work[row,column] : 0.0 for row in 1:3,column in 1:3]
    vectors=[row==permutation[column] ? 1.0 : 0.0 for row in 1:3,column in 1:3]
    maximum_diagonal=maximum(abs(triangular[axis,axis]) for axis in 1:3)
    for sweep in 1:100
        changed=false
        for (first,second) in ((2,1),(3,1),(3,2))
            threshold=max(floatmin(Float64),2.0*eps(Float64)*maximum_diagonal)
            max(abs(triangular[first,second]),abs(triangular[second,first]))>threshold || continue
            changed=true
            a=triangular[first,first];b=triangular[first,second]
            c=triangular[second,first];d=triangular[second,second]
            skew=c-b
            if abs(skew)<floatmin(Float64)
                polar_cosine=1.0;polar_sine=0.0
            else
                ratio=(a+d)/skew
                denominator=sqrt(1.0+ratio*ratio)
                polar_cosine=ratio/denominator;polar_sine=1.0/denominator
            end
            left_diagonal=polar_cosine*a+polar_sine*c
            right_diagonal=-polar_sine*b+polar_cosine*d
            symmetric=polar_cosine*b+polar_sine*d
            if abs(symmetric)<floatmin(Float64)
                cosine=1.0;sine=0.0
            else
                difference=(left_diagonal-right_diagonal)/(2.0*abs(symmetric))
                magnitude=1.0/(abs(difference)+sqrt(1.0+difference*difference))
                tangent=difference>0.0 ? magnitude : -magnitude
                cosine=1.0/sqrt(1.0+tangent*tangent)
                sine=-copysign(1.0,symmetric)*tangent*cosine
            end
            left_cosine=polar_cosine*cosine+polar_sine*sine
            left_sine=polar_sine*cosine-polar_cosine*sine
            for column in 1:3
                x=triangular[first,column];y=triangular[second,column]
                triangular[first,column]=left_cosine*x+left_sine*y
                triangular[second,column]=-left_sine*x+left_cosine*y
            end
            for row in 1:3
                x=triangular[row,first];y=triangular[row,second]
                triangular[row,first]=cosine*x-sine*y
                triangular[row,second]=sine*x+cosine*y
                x=vectors[row,first];y=vectors[row,second]
                vectors[row,first]=cosine*x-sine*y
                vectors[row,second]=sine*x+cosine*y
            end
            maximum_diagonal=max(maximum_diagonal,abs(triangular[first,first]),abs(triangular[second,second]))
        end
        if !changed
            values=[abs(triangular[axis,axis])*factor for axis in 1:3]
            order=sortperm(values;rev=true,alg=Base.Sort.MergeSort)
            return values[order],vectors[:,order]
        end
    end
    throw(ArgumentError("native mean Plane singular-vector iteration did not converge"))
end

function _gmsh_native_plane_small_svd_frame(points)
    center=zeros(3)
    for point in points,axis in 1:3
        center[axis]+=point[axis]
    end
    for axis in 1:3;center[axis]/=length(points);end
    matrix=[point[axis]-center[axis] for point in points,axis in 1:3]
    values,vectors=_gmsh_plane_right_singular_vectors(matrix)
    selected=values[1]<values[2] && values[1]<values[3] ? 1 :
        values[2]<values[1] && values[2]<values[3] ? 2 : 3
    normal=Mesh2D._gmsh_normalize3(ntuple(axis->vectors[axis,selected],3))
    reference=normal[1]==0.0 ? (1.0,0.0,0.0) :
        normal[2]==0.0 ? (0.0,1.0,0.0) : (0.0,0.0,1.0)
    first=Mesh2D._gmsh_normalize3(Mesh2D._gmsh_cross3(normal,reference))
    second=Mesh2D._gmsh_normalize3(Mesh2D._gmsh_cross3(first,normal))
    rhs=Mesh2D._gmsh_dot3(center,normal)
    dominant=abs(normal[1])>=abs(normal[2]) && abs(normal[1])>=abs(normal[3]) ? 1 :
        abs(normal[2])>=abs(normal[1]) && abs(normal[2])>=abs(normal[3]) ? 2 : 3
    origin=ntuple(axis->axis==dominant ? rhs/normal[dominant] : 0.0,3)
    all(isfinite,(origin...,first...,second...,normal...)) || throw(ArgumentError(
        "native mean Plane singular-vector fit is not finite"))
    return _GmshNativePlaneFrame(origin,first,second,normal)
end
