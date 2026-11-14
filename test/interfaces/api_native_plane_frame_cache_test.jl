module PlaneFrameCacheChecks
using Test, Tessella, Serialization
const A=Tessella.API
const M=Tessella.Model

# Preserve Float64 payloads and signed zeros when comparing retained and cached paths.
function exact_value(value)
    value isa Float64 && return (Float64,reinterpret(UInt64,value))
    if value isa AbstractVector || value isa Tuple
        return (typeof(value),map(exact_value,value))
    end
    return (typeof(value),value)
end
function outcome(f)
    try
        return (value=f(),error=nothing)
    catch error
        error isa InterruptException && rethrow()
        return (value=nothing,error=(typeof(error),sprint(showerror,error)))
    end
end
function snapshot(value)
    stream=IOBuffer();serialize(stream,value)
    return take!(stream)
end
function equal_outcomes(name,old,new;m,a,expected_error=false)
    before=(snapshot(m),snapshot(a))
    old_result=outcome(old)
    after_old=(snapshot(m),snapshot(a))
    new_result=outcome(new)
    after_new=(snapshot(m),snapshot(a))
    @testset "$name" begin
        @test before==after_old==after_new
        @test exact_value(old_result.value)==exact_value(new_result.value) &&
            old_result.error==new_result.error
        expected_error && @test old_result.error!==nothing
    end
    return new_result
end
function add_square!(m;offset=0,surface=1,coordinates=((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.)))
    for (index,point) in enumerate(coordinates)
        M.add_point!(m,point...;tag=offset+index,mesh_size=.5)
    end
    for (index,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        M.add_line!(m,offset+first,offset+last;tag=offset+index)
    end
    M.add_curve_loop!(m,collect(offset+1:offset+4);tag=offset+1)
    M.add_plane_surface!(m,[offset+1];tag=surface)
    return m
end
square()=add_square!(M.GeoModel())
function node(index,point,dimension,entity,parameter=nothing)
    return A._Dim01Node(UInt64(index),Tuple(Float64.(point)),(dimension,Int32(entity)),
        parameter===nothing ? nothing : Float64[parameter...])
end
function assembly(nodes)
    return A._Dim01Assembly(collect(nodes),A._Dim01Cell[],
        Dict(n.tag=>Int32(i) for (i,n) in enumerate(nodes)),UInt64(length(nodes)),UInt64(0),
        100000,100000,"plane frame control")
end
edge_old(m,a,i,s=Int32(1))=A._native_order_surface_uv(m,a,Int32(i),s)
face_old(m,a,i,s=Int32(1))=A._native_order_surface_face_uv(m,a,Int32(i),s)
edge_new(frames,m,a,i,s=Int32(1))=A._dim01_native_plane_surface_uv_cached!(frames,m,a,Int32(i),s)
face_new(frames,m,a,i,s=Int32(1))=A._dim01_native_plane_surface_face_uv_cached!(frames,m,a,Int32(i),s)

function endpoint_controls!()
    m=square();M._native_geo_kernel_admit_point!(m,1)
    a=assembly([node(1,(.25,.375,.125),0,1)])
    before=(snapshot(m),snapshot(a));frames=Dict{Int,Any}()
    edge=equal_outcomes("native Point EDGE uses CAD XYZ",()->edge_old(m,a,1),()->edge_new(frames,m,a,1);m,a)
    face=equal_outcomes("native Point FACE uses edited mesh XYZ",()->face_old(m,a,1),()->face_new(frames,m,a,1);m,a)
    @test edge.error===nothing && face.error===nothing
    @test exact_value(edge.value)!=exact_value(face.value)
    @test length(frames)==1
    frame=frames[1]
    edge_new(frames,m,a,1)
    @test frames[1]===frame
    @test before==(snapshot(m),snapshot(a))
    delete!(m.meshing.native_geo_kernel.points,1);frames=Dict{Int,Any}()
    result=equal_outcomes("non-native Point FACE delegates EDGE",()->face_old(m,a,1),()->face_new(frames,m,a,1);m,a)
    @test exact_value(result.value)==exact_value(outcome(()->edge_old(m,a,1)).value)
    # Saved Surface UV bypasses all validation and frame construction, including
    # deliberately unused trailing values and a missing target CAD entity.
    for uv in (Float64[.25,-0.0,NaN],Float64[NaN,Inf,-Inf])
        a=assembly([node(1,(NaN,Inf,-Inf),2,999,uv)]);frames=Dict{Int,Any}()
        result=equal_outcomes("Surface saved first-two UV $(repr(uv))",()->edge_old(m,a,1,Int32(999)),
            ()->edge_new(frames,m,a,1,Int32(999));m,a)
        @test exact_value(result.value)==exact_value(uv[1:2])
        @test isempty(frames)
    end
    for parameter in (Float64[.25],Float64[.25,NaN],Float64[.25,Inf],Float64[-.5],Float64[1.5],
            Float64[-0.0],Float64[NaN],Float64[Inf],Float64[-Inf],Float64[],nothing)
        a=assembly([node(1,(.7,.8,.9),1,1,parameter)]);frames=Dict{Int,Any}()
        result=equal_outcomes("Curve used first slice $(repr(parameter))",()->edge_old(m,a,1),
            ()->edge_new(frames,m,a,1);m,a)
        if parameter===nothing || isempty(parameter)
            @test result.value===nothing && isempty(frames)
        elseif !isfinite(first(parameter))
            @test result.error!==nothing && isempty(frames)
            @test startswith(result.error[2],"ArgumentError: model_reparametrize_on_surface:")
        else
            @test result.error===nothing && length(frames)==1
            expected=M.model_reparametrize_on_surface(m,1,1,parameter[1:1],1)
            @test exact_value(result.value)==exact_value(expected)
        end
    end
    for point in ((.125,.5,.25),(NaN,.5,.25),(.125,Inf,.25),(.125,.5,-Inf))
        a=assembly([node(1,point,0,999)]);frames=Dict{Int,Any}()
        result=equal_outcomes("physical fallback public validation $(repr(point))",()->edge_old(m,a,1),
            ()->edge_new(frames,m,a,1);m,a)
        if !all(isfinite,point)
            @test result.error!==nothing && isempty(frames)
            @test startswith(result.error[2],"ArgumentError: model_parametrization:")
        end
    end
end

function cache_and_error_order_controls!()
    m=square();M.add_point!(m,2.,2.,0.;tag=11);M.add_point!(m,3.,2.,0.;tag=12)
    M.add_line!(m,11,12;tag=11)
    m.curves[11]=(997,998) # Bad source Curve is independent of the Plane boundary.
    m.surfaces[1]=Int[] # Target metadata exists, but its frame is invalid.
    a=assembly([node(1,(0.,0.,0.),1,11,[.25])]);frames=Dict{Int,Any}()
    cold=equal_outcomes("cold invalid Plane precedes independent bad Curve",()->edge_old(m,a,1),
        ()->edge_new(frames,m,a,1);m,a,expected_error=true)
    @test cold.error==outcome(()->M._model_plane_frame(
        m,1,"model_reparametrize_on_surface")).error
    @test isempty(frames)
    a=assembly([node(1,(0.,0.,0.),1,11,[NaN])]);frames=Dict{Int,Any}()
    invalid=equal_outcomes("used parameter validation precedes cold invalid frame",()->edge_old(m,a,1),
        ()->edge_new(frames,m,a,1);m,a,expected_error=true)
    @test occursin("parametric coordinates",invalid.error[2]) && isempty(frames)
    # With an established frame, the same independent bad Curve must keep its
    # original caller and point index. CAD is unchanged throughout each call.
    m=square();M.add_point!(m,2.,2.,0.;tag=11);M.add_point!(m,3.,2.,0.;tag=12)
    M.add_line!(m,11,12;tag=11);m.curves[11]=(997,998)
    a=assembly([node(1,(0.,0.,0.),1,11,[.25])]);frames=Dict{Int,Any}()
    frames[1]=M._model_plane_frame(m,1,"model_reparametrize_on_surface")
    frame=frames[1]
    warm=equal_outcomes("warm valid Plane retains bad Curve evaluation error",()->edge_old(m,a,1),
        ()->edge_new(frames,m,a,1);m,a,expected_error=true)
    @test occursin("Line[11] references unknown Point[997]",warm.error[2])
    @test frames[1]===frame
    m=square();add_square!(m;offset=10,surface=2,
        coordinates=((0.,0.,1.),(1.,0.,2.),(1.,1.,2.),(0.,1.,1.)))
    a=assembly([node(1,(.1,.2,.3),0,999)]);frames=Dict{Int,Any}()
    for surface in Int32[1,2,1,2]
        equal_outcomes("two operation-local Plane keys $surface",()->edge_old(m,a,1,surface),
            ()->edge_new(frames,m,a,1,surface);m,a)
    end
    @test sort!(collect(keys(frames)))==[1,2]
    old_basis=exact_value(M._model_plane_point(frames[1],.25,.5,"control",1))
    for tag in 1:4
        p=m.points[tag];m.points[tag]=(p[1],p[2],p[3]+2.)
    end
    fresh=Dict{Int,Any}()
    equal_outcomes("fresh operation after CAD change recomputes Plane basis",()->edge_old(m,a,1),
        ()->edge_new(fresh,m,a,1);m,a)
    new_basis=exact_value(M._model_plane_point(fresh[1],.25,.5,"control",1))
    @test old_basis!=new_basis
    @test fresh!==frames && fresh[1]!==frames[1]
end

function exclusion_controls!()
    for kind in ("missing","implicit","discrete","OCC","ruled","unknown_kind")
        m=square();surface=Int32(1)
        if kind=="missing"
            surface=Int32(999)
        elseif kind=="implicit"
            delete!(m.surfaces,1)
        elseif kind=="discrete"
            m.discrete[(2,1)]=M.DiscreteEntity()
        elseif kind=="OCC"
            m.surface_geometry[1]=(occ=:unsupported_control,)
        elseif kind=="ruled"
            m.surface_types[1]=:ruled
        else
            m.surface_types[1]=:unsupported_control
        end
        a=assembly([node(1,(.2,.3,.4),0,1)]);frames=Dict{Int,Any}()
        @test !A._dim01_native_plane_frame_reusable(m,Int(surface))
        equal_outcomes("excluded $kind EDGE delegates retained helper",()->edge_old(m,a,1,surface),
            ()->edge_new(frames,m,a,1,surface);m,a)
        equal_outcomes("excluded $kind FACE delegates retained helper",()->face_old(m,a,1,surface),
            ()->face_new(frames,m,a,1,surface);m,a)
        @test isempty(frames)
    end
end

function support_controls!()
    m=square()
    a=assembly([node(1,(.2,.3,0.),0,1),node(2,(.8,.4,0.),1,1,nothing),node(3,(0.,1.,0.),0,4)])
    cell=A._Dim01Cell(2,(2,Int32(1)),Int32[1,2,3],UInt64(17),UInt64(17),Int32(3))
    frames=Dict{Int,Any}();point=(.5,.35,0.)
    before=(snapshot(m),snapshot(a))
    result=A._dim01_native_plane_support!(frames,m,a,cell,Int32[1,2,3],4,cell.owner,point)
    first=M.model_parametrization(m,2,1,collect(a.nodes[1].point))
    last=M.model_parametrization(m,2,1,collect(a.nodes[2].point))
    parameter=Float64[first[axis]+(last[axis]-first[axis])*.5 for axis in 1:2]
    placed=M._model_plane_point(frames[1],parameter[1],parameter[2],a.caller,1)
    @test exact_value(result)==exact_value((placed,parameter,Int32[1,2],true))
    cad_first=M.model_reparametrize_on_surface(m,0,1,Float64[],1)
    one_sided=Float64[cad_first[axis]+(last[axis]-cad_first[axis])*.5 for axis in 1:2]
    @test exact_value(parameter)!=exact_value(one_sided) # BOTH fallback, not one.
    @test before==(snapshot(m),snapshot(a))
    nodes=A._Dim01Node[]
    for tag in 1:4
        p=m.points[tag];M._native_geo_kernel_admit_point!(m,tag)
        push!(nodes,node(tag,(p[1]+.03125,p[2]-.015625,p[3]),0,tag))
    end
    for tag in 1:4
        p=M._model_curve_point(m,tag,.5,"control",1)
        push!(nodes,node(4+tag,(p[1]-.0078125,p[2]+.0234375,p[3]),1,tag,[.5]))
    end
    a=assembly(nodes);cell=A._Dim01Cell(16,(2,Int32(1)),collect(Int32,1:8),UInt64(18),UInt64(18),Int32(4))
    for unavailable in (false,true)
        if unavailable
            n=a.nodes[5];a.nodes[5]=A._Dim01Node(n.tag,n.point,n.owner,nothing)
        end
        frames=Dict{Int,Any}();point=(.5,.5,0.)
        before=(snapshot(m),snapshot(a))
        result=A._dim01_native_plane_support!(frames,m,a,cell,Int32[1,2,3,4],9,cell.owner,point)
        row=only(A._native_inner_rows(:qua,2))
        uv=Union{Nothing,Vector{Float64}}[face_old(m,a,node) for node in cell.nodes]
        if all(value->value!==nothing,uv)
            parameter=zeros(Float64,2)
            for (column,weight) in row
                parameter[1]+=weight*uv[column][1]
                parameter[2]+=weight*uv[column][2]
            end
        else
            physical=(0.,0.,0.)
            for (column,weight) in row
                p=a.nodes[cell.nodes[column]].point
                physical=(physical[1]+weight*p[1],physical[2]+weight*p[2],physical[3]+weight*p[3])
            end
            parameter=M.model_parametrization(m,2,1,collect(physical))
        end
        placed=M._model_plane_point(frames[1],parameter[1],parameter[2],a.caller,1)
        @test exact_value(result)==exact_value((placed,parameter,copy(cell.nodes),true))
        @test before==(snapshot(m),snapshot(a))
        @test cell.physical==Int32(4) && cell.tag==UInt64(18) && cell.source_tag==UInt64(18)
    end
end

function cached_sweep(m,a,n)
    frames=Dict{Int,Any}();total=0.
    for i in 1:n
        uv=edge_new(frames,m,a,1)
        total+=uv[1]+uv[2]
    end
    return total,length(frames)
end
function resource_controls!()
    allocations=Int[]
    for n in (32,64,128)
        m=square();m.curve_params[1]=Float64[i/(2n+5) for i in 0:2n+5]
        a=assembly([node(1,(.25,.5,0.),1,1,[.25])])
        equal_outcomes("resource n$n single endpoint full retained-helper bits",()->edge_old(m,a,1),
            ()->edge_new(Dict{Int,Any}(),m,a,1);m,a)
        before=(snapshot(m),snapshot(a))
        for _ in 1:2
            @test last(cached_sweep(m,a,n))==1
        end
        samples=Int[]
        for _ in 1:2
            GC.gc();push!(samples,@allocated cached_sweep(m,a,n))
        end
        @test before==(snapshot(m),snapshot(a))
        push!(allocations,maximum(samples))
    end
    for i in 2:length(allocations)
        @test allocations[i]<=2.15*allocations[i-1]+65536
    end
end

@testset "Operation-local plane frame caching" begin
    endpoint_controls!()
    cache_and_error_order_controls!()
    exclusion_controls!()
    support_controls!()
    resource_controls!()
end
end
