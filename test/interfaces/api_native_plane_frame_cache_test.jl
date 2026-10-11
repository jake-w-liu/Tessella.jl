module PlaneOperationFrameControlsV2
using Test, Tessella, SHA, TOML, Serialization, InteractiveUtils
const A=Tessella.API
const M=Tessella.Model
const HISTORICAL_SOURCE_SHA="e0416dd3573a0ffe1c778d72d7948ca5ac33664371ffa3896c56a177b2d66953"
const FRAME_BODY_SHA="13f35a5273bda8f2b953a09fe4ba3d87c6a935abea3ee6a777fd3f61f514291b"
const SOURCE_SHA=bytes2hex(sha256(read(joinpath(pkgdir(Tessella),"src","interfaces","APIGenerate01.jl"))))
const ROWS=Dict{String,Any}[]
const METHOD_ROWS=Dict{String,Any}[]
const RESOURCE_ROWS=Dict{String,Any}[]

function box_count(value)
    value isa Core.Box && return 1
    value isa GlobalRef && return Int(value.mod===Core && value.name===:Box)
    value isa Core.CodeInfo && return sum(box_count,value.code;init=0)
    value isa Expr && return sum(box_count,value.args;init=0)
    value isa Core.NewvarNode && return box_count(value.slot)
    value isa QuoteNode && return box_count(value.value)
    return value===Core.Box ? 1 : 0
end
function boxed_reader_positive(value)
    state=0
    reader=()->state
    state=value
    return reader
end
bits(x::Float64)=string(reinterpret(UInt64,x);base=16,pad=16)
function encoded(value)
    value===nothing && return Dict("kind"=>"nothing")
    value isa Float64 && return Dict("kind"=>"Float64","bits"=>bits(value))
    value isa Bool && return Dict("kind"=>"Bool","value"=>value)
    value isa Integer && return Dict("kind"=>string(typeof(value)),"value"=>string(value))
    value isa AbstractVector || value isa Tuple || return Dict("kind"=>string(typeof(value)),"repr"=>repr(value))
    return Dict("kind"=>string(typeof(value)),"values"=>[encoded(x) for x in value])
end
function outcome(f)
    try
        value=f()
        return (value=value,record=Dict("status"=>"return","result"=>encoded(value)))
    catch error
        error isa InterruptException && rethrow()
        return (value=nothing,record=Dict("status"=>"error","type"=>string(typeof(error)),
            "message"=>sprint(showerror,error)))
    end
end
function equal_outcomes(name,old,new;frames=nothing,m,a,expected_error=false)
    before=(model_digest(m),assembly_digest(a))
    old_result=outcome(old)
    after_old=(model_digest(m),assembly_digest(a))
    new_result=outcome(new)
    after_new=(model_digest(m),assembly_digest(a))
    @test before==after_old==after_new
    @test old_result.record==new_result.record
    expected_error && @test old_result.record["status"]=="error"
    push!(ROWS,Dict("case"=>name,"reference"=>old_result.record,"candidate"=>new_result.record,
        "cached_plane_keys"=>frames===nothing ? Int[] : sort!(collect(keys(frames))),
        "model_before_sha256"=>before[1],"model_after_reference_sha256"=>after_old[1],
        "model_after_candidate_sha256"=>after_new[1],"model_unchanged"=>before[1]==after_old[1]==after_new[1],
        "assembly_before_sha256"=>before[2],"assembly_after_reference_sha256"=>after_old[2],
        "assembly_after_candidate_sha256"=>after_new[2],"assembly_unchanged"=>before[2]==after_old[2]==after_new[2]))
    return new_result
end
function model_digest(m)
    stream=IOBuffer();serialize(stream,m)
    return bytes2hex(sha256(take!(stream)))
end
function assembly_digest(assembly)
    stream=IOBuffer();serialize(stream,assembly)
    return bytes2hex(sha256(take!(stream)))
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
    before=(model_digest(m),assembly_digest(a));frames=Dict{Int,Any}()
    edge=equal_outcomes("native Point EDGE uses CAD XYZ",()->edge_old(m,a,1),()->edge_new(frames,m,a,1);frames,m,a)
    face=equal_outcomes("native Point FACE uses edited mesh XYZ",()->face_old(m,a,1),()->face_new(frames,m,a,1);frames,m,a)
    @test edge.record["status"]==face.record["status"]=="return"
    @test edge.record["result"]!=face.record["result"]
    @test length(frames)==1
    frame=frames[1]
    edge_new(frames,m,a,1)
    @test frames[1]===frame
    @test before==(model_digest(m),assembly_digest(a))
    delete!(m.meshing.native_geo_kernel.points,1);frames=Dict{Int,Any}()
    result=equal_outcomes("non-native Point FACE delegates EDGE",()->face_old(m,a,1),()->face_new(frames,m,a,1);frames,m,a)
    @test result.record["result"]==outcome(()->edge_old(m,a,1)).record["result"]
    # Saved Surface UV bypasses all validation and frame construction, including
    # deliberately unused trailing values and a missing target CAD entity.
    for uv in (Float64[.25,-0.0,NaN],Float64[NaN,Inf,-Inf])
        a=assembly([node(1,(NaN,Inf,-Inf),2,999,uv)]);frames=Dict{Int,Any}()
        result=equal_outcomes("Surface saved first-two UV $(repr(uv))",()->edge_old(m,a,1,Int32(999)),
            ()->edge_new(frames,m,a,1,Int32(999));frames,m,a)
        @test result.record["result"]==encoded(uv[1:2])
        @test isempty(frames)
    end
    for parameter in (Float64[.25],Float64[.25,NaN],Float64[.25,Inf],Float64[-.5],Float64[1.5],
            Float64[-0.0],Float64[NaN],Float64[Inf],Float64[-Inf],Float64[],nothing)
        a=assembly([node(1,(.7,.8,.9),1,1,parameter)]);frames=Dict{Int,Any}()
        result=equal_outcomes("Curve used first slice $(repr(parameter))",()->edge_old(m,a,1),
            ()->edge_new(frames,m,a,1);frames,m,a)
        if parameter===nothing || isempty(parameter)
            @test result.value===nothing && isempty(frames)
        elseif !isfinite(first(parameter))
            @test result.record["status"]=="error" && isempty(frames)
            @test startswith(result.record["message"],"ArgumentError: model_reparametrize_on_surface:")
        else
            @test result.record["status"]=="return" && length(frames)==1
            expected=M.model_reparametrize_on_surface(m,1,1,parameter[1:1],1)
            @test result.record["result"]==encoded(expected)
        end
    end
    for point in ((.125,.5,.25),(NaN,.5,.25),(.125,Inf,.25),(.125,.5,-Inf))
        a=assembly([node(1,point,0,999)]);frames=Dict{Int,Any}()
        result=equal_outcomes("physical fallback public validation $(repr(point))",()->edge_old(m,a,1),
            ()->edge_new(frames,m,a,1);frames,m,a)
        if !all(isfinite,point)
            @test result.record["status"]=="error" && isempty(frames)
            @test startswith(result.record["message"],"ArgumentError: model_parametrization:")
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
        ()->edge_new(frames,m,a,1);frames,m,a,expected_error=true)
    @test occursin("Plane[1] has no boundary loops",cold.record["message"])
    @test isempty(frames)
    a=assembly([node(1,(0.,0.,0.),1,11,[NaN])]);frames=Dict{Int,Any}()
    invalid=equal_outcomes("used parameter validation precedes cold invalid frame",()->edge_old(m,a,1),
        ()->edge_new(frames,m,a,1);frames,m,a,expected_error=true)
    @test occursin("parametric coordinates",invalid.record["message"]) && isempty(frames)
    # With an established frame, the same independent bad Curve must keep its
    # original caller and point index. CAD is unchanged throughout each call.
    m=square();M.add_point!(m,2.,2.,0.;tag=11);M.add_point!(m,3.,2.,0.;tag=12)
    M.add_line!(m,11,12;tag=11);m.curves[11]=(997,998)
    a=assembly([node(1,(0.,0.,0.),1,11,[.25])]);frames=Dict{Int,Any}()
    frames[1]=M._model_plane_frame(m,1,"model_reparametrize_on_surface")
    frame=frames[1]
    warm=equal_outcomes("warm valid Plane retains bad Curve evaluation error",()->edge_old(m,a,1),
        ()->edge_new(frames,m,a,1);frames,m,a,expected_error=true)
    @test occursin("Line[11] references unknown Point[997]",warm.record["message"])
    @test frames[1]===frame
    m=square();add_square!(m;offset=10,surface=2,
        coordinates=((0.,0.,1.),(1.,0.,2.),(1.,1.,2.),(0.,1.,1.)))
    a=assembly([node(1,(.1,.2,.3),0,999)]);frames=Dict{Int,Any}()
    for surface in Int32[1,2,1,2]
        equal_outcomes("two operation-local Plane keys $surface",()->edge_old(m,a,1,surface),
            ()->edge_new(frames,m,a,1,surface);frames,m,a)
    end
    @test sort!(collect(keys(frames)))==[1,2]
    old_basis=encoded(M._model_plane_point(frames[1],.25,.5,"control",1))
    for tag in 1:4
        p=m.points[tag];m.points[tag]=(p[1],p[2],p[3]+2.)
    end
    fresh=Dict{Int,Any}()
    equal_outcomes("fresh operation after CAD change recomputes Plane basis",()->edge_old(m,a,1),
        ()->edge_new(fresh,m,a,1);frames=fresh,m,a)
    new_basis=encoded(M._model_plane_point(fresh[1],.25,.5,"control",1))
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
            ()->edge_new(frames,m,a,1,surface);frames,m,a)
        equal_outcomes("excluded $kind FACE delegates retained helper",()->face_old(m,a,1,surface),
            ()->face_new(frames,m,a,1,surface);frames,m,a)
        @test isempty(frames)
    end
end

function support_controls!()
    m=square()
    a=assembly([node(1,(.2,.3,0.),0,1),node(2,(.8,.4,0.),1,1,nothing),node(3,(0.,1.,0.),0,4)])
    cell=A._Dim01Cell(2,(2,Int32(1)),Int32[1,2,3],UInt64(17),UInt64(17),Int32(3))
    frames=Dict{Int,Any}();point=(.5,.35,0.)
    before=(model_digest(m),assembly_digest(a))
    result=A._dim01_native_plane_support!(frames,m,a,cell,Int32[1,2,3],4,cell.owner,point)
    first=M.model_parametrization(m,2,1,collect(a.nodes[1].point))
    last=M.model_parametrization(m,2,1,collect(a.nodes[2].point))
    parameter=Float64[first[axis]+(last[axis]-first[axis])*.5 for axis in 1:2]
    placed=M._model_plane_point(frames[1],parameter[1],parameter[2],a.caller,1)
    @test encoded(result)==encoded((placed,parameter,Int32[1,2],true))
    cad_first=M.model_reparametrize_on_surface(m,0,1,Float64[],1)
    one_sided=Float64[cad_first[axis]+(last[axis]-cad_first[axis])*.5 for axis in 1:2]
    @test encoded(parameter)!=encoded(one_sided) # BOTH fallback, not one.
    @test before==(model_digest(m),assembly_digest(a))
    push!(ROWS,Dict("case"=>"one unavailable EDGE endpoint forces BOTH physical fallbacks",
        "candidate"=>encoded(result),"reference"=>encoded((placed,parameter,Int32[1,2],true)),
        "one_sided_negative_control"=>encoded(one_sided),
        "model_unchanged"=>before[1]==model_digest(m),"assembly_unchanged"=>before[2]==assembly_digest(a)))
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
        before=(model_digest(m),assembly_digest(a))
        result=A._dim01_native_plane_support!(frames,m,a,cell,Int32[1,2,3,4],9,cell.owner,point)
        row=only(A._native_inner_rows(:qua,2))
        uv=[face_old(m,a,node) for node in cell.nodes]
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
        @test encoded(result)==encoded((placed,parameter,copy(cell.nodes),true))
        @test before==(model_digest(m),assembly_digest(a))
        @test cell.physical==Int32(4) && cell.tag==UInt64(18) && cell.source_tag==UInt64(18)
        push!(ROWS,Dict("case"=>"Quad8 ordered FACE $(unavailable ? "physical fallback" : "retained vertex UV")",
            "candidate"=>encoded(result),"reference"=>encoded((placed,parameter,copy(cell.nodes),true)),
            "legacy_carrier_physical"=>cell.physical,"public_element_tag"=>string(cell.tag),
            "model_unchanged"=>before[1]==model_digest(m),"assembly_unchanged"=>before[2]==assembly_digest(a)))
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
    for n in (32,64,128)
        m=square();m.curve_params[1]=Float64[i/(2n+5) for i in 0:2n+5]
        a=assembly([node(1,(.25,.5,0.),1,1,[.25])])
        equal_outcomes("resource n$n single endpoint full retained-helper bits",()->edge_old(m,a,1),
            ()->edge_new(Dict{Int,Any}(),m,a,1);m,a)
        before=(model_digest(m),assembly_digest(a))
        for _ in 1:2
            @test last(cached_sweep(m,a,n))==1
        end
        samples=Int[]
        for _ in 1:2
            GC.gc();push!(samples,@allocated cached_sweep(m,a,n))
        end
        @test before==(model_digest(m),assembly_digest(a))
        push!(RESOURCE_ROWS,Dict("endpoint_calls"=>n,"source_curve_parameters"=>2n+6,
            "allocation_bytes"=>maximum(samples),"post_warm_samples"=>samples,
            "charge"=>"one actual cold operation-local frame plus n actual cached CAD Curve endpoint calls; model/assembly validation outside charge"))
    end
    for i in 2:length(RESOURCE_ROWS)
        @test RESOURCE_ROWS[i]["allocation_bytes"]<=2.15*RESOURCE_ROWS[i-1]["allocation_bytes"]+65536
    end
    # This supplemental limit does not replace the unchanged original public
    # set_order(2) n64/128/256 cap in cumulative scope api_native_producer.
end

function actual_method_controls!(destination)
    m=square();a=assembly([node(1,(.2,.3,0.),1,1,[.25])]);frames=Dict{Int,Any}()
    cell=A._Dim01Cell(2,(2,Int32(1)),Int32[1,1,1],UInt64(1),UInt64(1),Int32(0))
    quad=A._Dim01Cell(16,(2,Int32(1)),collect(Int32,1:8),UInt64(2),UInt64(2),Int32(0))
    primary_view=@view quad.nodes[1:4]
    signatures=(
        (:_dim01_native_plane_frame_reusable,Tuple{typeof(m),Int}),
        (:_dim01_native_plane_parametrization_cached!,Tuple{typeof(frames),typeof(m),Int,NTuple{3,Float64}}),
        (:_dim01_native_plane_surface_uv_cached!,Tuple{typeof(frames),typeof(m),typeof(a),Int32,Int32}),
        (:_dim01_native_plane_surface_face_uv_cached!,Tuple{typeof(frames),typeof(m),typeof(a),Int32,Int32}),
        (:_dim01_native_plane_support!,Tuple{typeof(frames),typeof(m),typeof(a),typeof(cell),Vector{Int32},Int,Tuple{Int,Int32},NTuple{3,Float64}}),
        (:_dim01_native_plane_support!,Tuple{typeof(frames),typeof(m),typeof(a),typeof(quad),typeof(primary_view),Int,Tuple{Int,Int32},NTuple{3,Float64}}))
    root=realpath(pkgdir(Tessella));expected=realpath(joinpath(root,"src","interfaces","APIGenerate01.jl"))
    @test bytes2hex(sha256(read(expected)))==SOURCE_SHA
    source=replace(read(expected,String),"\r\n"=>"\n")
    before_body=split(source,"function _dim01_native_plane_frame_reusable(m,surface)")
    @test length(before_body)==2
    after_body=split(last(before_body),"function _dim01_finish!(")
    @test length(after_body)==2
    protected_body="function _dim01_native_plane_frame_reusable(m,surface)"*first(after_body)
    @test bytes2hex(sha256(codeunits(protected_body)))==FRAME_BODY_SHA
    open(joinpath(destination,"actual_helper_ir.txt"),"w") do stream
        for (name,signature) in signatures
            f=getfield(A,name);method=which(f,signature)
            @test realpath(String(method.file))==expected
            lowered=Base.uncompressed_ast(method)
            @test box_count(lowered)==0
            typed=code_typed(f,signature;optimize=false)
            @test !isempty(typed)
            @test all(pair->box_count(first(pair))==0,typed)
            row=Dict("name"=>String(name),"actual_method"=>string(method),"file"=>expected,"line"=>method.line,
                "source_sha256"=>SOURCE_SHA,"historical_source_sha256"=>HISTORICAL_SOURCE_SHA,
                "frame_body_sha256"=>FRAME_BODY_SHA,"signature"=>string(signature),"lowered_boxes"=>box_count(lowered),
                "typed_boxes"=>[box_count(first(pair)) for pair in typed],
                "lowered_sha256"=>bytes2hex(sha256(sprint(show,lowered))),
                "typed_sha256"=>bytes2hex(sha256(sprint(show,typed))))
            push!(METHOD_ROWS,row);println(stream,row)
            show(stream,MIME("text/plain"),lowered);println(stream)
            show(stream,MIME("text/plain"),typed);println(stream)
        end
        positive=which(boxed_reader_positive,Tuple{Int});lowered=Base.uncompressed_ast(positive)
        @test box_count(lowered)>0 && boxed_reader_positive(7)()==7
        push!(METHOD_ROWS,Dict("name"=>"REAL_BOXED_READER_POSITIVE","actual_method"=>string(positive),
            "file"=>String(positive.file),"line"=>positive.line,"lowered_boxes"=>box_count(lowered),
            "source_sha256"=>bytes2hex(sha256(read(String(positive.file)))),
            "lowered_sha256"=>bytes2hex(sha256(sprint(show,lowered)))))
        println(stream,"REAL BOXED READER POSITIVE ",positive);show(stream,MIME("text/plain"),lowered)
    end
    @test length(METHOD_ROWS)==7
end

function run!()
    destination=abspath(get(ENV,"TESSELLA_PLANE_FRAME_CONTROL_OUTPUT",joinpath(pkgdir(Tessella),"test","tmp","plane_frame_helper_controls_$(VERSION)")))
    mkpath(destination);success=false
    try
        @testset "Operation-scoped Plane frame endpoint/FACE/validation/delegation controls" begin
            endpoint_controls!()
            cache_and_error_order_controls!()
            exclusion_controls!()
            support_controls!()
            actual_method_controls!(destination)
            resource_controls!()
            @test length(ROWS)==46
            @test all(row->row["model_unchanged"] && row["assembly_unchanged"],ROWS)
        end
        success=true
    finally
        open(joinpath(destination,"plane_frame_helper_observed.toml"),"w") do stream
            TOML.print(stream,Dict("julia"=>string(VERSION),"success"=>success,"candidate_source_sha256"=>SOURCE_SHA,
                "historical_source_sha256"=>HISTORICAL_SOURCE_SHA,"frame_body_sha256"=>FRAME_BODY_SHA,
                "branch_rows"=>ROWS,"actual_methods"=>METHOD_ROWS,"resource_rows"=>RESOURCE_ROWS,
                "reference_kind"=>"unchanged actual original helper methods plus unchanged literal support stencil; complete native public fixtures remain in original33 scopes",
                "original_public_resource_cap_preserved"=>true,
                "source_model_and_assembly_unchanged_in_all_recorded_branches"=>
                    length(ROWS)==46 && all(row->row["model_unchanged"] && row["assembly_unchanged"],ROWS)))
        end
    end
    return nothing
end
end
PlaneOperationFrameControlsV2.run!()
