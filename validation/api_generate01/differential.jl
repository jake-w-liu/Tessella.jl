#!/usr/bin/env julia
# Public generation0/1 state, typed connectivity, classification and source tags.
using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Test
using Tessella

binding=get(ENV,"GMSH_JULIA_API","")
isfile(binding) || error("set GMSH_JULIA_API to the pinned Gmsh4.15.2 gmsh.jl")
include(binding)
gmsh.GMSH_API_VERSION=="4.15.2" || error("Gmsh4.15.2 binding required")
const api=Tessella.API
const CASE_COUNT=Ref(0)
const STAGE_COUNT=Ref(0)
const ALLOW_DUPLICATE_CELLS=Ref(false)
const LEGACY_SOURCE_TAG_BLOCKERS=Ref(0)

# Gmsh 4.15.2's generated Julia wrapper drops its CFunction after installing
# the callback. Keep the trampoline rooted until Gmsh has unregistered it.
function with_oracle_size_callback(f,callback)
    thunk(dim,tag,x,y,z,lc,data)=callback(dim,tag,x,y,z,lc)
    handle=@cfunction($thunk,Cdouble,
        (Cint,Cint,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cvoid}))
    GC.@preserve handle begin
        ierr=Ref{Cint}(0)
        ccall((:gmshModelMeshSetSizeCallback,gmsh.lib),Cvoid,
            (Ptr{Cvoid},Ptr{Cvoid},Ptr{Cint}),handle,C_NULL,ierr)
        ierr[]==0 || error(gmsh.logger.getLastError())
        try
            return f()
        finally
            gmsh.model.mesh.removeSizeCallback()
        end
    end
end

function fresh(name;order=1,only_empty=false,renumber=true)
    api.initialize()
    gmsh.clear();gmsh.model.add(name)
    for (key,value) in (("Mesh.ElementOrder",order),
                        ("Mesh.MeshOnlyEmpty",only_empty ? 1 : 0),
                        ("Mesh.Renumber",renumber ? 1 : 0),
                        ("Mesh.MeshOnlyVisible",0))
        api.option(key,value)
        gmsh.option.setNumber(key,value)
    end
    # The Gmsh option store survives clear(). Reset field-related controls.
    for (key,value) in (("Mesh.MeshSizeMin",0.),("Mesh.MeshSizeMax",1e22),
                        ("Mesh.MeshSizeFactor",1.),("Mesh.MeshSizeFromPoints",1),
                        ("Mesh.MeshSizeFromCurvature",0.),("Mesh.FlexibleTransfinite",0))
        api.option(key,value)
        gmsh.option.setNumber(key,value)
    end
    gmsh.option.setNumber("Mesh.SecondOrderLinear",0)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete",0)
    gmsh.model.mesh.removeSizeCallback()
end

function option(key,value)
    api.option(key,value);gmsh.option.setNumber(key,value)
end

function point(tag,x,y=0.,z=0.;size=0.2)
    api.model.add_point(x,y,z;tag,meshSize=size)
    gmsh.model.geo.addPoint(x,y,z,size,tag)
end

function line(tag,a,b;count=3,law="Progression",coefficient=1.)
    api.model.add_line(a,b;tag);gmsh.model.geo.addLine(a,b,tag)
    gmsh.model.geo.synchronize()
    if count!==nothing
        api.mesh.set_transfinite_curve(tag,count,law,coefficient)
        gmsh.model.mesh.setTransfiniteCurve(tag,count,law,coefficient)
    end
end

function discrete(dim,tag,boundary=Int[])
    api.model.add_discrete_entity(dim,tag,boundary)
    gmsh.model.addDiscreteEntity(dim,tag,boundary)
end

function nodes(dim,tag,tags,coordinates)
    api.mesh.add_nodes(dim,tag,tags,coordinates)
    gmsh.model.mesh.addNodes(dim,tag,tags,coordinates)
end

function elements(tag,msh,tags,connectivity)
    api.mesh.add_elements_by_type(tag,msh,tags,connectivity)
    gmsh.model.mesh.addElementsByType(tag,msh,tags,connectivity)
end

function raw_line(;native=true,quadratic=false)
    if native
        point(1,0.);point(2,1.);line(1,1,2)
    else
        discrete(0,1);discrete(0,2);discrete(1,1,[1,2])
    end
    nodes(0,1,[11],[0.,0,0]);nodes(0,2,[22],[1.,0,0])
    nodes(1,1,[33],[0.5,quadratic ? 0.2 : 0.,0])
    elements(1,15,[101],[11]);elements(2,15,[102],[22])
    quadratic ? elements(1,8,[201],[11,22,33]) :
        elements(1,1,[201,202],[11,33,33,22])
end

function native_cube(;point_attachment=false)
    api.model.add_box(0,0,0,1,1,1;tag=1)
    gmsh.model.occ.addBox(0,0,0,1,1,1,1);gmsh.model.occ.synchronize()
    if point_attachment
        nodes(0,1,[111],gmsh.model.getValue(0,1,Float64[]))
        elements(1,15,[777],[111])
    end
    for (_,tag) in api.model.get_entities(1)
        api.mesh.set_transfinite_curve(tag,2)
        gmsh.model.mesh.setTransfiniteCurve(tag,2)
    end
    for (_,tag) in api.model.get_entities(2)
        api.mesh.set_transfinite_surface(tag);api.mesh.set_recombine(2,tag)
        gmsh.model.mesh.setTransfiniteSurface(tag);gmsh.model.mesh.setRecombine(2,tag)
    end
    api.mesh.set_transfinite_volume(1);gmsh.model.mesh.setTransfiniteVolume(1)
    api.mesh.generate(3);gmsh.model.mesh.generate(3)
end

function snapshot(native)
    mesh=native ? api.mesh : gmsh.model.mesh
    node_query=native ? mesh.get_nodes : mesh.getNodes
    element_query=native ? mesh.get_elements : mesh.getElements
    entities=native ? api.model.get_entities() : gmsh.model.getEntities()
    tags,coords,param=node_query()
    types,element_tags,connectivity=element_query()
    owners=[(entity=Tuple(entity),nodes=node_query(entity...),
             elements=element_query(entity...)) for entity in entities]
    return (;tags,coords,param,types,element_tags,connectivity,owners)
end

function compare_stage(name;raw=false,compare_elements=true,coordinate_tolerance=5e-12)
    actual=snapshot(true);expected=snapshot(false)
    @testset "$name" begin
        @test actual.tags==expected.tags
        @test length(actual.coords)==length(expected.coords)
        if length(actual.coords)==length(expected.coords)
            @test raw ? actual.coords==expected.coords :
                all(isapprox.(actual.coords,expected.coords;atol=coordinate_tolerance,rtol=5e-12))
        end
        @test actual.param==expected.param
        @test actual.types==expected.types
        @test actual.element_tags==expected.element_tags
        @test actual.connectivity==expected.connectivity
        @test length(actual.owners)==length(expected.owners)
        if length(actual.owners)==length(expected.owners)
            for (left,right) in zip(actual.owners,expected.owners)
                @test left.entity==right.entity
                @test left.nodes[1]==right.nodes[1]
                @test length(left.nodes[2])==length(right.nodes[2])
                if length(left.nodes[2])==length(right.nodes[2])
                    @test raw ? left.nodes[2]==right.nodes[2] :
                        all(isapprox.(left.nodes[2],right.nodes[2];atol=coordinate_tolerance,rtol=5e-12))
                end
                @test length(left.nodes[3])==length(right.nodes[3])
                if length(left.nodes[3])==length(right.nodes[3])
                    @test all(isapprox.(left.nodes[3],right.nodes[3];atol=coordinate_tolerance,rtol=5e-12))
                end
                @test left.elements==right.elements
            end
        end
        @test api.mesh.get_max_node_tag()==gmsh.model.mesh.getMaxNodeTag()
        @test api.mesh.get_max_element_tag()==gmsh.model.mesh.getMaxElementTag()
        if compare_elements
            for tags in expected.element_tags,tag in tags
                @test api.mesh.get_element(tag)==gmsh.model.mesh.getElement(tag)
            end
        end
        for tag in expected.tags
            left=api.mesh.get_node(tag);right=gmsh.model.mesh.getNode(tag)
            @test left[3:4]==right[3:4]
            @test all(isapprox.(left[1],right[1];atol=coordinate_tolerance,rtol=5e-12))
            @test length(left[2])==length(right[2])
            if length(left[2])==length(right[2])
                @test all(isapprox.(left[2],right[2];atol=coordinate_tolerance,rtol=5e-12))
            end
        end
        @test allunique(actual.tags)
        @test Set(vcat(actual.connectivity...))⊆Set(actual.tags)
        if api.LAST_MESH[]===nothing
            @test isempty(actual.tags) && isempty(actual.types)
        else
            @test validate(api.mesh.get();reject_duplicate_cells=!ALLOW_DUPLICATE_CELLS[]).ok
        end
    end
    STAGE_COUNT[]+=1
    return actual,expected
end

function generate(dim,name;raw=false,coordinate_tolerance=5e-12)
    api.mesh.generate(dim);gmsh.model.mesh.generate(dim)
    compare_stage(name;raw,coordinate_tolerance)
end

function fixture(f,name;allow_duplicate_cells=false,kwargs...)
    ALLOW_DUPLICATE_CELLS[]=allow_duplicate_cells
    fresh(name;kwargs...)
    try
        @testset "$name" begin
            f(name)
        end
        CASE_COUNT[]+=1
    finally
        api.finalize()
    end
end

gmsh.initialize(String[],false)
try
    gmsh.option.setNumber("General.Terminal",0)
    gmsh.option.getString("General.Version")=="4.15.2" ||
        error("loaded Gmsh4.15.2 library required")
    @testset "API generation0/1 pinned differential" begin
        fixture("fresh_zero_then_two_independent_curves") do name
            for (tag,x) in ((1,0.),(2,1.),(3,0.),(4,1.),(5,0.5))
                point(tag,x)
            end
            line(1,1,2);line(2,3,4)
            # geo.synchronize restores GEO attributes on existing entities.
            for curve in (1,2);gmsh.model.mesh.setTransfiniteCurve(curve,3);end
            api.model.add_physical_group(1,[1];tag=7)
            gmsh.model.addPhysicalGroup(1,[1],7)
            api.model.add_physical_group(0,[5];tag=8)
            gmsh.model.addPhysicalGroup(0,[5],8)
            generate(0,name*"_zero")
            generate(1,name*"_one")
            for (dim,tag) in ((1,7),(0,8))
                left=api.mesh.get_nodes_for_physical_group(dim,tag)
                right=gmsh.model.mesh.getNodesForPhysicalGroup(dim,tag)
                @test left[1]==right[1]
                @test all(isapprox.(left[2],right[2];atol=5e-12,rtol=5e-12))
            end
            generate(0,name*"_repeat_zero")
        end
        for native in (true,false),explicit in (:none,:first,:both),renumber in (true,false)
            name="point_native$(native)_cells$(explicit)_renumber$(renumber)"
            fixture(name;renumber) do name
                native ? (point(1,0.);gmsh.model.geo.synchronize()) : discrete(0,1)
                nodes(0,1,[101,205],[0.,0,0,2.,0,0])
                explicit===:first && elements(1,15,[501],[101])
                explicit===:both && elements(1,15,[501,502],[101,205])
                generate(0,name*"_zero";raw=true)
                generate(1,name*"_one";raw=true)
            end
        end
        for dim in (0,1),renumber in (true,false)
            fixture("coincident_point_dim$(dim)_renumber$(renumber)";renumber) do name
                discrete(0,1);nodes(0,1,[901,902],zeros(6))
                generate(dim,name;raw=true)
            end
        end
        for only_empty in (false,true),renumber in (false,true)
            fixture("coincident_chains_empty$(only_empty)_renumber$(renumber)";
                    only_empty,renumber) do name
                discrete(1,1)
                nodes(1,1,[11,12,13,21,22,23],
                      [0.,0,0,.5,0,0,1,0,0,0,0,0,.5,0,0,1,0,0])
                elements(1,1,[101,102,103,104],[11,12,12,13,21,22,22,23])
                generate(1,name;raw=true)
            end
        end
        for only_empty in (false,true),renumber in (false,true),order in (1,2)
            fixture("raw_linear_empty$(only_empty)_renumber$(renumber)_order$(order)";
                    only_empty,renumber,order) do name
                raw_line();generate(1,name)
            end
        end
        for native in (false,true),dim in (0,1),order in (1,2)
            fixture("raw_P2_native$(native)_dim$(dim)_order$(order)";
                    only_empty=true,renumber=false,order) do name
                raw_line(;native,quadratic=true);generate(dim,name;raw=!native)
            end
        end
        for dim in (0,1),only_empty in (false,true),order in (1,2)
            fixture("discrete_all_dim$(dim)_empty$(only_empty)_order$(order)";
                    only_empty,renumber=false,order) do name
                for (dim,tag) in ((1,11),(2,22),(3,33));discrete(dim,tag);end
                nodes(3,33,[101,205,309,407],[0.,0,0,1,0,0,0,1,0,0,0,1])
                elements(11,1,[501],[101,205]);elements(22,2,[601],[101,205,309])
                elements(33,4,[701],[101,205,309,407])
                generate(dim,name;raw=true)
            end
        end
        fixture("cross_record_missing_node_atomicity";renumber=false,only_empty=true) do name
            discrete(0,1);discrete(0,2);discrete(1,1,[1,2])
            nodes(0,1,[11],[0.,0,0]);nodes(0,2,[22],[1.,0,0])
            nodes(1,1,[33],[.5,0,0]);elements(1,15,[101],[11])
            before=snapshot(true)
            @test_throws r"node 999" api.mesh.add_elements_by_type(1,1,[201],[11,999])
            @test snapshot(true)==before
            # Gmsh reports an error rather than partially installing a bad cell.
            @test_throws Exception gmsh.model.mesh.addElementsByType(1,1,[201],[11,999])
            elements(1,1,[201,202],[11,33,33,22]);generate(1,name;raw=true)
        end
        for only_empty in (false,true),order in (1,2)
            fixture("native_cube_3_to_1_empty$(only_empty)_order$(order)") do name
                native_cube()
                for (_,tag) in api.model.get_entities(1)
                    api.mesh.set_transfinite_curve(tag,3)
                    gmsh.model.mesh.setTransfiniteCurve(tag,3)
                end
                option("Mesh.MeshOnlyEmpty",only_empty ? 1 : 0)
                option("Mesh.ElementOrder",order)
                generate(1,name)
            end
        end
        fixture("element_order_option_independent_of_set_order") do name
            point(1,0.);point(2,1.);line(1,1,2)
            api.mesh.set_order(2);gmsh.model.mesh.setOrder(2)
            @test api.option("Mesh.ElementOrder")==gmsh.option.getNumber("Mesh.ElementOrder")==1
            generate(1,name*"_fresh_set_order2")
            option("Mesh.ElementOrder",2);generate(0,name*"_option2")
            api.mesh.set_order(1);gmsh.model.mesh.setOrder(1)
            @test api.option("Mesh.ElementOrder")==gmsh.option.getNumber("Mesh.ElementOrder")==2
            generate(0,name*"_after_set_order1")
            option("Mesh.ElementOrder",1);generate(0,name*"_option1")
        end
        for (law,coefficient,count) in (("Progression",2.,6),("Bump",2.,7),("Beta",1.2,7))
            fixture("graded_$(law)") do name
                point(1,0.);point(2,1.);line(1,1,2;count,law,coefficient)
                # Native Lines also use the sampled adaptive density primitive.
                # Keep this historical 3e-7 differential bound; independent
                # saved-coordinate unit checks enforce 2e-11 on these laws.
                generate(1,name;coordinate_tolerance=3e-7)
            end
        end
        fixture("periodic_master_grading") do name
            for (tag,x,y) in ((1,0.,0.),(2,1.,0.),(3,0.,2.),(4,1.,2.));point(tag,x,y);end
            line(1,1,2;count=5,coefficient=2.);line(2,3,4;count=8)
            gmsh.model.mesh.setTransfiniteCurve(1,5,"Progression",2.)
            transform=[1.,0,0,0,0,1.,0,2.,0,0,1.,0,0,0,0,1.]
            api.mesh.set_periodic(1,[2],[1],transform)
            gmsh.model.mesh.setPeriodic(1,[2],[1],transform)
            generate(1,name;coordinate_tolerance=3e-7)
            left=api.mesh.get_periodic_nodes(1,2);right=gmsh.model.mesh.getPeriodicNodes(1,2)
            @test left[1]==right[1]
            @test Dict(zip(left[2],left[3]))==Dict(zip(right[2],right[3]))
            @test collect(left[4])==right[4]
            generate(0,name*"_zero";coordinate_tolerance=3e-7)
        end
        fixture("uniform_field_clip_then_factor") do name
            point(1,0.);point(2,1.);line(1,1,2;count=nothing)
            option("Mesh.MeshSizeFromPoints",0);option("Mesh.MeshSizeMin",0.2)
            option("Mesh.MeshSizeMax",0.2);option("Mesh.MeshSizeFactor",2)
            native_field=api.mesh.field.add("MathEval");oracle_field=gmsh.model.mesh.field.add("MathEval")
            api.mesh.field.set_string(native_field,"F","0.1")
            gmsh.model.mesh.field.setString(oracle_field,"F","0.1")
            api.mesh.field.set_as_background_mesh(native_field)
            gmsh.model.mesh.field.setAsBackgroundMesh(oracle_field)
            native_lc=Float64[];oracle_lc=Float64[]
            api.mesh.set_size_callback((dim,tag,x,y,z,lc)->(push!(native_lc,lc);lc))
            with_oracle_size_callback((dim,tag,x,y,z,lc)->(push!(oracle_lc,lc);lc)) do
                GC.gc()
                generate(1,name)
                @test !isempty(native_lc) && !isempty(oracle_lc)
                @test Set(native_lc)==Set(oracle_lc)==Set([0.1])
            end
            api.mesh.remove_size_callback()
            option("Mesh.CharacteristicLengthFactor",1);generate(1,name*"_factor1")
        end
        fixture("native_visibility") do name
            point(1,0.);point(2,1.);line(1,1,2)
            api.model.set_visibility([(1,1)],0);gmsh.model.setVisibility([(1,1)],0)
            option("Mesh.MeshOnlyVisible",1);generate(1,name*"_hidden")
            api.model.set_visibility([(1,1)],1);gmsh.model.setVisibility([(1,1)],1)
            generate(1,name*"_visible")
        end
        fixture("sparse_native_queries_and_mutation_refresh";renumber=false) do name
            discrete(1,1);nodes(1,1,[500,100],[0.,0,0,1,0,0])
            elements(1,1,[700],[500,100]);generate(1,name)
            @test api.mesh.get_duplicate_nodes()==gmsh.model.mesh.getDuplicateNodes()==UInt64[]
            @test api.mesh.get_nodes_by_element_type(1)==gmsh.model.mesh.getNodesByElementType(1)
            @test api.mesh.get_element_edge_nodes(1)==gmsh.model.mesh.getElementEdgeNodes(1)
            @test api.mesh.get_keys(1,"Lagrange")==gmsh.model.mesh.getKeys(1,"Lagrange")
            @test api.mesh.get_keys_for_element(700,"Lagrange")==gmsh.model.mesh.getKeysForElement(700,"Lagrange")
            @test api.mesh.get_basis_functions_orientation(1,"H1Legendre1")==
                gmsh.model.mesh.getBasisFunctionsOrientation(1,"H1Legendre1")
            @test api.mesh.get_basis_functions_orientation_for_element(700,"H1Legendre1")==
                gmsh.model.mesh.getBasisFunctionsOrientationForElement(700,"H1Legendre1")
            @test api.mesh.get_element_by_coordinates(.5,0,0,1,true)==
                gmsh.model.mesh.getElementByCoordinates(.5,0,0,1,true)
            @test api.mesh.get_elements_by_coordinates(.5,0,0,1,true)==
                gmsh.model.mesh.getElementsByCoordinates(.5,0,0,1,true)
            @test api.mesh.get_local_coordinates_in_element(700,.25,0,0)==
                gmsh.model.mesh.getLocalCoordinatesInElement(700,.25,0,0)
            @test api.mesh.get_element_qualities([700],"minEdge")==
                gmsh.model.mesh.getElementQualities([700],"minEdge")
            @test api.mesh.get_jacobian(700,[0.,0,0])==gmsh.model.mesh.getJacobian(700,[0.,0,0])
            api.mesh.create_edges();gmsh.model.mesh.createEdges()
            @test api.mesh.get_all_edges()==gmsh.model.mesh.getAllEdges()
            @test api.mesh.get_edges([500,100,100,500])==gmsh.model.mesh.getEdges([500,100,100,500])
            api.mesh.reverse_elements([700]);gmsh.model.mesh.reverseElements([700])
            compare_stage(name*"_reverse";raw=true)
            generate(0,name*"_reverse_then_zero";raw=true)
            api.mesh.set_node(100,[2.,0,0]);gmsh.model.mesh.setNode(100,[2.,0,0],Float64[])
            compare_stage(name*"_set_node";raw=true)
            @test api.mesh.get_element_qualities([700],"minEdge")==
                gmsh.model.mesh.getElementQualities([700],"minEdge")
            generate(0,name*"_set_node_then_zero";raw=true)
            api.mesh.remove_elements(1,1,[700]);gmsh.model.mesh.removeElements(1,1,[700])
            compare_stage(name*"_remove";raw=true)
            api.mesh.clear();gmsh.model.mesh.clear()
            compare_stage(name*"_clear";raw=true)
            generate(0,name*"_clear_then_zero";raw=true)
        end
        for quadratic in (false,true),renumber in (false,true)
            fixture("sparse_line_refine_P2$(quadratic)_renumber$(renumber)";renumber=false) do name
                discrete(1,1)
                source_tags=quadratic ? [500,100,200] : [500,100]
                source_coords=quadratic ? [0.,0,0,1,0,0,.5,.2,0] : [0.,0,0,1,0,0]
                nodes(1,1,source_tags,source_coords)
                elements(1,quadratic ? 8 : 1,[700],source_tags)
                # refine itself applies the numbering option, so start from
                # the authoritative raw cells with their sparse source tags.
                api.mesh.generate(0)
                gmsh.model.mesh.generate(0)
                option("Mesh.Renumber",renumber ? 1 : 0)
                api.mesh.refine();gmsh.model.mesh.refine()
                compare_stage(name;raw=true)
            end
        end
        for reverse_tags in (false,true)
            fixture("duplicate_source_node_keeper$(reverse_tags)";renumber=false,allow_duplicate_cells=true) do name
                discrete(1,1)
                raw_tags=reverse_tags ? [100,500,501] : [500,100,501]
                nodes(1,1,raw_tags,[0.,0,0,0,0,0,1,0,0])
                elements(1,1,[701,702],[500,501,100,501]);generate(1,name)
                @test api.mesh.get_duplicate_nodes()==gmsh.model.mesh.getDuplicateNodes()
                api.mesh.remove_duplicate_nodes();gmsh.model.mesh.removeDuplicateNodes()
                compare_stage(name*"_removed";raw=true)
            end
            fixture("duplicate_source_element_keeper$(reverse_tags)";renumber=false,allow_duplicate_cells=true) do name
                discrete(1,1);nodes(1,1,[500,100],[0.,0,0,1,0,0])
                raw_tags=reverse_tags ? [300,700] : [700,300]
                elements(1,1,raw_tags,[500,100,500,100]);generate(1,name)
                api.mesh.remove_duplicate_elements();gmsh.model.mesh.removeDuplicateElements()
                compare_stage(name*"_removed";raw=true)
                @test validate(api.mesh.get()).ok
            end
        end
        for reverse_tags in (false,true)
            fixture("duplicate_cross_entity_keeper$(reverse_tags)";renumber=false) do name
                discrete(1,2);discrete(1,1)
                first,second=reverse_tags ? (100,500) : (500,100)
                nodes(1,1,[first,first+1],[0.,0,0,1,0,0])
                nodes(1,2,[second,second+1],[0.,0,0,2,0,0])
                elements(1,1,[701],[first,first+1]);elements(2,1,[702],[second,second+1])
                generate(1,name)
                api.mesh.remove_duplicate_nodes();gmsh.model.mesh.removeDuplicateNodes()
                compare_stage(name*"_removed";raw=true)
            end
        end
        renumber_cases=((:unknown,[999],[800]),(:unknown_and_matched,[999,500],[800,10]),
                        (:empty,Int[],Int[]),(:partial,[500],[800]),
                        (:collision_unmapped,[500],[100]),(:swap,[500,100],[100,500]),
                        (:duplicate_old,[500,500],[10,20]),
                        (:duplicate_old_high_previous,[500,500],[800,20]))
        for kind in (:node,:element),(label,old,new) in renumber_cases
            fixture("sparse_renumber_$(kind)_$(label)";renumber=false,allow_duplicate_cells=true) do name
                discrete(1,1);nodes(1,1,[500,100],[0.,0,0,1,0,0])
                elements(1,1,[500,100],[500,100,100,500]);generate(1,name)
                if kind===:node
                    api.mesh.renumber_nodes(old,new);gmsh.model.mesh.renumberNodes(old,new)
                else
                    api.mesh.renumber_elements(old,new);gmsh.model.mesh.renumberElements(old,new)
                end
                compare_stage(name*"_renumbered";raw=true)
                generate(0,name*"_then_zero";raw=true)
            end
        end
        for quad_first in (false,true)
            fixture("mixed_discrete_surface_family_order$(quad_first)") do name
                discrete(2,1);nodes(2,1,[101,205,309,407],[0.,0,0,1,0,0,1,1,0,0,1,0])
                for msh in (quad_first ? (3,2) : (2,3))
                    elements(1,msh,[msh==2 ? 700 : 300],msh==2 ? [101,205,309] : [101,205,309,407])
                end
                generate(0,name;raw=true)
            end
        end
        fixture("auto_additions_and_model_allocator_switch";renumber=false) do name
            discrete(1,1);nodes(1,1,[500,100],[0.,0,0,1,0,0])
            elements(1,1,[700],[500,100]);generate(0,name*"_a")
            native_a=api.model.get_current();oracle_a=gmsh.model.getCurrent()
            api.model.add("allocator_b");gmsh.model.add("allocator_b")
            discrete(1,1);nodes(1,1,[10,20],[0.,1,0,1,1,0])
            elements(1,1,[30],[10,20]);generate(0,name*"_b")
            api.model.set_current(native_a);gmsh.model.setCurrent(oracle_a)
            compare_stage(name*"_a_restored";raw=true)
            api.mesh.clear();gmsh.model.mesh.clear()
            nodes(1,1,Int[],[0.,0,0,1,0,0])
            elements(1,1,Int[],[501,502]);generate(0,name*"_auto_after_clear")
            api.model.set_current("allocator_b");gmsh.model.setCurrent("allocator_b")
            compare_stage(name*"_b_restored";raw=true)
        end
        fixture("selective_curve_clear_preserves_orphan_point_cells";renumber=false) do name
            for point in 1:3;discrete(0,point);end
            for curve in 1:2;discrete(1,curve,[curve,curve+1]);end
            for (tag,node,x) in ((1,11,0.),(2,22,1.),(3,33,2.))
                nodes(0,tag,[node],[x,0,0]);elements(tag,15,[100+tag],[node])
            end
            nodes(1,1,[44],[.5,0,0]);nodes(1,2,[55],[1.5,0,0])
            elements(1,1,[201,202],[11,44,44,22]);elements(2,1,[301,302],[22,55,55,33])
            generate(1,name)
            api.mesh.clear([(1,1)]);gmsh.model.mesh.clear([(1,1)])
            compare_stage(name*"_clear";raw=true)
            @test api.mesh.get_nodes(1,2,true)[1]==gmsh.model.mesh.getNodes(1,2,true)[1]
        end
        fixture("periodic_keys_full_quadratic_curve";order=2) do name
            for (tag,x,y) in ((1,0.,0.),(2,1.,0.),(3,0.,2.),(4,1.,2.));point(tag,x,y);end
            line(1,1,2);line(2,3,4)
            gmsh.model.mesh.setTransfiniteCurve(1,3)
            transform=[1.,0,0,0,0,1.,0,2.,0,0,1.,0,0,0,0,1.]
            api.mesh.set_periodic(1,[2],[1],transform)
            gmsh.model.mesh.setPeriodic(1,[2],[1],transform)
            generate(1,name)
            native=api.mesh.get_periodic_keys(8,"Lagrange",2)
            oracle=gmsh.model.mesh.getPeriodicKeys(8,"Lagrange",2)
            @test native[1:5]==oracle[1:5]
            @test all(isapprox.(native[6],oracle[6];atol=5e-12,rtol=5e-12))
            @test all(isapprox.(native[7],oracle[7];atol=5e-12,rtol=5e-12))
            no_coords=api.mesh.get_periodic_keys(8,"Lagrange",2,false)
            oracle_no_coords=gmsh.model.mesh.getPeriodicKeys(8,"Lagrange",2,false)
            @test no_coords==oracle_no_coords
            for include_high_order in (false,true)
                native_nodes=api.mesh.get_periodic_nodes(1,2,include_high_order)
                oracle_nodes=gmsh.model.mesh.getPeriodicNodes(1,2,include_high_order)
                @test native_nodes.master_entity==oracle_nodes[1]
                # Upstream iteration is over maps keyed by MVertex pointers;
                # compare the bijection independently of pointer allocation order.
                @test Dict(zip(native_nodes.slave_nodes,native_nodes.master_nodes))==
                    Dict(zip(oracle_nodes[2],oracle_nodes[3]))
                @test collect(native_nodes.affine)==oracle_nodes[4]
            end
        end
        fixture("sparse_visibility_reorder_and_homology";renumber=false) do name
            discrete(1,1);nodes(1,1,[500,100,300],[0.,0,0,.5,0,0,1.,0,0])
            elements(1,1,[700,300],[500,100,100,300])
            generate(0,name;raw=true)
            @test api.mesh.get_visibility([700,300,1])==gmsh.model.mesh.getVisibility([700,300,1])
            api.mesh.set_visibility([700,1],5);gmsh.model.mesh.setVisibility([700,1],5)
            @test api.mesh.get_visibility([700,300,1])==gmsh.model.mesh.getVisibility([700,300,1])
            api.mesh.reorder_elements(1,1,[1,0]);gmsh.model.mesh.reorderElements(1,1,[1,0])
            compare_stage(name*"_reorder";raw=true)
            @test api.mesh.get_visibility([700,300])==gmsh.model.mesh.getVisibility([700,300])
            api.mesh.add_homology_request("Homology",Int[],Int[],[0])
            gmsh.model.mesh.addHomologyRequest("Homology",Int[],Int[],[0])
            native=api.mesh.compute_homology();oracle=gmsh.model.mesh.computeHomology()
            @test native==oracle
            @test length(api.mesh.get_elements_by_type(15)[1])==
                length(gmsh.model.mesh.getElementsByType(15)[1])==1
            @test Set(api.mesh.get_elements_by_type(15)[2])⊆Set(api.mesh.get_nodes()[1])
        end
        fixture("point_and_line_endpoint_locator_priority") do name
            point(1,0.);point(2,1.);line(1,1,2)
            generate(1,name)
            for x in (0.,1.),dimension in (-1,0,1),strict in (false,true)
                native=api.mesh.get_element_by_coordinates(x,0,0,dimension,strict)
                oracle=gmsh.model.mesh.getElementByCoordinates(x,0,0,dimension,strict)
                @test native[1:3]==oracle[1:3]
                @test collect(native[4:6])≈collect(oracle[4:6]) atol=1e-10 rtol=1e-10
                @test api.mesh.get_elements_by_coordinates(x,0,0,dimension,strict)==
                    gmsh.model.mesh.getElementsByCoordinates(x,0,0,dimension,strict)
            end
        end
        for action in (:set_node,:reverse,:reverse_elements,:affine,:add_node,:add_element,
                       :renumber_node,:renumber_element,:remove_other_element,:remove_element)
            fixture("tagged_visibility_$(action)";renumber=false) do name
                discrete(1,1);nodes(1,1,[500,100,300],[0.,0,0,.5,0,0,1.,0,0])
                elements(1,1,[700,300],[500,100,100,300]);generate(0,name;raw=true)
                api.mesh.set_visibility([700],5);gmsh.model.mesh.setVisibility([700],5)
                tag=700
                if action===:set_node
                    api.mesh.set_node(100,[.6,0,0]);gmsh.model.mesh.setNode(100,[.6,0,0],Float64[])
                elseif action===:reverse
                    api.mesh.reverse([(1,1)]);gmsh.model.mesh.reverse([(1,1)])
                elseif action===:reverse_elements
                    api.mesh.reverse_elements([700]);gmsh.model.mesh.reverseElements([700])
                elseif action===:affine
                    transform=[1.,0,0,0,0,1.,0,1.,0,0,1.,0,0,0,0,1.]
                    api.mesh.affine_transform(transform);gmsh.model.mesh.affineTransform(transform)
                elseif action===:add_node
                    nodes(1,1,[800],[2.,0,0])
                elseif action===:add_element
                    elements(1,1,[900],[500,300])
                elseif action===:renumber_node
                    api.mesh.renumber_nodes([500],[800]);gmsh.model.mesh.renumberNodes([500],[800])
                elseif action===:renumber_element
                    api.mesh.renumber_elements([700],[800]);gmsh.model.mesh.renumberElements([700],[800]);tag=800
                elseif action===:remove_other_element
                    api.mesh.remove_elements(1,1,[300]);gmsh.model.mesh.removeElements(1,1,[300])
                else
                    api.mesh.remove_elements(1,1,[700]);gmsh.model.mesh.removeElements(1,1,[700])
                end
                compare_stage(name*"_after_edit";raw=true)
                @test api.mesh.get_visibility([tag])==gmsh.model.mesh.getVisibility([tag])
            end
        end
        fixture("native_quadratic_same_order_visibility";renumber=false,order=2) do name
            point(1,0.);point(2,1.);line(1,1,2)
            generate(1,name)
            tag=first(gmsh.model.mesh.getElementsByType(8)[1])
            api.mesh.set_visibility([tag],5);gmsh.model.mesh.setVisibility([tag],5)
            api.mesh.set_order(2);gmsh.model.mesh.setOrder(2)
            compare_stage(name*"_same_order")
            @test api.mesh.get_visibility([tag])==gmsh.model.mesh.getVisibility([tag])==Int32[5]
            api.mesh.set_order(1);gmsh.model.mesh.setOrder(1)
            compare_stage(name*"_changed_order")
            @test api.mesh.get_visibility([tag])==gmsh.model.mesh.getVisibility([tag])==Int32[0]
            current=gmsh.model.mesh.getElementsByType(1)[1]
            @test api.mesh.get_visibility(current)==gmsh.model.mesh.getVisibility(current)==Int32[1,1]
        end
        fixture("discrete_quadratic_generate_zero_visibility";renumber=false) do name
            raw_line(;native=false,quadratic=true);generate(0,name;raw=true)
            api.mesh.set_visibility([201],5);gmsh.model.mesh.setVisibility([201],5)
            generate(0,name*"_same_cells";raw=true)
            @test api.mesh.get_visibility([201])==gmsh.model.mesh.getVisibility([201])==Int32[5]
        end
        fixture("legacy_higher_cache_sparse_attachment_blocker";renumber=false) do name
            native_cube(;point_attachment=true)
            for (_,tag) in api.model.get_entities(1)
                api.mesh.set_transfinite_curve(tag,3);gmsh.model.mesh.setTransfiniteCurve(tag,3)
            end
            model=api.CURRENT[];cache=api.LAST_MESH[];class=api.LAST_MESH_CLASS[]
            before=snapshot(true);parameters=deepcopy(model.curve_params)
            counters=(api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())
            @test_throws r"untracked source tag allocations" api.mesh.generate(1)
            @test api.CURRENT[]===model && api.LAST_MESH[]===cache && api.LAST_MESH_CLASS[]===class
            @test snapshot(true)==before
            @test model.curve_params==parameters
            @test (api.mesh.get_max_node_tag(),api.mesh.get_max_element_tag())==counters
            gmsh.model.mesh.generate(1)
            @test gmsh.model.mesh.getNodes()[1]==UInt64.(111:130)
            @test (gmsh.model.mesh.getMaxNodeTag(),gmsh.model.mesh.getMaxElementTag())==(130,846)
            @test gmsh.model.mesh.getElementsByType(15,1)==(UInt64[777],UInt64[111])
            LEGACY_SOURCE_TAG_BLOCKERS[]+=1
        end
        for quadratic in (false,true)
            fixture("raw_only_immediate_order_$(quadratic ? "lower" : "elevate")";renumber=false) do name
                raw_line(;native=false,quadratic)
                @test api.LAST_MESH[]===nothing
                order=quadratic ? 1 : 2
                api.mesh.set_order(order);gmsh.model.mesh.setOrder(order)
                compare_stage(name;raw=true)
                @test api.option("Mesh.ElementOrder")==gmsh.option.getNumber("Mesh.ElementOrder")==1
            end
        end
        fixture("periodic_p2_compaction_and_half_clear") do name
            for (tag,x,y) in ((1,0.,0.),(2,1.,0.),(3,0.,2.),(4,1.,2.));point(tag,x,y);end
            line(1,1,2);line(2,3,4);gmsh.model.mesh.setTransfiniteCurve(1,3)
            affine=[1.,0,0,0,0,1.,0,2.,0,0,1.,0,0,0,0,1.]
            api.mesh.set_periodic(1,[2],[1],affine);gmsh.model.mesh.setPeriodic(1,[2],[1],affine)
            generate(1,name)
            for order in (2,1)
                api.mesh.set_order(order);gmsh.model.mesh.setOrder(order)
                compare_stage(name*"_order$(order)")
                native=api.mesh.get_periodic_nodes(1,2,true)
                oracle=gmsh.model.mesh.getPeriodicNodes(1,2,true)
                @test native.master_entity==oracle[1]
                @test Dict(zip(native.slave_nodes,native.master_nodes))==Dict(zip(oracle[2],oracle[3]))
                @test collect(native.affine)==oracle[4]
            end
            before=snapshot(true);model=api.CURRENT[];cache=api.LAST_MESH[]
            @test_throws ArgumentError api.mesh.clear([(1,2),(1,999)])
            @test_throws Exception gmsh.model.mesh.clear([(1,2),(1,999)])
            @test snapshot(true)==before && api.CURRENT[]===model && api.LAST_MESH[]===cache
            compare_stage(name*"_invalid_scope_unchanged")
            api.mesh.clear([(1,2)]);gmsh.model.mesh.clear([(1,2)])
            compare_stage(name*"_half_clear")
            native=api.mesh.get_periodic_nodes(1,2,true)
            oracle=gmsh.model.mesh.getPeriodicNodes(1,2,true)
            @test Dict(zip(native.slave_nodes,native.master_nodes))==Dict(zip(oracle[2],oracle[3]))==Dict()
            @test all(link->isempty(link.slave_nodes) && isempty(link.master_nodes),api.mesh.get().periodic_links)
        end
        fixture("hidden_unmeshed_circle_straight_line_refine") do name
            source="Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={2,0,0,1};" *
                "Point(4)={2,1,0,1};Point(5)={3,1,0,1};Line(1)={1,2};" *
                "Circle(2)={3,4,5};Transfinite Curve{1}=3;"
            mktempdir() do directory
                path=joinpath(directory,"hidden_refine.geo")
                write(path,source);api.open_geo!(path)
                for (tag,x,y) in ((1,0.,0.),(2,1.,0.),(3,2.,0.),(4,2.,1.),(5,3.,1.))
                    gmsh.model.geo.addPoint(x,y,0,1,tag)
                end
                gmsh.model.geo.addLine(1,2,1);gmsh.model.geo.addCircleArc(3,4,5,2)
                gmsh.model.geo.synchronize();gmsh.model.mesh.setTransfiniteCurve(1,3)
                api.model.set_visibility([(1,2)],0);gmsh.model.setVisibility([(1,2)],0)
                option("Mesh.MeshOnlyVisible",1)
                generate(1,name)
                api.mesh.refine();gmsh.model.mesh.refine()
                compare_stage(name*"_refine")
            end
        end
    end
    println("API_GENERATE01_DIFFERENTIAL_OK gmsh=4.15.2 cases=$(CASE_COUNT[]) stages=$(STAGE_COUNT[]) legacy_source_tag_blockers=$(LEGACY_SOURCE_TAG_BLOCKERS[])")
finally
    api.finalize();gmsh.finalize()
end
