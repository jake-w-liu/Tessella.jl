using Test,Tessella,SHA
const _MBO=Tessella.API
include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_orientation_public_v1.jl"))

function _mbo_cube(;occ=false,extra=true)
    _MBO.initialize();_MBO.option("Mesh.Renumber",0)
    if occ
        _MBO.model.add_box(0,0,0,1,1,1;tag=1)
    else
        for (tag,point) in enumerate(((0,0,0),(1,0,0),(1,1,0),(0,1,0),
                                      (0,0,1),(1,0,1),(1,1,1),(0,1,1)))
            _MBO.model.add_point(point...;tag)
        end
        for (tag,edge) in enumerate(((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
                                     (7,8),(8,5),(1,5),(2,6),(3,7),(4,8)))
            _MBO.model.add_line(edge...;tag)
        end
        loops=([1,2,3,4],[5,6,7,8],[1,10,-5,-9],
               [2,11,-6,-10],[3,12,-7,-11],[4,9,-8,-12])
        for (tag,curves) in enumerate(loops)
            _MBO.model.add_curve_loop(curves;tag)
            _MBO.model.add_plane_surface([tag];tag)
        end
        _MBO.model.add_surface_loop([-1,2,3,4,5,6];tag=1)
        _MBO.model.add_volume([1];tag=1)
    end
    if extra
        first=occ ? 9 : 21;curve=occ ? 13 : 21
        for (offset,point) in enumerate(((3,0,0),(4,0,0),(4,1,0),(3,1,0)))
            _MBO.model.add_point(point...;tag=first+offset-1)
        end
        for offset in 0:3
            _MBO.model.add_line(first+offset,first+mod(offset+1,4);tag=curve+offset)
        end
        _MBO.model.add_curve_loop(collect(curve:curve+3);tag=7)
        _MBO.model.add_plane_surface([7];tag=7)
    end
end

function _mbo_load(case;quad_first=false)
    factory=case["factory"]
    if factory=="discrete"
        _MBO.initialize();_MBO.option("Mesh.Renumber",0)
        for surface in 1:6;_MBO.model.add_discrete_entity(2,surface);end
        _MBO.model.add_discrete_entity(3,1,[1,2,3,4,5,6])
    else
        _mbo_cube(;occ=factory=="OpenCASCADE")
        if factory=="periodic"
            _MBO.mesh.set_periodic(2,[2],[1],Float64[1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1])
        end
        for surface in 1:7
            reverse=surface==7 || factory=="periodic" && surface>1 ||
                get(case,"pattern","")=="all_inward" &&
                    (factory=="OpenCASCADE" || surface>1)
            reverse && _MBO.mesh.set_reverse(2,surface,true)
        end
    end
    by_owner=Dict{Tuple{Int,Int},Vector{Any}}()
    for node in case["nodes"]
        push!(get!(()->Any[],by_owner,(node["dim"],node["owner"])),node)
    end
    for ((dimension,tag),nodes) in sort!(collect(by_owner);by=first)
        tags=[node["tag"] for node in nodes]
        coordinates=vcat((node["coords"] for node in nodes)...)
        parameters=vcat((node["params"] for node in nodes)...)
        _MBO.mesh.add_nodes(dimension,tag,tags,coordinates,parameters)
    end
    quad_first && _MBO.mesh.add_elements_by_type(1,3,[9001],[1,4,3,2])
    cells=vcat(case["before"],get(case,"volume_cells",[]))
    for cell in cells
        _MBO.mesh.add_elements_by_type(cell["owner"],cell["type"],
            [cell["tag"]],cell["nodes"])
    end
end

function _mbo_assert_cells(cells)
    for cell in cells
        found=_MBO.mesh.get_element(cell["tag"])
        @test found[1]==cell["type"]
        @test found[2]==cell["nodes"]
        @test found[4]==cell["owner"]
    end
end

@testset "Actual public retained Surface orientation" begin
    for (index,case) in enumerate(ORIENTATION_PUBLIC_FIXTURES)
        @testset "primary case $index $(case["factory"]) order$(case["order"])" begin
            try
                _mbo_load(case)
                _mbo_assert_cells(case["before"])
                nodes=Dict(node["tag"]=>_MBO.mesh.get_node(node["tag"]) for node in case["nodes"])
                maximum=(_MBO.NODE_TAG_MAX[],_MBO.ELEMENT_TAG_MAX[])
                periodic=deepcopy(_MBO.CURRENT[].periodic)
                labels=[cell["tag"] for cell in case["before"]]
                _MBO.mesh.set_visibility(labels,0)
                visibility=copy(_MBO.ELEMENT_VISIBILITY[])
                _MBO.mesh.set_outward_orientation(1)
                _mbo_assert_cells(case["after"])
                @test Dict(tag=>_MBO.mesh.get_node(tag) for tag in keys(nodes))==nodes
                @test (_MBO.NODE_TAG_MAX[],_MBO.ELEMENT_TAG_MAX[])==maximum
                @test _MBO.CURRENT[].periodic==periodic
                @test _MBO.ELEMENT_VISIBILITY[]==visibility
                for cell in get(case,"volume_cells",[])
                    @test _MBO.mesh.get_element(cell["tag"])[2]==cell["nodes"]
                end
                _MBO.mesh.set_outward_orientation(1)
                _mbo_assert_cells(get(case,"second",case["after"]))
                @test Dict(tag=>_MBO.mesh.get_node(tag) for tag in keys(nodes))==nodes
            finally
                _MBO.finalize()
            end
        end
    end
end

@testset "Native Box face sense and parameter frame" begin
    try
        _MBO.initialize();_MBO.model.add_box(2,-3,5,2,3,4;tag=1)
        normals=([-1.,0,0],[1.,0,0],[0.,-1,0],[0.,1,0],[0.,0,-1],[0.,0,1])
        points=([2.,-2.25,8.],[4.,-2.25,8.],[3.5,-3.,6.],[3.5,0.,6.],[2.5,-.75,5.],[2.5,-.75,9.])
        parameters=([3.,-.75],[3.,-.75],[1.,1.5],[1.,1.5],[.5,2.25],[.5,2.25])
        for surface in 1:6
            @test _MBO.model.get_normal(surface,parameters[surface])==normals[surface]
            @test _MBO.model.get_parametrization(2,surface,points[surface])≈parameters[surface] atol=2e-11
            @test _MBO.model.get_value(2,surface,parameters[surface])≈points[surface] atol=2e-11
        end
        _MBO.mesh.set_outward_orientation(1)
        @test all(!get(_MBO.CURRENT[].meshing.reverse,(2,s),false) for s in 1:6)
        @test length(_MBO.CURRENT[].meshing.orientation_stl)==6
        @test isempty(_MBO.CURRENT[].meshing.outward_orientation)
    finally
        _MBO.finalize()
    end
end

@testset "Edited native P2 support row preservation" begin
    try
        _mbo_cube(;extra=false)
        for curve in 1:12;_MBO.mesh.set_transfinite_curve(curve,2);end
        for surface in 1:6;_MBO.mesh.set_transfinite_surface(surface);end
        _MBO.mesh.generate(2);_MBO.mesh.set_order(2);_MBO.mesh.renumber_nodes()
        @test _MBO.LAST_MESH[] isa Tessella.MixedMesh
        @test _MBO.LAST_MESH_HIGH_ORDER[]===nothing
        labels=_MBO.mesh.get_elements(2,1)[2][1]
        element=_MBO.mesh.get_element(first(labels))
        node=element[2][5]
        old=_MBO.mesh.get_node(node)
        coordinates=copy(old[1]);coordinates[3]=.07
        _MBO.mesh.set_node(node,coordinates,old[2])
        before=_MBO.LAST_MESH[]
        cells=Dict(label=>_MBO.mesh.get_element(label) for labels in
            _MBO.mesh.get_elements(2,-1)[2] for label in labels)
        points=copy(before.coords)
        node_labels=copy(_MBO.LAST_MESH_CLASS[].public_tags.node_tags)
        parameters=deepcopy(_MBO.LAST_MESH_CLASS[].node_parametric)
        _MBO.mesh.set_outward_orientation(1)
        after=_MBO.LAST_MESH[]
        @test after.coords==points
        @test _MBO.LAST_MESH_CLASS[].public_tags.node_tags==node_labels
        @test _MBO.LAST_MESH_CLASS[].node_parametric==parameters
        @test _MBO.mesh.get_node(node)[1]==coordinates
        for (label,cell) in cells
            actual=_MBO.mesh.get_element(label)
            @test actual[1]==cell[1]==9
            @test actual[3:4]==cell[3:4]
            @test actual[2]==(cell[4]==1 ? cell[2][[1,3,2,6,5,4]] : cell[2])
        end
    finally
        _MBO.finalize()
    end
end

@testset "Actual same-Surface Triangle-before-Quadrangle storage" begin
    case=first(ORIENTATION_PUBLIC_FIXTURES)
    for quad_first in (true,false)
        try
            _mbo_load(case;quad_first)
            quad_first || _MBO.mesh.add_elements_by_type(1,3,[9001],[1,4,3,2])
            @test _MBO.mesh.get_element(9001)[2]==[1,4,3,2]
            nodes=Dict(node["tag"]=>_MBO.mesh.get_node(node["tag"]) for node in case["nodes"])
            _MBO.mesh.set_outward_orientation(1)
            _mbo_assert_cells(case["after"])
            @test _MBO.mesh.get_element(9001)[2]==[1,2,3,4]
            @test Dict(tag=>_MBO.mesh.get_node(tag) for tag in keys(nodes))==nodes
        finally
            _MBO.finalize()
        end
    end
end

@testset "Net-zero orientation retains explicit carrier identity" begin
    carrier=Tessella.Mesh([0. 1 0;0 0 1;0 0 0];tris=reshape(Int32[1,2,3],3,1))
    actual=Tessella.simplex_to_mixed(carrier)
    changes=Dict(1=>Int32[1,2,3])
    replacement=Tessella.Model._model_orientation_mixed_replacement(actual,changes)
    @test replacement===actual
    @test Tessella.Model._model_orientation_carrier(actual,replacement,carrier,Int32[1,2,3],changes)===carrier
    reversed=Tessella.Model._model_orientation_mixed_replacement(actual,Dict(1=>Int32[1,3,2]))
    rebound=Tessella.Model._model_orientation_carrier(actual,reversed,carrier,Int32[1,2,3],Dict(1=>Int32[1,3,2]))
    @test rebound.tris==reshape(Int32[1,3,2],3,1)
    @test rebound.coords==carrier.coords
end

include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_orientation_box_edges_v1.jl"))
@testset "Native Box actual length domains and face PCurves" begin
    for row in ORIENTATION_BOX_EDGE_FIXTURES
        try
            _MBO.initialize();_MBO.model.add_box(row["origin"]...,row["extents"]...;tag=1)
            for edge in row["edges"]
                tag=edge["tag"];parameter=edge["parameter"]
                lower,upper=_MBO.model.get_parametrization_bounds(1,tag)
                @test [only(lower),only(upper)]≈edge["bounds"] atol=2e-11
                @test _MBO.model.get_value(1,tag,[parameter])≈edge["value"] atol=2e-11
                @test _MBO.model.get_parametrization(1,tag,edge["value"])≈edge["inverse"] atol=2e-11
                for face in edge["surfaces"]
                    surface=face["surface"]
                    uv=_MBO.model.reparametrize_on_surface(1,tag,[parameter],surface)
                    @test uv≈face["uv"] atol=2e-11
                    @test _MBO.model.get_value(2,surface,uv)≈face["value"] atol=2e-11
                end
            end
        finally
            _MBO.finalize()
        end
    end
end

@testset "Pinned Box endpoint collapse rejects before entity allocation" begin
    try
        _MBO.initialize();_MBO.model.add_point(0,0,0;tag=21)
        model=_MBO.CURRENT[]
        before=Dict(name=>deepcopy(getfield(model,name)) for name in fieldnames(typeof(model)))
        @test_throws ArgumentError _MBO.model.add_box(1e16,0,0,1,1,1;tag=1)
        @test _MBO.CURRENT[]===model
        @test all(name->getfield(model,name)==before[name],fieldnames(typeof(model)))
        @test _MBO.model.get_entities()==[(0,21)]
    finally
        _MBO.finalize()
    end
    # Actual primary admits the nearby representable widths. The rejected
    # OCC1e-9 absolute-tolerance case remains outside this exact collapse fix.
    for (origin,width) in ((1e16,4.0),(1e15,1.0))
        try
            _MBO.initialize();_MBO.model.add_box(origin,0,0,width,1,1;tag=1)
            @test length(_MBO.model.get_entities())==27
            for curve in 9:12
                lower,upper=_MBO.model.get_parametrization_bounds(1,curve)
                @test [only(lower),only(upper)]==[0.0,width]
                @test all(isfinite,_MBO.model.get_value(1,curve,[width/2]))
            end
        finally
            _MBO.finalize()
        end
    end
end

println("ORIENTATION_PUBLIC_FOCUS_FINISHED")

using Test,Tessella,Logging
const _MC_API=Tessella.API
include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_orientation_missing_curve_v1.jl"))
function _mc_cube()
    _MC_API.initialize();_MC_API.option("Mesh.Renumber",0)
    for (tag,point) in enumerate(((0,0,0),(1,0,0),(1,1,0),(0,1,0),
                                 (0,0,1),(1,0,1),(1,1,1),(0,1,1)))
        _MC_API.model.add_point(point...;tag)
    end
    for (tag,edge) in enumerate(((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
                                (7,8),(8,5),(1,5),(2,6),(3,7),(4,8)))
        _MC_API.model.add_line(edge...;tag)
    end
    for (tag,curves) in enumerate(([1,2,3,4],[5,6,7,8],[1,10,-5,-9],
                                  [2,11,-6,-10],[3,12,-7,-11],[4,9,-8,-12]))
        _MC_API.model.add_curve_loop(curves;tag)
        _MC_API.model.add_plane_surface([tag];tag)
    end
    _MC_API.model.add_surface_loop([-1,2,3,4,5,6];tag=1)
    _MC_API.model.add_volume([1];tag=1)
end
@testset "Original Curve vertex U storage controls barycenter orientation" begin
    for row in ORIENTATION_MISSING_CURVE_FIXTURES
        @testset "P$(row["order"]) storedU=$(row["curve_parameters_supplied"])" begin
            try
                _mc_cube()
                nodes=row["nodes"]
                _MC_API.mesh.add_nodes(1,1,[node["tag"] for node in nodes],
                    vcat((node["coords"] for node in nodes)...),
                    vcat((node["params"] for node in nodes)...))
                for cell in row["before"]
                    _MC_API.mesh.add_elements_by_type(cell["owner"],cell["type"],
                        [cell["tag"]],cell["nodes"])
                end
                maximum=(_MC_API.NODE_TAG_MAX[],_MC_API.ELEMENT_TAG_MAX[])
                logger=Test.TestLogger()
                with_logger(logger) do
                    _MC_API.mesh.set_outward_orientation(1)
                end
                @test [string(log.message) for log in logger.logs]==
                    [replace(message,"Warning: "=>"") for message in row["warnings_first"]]
                for cell in row["after"]
                    found=_MC_API.mesh.get_element(cell["tag"])
                    @test found[1]==cell["type"]
                    @test found[2]==cell["nodes"]
                    @test found[4]==cell["owner"]
                end
                after=_MC_API.mesh.get_elements(2,-1)
                repeated=Test.TestLogger()
                with_logger(repeated) do
                    _MC_API.mesh.set_outward_orientation(1)
                end
                @test [string(log.message) for log in repeated.logs]==
                    [replace(message,"Warning: "=>"") for message in row["warnings_second"]]
                @test _MC_API.mesh.get_elements(2,-1)==after
                @test (_MC_API.NODE_TAG_MAX[],_MC_API.ELEMENT_TAG_MAX[])==maximum
                for node in nodes
                    found=_MC_API.mesh.get_node(node["tag"])
                    @test found[1]==node["coords"]
                    @test found[2]==node["params"]
                    @test found[3:4]==(1,1)
                end
            finally
                _MC_API.finalize()
            end
        end
    end
end
println("ORIENTATION_MISSING_CURVE_FOCUS_FINISHED")
