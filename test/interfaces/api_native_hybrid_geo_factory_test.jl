module NativeHybridGeoFactoryChecks
using Test,TOML,SHA,Tessella
using Serialization
const M=Tessella.Model
const E=Tessella.Elements
const G=Tessella.GeoExec
const I=Tessella.IO
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const native=TOML.parsefile(joinpath(ROOT,"test/artifacts/native_hybrid_geo_import_4152_v16.toml"))
hexpoint(xs)=Float64[Core.eval(@__MODULE__,Meta.parse(x)) for x in xs]

function box_model()
    m=M.GeoModel()
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (tag,p) in enumerate(points);M.add_point!(m,p...;tag,mesh_size=.4);end
    curves=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (tag,ends) in enumerate(curves);M.add_line!(m,ends...;tag);end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),(2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (tag,curves) in enumerate(loops)
        M.add_curve_loop!(m,collect(curves);tag);M.add_plane_surface!(m,[tag];tag)
    end
    M.add_surface_loop!(m,collect(1:6);tag=1);M.add_volume!(m,[1];tag=1)
    return m
end

function pending_input(case)
    m=box_model();context=I._GeoNumericContext()
    done=case["warm"] ? Set{Tuple{Int,Int32}}((dim,Int32(tag))
        for (dim,tag) in M.model_entities(m) if dim in (1,2)) : Set{Tuple{Int,Int32}}()
    authority=M._ModelMeshAuthority(nothing,Dict{Tuple{Int,Int32},M._ModelMeshProduct}(),done,UInt64[])
    m.meshing.mesh_authority=authority
    context.mesh_node_max=UInt64(case["initial_node_max"])
    context.mesh_element_max=UInt64(case["initial_element_max"])
    return m,context,authority,done
end

function publication(context,m)
    stream=IOBuffer()
    serialize(stream,(context.mesh,context.mesh_parts,context.mesh_part_rows,
        context.mesh_geometry_parts,context.mesh_node_owner,
        context.mesh_node_parametric,context.mesh_node_max,context.mesh_element_max,
        context.native_element_lookup_receipt,m.meshing.mesh_authority,
        m.meshing.imported_entities,m.meshing.imported_ghost_elements))
    return bytes2hex(sha256(take!(stream)))
end

function check_products(actual,products)
    data=actual.entity_data
    @test data.external_node_tags==UInt64[uid for product in products for uid in product["node_tags"]]
    for product in products
        owner=(Int(product["owner"][1]),Int32(product["owner"][2]))
        rows=findall(==(owner),data.node_entities)
        @test data.external_node_tags[rows]==UInt64.(product["node_tags"])
        @test isequal(actual.coords[:,rows],isempty(rows) ? zeros(3,0) : reduce(hcat,hexpoint.(product["xyz_hex"])))
        @test isequal(reduce(vcat,(data.node_parametric[row]===nothing ? Float64[] : data.node_parametric[row]
            for row in rows);init=Float64[]),hexpoint(product["uv_hex"]))
        observed=NamedTuple[]
        for (block,owners,labels) in zip(actual.blocks,data.block_entities,data.external_element_tags)
            E.msh_dimension(block.msh)==owner[1] || continue
            columns=findall(==(owner[2]),owners);isempty(columns) && continue
            push!(observed,(msh=Int(block.msh),labels=labels[columns],
                nodes=[Tuple(data.external_node_tags[block.nodes[:,column]]) for column in columns]))
        end
        @test observed==[(msh=Int(cell["msh"]),labels=UInt64.(cell["labels"]),
            nodes=Tuple.(Vector{UInt64}.(cell["nodes"]))) for cell in get(product,"cell",[])]
    end
end

@testset "GEO imports own factory events and preserve completed lower stages" begin
    @test bytes2hex(sha256(read(joinpath(ROOT,native["primary"]))))==native["primary_sha256"]
    @test length(native["case"])==9
    for case in native["case"]
        println("GEO_FACTORY_IMPORT_CASE ",case["name"]);flush(stdout)
        m,context,previous,done=pending_input(case)
        path=joinpath(ROOT,"test/data/native_hybrid_import_events_4152_v1",case["name"]*".msh")
        saved=publication(context,m)
        @test G._geo_native_lower_factory_scope(m,context,path)
        @test isequal(publication(context,m),saved)
        @test G._geo_import_mesh_file!(m,context,path,"native GEO factory")===nothing
        authority=m.meshing.mesh_authority
        @test authority!==previous && authority.done==done && authority.done!==previous.done
        @test previous.mesh===nothing && previous.done==done
        @test M._model_mesh_status(m,authority)==(case["warm"] ? 2 : 0)
        @test (context.mesh_node_max,context.mesh_element_max)==
            (UInt64(case["import_node_max"]),UInt64(case["import_element_max"]))
        check_products(authority.mesh,case["imported"])
        @test authority.carrier===context.mesh && authority.carrier_rows==collect(Int32,1:size(authority.mesh.coords,2))
        @test isequal(context.mesh.coords,authority.mesh.coords)
        @test context.mesh_geometry_parts.carrier===context.mesh
        @test context.mesh_geometry_parts.node_tags==authority.mesh.entity_data.external_node_tags
        @test length(context.mesh_parts)==length(authority.products)
        for owner in keys(authority.surface_constructor_inputs)
            actual=M._model_mesh_product(authority,2,Int(owner[2]))
            @test M._model_surface_constructor_binding(actual,authority.surface_constructor_inputs[owner],"GEO native")===nothing
        end
        lookup=context.native_element_lookup_receipt
        @test lookup.factory_owner.entity_data.external_node_tags==authority.mesh.entity_data.external_node_tags
        @test M._model_mesh_factory_lookup_binding(lookup,authority.mesh,"GEO native")===nothing
        if occursin("zero",case["name"])
            @test haskey(lookup.node_targets,UInt64(0)) && haskey(lookup.element_targets,UInt64(0))
            event=only(filter(row->row[3]==0,lookup.factory_identities.factory_events.node_events))
            expected=only(filter(product->product["owner"]==[event[1],event[2]],case["imported"]))
            @test lookup.node_targets[UInt64(0)]==UInt64(first(expected["node_tags"]))
            event=only(filter(row->row[3]==0,lookup.factory_identities.factory_events.cell_events))
            expected=only(filter(product->product["owner"]==[E.msh_dimension(event[1]),event[2]],case["imported"]))
            cells=only(filter(cell->cell["msh"]==event[1],expected["cell"]))
            position=occursin("within_family_reverse",case["name"]) ? lastindex(cells["labels"]) : firstindex(cells["labels"])
            @test lookup.element_targets[UInt64(0)]==UInt64(cells["labels"][position])
        end
        saved=publication(context,m)
        @test !G._geo_native_lower_factory_scope(m,context,path)
        @test !G._geo_import_cleared_native_factory!(m,context,path,"native GEO live append")
        @test isequal(publication(context,m),saved)
        clocks=(context.mesh_node_max,context.mesh_element_max)
        @test G._geo_exec_delete!(m,false,"Meshes;",context,nothing)===nothing
        cleared=m.meshing.mesh_authority
        @test cleared.mesh===nothing && cleared.done==done && cleared.done!==authority.done
        @test isempty(cleared.surface_constructor_inputs) && isempty(cleared.products)
        @test context.mesh===nothing && isempty(context.mesh_parts) && isempty(context.mesh_part_rows)
        @test context.mesh_geometry_parts===nothing && isempty(context.mesh_node_parametric)
        @test context.native_element_lookup_receipt===nothing
        @test (context.mesh_node_max,context.mesh_element_max)==clocks
    end
end

@testset "GEO factory refuses unsupported or invalid histories before publication" begin
    case=first(native["case"])
    path=joinpath(ROOT,"test/data/native_hybrid_import_events_4152_v1","ordinary.msh")
    m,context,previous,done=pending_input(case)
    source=read(path,String)
    mktempdir() do directory
        for (name,bytes) in (("legacy.msh",replace(source,"4.1 0 8"=>"2.2 0 8";count=1)),
                ("binary.msh",replace(source,"4.1 0 8"=>"4.1 1 8";count=1)),
                ("ancillary.msh",source*"\n\$NodeData\n\$EndNodeData\n"))
            file=joinpath(directory,name);write(file,bytes);saved=publication(context,m)
            @test !G._geo_native_lower_factory_scope(m,context,file)
            @test !G._geo_import_cleared_native_factory!(m,context,file,"native unsupported")
            @test isequal(publication(context,m),saved)
        end
    end
    delete!(m.surfaces,1);saved=publication(context,m)
    @test !G._geo_native_lower_factory_scope(m,context,path)
    @test !G._geo_import_cleared_native_factory!(m,context,path,"foreign CAD")
    @test isequal(publication(context,m),saved)
    for clock in (UInt64(20),typemax(UInt64))
        m,context,previous,done=pending_input(case);context.mesh_node_max=clock
        path=joinpath(ROOT,"test/data/native_hybrid_import_events_4152_v1","zero.msh")
        saved=publication(context,m)
        @test G._geo_native_lower_factory_scope(m,context,path)
        @test_throws ArgumentError G._geo_import_cleared_native_factory!(m,context,path,"colliding native factory")
        @test isequal(publication(context,m),saved)
    end
end
end
