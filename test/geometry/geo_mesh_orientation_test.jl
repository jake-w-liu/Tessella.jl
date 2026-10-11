using Test,Tessella
using Tessella.MeshTypes: nnodes
include(joinpath(@__DIR__,"..","..","validation","native_geometry","native_orientation_geo_v1.jl"))

function _mgo_execute(source)
    mktempdir() do directory
        path=joinpath(directory,"orientation.geo")
        write(path,source)
        return Tessella.GeoExec.execute_geo(path)
    end
end

function _mgo_assert_signs(execution,signs,factory)
    outward=factory=="OpenCASCADE" ?
        ((-1.,0,0),(1.,0,0),(0.,-1,0),(0.,1,0),(0.,0,-1),(0.,0,1)) :
        ((0.,0,-1),(0.,0,1),(0.,-1,0),(1.,0,0),(0.,1,0),(-1.,0,0))
    found=Set{Int}()
    for (dimension,tag,part) in execution.mesh_parts
        dimension==2 || continue
        cells=part isa Tessella.Mesh ? [(2,part.tris)] :
            [(Int(block.msh),block.nodes) for block in part.blocks]
        for (msh,nodes) in cells
            Tessella.msh_dimension(msh)==2 || continue
            for column in axes(nodes,2)
                a=part.coords[:,nodes[1,column]]
                b=part.coords[:,nodes[2,column]]
                c=part.coords[:,nodes[3,column]]
                u=b-a;v=c-a
                normal=[u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1]]
                normal./=sqrt(sum(abs2,normal))
                if tag==7
                    @test normal[3]≈-1 atol=2e-11
                else
                    @test sum(normal.*collect(outward[tag]))≈signs[string(tag)] atol=2e-11
                    push!(found,tag)
                end
            end
        end
    end
    @test found==Set(1:6)
end

@testset "Actual GEO orientation dispatcher and retained STL lifecycle" begin
    for (index,row) in enumerate(ORIENTATION_GEO_FIXTURES)
        @testset "independent GEO row $index" begin
            prefix=row["prefix"]
            command=row["stage"]=="before_mesh" ?
                "ReorientMesh Volume{1};Mesh 2;" : "Mesh 2;ReorientMesh Volume{1};"
            execution=_mgo_execute(prefix*command)
            signs=row["stage"]=="before_mesh" && row["factory"]=="built_in" ?
                row["expected_initial_outward_signs"] : row["expected_materialized_outward_signs"]
            _mgo_assert_signs(execution,signs,row["factory"])
            @test execution.mesh!==nothing
            source_flags=row["initial_reverse_pattern"]=="all_inward" ?
                row["factory"]=="built_in" ? Set(2:6) : Set(1:6) : Set{Int}()
            if row["factory"]=="built_in"
                @test all(get(execution.model.meshing.reverse,(2,s),false)==(s in source_flags) for s in 1:6)
            end
            repeated=_mgo_execute(prefix*command*"ReorientMesh Volume{1};ReorientMesh Volume{1};")
            _mgo_assert_signs(repeated,row["expected_materialized_outward_signs"],row["factory"])
            regeneration=_mgo_execute(prefix*command*"ReorientMesh Volume{1};Delete Meshes;Mesh 2;")
            _mgo_assert_signs(regeneration,row["expected_materialized_outward_signs"],row["factory"])
            if row["factory"]=="built_in" && (row["stage"]=="before_mesh" || row["family"]=="quad")
                @test any(message->occursin("No valid STL triangulation",message),execution.warnings)
            end
        end
    end
    prefix=first(ORIENTATION_GEO_FIXTURES)["prefix"]
    reference=_mgo_execute(prefix*"Mesh 2;")
    for selector in ("99","")
        selected=_mgo_execute(prefix*"Mesh 2;ReorientMesh Volume{"*selector*"};")
        @test Tessella.mesh_crc(selected.mesh)==Tessella.mesh_crc(reference.mesh)
        @test isempty(selected.model.meshing.orientation_stl)
    end
end

println("ORIENTATION_GEO_FOCUS_FINISHED")
