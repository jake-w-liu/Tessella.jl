module ActualVolumePLCRegression
using SHA,Test,Tessella
using Tessella.Elements: MixedMesh,msh_dimension
const DATA=normpath(joinpath(@__DIR__,"..","artifacts","actual_volume_plc"))
const GEO=joinpath(DATA,"box_embedded_sheet.geo")
const PRIMARY=joinpath(DATA,"sdk4152_exact_sheet.json")
const CASES=Dict{String,Any}[]
@assert bytes2hex(sha256(read(GEO)))=="79108e0f08686fb3fcfde24d01251c96a430496d72781d9db30f97d982553b79"
@assert bytes2hex(sha256(read(PRIMARY)))=="39643c3330c637e910b957a0bd59709cd3c2e1d6d643fc7c4c7418effb24e933"
@assert bytes2hex(sha256(read(joinpath(DATA,"plain_box.geo"))))=="2c7a2cf7067759fcc91e0744f2e81b402a11bc9a1212c20791dbd5b5e01dbe6b"
@assert bytes2hex(sha256(read(joinpath(DATA,"sdk4152_conflicting_line.json"))))=="3e0260b0cb295812138886e8fb2bb1bc7a51d7f246799b48224321a9cc3c8add"
@assert bytes2hex(sha256(read(joinpath(DATA,"sdk4152_staged_sheet.json"))))=="dbdd729bb68414585fcffbe30f9f7334244e9c0ed27cbc5639777b5bfabdca88"
const A=Tessella.API
bits(values)=[string(reinterpret(UInt64,Float64(value));base=16,pad=16) for value in values]
function science(actual)
    @test actual isa MixedMesh && actual.entity_data!==nothing
    data=actual.entity_data
    @test allunique(data.external_node_tags)
    @test all(>(0),data.external_node_tags)
    edges=Set{Tuple{Int32,Int32}}();faces=Dict{NTuple{3,Int32},Int}()
    volume=zero(Rational{BigInt});tet_count=0;all_positive=true
    outer=NTuple{3,Int32}[];sheet=NTuple{3,Int32}[];lines=Tuple{Int32,Int32}[]
    for block in actual.blocks
        block.msh==4 || continue
        for tet in eachcol(block.nodes)
            tet_count+=1
            p,q,r,s=[Rational{BigInt}.(actual.coords[:,row]) for row in tet]
            u=q-p;v=r-p;w=s-p
            determinant=u[1]*(v[2]*w[3]-v[3]*w[2])-u[2]*(v[1]*w[3]-v[3]*w[1])+u[3]*(v[1]*w[2]-v[2]*w[1])
            @test determinant>0;all_positive &= determinant>0;volume+=determinant/6
            for i in 1:4,j in i+1:4;push!(edges,minmax(tet[i],tet[j]));end
            for slots in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
                face=Tuple(sort!(collect(tet[collect(slots)])))
                faces[face]=get(faces,face,0)+1
            end
        end
    end
    @test tet_count>0
    @test volume==1//1
    @test all(count->count in (1,2),values(faces))
    for (index,block) in enumerate(actual.blocks)
        dimension=msh_dimension(block.msh)
        for (column,cell) in enumerate(eachcol(block.nodes))
            entity=Int(data.block_entities[index][column])
            if block.msh==1 && entity in (13,14,15)
                push!(lines,minmax(cell...))
            elseif block.msh==2 && entity==7
                push!(sheet,Tuple(sort!(collect(cell))))
            elseif dimension==2 && entity in 1:6
                @test block.msh==2
                push!(outer,Tuple(sort!(collect(cell))))
            end
        end
    end
    @test length(lines)==6 && allunique(lines)
    @test all(edge->edge in edges,lines)
    @test length(sheet)==4 && allunique(sheet)
    @test all(face->get(faces,face,0)==2,sheet)
    @test allunique(outer)
    @test Set(outer)==Set(face for (face,count) in faces if count==1)
    return Dict("tet_count"=>tet_count,"volume_exact"=>string(volume),"all_exact_positive"=>all_positive,
        "strict_curve_pairs"=>[collect(edge) for edge in lines],"strict_sheet_faces"=>[collect(face) for face in sheet],
        "outer_shell_exact"=>Set(outer)==Set(face for (face,count) in faces if count==1),
        "line_pairs_exact"=>all(edge->edge in edges,lines),"sheet_faces_exact"=>all(face->get(faces,face,0)==2,sheet),
        "tet_cells"=>[collect(cell) for block in actual.blocks if block.msh==4 for cell in eachcol(block.nodes)],"nodes"=>length(data.external_node_tags),
        "node_tags"=>collect(data.external_node_tags),"xyz_bits"=>bits(vec(actual.coords)))
end
function lower_snapshot()
    records=Dict{String,Any}()
    for (dimension,tag) in ((1,13),(1,14),(1,15),(2,7))
        ids,xyz,params=A.mesh.get_nodes(dimension,tag)
        types,labels,cells=A.mesh.get_elements(dimension,tag)
        records["$dimension:$tag"]=Dict("nodes"=>collect(ids),"xyz_bits"=>bits(xyz),"params_bits"=>bits(params),
            "types"=>collect(types),"labels"=>[collect(v) for v in labels],"cells"=>[collect(v) for v in cells])
    end
    return records
end

const BOX=joinpath(DATA,"plain_box.geo")
function cube_science(mesh)
    mixed=mesh isa MixedMesh
    catalogue=mixed ? [(Int(block.msh),block.nodes,index) for (index,block) in enumerate(mesh.blocks)] : [(4,mesh.tets,1)]
    faces=Dict{NTuple{3,Int32},Int}();edges=Set{Tuple{Int32,Int32}}();volume=zero(Rational{BigInt});count=0
    for (msh,cells,_) in catalogue
        msh==4 || continue
        for cell in eachcol(cells)
            count+=1;p,q,r,t=[Rational{BigInt}.(mesh.coords[:,i]) for i in cell]
            u=q-p;v=r-p;w=t-p
            det=u[1]*(v[2]*w[3]-v[3]*w[2])-u[2]*(v[1]*w[3]-v[3]*w[1])+u[3]*(v[1]*w[2]-v[2]*w[1])
            @test det>0;volume+=det/6
            for i in 1:4,j in i+1:4;push!(edges,minmax(cell[i],cell[j]));end
            for skip in 1:4
                key=Tuple(sort!(Int32[cell[i] for i in 1:4 if i!=skip]));faces[key]=get(faces,key,0)+1
            end
        end
    end
    @test count>0;@test volume==1//1;@test all(n->n in (1,2),values(faces))
    if mixed
        outer=Set{NTuple{3,Int32}}()
        for (msh,cells,_) in catalogue,cell in eachcol(cells)
            if msh==2;push!(outer,Tuple(sort!(collect(cell))))
            elseif msh==1;@test minmax(cell...) in edges
            end
        end
        @test outer==Set(f for (f,n) in faces if n==1)
    end
    return Dict("tet_count"=>count,"exact_volume"=>string(volume),"tet_cells"=>[collect(c) for (msh,cells,_) in catalogue if msh==4 for c in eachcol(cells)],"xyz_bits"=>bits(vec(mesh.coords)))
end

function public_snapshot()
    records=Dict{String,Any}()
    for (dimension,tag) in A.model.get_entities()
        ids,xyz,params=A.mesh.get_nodes(dimension,tag)
        types,labels,cells=A.mesh.get_elements(dimension,tag)
        records["$dimension:$tag"]=(collect(ids),bits(xyz),bits(params),collect(types),[collect(v) for v in labels],[collect(v) for v in cells])
    end
    return (records,A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag())
end

function lower_signatures(mesh)
    data=mesh.entity_data
    return sort!([(data.node_entities[i],Tuple(bits(mesh.coords[:,i])),data.node_parametric[i]===nothing ? () : Tuple(bits(data.node_parametric[i]))) for i in eachindex(data.node_entities) if data.node_entities[i][1]<3])
end


try
    @testset "Current embedded-sheet exact actual input" begin
        @testset "Exact original GEO production" begin
            entry=Dict{String,Any}("route"=>"GEO original");push!(CASES,entry)
            try
                result=Tessella.execute_geo(GEO;mesh_dim=3)
                entry["science"]=science(result.model.meshing.mesh_authority.mesh)
                println("PLC_CASE_COMPLETE route=GEO tets=",entry["science"]["tet_count"]);flush(stdout)
            catch error
                entry["error"]=sprint(showerror,error);rethrow()
            end
        end
        @testset "API generated lower inputs" begin
            A.initialize()
            entry=Dict{String,Any}("route"=>"API lower then Volume");push!(CASES,entry)
            try
                A.open_geo!(GEO;mesh_dim=0);A.mesh.generate(1);A.mesh.generate(2)
                before=lower_snapshot();entry["before"]=before
                A.mesh.generate(3)
                after=lower_snapshot();entry["after"]=after
                for key in keys(before),field in ("nodes","xyz_bits","params_bits","types","labels")
                    @test before[key][field]==after[key][field]
                end
                for key in keys(before)
                    if startswith(key,"1:")
                        @test before[key]["cells"]==after[key]["cells"]
                    else
                        # SDK4.15.2 reorders Sheet Triangle columns/cycles on Mesh3.
                        # Preserve the exact public node support set of each cell.
                        triangle_keys(record)=sort!([Tuple(sort!(collect(flat[i:i+2])))
                            for flat in record["cells"] for i in 1:3:length(flat)])
                        @test triangle_keys(before[key])==triangle_keys(after[key])
                    end
                end
                entry["science"]=science(A.LAST_MESH[])
                println("PLC_CASE_COMPLETE route=API tets=",entry["science"]["tet_count"]);flush(stdout)
            catch error
                entry["error"]=sprint(showerror,error);rethrow()
            finally
                A.finalize()
            end
        end
        @testset "Ordinary Box actual PLC and OnlyEmpty" begin
            for onlyempty in (0,1)
                A.initialize()
                try
                    A.option("Mesh.MeshOnlyEmpty",onlyempty)
                    A.open_geo!(BOX;mesh_dim=0);A.mesh.generate(2)
                    before=lower_signatures(A.LAST_MESH[])
                    for iteration in 1:2
                        A.mesh.generate(3)
                        @test lower_signatures(A.LAST_MESH[])==before
                        entry=Dict("route"=>"plainBox OnlyEmpty$onlyempty iteration$iteration","science"=>cube_science(A.LAST_MESH[]))
                        push!(CASES,entry)
                        println("PLC_CASE_COMPLETE route=Box onlyempty=",onlyempty," iteration=",iteration);flush(stdout)
                    end
                finally
                    A.finalize()
                end
            end
        end
        @testset "Actual PLC explicit optimizer settings" begin
            for optimize in (0,1)
                A.initialize()
                try
                    A.option("Mesh.Optimize",optimize)
                    optimize==1 && A.option("Mesh.OptimizeThreshold",0.99)
                    A.open_geo!(GEO;mesh_dim=0);A.mesh.generate(2)
                    before=lower_snapshot()
                    A.mesh.generate(3)
                    after=lower_snapshot()
                    for key in keys(before),field in ("nodes","xyz_bits","params_bits","types","labels")
                        @test before[key][field]==after[key][field]
                    end
                    science(A.LAST_MESH[])
                finally
                    A.finalize()
                end
            end
        end
        @testset "Actual UID/coordinate input validation" begin
            A.initialize()
            try
                A.open_geo!(GEO;mesh_dim=0);A.mesh.generate(2)
                before=public_snapshot()
                model=A._with_model() do current
                    deepcopy(current)
                end
                authority=model.meshing.mesh_authority
                geometry=Tessella.Model._model_explicit_volume_geometry(model,1,"actualPLC regression";_mesh_input=authority)
                input=Tessella.Model._model_native_initial_embedded_input(model,1,geometry,authority,"actualPLC regression")
                duplicate=copy(input.native_node_tags);duplicate[2]=duplicate[1]
                @test_throws ArgumentError Tessella.Model._model_native_volume_input_constraints(model,1,merge(input,(native_node_tags=duplicate,)),authority,"actualPLC regression")
                missing=copy(input.native_node_tags);missing[1]=maximum(missing)+UInt64(1)
                @test_throws ArgumentError Tessella.Model._model_native_volume_input_constraints(model,1,merge(input,(native_node_tags=missing,)),authority,"actualPLC regression")
                coordinates=copy(input.surface.coords);coordinates[1,1]=nextfloat(coordinates[1,1])
                altered=Tessella.Elements.Mesh(coordinates;tris=input.surface.tris)
                @test_throws ArgumentError Tessella.Model._model_native_volume_input_constraints(model,1,merge(input,(surface=altered,)),authority,"actualPLC regression")
                @test public_snapshot()==before
            finally
                A.finalize()
            end
        end
        @testset "Conflicting actual Line rejects atomically" begin
            A.initialize()
            try
                A.open_geo!(GEO;mesh_dim=0);A.mesh.generate(2)
                first_node=first(A.mesh.get_nodes(0,10)[1]);last_node=first(A.mesh.get_nodes(0,12)[1])
                A.mesh.add_elements_by_type(13,1,UInt64[9001],UInt64[first_node,last_node])
                A.option("Mesh.MeshOnlyEmpty",1)
                before=public_snapshot()
                # SDK4.15.2 accepts this long Line but leaves it outside the Tet
                # edge complex. Tessella's FEM validity guard rejects publication;
                # this is a documented stricter-validity behavior, not SDK error parity.
                @test_throws ArgumentError A.mesh.generate(3)
                @test public_snapshot()==before
            finally
                A.finalize()
            end
        end
        @testset "Legacy no-input initial Box" begin
            result=Tessella.execute_geo(BOX;mesh_dim=0)
            result.model.meshing.algorithm_3d=3
            legacy=Tessella.Model.mesh_model_volume(result.model,1)
            push!(CASES,Dict("route"=>"legacy no-input Algorithm3 Box","science"=>cube_science(legacy)))
        end
    end
finally
    # Test captures have local lifetime and do not accumulate across the suite.
    empty!(CASES)
end
end # module
