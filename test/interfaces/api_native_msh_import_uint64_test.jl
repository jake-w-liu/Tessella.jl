module APINativeMshImportUInt64Tests
using Tessella,Test
const E=Tessella.Elements
const G=Tessella.GeoExec
const M=Tessella.Model
const A=Tessella.API

function import_source(path;command="Merge",tail="")
    mktempdir() do directory
        geo=joinpath(directory,"import.geo")
        write(geo,"$command \"$(replace(abspath(path),'\\'=>'/'))\";\n"*tail)
        return A.open_geo!(geo)
    end
end
function with_session(f)
    A.initialize()
    try f() finally A.finalize() end
end

@testset "fixed-cell corner projection retains columns and constant cells" begin
    for (msh,width,coords,nodes,property) in (
        (8,2,Float64[0 1 2 3 .5 2.5;zeros(2,6)],Int32[1 3;2 4;5 6],:segs),
        (9,3,Float64[0 1 0 2 3 2 .5 .5 0 2.5 2.5 2;
                      0 0 1 0 0 1 0 .5 .5 0 .5 .5;zeros(1,12)],
            Int32[1 4;2 5;3 6;7 10;8 11;9 12],:tris),
        (11,4,Float64[0 1 0 0 2 3 2 2 .5 .5 0 0 .5 0 2.5 2.5 2 2 2.5 2;
                       0 0 1 0 0 0 1 0 0 .5 .5 0 0 .5 0 .5 .5 0 0 .5;
                       0 0 0 1 0 0 0 1 0 0 0 .5 .5 .5 0 0 0 .5 .5 .5],
            Int32[1 5;2 6;3 7;4 8;9 15;10 16;11 17;12 18;13 19;14 20],:tets))
        mixed=E.MixedMesh(coords,[E.ElementBlock(msh,nodes,Int32[0,0])])
        folded=G._geo_mixed_simplex_mesh(mixed,"two-column corner projection")
        @test getproperty(folded,property)==nodes[1:width,:]
    end
    for msh in (84,85,87)
        mixed=E.MixedMesh(zeros(3,1),[E.ElementBlock(msh,ones(Int32,1,1),Int32[0])])
        folded=G._geo_mixed_simplex_mesh(mixed,"constant-cell import")
        @test isempty(folded.segs) && isempty(folded.tris) && isempty(folded.tets)
    end
end

@testset "strict UInt64 record constructors" begin
    @test !applicable(M._record_append_element!,M.DiscreteEntity(),1,Int32(7),Int32[1,2])
    entity=E.MixedEntity(1,7,(0.,0.,0.,1.,0.,0.))
    for value in (true,-1,BigInt(typemax(UInt64))+1,1.0)
        @test_throws ArgumentError E.MixedEntityData(Dict((1,7)=>entity);
            node_entities=[(1,7)],node_parametric=[nothing],external_node_tags=[value])
    end
end

@testset "UInt64 import records retain linear allocation growth" begin
    allocated=Int[]
    for count in (256,512,1024)
        nodes=reshape(collect(Int32,1:count),2,:)
        coords=zeros(3,count);coords[1,:].=1:count
        mixed=E.MixedMesh(coords,[E.ElementBlock(1,nodes,zeros(Int32,count÷2))])
        function ingest()
            record=M.DiscreteEntity()
            G._geo_append_discrete_import_nodes!(record,mixed,collect(Int32,1:count),nothing)
            return record
        end
        ingest();GC.gc()
        push!(allocated,@allocated ingest())
    end
    @test allocated[2]<=2.8allocated[1]+65_536
    @test allocated[3]<=2.8allocated[2]+65_536
    println("IO_WIDTH_BATCH_NODE_ALLOCATIONS ",allocated)
    public_allocations=Dict{String,Vector{Int}}()
    for shape in ("one_curve","many_points")
        values=Int[]
        mktempdir() do directory
            for count in (128,256,512)
                path=joinpath(directory,"growth-$shape-$count.msh")
                open(path,"w") do stream
                    println(stream,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n",count)
                    for row in 1:count
                        println(stream,UInt64(2)^63+UInt64(row)," ",row," 0 0")
                    end
                    println(stream,"\$EndNodes\n\$Elements\n",shape=="many_points" ? count : count÷2)
                    for cell in 1:(shape=="many_points" ? count : count÷2)
                        uid=UInt64(2)^63+UInt64(cell)
                        if shape=="many_points"
                            println(stream,uid," 15 2 0 ",cell," ",UInt64(2)^63+UInt64(cell))
                        else
                            println(stream,uid," 1 2 0 7 ",UInt64(2)^63+UInt64(2cell-1),
                                " ",UInt64(2)^63+UInt64(2cell))
                        end
                    end
                    println(stream,"\$EndElements")
                end
                function public_import()
                    with_session() do
                        import_source(path)
                        @assert length(A.mesh.get_nodes()[1])==count
                        @assert A.mesh.get_max_node_tag()==UInt64(2)^63+UInt64(count)
                    end
                end
                with_session() do
                    import_source(path)
                    expected_nodes=UInt64(2)^63 .+ UInt64.(1:count)
                    expected_cells=shape=="many_points" ? count : count÷2
                    kind=shape=="many_points" ? 15 : 1
                    @test A.mesh.get_nodes()[1]==expected_nodes
                    @test A.mesh.get_elements()==(Int32[kind],
                        [expected_nodes[1:expected_cells]],[expected_nodes])
                    @test A.mesh.get_max_element_tag()==UInt64(2)^63+UInt64(expected_cells)
                    for cell in (1,expected_cells)
                        connected=shape=="many_points" ? expected_nodes[cell:cell] : expected_nodes[2cell-1:2cell]
                        owner=shape=="many_points" ? cell : 7
                        dimension=shape=="many_points" ? 0 : 1
                        @test A.mesh.get_element(expected_nodes[cell])==
                            (Int32(kind),connected,Int32(dimension),Int32(owner))
                    end
                end
                public_import();GC.gc()
                push!(values,@allocated public_import())
            end
        end
        @test values[2]<=3.1values[1]+262_144
        @test values[3]<=3.1values[2]+262_144
        public_allocations[shape]=values
    end
    println("IO_WIDTH_FULL_PUBLIC_IMPORT_ALLOCATIONS ",public_allocations)
end
end
