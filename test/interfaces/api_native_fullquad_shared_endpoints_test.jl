module APINativeFullQuadSharedEndpointsTests
using Test,Tessella
if Sys.iswindows() && Sys.ARCH===:x86_64
const A=Tessella.API
call(f,args...;kwargs...)=Base.invokelatest(f,args...;kwargs...)
function shared_endpoints_geo()
    geo=joinpath(tempdir(),"fullquad_shared_endpoints.geo")
    # Two distinct Curves span the same two vertices: their paired
    # super-edges would share the (5,6) endpoint identity. The full-quad
    # coarse view cannot coarsen either Curve onto the same triangulation
    # edge, so both stay uncoarsened; subdivision still births one
    # Curve-owned middle per literal cell edge and the committed dim-1
    # products grow to match.
    open(geo,"w") do io
        print(io,"""
h=10;
Point(1)={0,0,0,h}; Point(2)={1,0,0,h}; Point(3)={1,1,0,h}; Point(4)={0,1,0,h};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
Point(5)={0.3,0.5,0,h}; Point(6)={0.7,0.5,0,h};
Point(7)={0.5,0.8,0,h}; Point(8)={0.5,0.2,0,h};
Spline(10)={5,7,6}; Spline(11)={6,8,5};
Curve{10,11} In Surface{1};
Transfinite Curve{1,2,3,4,10,11}=3;
Recombine Surface{1};
Mesh.RecombinationAlgorithm=2;
Mesh.Algorithm=3;
Mesh.MeshSizeMin=10; Mesh.MeshSizeMax=10;
Mesh.Smoothing=0; Mesh.Renumber=0;
Mesh.RecombineNodeRepositioning=0; Mesh.RecombineOptimizeTopology=0;
""")
    end
    return geo
end
function digon_face_geo()
    path=joinpath(tempdir(),"fullquad_digon_face.geo")
    open(path,"w") do io
        print(io,"""
h=10;
Point(1)={0,0,0,h}; Point(2)={2,0,0,h};
Point(3)={1,0.7,0,h}; Point(4)={1,-0.7,0,h};
Spline(1)={1,3,2}; Spline(2)={2,4,1};
Curve Loop(1)={1,2}; Plane Surface(1)={1};
Transfinite Curve{1,2}=3;
Recombine Surface{1};
Mesh.RecombinationAlgorithm=2;
Mesh.Algorithm=3;
Mesh.MeshSizeMin=10; Mesh.MeshSizeMax=10;
Mesh.Smoothing=0; Mesh.Renumber=0;
Mesh.RecombineNodeRepositioning=0; Mesh.RecombineOptimizeTopology=0;
""")
    end
    return path
end
function dim2_edges(elements)
    types,_,cells=elements
    edges=Set{Tuple{UInt64,UInt64}}()
    for (et,cc) in zip(types,cells)
        et!=3 && continue
        for column in 1:4:length(cc)
            cell=cc[column:column+3]
            for index in 1:4
                push!(edges,minmax(cell[index],cell[mod1(index+1,4)]))
            end
        end
    end
    return edges
end
@testset "native full-quad shared-endpoint Curve pairings" begin
@testset "embedded Curve digon keeps both chains as dim-2 edges" begin
    geo=shared_endpoints_geo()
    for algorithm in (2,3)
        call(A.initialize)
        call(A.mesh.generate,0)
        try
            text=replace(read(geo,String),"Mesh.RecombinationAlgorithm=2"=>
                "Mesh.RecombinationAlgorithm=$algorithm")
            path=tempname()*".geo";write(path,text)
            call(A.open_geo!,path;mesh_dim=2)
            elements=call(A.mesh.get_elements)
            types,tags,cells=elements
            @test Int32(3) in types
            edges=dim2_edges(elements)
            li=findfirst(==(Int32(1)),types)
            line_cells=[(cells[li][k],cells[li][k+1])
                        for k in 1:2:length(cells[li])]
            # The boundary square keeps its halved 2-cell products; the two
            # embedded splines grew to 4-cell chains whose cells all appear
            # as surface edges.
            @test length(line_cells)==16
            @test all(cell->minmax(cell...) in edges,line_cells)
            node_owner=Dict{UInt64,Tuple}()
            ntags,=call(A.mesh.get_nodes)
            for uid in ntags
                _,_,dim,own=call(A.mesh.get_node,uid)
                node_owner[uid]=(dim,own)
            end
            for curve in (10,11)
                chain=[cell for cell in line_cells if
                    node_owner[cell[1]]==(1,curve) || node_owner[cell[2]]==(1,curve)]
                @test length(chain)==4
            end
            @test call(A.mesh.get_elements,2,1)[1]==[3]
        finally
            call(A.finalize)
        end
    end
end
@testset "boundary digon Face pairs two shared-endpoint Curves" begin
    path=digon_face_geo()
    call(A.initialize)
    call(A.mesh.generate,0)
    try
        call(A.open_geo!,path;mesh_dim=2)
        elements=call(A.mesh.get_elements)
        types,tags,cells=elements
        @test Int32(3) in types
        @test !(Int32(2) in types)
        edges=dim2_edges(elements)
        li=findfirst(==(Int32(1)),types)
        line_cells=[(cells[li][k],cells[li][k+1])
                    for k in 1:2:length(cells[li])]
        @test all(cell->minmax(cell...) in edges,line_cells)
    finally
        call(A.finalize)
    end
end
end
end
end
