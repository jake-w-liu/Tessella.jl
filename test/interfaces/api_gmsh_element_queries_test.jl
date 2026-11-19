# Element queries on meshes loaded from MSH files answer what Gmsh 4.15.2 answers.
# The test/fixtures/gmsh_4152_api_*.msh files were written by Gmsh 4.15.2 and
# every expected value below is what `gmsh` returned after reloading the file.
module GmshElementQueryParityChecks
using Test, Tessella
const A=Tessella.API
const DATA=joinpath(@__DIR__,"..","fixtures")

function load(name)
    A.finalize();A.initialize();A.option("General.Terminal",0)
    A.open(joinpath(DATA,"gmsh_4152_api_"*name*".msh"))
    return nothing
end
quality(tags,name)=A.mesh.get_element_qualities(tags,name)
sorted_by_tag(tags,nodes,width)=(perm=sortperm(tags);
    (Int.(tags[perm]),Int[nodes[width*(i-1)+k] for i in perm for k in 1:width]))

@testset "element qualities of every family" begin
    try
        load("square_o1")
        # eta: one minus the deviation of the smallest angle from 60 degrees
        # on a triangle, zero on a line or a point.
        @test quality([9,10],"eta")≈[0.75,0.75] rtol=1e-14
        @test quality(5:8,"eta")==zeros(4) && quality(1:4,"eta")==zeros(4)
        # A line has no gradient-error measure: Gmsh answers DBL_MAX.
        @test quality(5:8,"minSIGE")==fill(floatmax(Float64),4)
        # The MElement defaults of a point element.
        @test quality(1:4,"minSIGE")==fill(floatmax(Float64),4)
        @test quality(1:4,"angleShape")==ones(4)
        @test quality(1:4,"gamma")==zeros(4) && quality(1:4,"innerRadius")==zeros(4)
        @test quality(1:4,"outerRadius")==zeros(4)
        @test_throws ArgumentError quality(1:4,"minSJ")
        load("tets")
        @test quality(33:38,"eta")≈[0.8399473665965821,0.6872296635790217,
            0.7559526299369239,0.7559526299369239,0.6872296635790217,0.8399473665965821] rtol=1e-14
        # Tags in any order, repeated, and across element blocks (tetrahedra,
        # a triangle, a point).
        mixed_tags=[38,33,33,21,36,1,37]
        @test quality(mixed_tags,"volume")≈[1/6,1/6,1/6,0.5,1/6,1.0,1/6] rtol=1e-14
        @test quality(mixed_tags,"gamma")≈[0.7320508075688774,0.7320508075688774,
            0.7320508075688774,0.82842712474619,0.7174389352143008,0.0,0.6229863075268959] rtol=1e-13
        @test quality(mixed_tags,"eta")≈[0.8399473665965821,0.8399473665965821,
            0.8399473665965821,0.75,0.7559526299369239,0.0,0.6872296635790217] rtol=1e-14
        load("pyramid")
        @test quality([1],"eta")==[0.0]
        load("prism")
        @test quality([21],"eta")==[0.0]
        # _computeIGE has its own prism formula.
        @test quality([21],"minSIGE")≈[0.9754431360377266] rtol=1e-13
    finally
        A.finalize()
    end
end

@testset "element qualities above order two and on curved elements" begin
    try
        load("square_o3")
        # Cubic triangles take the Jacobian-based paths and keep the radii of
        # their primary vertices.
        @test quality([9,10],"minSICN")≈[0.8660254037844383,0.8660254037844335] rtol=1e-12
        @test quality([9,10],"minSIGE")≈[0.9292312334115675,0.9292312334115649] rtol=1e-12
        @test quality([9,10],"minSJ")≈[0.9999999999999949,0.9999999999999731] rtol=1e-12
        @test quality([9,10],"volume")≈[0.5,0.5] rtol=1e-13
        @test quality([9,10],"eta")≈[0.75,0.75] rtol=1e-14
        @test quality([9,10],"angleShape")≈fill(0.0029285737614364003,2) rtol=1e-9
        @test quality([9,10],"innerRadius")≈fill(0.29289321881345237,2) rtol=1e-14
        @test quality([9,10],"outerRadius")≈fill(0.7071067811865478,2) rtol=1e-14
        # Cubic lines: the chord measures and the 1-D defaults.
        @test quality(5:8,"minSJ")≈ones(4) rtol=1e-14
        @test quality(5:8,"volume")≈ones(4) rtol=1e-15
        @test quality(5:8,"innerRadius")≈fill(0.5,4) rtol=1e-15
        @test quality(5:8,"minSICN")==zeros(4) && quality(5:8,"angleShape")==ones(4)
        @test quality(5:8,"minSIGE")==fill(floatmax(Float64),4)
        load("disk_o2")
        # A 9-node quadrangle bent onto the unit circle: the angle measure
        # samples its own Jacobian, not the square of its corners.
        @test quality([6],"angleShape")≈[0.0038839754673125285] rtol=1e-9
        @test quality([6],"eta")==[1.0]
        @test quality([6],"minSJ")≈[0.31370849898475905] rtol=1e-12
        @test quality([6],"minSICN")≈[0.18603456092608572] rtol=1e-12
        @test quality(2:5,"minSJ")≈ones(4) rtol=1e-14
        @test quality(2:5,"volume")≈fill(1.4142135623730951,4) rtol=1e-15
        @test quality(2:5,"eta")==zeros(4)
    finally
        A.finalize()
    end
end

@testset "orientation queries rank the vertices of any family" begin
    try
        load("pyramid")
        @test A.mesh.get_number_of_orientations(7,"H1Legendre1")==120
        @test A.mesh.get_basis_functions_orientation(7,"H1Legendre1")==[32]
        @test A.mesh.get_basis_functions_orientation(7,"Lagrange")==[0]
        @test A.mesh.get_basis_functions_orientation_for_element(1,"HcurlLegendre1")==32
        load("square_o1")
        @test A.mesh.get_number_of_orientations(15,"H1Legendre2")==1
        @test A.mesh.get_basis_functions_orientation(15,"HcurlLegendre1")==zeros(4)
        # HierarchicalBasisH1Point reports (1, 0, 0) as its gradient.
        components,values,orientations=A.mesh.get_basis_functions(15,[0.,0.,0.],"GradH1Legendre2")
        @test components==3 && values==[1.0,0.0,0.0] && orientations==1
        # An unknown family name answers 0; "Quadrangle" is the family.
        @test A.mesh.get_element_type("Quadrilateral",1)==0
        @test A.mesh.get_element_type("Quadrangle",1)==3
    finally
        A.finalize()
    end
end

@testset "edges and faces of every order, in Gmsh's numbering" begin
    try
        load("square_o3")
        # Hierarchical keys on cubic triangles number the edges and faces of
        # the primary vertices.
        types,keys,_=A.mesh.get_keys(21,"H1Legendre3",-1,false)
        @test types==[0,0,0,1,2,1,2,1,2,3,0,0,0,1,2,1,2,1,2,3]
        @test keys==[1,2,4,1,1,2,2,3,3,1,4,2,3,2,2,4,4,5,5,2]
        @test sorted_by_tag(A.mesh.get_all_edges()...,2)==(collect(1:5),[1,2,2,4,4,1,2,3,3,4])
        @test sorted_by_tag(A.mesh.get_all_faces(3)...,3)==([1,2],[1,2,4,4,2,3])
        A.mesh.create_edges()
        @test sorted_by_tag(A.mesh.get_all_edges()...,2)==(collect(1:5),[1,2,2,4,4,1,2,3,3,4])
        load("cube_o2")
        # A 20-node hexahedron; faces registered by the key query keep the
        # vertex order of MElement::getFaceSolin.
        types,keys,_=A.mesh.get_keys(17,"H1Legendre2",-1,false)
        @test keys==[3,1,2,4,7,5,6,8,1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,27]
        @test sorted_by_tag(A.mesh.get_all_faces(4)...,4)==(collect(1:6),
            [3,1,4,2,3,1,7,5,3,4,7,8,1,2,5,6,4,2,8,6,7,5,8,6])
        @test sorted_by_tag(A.mesh.get_all_edges()...,2)==(collect(1:12),
            [3,1,3,4,3,7,1,2,1,5,2,4,2,6,4,8,7,5,7,8,5,6,6,8])
        load("tets")
        A.mesh.get_keys(4,"H1Legendre3",-1,false)
        faces=sorted_by_tag(A.mesh.get_all_faces(3)...,3)
        @test faces[1]==1:18 && faces[2][1:12]==[1,2,4,1,2,6,1,4,6,2,4,6]
        load("prism")
        # At order two a prism owns functions on its quadrangles only: the
        # triangles are not registered.
        types,keys,_=A.mesh.get_keys(6,"H1Legendre2",-1,false)
        @test keys==[1,2,3,4,5,6,1,2,3,4,5,6,7,8,9,1,2,3]
        @test sorted_by_tag(A.mesh.get_all_faces(4)...,4)==([1,2,3],[1,2,4,5,1,3,4,6,2,3,5,6])
        @test isempty(A.mesh.get_all_faces(3)[1])
        load("disk_o2")
        # createEdges walks the entities by dimension: the curve's lines first.
        A.mesh.create_edges()
        @test sorted_by_tag(A.mesh.get_all_edges()...,2)==([1,2,3,4],[1,2,2,3,3,4,4,1])
    finally
        A.finalize()
    end
end

end # module
