module VTKLagrangePyramidChecks
using Test,Tessella,LinearAlgebra
const A=Tessella.API
const E=Tessella.Elements

# VTK ships no vtkLagrangePyramid class, so VTK_LAGRANGE_PYRAMID (74) has no
# authored node ordering. Tessella writes the Lagrange template over the
# vtkPyramid/vtkTriQuadraticPyramid topology (public VTK tables): five corner
# vertices (counterclockwise base, apex last), the eight edges in VTK order
# ((0,1),(1,2),(2,3),(3,0),(0,4),(1,4),(2,4),(3,4)), the four lateral
# triangular faces in VTK face order ((0,1,4),(1,2,4),(2,3,4),(3,0,4)) each
# carrying the interior of a Lagrange triangle of degree p-3 in its vertex
# frame, the quadrilateral base interior in flat i-fastest rows, then the
# interior nodes. On the integer lattice {(i,j,k):0<=i,j<=k<=p} a node
# (i,j,k) sits at reference position ((2i-k)/p,(2j-k)/p,1-k/p). The checks
# below verify the convention semantically: exported node coordinates are
# compared against affine combinations of the exported corner positions, not
# against the generator.

const VTK_PYR_EDGES=((1,2),(2,3),(3,4),(4,1),(1,5),(2,5),(3,5),(4,5))
const VTK_PYR_TRIFACES=((1,2,5),(2,3,5),(3,4,5),(4,1,5))

# Reduced-triangle lattice in VTK order (vertices, then edges, then its own
# interior) — the within-face node sequence of a Lagrange triangle.
function reduced_triangle_lattice(q)
    q==0 && return NTuple{2,Int}[(0,0)]
    nodes=NTuple{2,Int}[(0,0),(q,0),(0,q)]
    for i in 1:q-1;push!(nodes,(i,0));end
    for i in 1:q-1;push!(nodes,(q-i,i));end
    for i in 1:q-1;push!(nodes,(0,q-i));end
    q>=3 && append!(nodes,((x+1,y+1) for (x,y) in reduced_triangle_lattice(q-3)))
    nodes
end

# Expected reference coordinates (3 x n) under the documented ordering, built
# from the rule above for an independent full-sequence comparison.
function expected_lattice(p)
    p==0 && return NTuple{3,Int}[(0,0,0)]
    vertices=[(0,0,p),(p,0,p),(p,p,p),(0,p,p),(0,0,0)]
    out=NTuple{3,Int}[vertices...]
    for (a,b) in VTK_PYR_EDGES,i in 1:p-1
        va=vertices[a];vb=vertices[b]
        push!(out,ntuple(d->va[d]+(vb[d]-va[d])*i÷p,3))
    end
    if p>=3
        for (a,b,c) in VTK_PYR_TRIFACES,(x,y) in reduced_triangle_lattice(p-3)
            u=x+1;v=y+1
            push!(out,ntuple(d->(vertices[a][d]*(p-u-v)+vertices[b][d]*u+
                vertices[c][d]*v)÷p,3))
        end
    end
    p>=2 && append!(out,((i,j,p) for j in 1:p-1 for i in 1:p-1))
    p>=3 && append!(out,((i+1,j+1,k+2) for (i,j,k) in expected_lattice(p-3)))
    out
end

refnode(p,(i,j,k))=[(2i-k)/p,(2j-k)/p,1-k/p]

function read_vtk(path)
    lines=readlines(path)
    modern=first(lines)=="# vtk DataFile Version 5.1"
    (modern || first(lines)=="# vtk DataFile Version 2.0") && lines[3]=="ASCII" ||
        error("invalid VTK header")
    tokens=split(join(lines[4:end],' '));index=Ref(1)
    take()=let value=tokens[index[]];index[]+=1;value;end
    integer()=parse(Int,take())
    take()=="DATASET" && take()=="UNSTRUCTURED_GRID" || error("invalid VTK dataset")
    points=zeros(3,0);cells=Vector{Int}[];types=Int[]
    fields=Dict{Tuple{String,String},Any}();association="global"
    while index[]<=length(tokens)
        key=take()
        if key=="POINTS"
            count=integer();take()=="double" || error("invalid point type")
            points=reshape([parse(Float64,take()) for _ in 1:3count],3,count)
        elseif key=="CELLS"
            count=integer();extent=integer()
            if modern
                take()=="OFFSETS" && take()=="vtktypeint64" || error("invalid cell offsets")
                offsets=[integer() for _ in 1:count]
                count>=1 && first(offsets)==0 && last(offsets)==extent && issorted(offsets) ||
                    error("invalid cell extent")
                take()=="CONNECTIVITY" && take()=="vtktypeint64" || error("invalid cell connectivity")
                nodes=[integer() for _ in 1:extent]
                for cell in 1:count-1;push!(cells,nodes[offsets[cell]+1:offsets[cell+1]]);end
            else
                consumed=0
                for _ in 1:count
                    width=integer();push!(cells,[integer() for _ in 1:width]);consumed+=width+1
                end
                consumed==extent || error("invalid cell extent")
            end
        elseif key=="CELL_TYPES"
            count=integer();types=[integer() for _ in 1:count]
        elseif key in ("POINT_DATA","CELL_DATA")
            association=key;count=integer()
            count==(key=="POINT_DATA" ? size(points,2) : length(cells)) ||
                error("invalid data count")
        elseif key=="FIELD"
            take();count=integer()
            for _ in 1:count
                name=take();components=integer();tuples=integer();type=take()
                T=type=="double" ? Float64 : type=="vtktypeuint64" ? UInt64 : Int
                fields[(association,name)]=reshape(
                    [parse(T,take()) for _ in 1:components*tuples],components,tuples)
            end
        else
            error("unexpected VTK token $key")
        end
    end
    (;points,cells,types,fields,modern)
end

world(reference)=[2. 1. .5;.5 3. 1.;1. .5 4.]*reference .+ [7.,-2.,5.]

# Build one complete pyramid element on a discrete volume entity. Each node
# has a distinctive UID (tag = base+ordinal reversed) and carries its Gmsh
# ordinal triply as a distinctive nodal pattern, so a misplaced slot shows up
# in the position, the UID and the value columns alike.
function fixture(msh,name)
    A.model.add(name)
    _,dimension,_,count,local_coordinates,_=A.mesh.get_element_properties(msh)
    reference=vcat(reshape(local_coordinates,dimension,count),zeros(3-dimension,count))
    coordinates=world(reference)
    base=(UInt64(1)<<40)+UInt64(1000msh)
    tags=base.+UInt64.(count:-1:1)
    A.model.add_discrete_entity(dimension,7)
    A.mesh.add_nodes(dimension,7,tags,vec(coordinates))
    A.mesh.add_elements_by_type(7,msh,[5000+msh],tags)
    _,actual,_=A.mesh.get_elements(dimension,7);cell=only(only(actual))
    values=[[Float64(node),-Float64(node^2),Float64(node)*.5] for node in 1:count]
    (;name,tags,cell,coordinates,reference,values)
end

function check_pyramid_export(parsed,model,p)
    slots=only(parsed.cells).+1
    expected=world(reduce(hcat,(refnode(p,n) for n in expected_lattice(p))))
    uids=vec(parsed.fields[("POINT_DATA","node_uid")])
    positions=Dict(uid=>index for (index,uid) in enumerate(model.tags))
    @test parsed.modern && parsed.types==[74] && length(slots)==length(model.tags)
    @test isapprox(parsed.points[:,slots],expected;atol=1e-12)
    @test Set(uids)==Set(model.tags) && all(
        isapprox(parsed.points[:,slot],model.coordinates[:,positions[uids[slot]]];atol=1e-12)
        for slot in slots)
    @test parsed.fields[("POINT_DATA","step_0")][:,slots]==
        hcat((model.values[positions[uids[slot]]] for slot in slots)...)
    @test all(==(1),parsed.fields[("POINT_DATA","present_0")])
end

# Deep-debug ordering check: read every exported slot position back against
# the corners and the lattice classes, all in the written coordinates.
function check_slot_semantics(parsed,p)
    slots=only(parsed.cells).+1
    corners=parsed.points[:,slots[1:5]]
    n=5
    for (a,b) in VTK_PYR_EDGES,i in 1:p-1
        n+=1
        @test parsed.points[:,slots[n]]≈corners[:,a]+(i/p)*(corners[:,b]-corners[:,a]) atol=1e-12
    end
    if p>=3
        for (a,b,c) in VTK_PYR_TRIFACES,(x,y) in reduced_triangle_lattice(p-3)
            u=x+1;v=y+1;n+=1
            @test parsed.points[:,slots[n]]≈
                ((p-u-v)*corners[:,a]+u*corners[:,b]+v*corners[:,c])/p atol=1e-12
        end
    end
    for j in 1:p-1,i in 1:p-1
        n+=1
        @test parsed.points[:,slots[n]]≈
            corners[:,1]+(i/p)*(corners[:,2]-corners[:,1])+(j/p)*(corners[:,4]-corners[:,1]) atol=1e-12
    end
    d1=corners[:,2]-corners[:,1];d2=corners[:,4]-corners[:,1]
    frame=[corners[:,5]-corners[:,1] d1 d2]
    for slot in n+1:length(slots)
        # A volume-interior node is apex*height + base_lattice_point*(1-height)
        # with 0<height<1 and the base point strictly inside the base
        # parallelogram. Writing point = corner1 + l*(apex-c1) + s*d1 + t*d2
        # gives s/(1-l),t/(1-l) in (0,1) exactly when that is true.
        l,s,t=frame\(parsed.points[:,slots[slot]]-corners[:,1])
        @test 0<l<1
        @test -1e-9<s/(1-l)<1+1e-9 && -1e-9<t/(1-l)<1+1e-9
    end
    @test length(slots)==5+8(p-1)+
        4*(p>=3 ? length(reduced_triangle_lattice(p-3)) : 0)+
        (p-1)^2+(p>=3 ? length(expected_lattice(p-3)) : 0)
end

@testset "complete pyramids export as VTK_LAGRANGE_PYRAMID 74" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            for (msh,p) in ((14,2),(118,3),(119,4),(124,9))
                model=fixture(msh,"lagrange-pyr-$msh")
                for kind in ("NodeData","ElementNodeData")
                    view=A.view.add(kind)
                    tags=kind=="NodeData" ? model.tags : [model.cell]
                    data=kind=="NodeData" ? model.values : [vcat(model.values...)]
                    A.view.add_model_data(view,0,model.name,kind,tags,data,0.5,3)
                    path=joinpath(directory,"pyr$msh-$kind.vtk")
                    A.view.write(view,path)
                    parsed=read_vtk(path)
                    check_pyramid_export(parsed,model,p)
                    kind=="NodeData" && check_slot_semantics(parsed,p)
                end
            end
        end
    finally
        A.finalize()
    end
end

@testset "unsupported pyramid cells keep their paths" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir() do directory
            # Serendipity Pyramid13 stays on VTK_QUADRATIC_PYRAMID 27.
            model=fixture(19,"pyr-serendipity")
            view=A.view.add("serendipity")
            A.view.add_model_data(view,0,model.name,"NodeData",model.tags,model.values,0.,3)
            path=joinpath(directory,"s19.vtk");A.view.write(view,path)
            parsed=read_vtk(path)
            @test parsed.types==[27] && length(only(parsed.cells))==13
            # Linear Pyramid5 still writes as plain VTK_PYRAMID 14.
            model=fixture(7,"pyr-linear")
            view=A.view.add("linear")
            A.view.add_model_data(view,0,model.name,"NodeData",model.tags,model.values,0.,3)
            path=joinpath(directory,"p1.vtk");A.view.write(view,path)
            parsed=read_vtk(path)
            @test parsed.types==[14] && length(only(parsed.cells))==5
            # P0 and incomplete pyramids still reject before any output.
            for msh in (132,125,127,131)
                model=fixture(msh,"pyr-unsupported-$msh")
                view=A.view.add("unsupported")
                A.view.add_model_data(view,0,model.name,"NodeData",model.tags,model.values,0.,3)
                path=joinpath(directory,"bad$msh.vtk");write(path,"sentinel")
                @test_throws ArgumentError A.view.write(view,path)
                @test read(path,String)=="sentinel"
            end
        end
    finally
        A.finalize()
    end
end

# The Bergot rational basis u^i v^j w^k/(1-w)^min(i,j), max(i,j)+k<=p, must
# interpolate the Lagrange nodes: the exported coefficients solve V*C=I where
# V[r,c] is monomial r at reference node c — rebuilt here from lagrange_nodes.
function pyramid_vandermonde(msh)
    spec=E.msh_spec(msh);p=spec.order
    lattice=NTuple{3,Int}[(i,j,k) for i in 0:p for j in 0:p
        for k in 0:p-max(i,j) if !(spec.serendipity && i==p && j==p)]
    nodes=E.lagrange_nodes(msh);n=length(lattice)
    matrix=Matrix{Float64}(undef,n,n)
    for row in 1:n,column in 1:n
        i,j,k=lattice[row];u,v,w=nodes[:,column]
        matrix[row,column]=w==1.0 && max(i,j)>0 ? 0.0 :
            u^i*v^j*w^k/(1.0-w)^min(i,j)
    end
    return matrix
end

@testset "interpolation basis covers complete pyramids of every order" begin
    # The Vandermonde conditioning grows quickly with p; residuals of ~1e-5 at
    # n=385 are floating-point accumulation in the solve, not a basis defect.
    for (msh,atol) in ((14,1e-10),(118,1e-10),(119,1e-9),(120,1e-8),(124,1e-4))
        basis=A._native_view_writer_pyramid_basis(msh)
        spec=E.msh_spec(msh)
        @test basis.family==6 && size(basis.coefficients)==(spec.nnodes,spec.nnodes) &&
            size(basis.exponents)==(spec.nnodes,3)
        @test norm(basis.coefficients*pyramid_vandermonde(msh)-I)<atol
    end
    basis19=A._native_view_writer_pyramid_basis(19)
    @test basis19.family==6 && size(basis19.coefficients,1)==13
    for msh in (125,127,131,7,132)
        @test_throws ArgumentError A._native_view_writer_pyramid_basis(msh)
    end
end

@testset "msh view export writes the scheme for a cubic pyramid" begin
    A.initialize()
    try
        A.option("General.Terminal",0);A.option("PostProcessing.SaveInterpolationMatrices",1)
        model=fixture(118,"pyr-p3-scheme")
        view=A.view.add("scheme")
        A.view.add_model_data(view,0,model.name,"NodeData",model.tags,model.values,0.,3)
        mktempdir() do directory
            path=joinpath(directory,"scheme.msh");A.view.write(view,path)
            text=replace(read(path,String),"\r\n"=>"\n")
            @test occursin("\$InterpolationScheme",text)
            lines=split(text,'\n')
            row=findfirst(==("6"),lines)
            @test row!==nothing
            dims=split(lines[row+2],' ')
            @test parse(Int,dims[1])==30 && parse(Int,dims[2])==30
        end
    finally
        A.finalize()
    end
end
end
