module MeshAdaptFunctionalChecks
using Test,Tessella,LinearAlgebra
const A=Tessella.API
const K=Tessella.Model

function node_map()
    tags,xyz,_=A.mesh.get_nodes()
    Dict(uid=>Tuple(xyz[3i-2:3i]) for (i,uid) in enumerate(tags))
end
function face_records(records)
    kinds,labels,connected=records
    Dict(uid=>(kind,Tuple(cell)) for (kind,uids,nodes) in zip(kinds,labels,connected)
        for (uid,cell) in zip(uids,eachcol(reshape(nodes,Tessella.Elements.msh_num_nodes(kind),:))))
end
function surface_edges(face)
    kinds,tags,nodes=A.mesh.get_elements(2,face)
    incidence=Dict{Tuple{UInt64,UInt64},Int}()
    cells=Vector{Vector{UInt64}}()
    for (kind,labels,connected) in zip(kinds,tags,nodes)
        @test kind in (2,3)
        width=kind==2 ? 3 : 4
        for cell in eachcol(reshape(connected,width,length(labels)))
            push!(cells,collect(cell))
            for i in 1:width
                edge=minmax(cell[i],cell[mod1(i+1,width)])
                incidence[edge]=get(incidence,edge,0)+1
            end
        end
    end
    return incidence,cells
end
function patch(warp,algorithm;recombine=false,embedded=false,bounded=true)
    A.initialize()
    try
        for (name,value) in (("General.Terminal",0),("Mesh.Renumber",0),
                ("Mesh.Algorithm",algorithm),("Mesh.Smoothing",0))
            A.option(name,value)
        end
        if bounded
            A.option("Mesh.MeshSizeMin",.35);A.option("Mesh.MeshSizeMax",.35)
        end
        if warp==0
            A.model.add("local edge adaptation")
            points=[A.model.add_point(p...;meshSize=.35) for p in
                ((0.,0.,0.),(1.,0.,0.),(1.,1.,warp),(0.,1.,0.))]
            curves=[A.model.add_line(points[i],points[mod1(i+1,4)]) for i in 1:4]
            for curve in curves;A.mesh.set_transfinite_curve(curve,5);end
            loop=A.model.add_curve_loop(curves)
            face=A.model.add_plane_surface([loop])
        else
            # Ruled Surface construction is exposed by the public GEO
            # frontend rather than the model namespace's plane-only helper.
            path=joinpath(mktempdir(),"ruled.geo")
            io=open(path,"w")
            write(io,"""
Point(1)={0,0,0,.35}; Point(2)={1,0,0,.35};
Point(3)={1,1,$warp,.35}; Point(4)={0,1,0,.35};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Curve Loop(1)={1,2,3,4}; Surface(1)={1};
Transfinite Curve{1,2,3,4}=5;
""")
            close(io);A.open_geo!(path)
            curves=collect(1:4);face=1
        end
        if recombine
            A.option("Mesh.RecombinationAlgorithm",0)
            A.mesh.set_recombine(2,face)
        end
        embedded_curve=0
        if embedded
            p=A.model.add_point(.25,.5,0.;meshSize=.35)
            q=A.model.add_point(.75,.5,0.;meshSize=.35)
            embedded_curve=A.model.add_line(p,q)
            A.mesh.set_transfinite_curve(embedded_curve,3)
            A.mesh.embed(1,[embedded_curve],2,face)
        end
        A.mesh.generate(1)
        lower=node_map()
        inherited_bound=maximum(curves) do curve
            _,_,connected=A.mesh.get_elements(1,curve)
            maximum(eachcol(reshape(only(connected),2,:))) do edge
                norm(collect(lower[edge[1]]).-collect(lower[edge[2]]))
            end
        end
        A.mesh.generate(2)
        points_by_uid=node_map()
        incidence,cells=surface_edges(face)
        @test all(points_by_uid[uid]==p for (uid,p) in lower)
        @test all(in((1,2)),values(incidence))
        @test !isempty(cells)
        @test length(points_by_uid)>length(lower)
        @test all(values(points_by_uid)) do p
            0<=p[1]<=1 && 0<=p[2]<=1 && abs(p[3]-warp*p[1]*p[2])<=1e-10
        end
        projected_area=sum(cells) do cell
            sum(begin
                p=points_by_uid[cell[i]];q=points_by_uid[cell[mod1(i+1,length(cell))]]
                p[1]*q[2]-q[1]*p[2]
            end for i in eachindex(cell))/2
        end
        @test projected_area≈1. atol=1e-10
        @test length(points_by_uid)-length(incidence)+length(cells)==1
        for curve in curves
            _,_,nodes=A.mesh.get_elements(1,curve)
            @test all(edge->get(incidence,minmax(edge...),0)==1,
                eachcol(reshape(only(nodes),2,:)))
        end
        if embedded
            _,_,nodes=A.mesh.get_elements(1,embedded_curve)
            @test all(edge->get(incidence,minmax(edge...),0)==2,
                eachcol(reshape(only(nodes),2,:)))
        end
        if !recombine
            limit=bounded ? .35 : inherited_bound
            @test all(keys(incidence)) do (a,b)
                norm(collect(points_by_uid[a]).-collect(points_by_uid[b]))<=limit+64eps(Float64)*max(1.,limit)
            end
        end
        @test A.option("Mesh.Algorithm")==algorithm
        return length(points_by_uid),length(cells)
    finally
        A.finalize()
    end
end

function occ_sector()
    A.initialize()
    try
        path=joinpath(mktempdir(),"sector.geo")
        write(path,"""
SetFactory("OpenCASCADE");
General.Terminal=0; Mesh.Renumber=0; Mesh.Algorithm=1;
Mesh.MeshSizeMin=.2; Mesh.MeshSizeMax=.2; Mesh.Smoothing=0;
Cylinder(1)={0,0,0,0,0,.4,.3,Pi/2};
""")
        A.open_geo!(path)
        A.mesh.generate(1)
        lower=node_map()
        A.mesh.generate(2)
        points=node_map()
        @test all(points[uid]==point for (uid,point) in lower)
        face=first(tag for (tag,kind) in A.CURRENT[].surface_types if kind===:cylinder)
        tags,_,_=A.mesh.get_nodes(2,face,true)
        @test !isempty(tags)
        @test all(uid->abs(points[uid][1]^2+points[uid][2]^2-.09)<1e-11,tags)
        @test all(uid->-1e-12<=points[uid][3]<=.4+1e-12,tags)
        incidence,cells=surface_edges(face)
        @test all(in((1,2)),values(incidence))
        @test length(tags)-length(incidence)+length(cells)==1
        @test all(keys(incidence)) do (a,b)
            norm(collect(points[a]).-collect(points[b]))<=.2+1e-12
        end
    finally
        A.finalize()
    end
end

function primitive_volume()
    A.initialize()
    try
        for (name,value) in (("General.Terminal",0),("Mesh.Renumber",0),
                ("Mesh.Algorithm",1),("Mesh.MeshSizeMin",.3),
                ("Mesh.MeshSizeMax",.3),("Mesh.Smoothing",0))
            A.option(name,value)
        end
        volume=A.model.add_cylinder(0.,0.,0.,0.,0.,.4,.3)
        A.mesh.generate(2)
        lower=node_map()
        faces=Dict(face=>A.mesh.get_elements(2,face) for (_,face) in A.model.get_entities(2))
        boundary=Set(Tuple(sort!(collect(cell))) for (_,_,nodes) in values(faces)
            for cell in eachcol(reshape(only(nodes),3,:)))
        A.option("Mesh.MeshOnlyEmpty",1)
        A.mesh.generate(3)
        points=node_map()
        @test all(points[uid]==point for (uid,point) in lower)
        # Generation's first-order constructors renew cell labels even when
        # MeshOnlyEmpty retains the geometric face. The oriented cells and
        # node identities survive; current cell labels bind those actual cells.
        @test all(sort!(collect(values(face_records(A.mesh.get_elements(2,face)))))==
                  sort!(collect(values(face_records(records)))) for (face,records) in faces)
        @test all(keys(faces)) do face
            all(face_records(A.mesh.get_elements(2,face))) do (uid,(kind,cell))
                A.mesh.get_element(uid)==(kind,collect(cell),2,face)
            end
        end
        @test Tessella.validate(A.mesh.get()).ok
        kinds,labels,nodes=A.mesh.get_elements(3,volume)
        @test kinds==Int32[4] && !isempty(only(labels))
        incidence=Dict{NTuple{3,UInt64},Int}()
        volumes=Float64[]
        for cell in eachcol(reshape(only(nodes),4,:))
            push!(volumes,Tessella.MeshTypes.tet_signed_volume((points[uid] for uid in cell)...))
            for slots in ((1,2,3),(1,4,2),(1,3,4),(2,4,3))
                face=Tuple(sort!(UInt64[cell[slot] for slot in slots]))
                incidence[face]=get(incidence,face,0)+1
            end
        end
        @test all(>(0),volumes)
        @test 0<sum(volumes)<=π*.3^2*.4
        @test Set(face for (face,count) in incidence if count==1)==boundary
        @test all(p->hypot(p[1],p[2])<=.3+1e-10 && -1e-10<=p[3]<=.4+1e-10,values(points))
    finally
        A.finalize()
    end
end

@testset "MeshAdapt owns local refinement on planar and ruled surfaces" begin
    @testset "planar selector 1" begin;patch(0.,1);end
    @testset "warped ruled Auto selects local adaptation" begin;patch(.2,2);end
    @testset "recombined embedded Curve" begin;patch(0.,1;recombine=true,embedded=true);end
    @testset "boundary sizes propagate without a global field" begin;patch(0.,1;bounded=false);end
    @testset "OCC chart retains its actual lower products" begin;occ_sector();end
    @testset "primitive volume consumes the retained adapted faces" begin;primitive_volume();end
end
end
