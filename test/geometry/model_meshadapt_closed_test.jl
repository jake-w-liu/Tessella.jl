module MeshAdaptClosedChecks
using Test, Tessella, LinearAlgebra
const A=Tessella.API
const M=Tessella.Model

function nodes()
    tags,xyz,_=A.mesh.get_nodes()
    Dict(uid=>Tuple(xyz[3i-2:3i]) for (i,uid) in enumerate(tags))
end

function closed_primitive(kind,h;algorithm=1,volume=false,extend=0,edit_parameter=false)
    A.initialize()
    try
        for (name,value) in (("General.Terminal",0),("Mesh.Renumber",0),
                ("Mesh.Algorithm",algorithm),("Mesh.Smoothing",0),
                ("Mesh.MeshSizeMin",h),("Mesh.MeshSizeMax",h),
                ("Mesh.MeshSizeExtendFromBoundary",extend))
            A.option(name,value)
        end
        solid=if kind===:sphere
            A.model.occ.add_sphere(0.,0.,0.,.3)
        elseif kind===:cone
            A.model.occ.add_cone(0.,0.,0.,0.,0.,.4,.3,0.)
        elseif kind===:box
            A.model.add_box(0.,0.,0.,.4,.4,.4)
        else
            A.model.occ.add_torus(0.,0.,0.,.4,.15)
        end
        A.model.occ.synchronize()
        face=kind===:box ? last(first(A.model.get_entities(2))) :
            first(t for (t,k) in A.CURRENT[].surface_types if k===kind)
        @test A.option("Mesh.Algorithm")==algorithm
        selected=deepcopy(A.CURRENT[])
        selected.meshing.algorithm_2d=algorithm
        @test M._normal_surface_mesh_algorithm(selected,face,"test")==1
        A.mesh.generate(1)
        lower=nodes()
        curves=Dict(tag=>(records[1],records[3]) for (_,tag) in A.model.get_entities(1)
            for records in (A.mesh.get_elements(1,tag),))
        A.mesh.generate(2)
        points=nodes()
        @test all(points[uid]==p for (uid,p) in lower)
        labels,connected=A.mesh.get_elements_by_type(2)
        cells=collect(eachcol(reshape(connected,3,length(labels))))
        @test !isempty(cells)
        incidence=Dict{Tuple{UInt64,UInt64},Int}()
        for cell in cells,i in 1:3
            edge=minmax(cell[i],cell[mod1(i+1,3)])
            incidence[edge]=get(incidence,edge,0)+1
        end
        @test all(==(2),values(incidence))
        used=Set(uid for cell in cells for uid in cell)
        @test length(used)-length(incidence)+length(cells)==(kind===:torus ? 0 : 2)
        @test all(keys(incidence)) do (a,b)
            norm(collect(points[a]).-collect(points[b]))<=h+128eps(Float64)
        end
        @test all(cells) do cell
            p,q,r=(collect(points[uid]) for uid in cell)
            norm(cross(q-p,r-p))>0
        end
        @test all(used) do uid
            x,y,z=points[uid];rho=hypot(x,y)
            if kind===:sphere
                isapprox(hypot(x,y,z),.3;atol=1e-11)
            elseif kind===:cone
                -1e-12<=z<=.4+1e-12 &&
                    (isapprox(rho,.3*(1-z/.4);atol=1e-11) ||
                     (abs(z)<=1e-12 && rho<=.3+1e-12))
            elseif kind===:box
                all(v->-1e-12<=v<=.4+1e-12,(x,y,z)) &&
                    any(v->abs(v)<=1e-12 || abs(v-.4)<=1e-12,(x,y,z))
            else
                isapprox(hypot(rho-.4,z),.15;atol=1e-11)
            end
        end
        boundary=Set(Tuple(sort!(collect(cell))) for cell in cells)
        area(cells,points)=sum(cells) do cell
            p,q,r=(collect(points[uid]) for uid in cell)
            norm(cross(q-p,r-p))/2
        end
        old_area=area(cells,points)
        if volume
            if edit_parameter
                uid=first(A.mesh.get_nodes(2,face,false)[1])
                xyz,uv,_,_=A.mesh.get_node(uid)
                uv=copy(uv);uv[1]+=2pi
                A.mesh.set_node(uid,xyz,uv)
                @test !haskey(A.CURRENT[].meshing.mesh_authority.surface_constructor_inputs,
                    (2,Int32(face)))
            end
            A.option("Mesh.MeshOnlyEmpty",1)
            A.mesh.generate(3)
            after=nodes()
            @test all(after[uid]==p for (uid,p) in points)
            @test all(curves) do (curve,(kinds,connected))
                current=A.mesh.get_elements(1,curve)
                current[1]==kinds && current[3]==connected
            end
            @test Tessella.validate(A.mesh.get()).ok
            kinds,tags,connected=A.mesh.get_elements(3,solid)
            @test kinds==Int32[4] && !isempty(only(tags))
            face_counts=Dict{NTuple{3,UInt64},Int}()
            volumes=Float64[]
            for cell in eachcol(reshape(only(connected),4,:))
                push!(volumes,Tessella.MeshTypes.tet_signed_volume((after[uid] for uid in cell)...))
                for slots in ((1,2,3),(1,4,2),(1,3,4),(2,4,3))
                    key=Tuple(sort!(UInt64[cell[i] for i in slots]))
                    face_counts[key]=get(face_counts,key,0)+1
                end
            end
            @test all(>(0),volumes)
            @test sum(volumes)>0
            @test all(in((1,2)),values(face_counts))
            # Volume recovery can retriangulate a coplanar Surface region.
            # It must publish those actual cells with valid label bindings;
            # the retained node coordinates, Curve edges and surface area
            # remain fixed even when a flat cap's diagonals change.
            published_labels,published_nodes=A.mesh.get_elements_by_type(2)
            published=collect(eachcol(reshape(published_nodes,3,length(published_labels))))
            @test Set(key for (key,count) in face_counts if count==1)==
                Set(Tuple(sort!(collect(cell))) for cell in published)
            if kind===:sphere
                @test Set(Tuple(sort!(collect(cell))) for cell in published)==boundary
            end
            @test area(published,after)≈old_area atol=1e-11 rtol=1e-11
            @test all(A.model.get_entities(2)) do (_,surface)
                kinds,labels,connected=A.mesh.get_elements(2,surface)
                all(zip(kinds,labels,connected)) do (kind,uids,flat)
                    all(zip(uids,eachcol(reshape(flat,3,:)))) do (uid,cell)
                        A.mesh.get_element(uid)==(kind,collect(cell),2,surface)
                    end
                end
            end
            @test all(A.model.get_entities(2)) do (_,surface)
                M._model_mesh_product(A.CURRENT[].meshing.mesh_authority,2,surface;
                    producer_order=true,first_order=true,caller="actual recovered source")!==nothing
            end
        end
        return length(cells)
    finally
        A.finalize()
    end
end

function second_order_torus()
    A.initialize()
    try
        for (name,value) in (("General.Terminal",0),("Mesh.Renumber",0),
                ("Mesh.Algorithm",1),("Mesh.Smoothing",0),
                ("Mesh.MeshSizeMin",.5),("Mesh.MeshSizeMax",.5),
                ("Mesh.MeshSizeExtendFromBoundary",0))
            A.option(name,value)
        end
        solid=A.model.occ.add_torus(0.,0.,0.,.4,.15)
        A.model.occ.synchronize();A.mesh.generate(2)
        before=nodes()
        A.option("Mesh.MeshOnlyEmpty",1);A.option("Mesh.ElementOrder",2)
        A.option("Mesh.SecondOrderLinear",1)
        for iteration in 1:2
            iteration==1 || A.mesh.clear([(3,solid)])
            A.mesh.generate(3)
            @test Tessella.validate(A.mesh.get()).ok
            after=nodes()
            @test all(after[uid]==point for (uid,point) in before)
            @test A.mesh.get_elements(3,solid)[1]==Int32[11]
            # Recovery-born primary boundary vertices must belong to the
            # actual source vector after elevation, so clearing the Region
            # can consume the retained Face again through its linear view.
            model=A.CURRENT[]
            face=first(t for (t,kind) in model.surface_types if kind===:torus)
            source=M._model_mesh_product(model.meshing.mesh_authority,2,face;
                producer_order=true,first_order=true,caller="retained P2 torus")
            @test source!==nothing && any(block->!isempty(block.nodes),source.blocks)
        end
    finally
        A.finalize()
    end
end

@testset "Surface quotient vertex links" begin
    tetra=NTuple{3,Int32}[(1,2,3),(1,4,2),(1,3,4),(2,4,3)]
    coordinates=[(0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.),
        (-1.,0.,0.),(0.,-1.,0.),(0.,0.,-1.)]
    empty_rim=Set{NTuple{2,Int32}}()
    @test M._metric_surface_quotient_ok(tetra,coordinates,empty_rim)
    pinched=vcat(tetra,NTuple{3,Int32}[(1,5,6),(1,7,5),(1,6,7),(5,7,6)])
    @test !M._metric_surface_quotient_ok(pinched,coordinates,empty_rim)
    @test !M._metric_surface_quotient_ok(vcat(tetra,[tetra[1]]),coordinates,empty_rim)
end

@testset "MeshAdapt closed periodic and pole charts" begin
    counts=Dict{Symbol,Int}()
    for kind in (:sphere,:cone,:torus)
        @testset "$kind" begin
            coarse=closed_primitive(kind,.5;volume=true)
            fine=closed_primitive(kind,.2)
            @test fine>coarse
            counts[kind]=coarse
        end
    end
    @test closed_primitive(:sphere,.5;algorithm=2)==counts[:sphere]
    # Pole-row preimages have zero physical span, which supplies no initial
    # inherited size. Ordinary positive incident edges still extend sizes.
    closed_primitive(:cone,.5;extend=1)
end

@testset "MeshAdapt Box consumes its actual faces" begin
    closed_primitive(:box,.3;volume=true)
end

@testset "MeshAdapt retained second-order recovery" begin
    second_order_torus()
end

@testset "MeshAdapt adopts an edited actual source" begin
    closed_primitive(:torus,.5;volume=true,edit_parameter=true)
end
end
