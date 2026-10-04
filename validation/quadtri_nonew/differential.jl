#!/usr/bin/env julia
# Gmsh supplies an independent, allocation-dependent admissible mesh. Production
# never invokes it, and this gate never pins a random pointer-selected topology.
using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.Elements: MixedMesh, ElementBlock, msh_spec
using Tessella.MeshTypes: nnodes

include(joinpath(@__DIR__,"..","..","test","geometry","quadtri_nonew_certificates.jl"))
const QTNN=QuadTriNoNewCertificates
include(joinpath(@__DIR__,"..","..","test","geometry","quadtri_nonew_triangle_certificates.jl"))
const QTNT=QuadTriNoNewTriangleCertificates
binding=get(ENV,"GMSH_JULIA_API","")
isfile(binding) || error("set GMSH_JULIA_API to pinned Gmsh4.15.2 gmsh.jl")
include(binding)
gmsh.GMSH_API_VERSION=="4.15.2" || error("Gmsh4.15.2 binding required")

function oracle_volume(entity=1)
    node_tags,xyz,_=gmsh.model.mesh.getNodes()
    coordinates=reshape(xyz,3,:)
    source=Dict(tag=>index for (index,tag) in enumerate(node_tags))
    types,tags,nodes=gmsh.model.mesh.getElements(3,entity)
    used=sort!(unique(vcat(nodes...)))
    remap=Dict(tag=>Int32(index) for (index,tag) in enumerate(used))
    blocks=ElementBlock[]
    for (type,cell_nodes) in zip(types,nodes)
        width=gmsh.model.mesh.getElementProperties(type)[4]
        connectivity=reshape(Int32[remap[tag] for tag in cell_nodes],width,:)
        push!(blocks,ElementBlock(type,connectivity))
    end
    return MixedMesh(coordinates[:,[source[tag] for tag in used]],blocks),vcat(tags...)
end

function native_api(path,action)
    api=Tessella.API
    api.initialize()
    try
        api.open_geo!(path;mesh_dim=0)
        return action(api)
    finally
        api.finalize()
    end
end

# Read the actual public cell product, including the legacy simplex P2 overlay.
# mesh.get() deliberately continues to expose the primary Mesh in that route.
function triangle_api_volume(api,entity=1)
    tags,xyz,_=api.mesh.get_nodes();coordinates=reshape(xyz,3,:)
    positions=Dict(tag=>index for (index,tag) in enumerate(tags))
    types,_,families=api.mesh.get_elements(3,entity)
    used=sort!(unique(vcat(families...)))
    remap=Dict(tag=>Int32(index) for (index,tag) in enumerate(used))
    blocks=ElementBlock[]
    for (type,nodes) in zip(types,families)
        width=msh_spec(type).nnodes
        push!(blocks,ElementBlock(type,reshape(Int32[remap[tag] for tag in nodes],width,:)))
    end
    return MixedMesh(coordinates[:,[positions[tag] for tag in used]],blocks)
end

function quadratic_volume_from_jacobians(mesh,interface;expected=1.)
    total=0.
    for block in mesh.blocks
        if interface===gmsh.model.mesh
            locations,weights=interface.getIntegrationPoints(block.msh,"Gauss4")
            _,determinants,_=interface.getJacobians(block.msh,locations)
        else
            locations,weights=interface.get_integration_points(block.msh,"Gauss4")
            _,determinants,_=interface.get_jacobians(block.msh,locations)
        end
        all(d->isfinite(d)&&d>0,determinants) || error("nonpositive P2 quadrature Jacobian")
        length(determinants)==length(weights)*size(block.nodes,2) || error("P2 quadrature shape")
        total+=sum(determinants .* repeat(weights,size(block.nodes,2)))
    end
    abs(total-expected)<=2e-11 || error("P2 integrated swept volume $total differs from $expected")
    return total
end

function oracle_surface_faces(volume)
    lookup=Dict(QTNN.point(volume,index)=>index for index in 1:nnodes(volume))
    tags,xyz,_=gmsh.model.mesh.getNodes()
    coordinates=reshape(xyz,3,:)
    points=Dict(tag=>Tuple(coordinates[:,index]) for (index,tag) in enumerate(tags))
    result=Set{Tuple}()
    types,_,families=gmsh.model.mesh.getElements(2)
    for (type,nodes) in zip(types,families)
        width=gmsh.model.mesh.getElementProperties(type)[4]
        for cell in eachcol(reshape(nodes,width,:))
            push!(result,QTNN.key(lookup[points[node]] for node in cell))
        end
    end
    return result
end

function same_grid_and_centroid(native,oracle,fixture)
    nc=QTNN.certify_fixture(native,fixture)
    oc=QTNN.certify_fixture(oracle,fixture)
    nnodes(native)==nnodes(oracle) || error("$(fixture.name): node counts differ")
    for i in eachindex(nc.grid_ids)
        p=QTNN.point(native,nc.grid_ids[i]);q=QTNN.point(oracle,oc.grid_ids[i])
        maximum(abs.(p.-q))<=2e-11 || error("$(fixture.name): swept coordinates differ")
    end
    if fixture.laterals
        p=QTNN.point(native,only(nc.extras));q=QTNN.point(oracle,only(oc.extras))
        maximum(abs.(p.-q))<=2e-11 || error("$(fixture.name): final centroid differs")
        nc.family_counts==oc.family_counts || error("$(fixture.name): recombined families differ")
    end
    # Free lateral diagonals and even family counts can legitimately differ.
    QTNN.surface_faces(native_execution[].mesh_parts,native)==Set(keys(nc.boundary)) ||
        error("$(fixture.name): native volume/surface complex differs")
    oracle_surface_faces(oracle)==Set(keys(oc.boundary)) ||
        error("$(fixture.name): oracle volume/surface complex differs")
    return nc,oc
end

const native_execution=Ref{Any}(nothing)
function initialize_oracle()
    gmsh.initialize(["-nopopup"],false,false)
    gmsh.option.setNumber("General.Terminal",0)
    gmsh.option.setNumber("Mesh.ElementOrder",1)
    gmsh.option.setNumber("Mesh.SecondOrderLinear",0)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete",0)
end
initialize_oracle()
oracle_only=get(ENV,"QUADTRI_NONEW_ORACLE_ONLY","")=="1"
completed=Ref(0);oracle_samples=Ref(0);quadratic_cases=Ref(0)
height_errors=Ref(0);isolated_cases=Ref(0)
native_helical_cases=Ref(0);helical_oracle_gaps=Ref(0)
triangle_cases=Ref(0);triangle_quadratic_cases=Ref(0);triangle_isolated_cases=Ref(0)
triangle_parameter_provenance_gaps=Ref(0)
try
    startswith(gmsh.option.getString("General.Version"),"4.15.2") ||
        error("Gmsh4.15.2 runtime required")
    selected=get(ENV,"QUADTRI_NONEW_CASE","")
    for fixture in QTNN.fixtures()
        isempty(selected) || occursin(selected,fixture.name) || continue
        mktempdir() do directory
            path=joinpath(directory,"$(fixture.name).geo")
            write(path,fixture.source*"Mesh 3;\n")
            # Repeats expose allocator choices; validity does not depend on
            # seeing every possible outcome during this bounded gate.
            repeat_count=startswith(fixture.name,"uniform") ? 3 : 1
            native=oracle_only ? nothing : QTNN.execute(fixture.source)
            native_execution[]=native
            nvolume=oracle_only ? nothing : geo_entity_mesh(native,3,1)
            outcomes=Set{Tuple}()
            for repeat in 1:repeat_count
                gmsh.clear();gmsh.open(path)
                oracle,element_tags=oracle_volume()
                oc=QTNN.certify_fixture(oracle,fixture)
                oracle_surface_faces(oracle)==Set(keys(oc.boundary)) ||
                    error("$(fixture.name): oracle surface boundary mismatch")
                qualities=gmsh.model.mesh.getElementQualities(element_tags,"minDetJac")
                all(q->isfinite(q)&&q>0,qualities) || error("$(fixture.name): invalid oracle Jacobian")
                oracle_only || same_grid_and_centroid(nvolume,oracle,fixture)
                push!(outcomes,Tuple(sort!(collect(oc.family_counts))))
                oracle_samples[]+=1
            end
            if !oracle_only
                validate(native.mesh).ok || error("$(fixture.name): invalid native global mesh")
                again=QTNN.execute(fixture.source)
                mixed_crc(again.mesh)==mixed_crc(native.mesh) ||
                    error("$(fixture.name): nondeterministic native CRC")
                native_api(path,api->begin
                    actual=api.mesh.generate(3)
                    QTNN.certify_fixture(actual,fixture)
                    QTNN.typed_signature(actual)==QTNN.typed_signature(nvolume) ||
                        error("$(fixture.name): API/GEO volume differs")
                    isempty(api.mesh.get_elements(2)[1]) ||
                        error("legacy API3 cache lower-cell contract changed")
                    length(api.mesh.get_nodes(3,1,true)[1])==nnodes(actual) ||
                        error("$(fixture.name): API volume closure nodes differ")
                    sum(length,api.mesh.get_elements(3,1)[2])==
                        QTNN.certify_complex(actual).ncell || error("API3 cell query differs")
                end)
            end
            completed[]+=1
            println("QUADTRI_NONEW_CASE_OK name=$(fixture.name) oracle_outcomes=$(collect(outcomes))")
        end
    end
    if isempty(selected)
        for fixture in filter(f->startswith(f.name,"uniform"),QTNN.fixtures())
            mktempdir() do directory
                path=joinpath(directory,"quadratic.geo");write(path,fixture.source)
                gmsh.clear();gmsh.open(path);gmsh.model.mesh.generate(3)
                linear,_=oracle_volume()
                gmsh.model.mesh.setOrder(2)
                elevated,tags=oracle_volume()
                QTNN.certify_quadratic(elevated,linear)
                quadratic_volume_from_jacobians(elevated,gmsh.model.mesh)
                all(>(0),gmsh.model.mesh.getElementQualities(tags,"minDetJac")) ||
                    error("oracle P2 quality")
                if !oracle_only
                    native_api(path,api->begin
                        first=api.mesh.generate(3)
                        api.mesh.set_order(2);quadratic=api.mesh.get()
                        QTNN.certify_quadratic(quadratic,first)
                        quadratic_volume_from_jacobians(quadratic,api.mesh)
                        api.mesh.set_order(1)
                        QTNN.typed_signature(api.mesh.get())==QTNN.typed_signature(first) ||
                            error("P2/P1 roundtrip cell set changed")
                    end)
                end
                quadratic_cases[]+=1;oracle_samples[]+=1
                println("QUADTRI_NONEW_QUADRATIC_OK name=$(fixture.name)")
            end
        end
        for coincident in (false,true)
            mktempdir() do directory
                path=joinpath(directory,"isolated.geo");write(path,QTNN.paired_source())
                gmsh.clear();gmsh.open(path)
                if coincident
                    for (_,tag) in gmsh.model.getEntities(0)
                        x,y,z=gmsh.model.getValue(0,tag,Float64[])
                        x>=3 || continue
                        gmsh.model.setCoordinates(tag,x-3,y,z)
                    end
                end
                gmsh.model.mesh.generate(3)
                entities=last.(gmsh.model.getEntities(3))
                length(entities)==2 || error("isolated oracle region count")
                volumes=[oracle_volume(tag)[1] for tag in entities]
                fixtures=(QTNN.fixture("uniform3_false","Layers{3}",[0.,1/3,2/3,1.],false),
                          QTNN.fixture("uniform3_true","Layers{3}",[0.,1/3,2/3,1.],true))
                for (i,volume) in enumerate(volumes)
                    shifted=i==2 && !coincident ? MixedMesh(volume.coords.-[3.,0,0],volume.blocks) : volume
                    QTNN.certify_fixture(shifted,fixtures[i])
                end
                tags=gmsh.model.mesh.getNodes()[1]
                length(tags)==33 || error("isolated oracle nodes $(length(tags)) differ from 33")
                a=gmsh.model.mesh.getNodes(3,entities[1],true)[1]
                b=gmsh.model.mesh.getNodes(3,entities[2],true)[1]
                isempty(intersect(a,b)) || error("independent oracle topology welded")
                if !oracle_only
                    native_api(path,api->begin
                        if coincident
                            for (tag,p) in collect(api.CURRENT[].points)
                                p[1]>=3 || continue
                                api.model.set_coordinates(tag,p[1]-3,p[2],p[3])
                            end
                        end
                        mesh=api.mesh.generate(3)
                        nnodes(mesh)==33 || error("isolated native nodes differ")
                        a=api.mesh.get_nodes(3,entities[1],true)[1]
                        b=api.mesh.get_nodes(3,entities[2],true)[1]
                        length(a)==16 && length(b)==17 && isempty(intersect(a,b)) ||
                            error("isolated native ownership differs")
                        validate(mesh).ok || error("isolated native mesh invalid")
                    end)
                end
                isolated_cases[]+=1;oracle_samples[]+=1
                println("QUADTRI_NONEW_ISOLATED_OK coincident=$coincident")
            end
        end
        # Large-angle Gmsh extrusion takes a different route from this bounded
        # native planner. Record actual outcomes without claiming cell parity.
        for (pitch,laterals) in ((-3.,false),(-3.,true),(3.,false),(3.,true),
                                 (-.1,true),(.1,true))
            fixture=QTNN.helical_fixture("helical_oracle_gap",pitch,laterals)
            mktempdir() do directory
                path=joinpath(directory,"helix.geo");write(path,fixture.source)
                # After a fatal FindDiagonalEdgeIndices error, clear() alone
                # leaves upstream parser state unusable for later fixtures.
                gmsh.finalize();initialize_oracle();gmsh.open(path)
                failure=try gmsh.model.mesh.generate(3);nothing catch e;e end
                source_types,source_tags,_=gmsh.model.mesh.getElements(2,1)
                source_types==[3] && length(only(source_tags))==1 ||
                    error("helical oracle source quad changed")
                types,tags,_=gmsh.model.mesh.getElements(3,1)
                if !laterals
                    failure!==nothing && occursin("FindDiagonalEdgeIndices",sprint(showerror,failure)) ||
                        error("unexpected free helical oracle outcome")
                    isempty(types) || error("failed free helical oracle emitted volume cells")
                    println("QUADTRI_NONEW_HELIX_ORACLE_GAP pitch=$pitch laterals=$laterals kind=generation_error volume_cells=0")
                else
                    failure===nothing || error("unexpected recombined helical oracle error: $failure")
                    oracle,element_tags=oracle_volume()
                    all(>(0),gmsh.model.mesh.getElementQualities(element_tags,"minDetJac")) ||
                        error("large-angle oracle local quality changed")
                    # A positive sampled quality alone is not a complete cell,
                    # face-complex or global disjointness certificate.
                    local_result=try
                        QTNN.certify_complex(oracle);"certified"
                    catch e
                        sprint(showerror,e)
                    end
                    parity_failure=try QTNN.certify_fixture(oracle,fixture);nothing catch e;e end
                    parity_failure!==nothing || error("large-angle oracle gap disappeared; reassess parity")
                    families=Tuple((Int(type),length(family)) for (type,family) in zip(types,tags))
                    println("QUADTRI_NONEW_HELIX_ORACLE_GAP pitch=$pitch laterals=$laterals kind=different_cell_route volume_nodes=$(nnodes(oracle)) families=$families local_complex=$(repr(local_result))")
                end
                if !oracle_only
                    if abs(pitch)==3
                        native=QTNN.execute(fixture.source)
                        nvolume=geo_entity_mesh(native,3,1)
                        QTNN.certify_fixture(nvolume,fixture)
                        native_api(path,api->begin
                            actual=api.mesh.generate(3)
                            QTNN.certify_fixture(actual,fixture)
                            QTNN.typed_signature(actual)==QTNN.typed_signature(nvolume) ||
                                error("certified native helical API/GEO cells differ")
                        end)
                        native_helical_cases[]+=1
                    else
                        failure=try QTNN.execute(fixture.source);nothing catch e;e end
                        failure isa ArgumentError && occursin("global P1 cell-hull disjointness is not certified",sprint(showerror,failure)) ||
                            error("native overlapping helix was not rejected precisely")
                    end
                end
                helical_oracle_gaps[]+=1;oracle_samples[]+=1
            end
        end
        # Literal layer heights do not move the CAD cap, which stays at u=1.
        # These upstream errors are recorded separately from valid fixtures.
        for height in (.5,2.),laterals in (false,true)
            fixture=QTNN.fixture("terminal_height","Layers{{1},{$height}}",[0.,height],laterals)
            mktempdir() do directory
                path=joinpath(directory,"invalid_height.geo");write(path,fixture.source)
                gmsh.finalize();initialize_oracle();gmsh.open(path)
                failure=try gmsh.model.mesh.generate(3);nothing catch e;e end
                failure!==nothing && occursin("Could not find extruded node",sprint(showerror,failure)) ||
                    error("unexpected non-unit height oracle behavior")
                isempty(gmsh.model.mesh.getElements(3,1)[1]) || error("invalid height created volume cells")
                if !oracle_only
                    failure=try QTNN.execute(fixture.source);nothing catch e;e end
                    failure isa ArgumentError && occursin("normalized final layer height of 1.0",sprint(showerror,failure)) ||
                        error("native non-unit height was not rejected precisely")
                end
                height_errors[]+=1;oracle_samples[]+=1
            end
        end
    end
    for fixture in QTNT.fixtures()
        isempty(selected) || occursin(selected,fixture.name) || continue
        mktempdir() do directory
            path=joinpath(directory,"$(fixture.name).geo")
            write(path,fixture.source*"Mesh 3;\n")
            gmsh.clear();gmsh.open(path)
            oracle,element_tags=oracle_volume()
            oc=QTNT.certify(oracle,fixture)
            source_types,_,source_nodes=gmsh.model.mesh.getElements(2,1)
            source_types==[2] && length(only(source_nodes))==3 || error("triangle oracle source cell changed")
            source_points=Tuple(Tuple(gmsh.model.mesh.getNode(node)[1]) for node in only(source_nodes))
            q=map(QTNT.exact_point,source_points)
            sign(QTNT.determinant(q[2].-q[1],q[3].-q[1],(0,0,1)))==fixture.winding ||
                error("triangle oracle source winding changed")
            oracle_surface_faces(oracle)==Set(keys(oc.boundary)) ||
                error("$(fixture.name): triangle oracle surface boundary mismatch")
            all(>(0),gmsh.model.mesh.getElementQualities(element_tags,"minDetJac")) ||
                error("$(fixture.name): triangle oracle Jacobian quality")
            if !oracle_only
                execution=QTNN.execute(fixture.source)
                native=geo_entity_mesh(execution,3,1);nc=QTNT.certify(native,fixture)
                QTNT.certify_cap_winding(execution.mesh_parts,native,fixture,nc)
                QTNN.surface_faces(execution.mesh_parts,native)==Set(keys(nc.boundary)) ||
                    error("$(fixture.name): native triangle surface boundary mismatch")
                nnodes(native)==nnodes(oracle) || error("triangle volume node counts differ")
                for i in eachindex(nc.grid_ids)
                    p=QTNN.point(native,nc.grid_ids[i]);q=QTNN.point(oracle,oc.grid_ids[i])
                    maximum(abs.(p.-q))<=2e-11 || error("triangle swept coordinates differ")
                end
                nc.family_counts==oc.family_counts || error("triangle volume family counts differ")
                QTNT.crc(QTNN.execute(fixture.source).mesh)==QTNT.crc(execution.mesh) ||
                    error("nondeterministic native triangle CRC")
                native_api(path,api->begin
                    actual=api.mesh.generate(3)
                    QTNT.certify(actual,fixture)
                    QTNN.typed_signature(actual)==QTNN.typed_signature(native) ||
                        error("triangle API/GEO volume differs")
                    isempty(api.mesh.get_elements(2)[1]) || error("legacy API3 lower-cell contract changed")
                    length(api.mesh.get_nodes(3,1,true)[1])==nnodes(actual) ||
                        error("triangle API volume closure differs")
                    isempty(api.mesh.get_nodes(3,1,false)[1]) || error("triangle API invented body nodes")
                end)
            end
            triangle_cases[]+=1;oracle_samples[]+=1
            println("QUADTRI_NONEW_TRIANGLE_OK name=$(fixture.name) families=$(oc.family_counts)")
        end
    end
    if isempty(selected) || selected=="triangle"
        for laterals in (false,true),winding in (-1,1)
            fixture=QTNT.fixture("triangle_p2","Layers{3}",[0.,1/3,2/3,1.],laterals;winding)
            mktempdir() do directory
                path=joinpath(directory,"triangle_p2.geo");write(path,fixture.source)
                gmsh.clear();gmsh.open(path);gmsh.model.mesh.generate(3)
                linear,_=oracle_volume();QTNT.certify(linear,fixture)
                gmsh.model.mesh.setOrder(2);quadratic,tags=oracle_volume()
                QTNN.certify_quadratic(quadratic,linear)
                quadratic_volume_from_jacobians(quadratic,gmsh.model.mesh;expected=.5)
                all(>(0),gmsh.model.mesh.getElementQualities(tags,"minDetJac")) ||
                    error("triangle oracle P2 quality")
                empty_face_parameters=0
                oracle_lateral_count=0
                for (_,surface) in gmsh.model.getEntities(2)
                    owned,_,parameters=gmsh.model.mesh.getNodes(2,surface,false,true)
                    isempty(owned) && continue
                    oracle_lateral_count+=1
                    length(owned)==5 && length(parameters)==(laterals ? 4 : 10) ||
                        error("oracle lateral ownership/parameter dimensions changed")
                    boundary_nodes,_,boundary_parameters=gmsh.model.mesh.getNodes(2,surface,true,true)
                    length(boundary_nodes)==21 && length(boundary_parameters)==(laterals ? 36 : 42) ||
                        error("oracle partial boundary parameter provenance changed")
                    for node in owned
                        p,uv,dim,owner=gmsh.model.mesh.getNode(node)
                        dim==2 && owner==surface || error("oracle P2 lateral ownership")
                        if isempty(uv)
                            empty_face_parameters+=1
                            # The stored absence is distinct from the complete
                            # geometric inverse available on the same surface.
                            computed=gmsh.model.getParametrization(2,surface,p)
                            maximum(abs.(gmsh.model.getValue(2,surface,computed).-p))<=2e-11 ||
                                error("oracle geometric inverse roundtrip")
                        end
                    end
                end
                oracle_lateral_count==3 || error("oracle P2 lateral ownership missing")
                empty_face_parameters==(laterals ? 9 : 0) || error("oracle face-center parameter provenance")
                if laterals
                    triangle_parameter_provenance_gaps[]+=1
                    println("QUADTRI_NONEW_TRIANGLE_PARAMETER_PROVENANCE_GAP winding=$winding oracle_empty_face_nodes=$empty_face_parameters native_contract=computed_uv")
                end
                if !oracle_only
                    native_api(path,api->begin
                        api.mesh.generate(3)
                        first=triangle_api_volume(api)
                        api.mesh.set_order(2);actual=triangle_api_volume(api)
                        QTNN.certify_quadratic(actual,first)
                        quadratic_volume_from_jacobians(actual,api.mesh;expected=.5)
                        node_tags=api.mesh.get_nodes()[1]
                        for node in node_tags
                            p,uv,dim,entity=api.mesh.get_node(node)
                            dim==2 || continue
                            length(uv)==2 || error("native computed lateral UV missing")
                            maximum(abs.(api.model.get_value(2,entity,uv).-p))<=2e-11 ||
                                error("native computed lateral UV does not reproduce actual node")
                        end
                        lateral_count=0
                        for (_,surface) in api.model.get_entities(2)
                            # The source and cap have no interior P2 nodes.
                            # All three actual laterals own five supports.
                            owned,_,parameters=api.mesh.get_nodes(2,surface,false,true)
                            isempty(owned) && continue
                            lateral_count+=1
                            length(owned)==5 && length(parameters)==10 ||
                                error("native computed lateral parameter dimensions")
                            boundary_nodes,coordinates,uv=api.mesh.get_nodes(2,surface,true,true)
                            length(boundary_nodes)==21 && length(uv)==42 ||
                                error("native computed boundary parameter dimensions")
                            maximum(abs.(api.model.get_value(2,surface,uv).-coordinates))<=2e-11 ||
                                error("native computed boundary UV roundtrip")
                        end
                        lateral_count==3 || error("native P2 lateral ownership missing")
                        api.mesh.set_order(1)
                        QTNN.typed_signature(triangle_api_volume(api))==QTNN.typed_signature(first) ||
                            error("triangle P2/P1 actual cell roundtrip differs")
                    end)
                end
                triangle_quadratic_cases[]+=1;oracle_samples[]+=1
                println("QUADTRI_NONEW_TRIANGLE_P2_OK laterals=$laterals winding=$winding")
            end
        end
        for coincident in (false,true)
            mktempdir() do directory
                path=joinpath(directory,"triangle_isolated.geo");write(path,QTNT.paired_source())
                gmsh.clear();gmsh.open(path)
                if coincident
                    for (_,tag) in gmsh.model.getEntities(0)
                        x,y,z=gmsh.model.getValue(0,tag,Float64[])
                        x>=3 || continue
                        gmsh.model.setCoordinates(tag,x-3,y,z)
                    end
                end
                gmsh.model.mesh.generate(3)
                entities=last.(gmsh.model.getEntities(3));length(entities)==2 || error("triangle region count")
                for (index,entity) in enumerate(entities)
                    volume=oracle_volume(entity)[1]
                    if index==2 && !coincident
                        volume=MixedMesh(volume.coords.-[3.,0,0],volume.blocks)
                    end
                    f=QTNT.fixture("triangle_isolated","Layers{3}",[0.,1/3,2/3,1.],index==2)
                    QTNT.certify(volume,f)
                end
                a=gmsh.model.mesh.getNodes(3,entities[1],true)[1]
                b=gmsh.model.mesh.getNodes(3,entities[2],true)[1]
                length(a)==length(b)==12 && isempty(intersect(a,b)) || error("triangle oracle identities welded")
                length(gmsh.model.mesh.getNodes()[1])==24 || error("triangle oracle global node count")
                if !oracle_only
                    native_api(path,api->begin
                        if coincident
                            for (tag,p) in collect(api.CURRENT[].points)
                                p[1]>=3 || continue
                                api.model.set_coordinates(tag,p[1]-3,p[2],p[3])
                            end
                        end
                        generated=api.mesh.generate(3)
                        a=api.mesh.get_nodes(3,entities[1],true)[1]
                        b=api.mesh.get_nodes(3,entities[2],true)[1]
                        length(a)==length(b)==12 && isempty(intersect(a,b)) || error("triangle native identities welded")
                        nnodes(generated)==24 || error("triangle native global node count")
                    end)
                end
                triangle_isolated_cases[]+=1;oracle_samples[]+=1
                println("QUADTRI_NONEW_TRIANGLE_ISOLATED_OK coincident=$coincident")
            end
        end
    end
    println("QUADTRI_NONEW_DIFFERENTIAL_OK gmsh=4.15.2 p1_cases=$(completed[]) p2_cases=$(quadratic_cases[]) native_helical_cases=$(native_helical_cases[]) helical_oracle_gaps=$(helical_oracle_gaps[]) height_errors=$(height_errors[]) isolated_cases=$(isolated_cases[]) triangle_cases=$(triangle_cases[]) triangle_p2_cases=$(triangle_quadratic_cases[]) triangle_isolated_cases=$(triangle_isolated_cases[]) triangle_parameter_provenance_gaps=$(triangle_parameter_provenance_gaps[]) oracle_samples=$(oracle_samples[]) oracle_only=$oracle_only")
finally
    gmsh.finalize()
end
