#!/usr/bin/env julia
# P6: Tessella translation/rotation-periodic node correspondence vs Gmsh 4.15.2.

using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.MeshTypes: Mesh,validate,mesh_crc

function find_gmsh_api()
    explicit=get(ENV,"GMSH_JULIA_API","")
    !isempty(explicit) && isfile(explicit) && return explicit
    executable=Sys.which("gmsh")
    executable===nothing && error("gmsh is not on PATH")
    prefix=dirname(dirname(realpath(executable)))
    for path in (joinpath(prefix,"lib","gmsh.jl"),
                 "/opt/homebrew/opt/gmsh/lib/gmsh.jl")
        isfile(path) && return path
    end
    error("could not locate gmsh.jl")
end

include(find_gmsh_api())
const GEO=joinpath(@__DIR__,"periodic_translation.geo")
const NATIVE_GEO=joinpath(@__DIR__,"periodic_native.geo")
const TWO_DIRECTION_GEO=joinpath(@__DIR__,"periodic_two_direction.geo")
gmsh.initialize(["gmsh","-v","0"])
try
    startswith(gmsh.GMSH_API_VERSION,"4.15.2") || error(
        "periodic differential requires Gmsh 4.15.2, got $(gmsh.GMSH_API_VERSION)")
    gmsh.open(GEO)
    gmsh.model.mesh.generate(2)
    master_entity,slave_tags,master_tags,affine=
        gmsh.model.mesh.getPeriodicNodes(1,2)
    master_entity==4 || error("Gmsh periodic master curve is $master_entity, expected 4")
    length(slave_tags)==length(master_tags)==5 || error(
        "Gmsh periodic curve pair count is $(length(slave_tags)), expected 5")
    expected_affine=[1.0,0.0,0.0,1.0,
                     0.0,1.0,0.0,0.0,
                     0.0,0.0,1.0,0.0,
                     0.0,0.0,0.0,1.0]
    length(affine)==16 || error("Gmsh periodic affine transform is not 4×4")
    maximum(abs.(affine.-expected_affine))<=1e-14 || error(
        "Gmsh periodic affine transform is $affine")

    node_tags,node_coordinates,_=gmsh.model.mesh.getNodes()
    coordinates=Dict{Int,NTuple{3,Float64}}()
    for (i,tag) in enumerate(node_tags)
        coordinates[Int(tag)]=(node_coordinates[3i-2],node_coordinates[3i-1],
                               node_coordinates[3i])
    end
    all(tag->haskey(coordinates,Int(tag)),master_tags) || error(
        "Gmsh omitted a periodic master node from getNodes")
    all(tag->haskey(coordinates,Int(tag)),slave_tags) || error(
        "Gmsh omitted a periodic slave node from getNodes")
    order=sortperm(eachindex(master_tags);by=i->coordinates[Int(master_tags[i])][2])
    n=length(order)
    tessella_coordinates=Matrix{Float64}(undef,3,2n)
    max_gmsh_error=0.0
    for (column,pair_index) in enumerate(order)
        master=coordinates[Int(master_tags[pair_index])]
        slave=coordinates[Int(slave_tags[pair_index])]
        expected=(master[1]+1.0,master[2],master[3])
        max_gmsh_error=max(max_gmsh_error,
            hypot(slave[1]-expected[1],slave[2]-expected[2],slave[3]-expected[3]))
        tessella_coordinates[:,column].=master
        tessella_coordinates[:,n+column].=slave
    end
    max_gmsh_error<=1e-11 || error(
        "Gmsh periodic node correspondence error is $max_gmsh_error")

    # Execute the persistent native model path from `.geo` against the same five
    # physical node pairs. Tessella's 0.5 characteristic length and Gmsh's
    # explicit Transfinite count both produce quarter-point subdivisions here.
    native_execution=execute_geo(NATIVE_GEO;mesh_dim=2)
    native_model=native_execution.model
    native_mesh=native_execution.mesh
    native_mesh===nothing && error("Tessella native periodic `.geo` did not mesh")
    validate(native_mesh).ok || error(
        "Tessella native periodic `.geo` mesh is invalid")
    native_crc=mesh_crc(native_mesh).sha
    native_crc=="3511d556ca0894daa79152eaf56abc6961024a72fa4f7e94f3357a7aa3cf0ff5" ||
        error("Tessella native periodic `.geo` mesh CRC changed to $native_crc")
    native_constraint=only(model_periodic_constraints(native_model))
    native_constraint.reversed || error(
        "Tessella did not retain the reversed native `.geo` curve orientation")
    native_mapping=model_periodic_nodes(native_model,native_mesh,1,2)
    native_mapping.master_entity==4 || error(
        "Tessella native periodic master curve is $(native_mapping.master_entity), expected 4")
    native_mapping.affine==Tuple(expected_affine) || error(
        "Tessella changed the native periodic affine transform")
    length(native_mapping.slave_nodes)==length(native_mapping.master_nodes)==n ||
        error("Tessella native periodic curve pair count is not $n")
    native_order=sortperm(eachindex(native_mapping.master_nodes);
                          by=i->native_mesh.coords[2,native_mapping.master_nodes[i]])
    max_native_difference=0.0
    for (native_index,gmsh_index) in zip(native_order,order)
        native_master=Tuple(native_mesh.coords[:,native_mapping.master_nodes[native_index]])
        native_slave=Tuple(native_mesh.coords[:,native_mapping.slave_nodes[native_index]])
        gmsh_master=coordinates[Int(master_tags[gmsh_index])]
        gmsh_slave=coordinates[Int(slave_tags[gmsh_index])]
        max_native_difference=max(
            max_native_difference,
            hypot((native_master.-gmsh_master)...),
            hypot((native_slave.-gmsh_slave)...))
        native_slave==(native_master[1]+1.0,native_master[2],native_master[3]) ||
            error("Tessella native periodic pair was not snapped exactly")
    end
    max_native_difference<=1e-11 || error(
        "Tessella/Gmsh native periodic pair difference is $max_native_difference")

    # Project the native geometry ownership and all three endpoint/curve links
    # to each supported MSH mode, then require Gmsh to recover the live curve
    # relation. MSH2 parsing needs Mesh.IgnorePeriodicity disabled explicitly.
    native_projection=model_to_mixed(native_model,native_mesh,1)
    validate(native_projection).ok || error(
        "Tessella native periodic MSH projection is invalid")
    projected_crcs=Dict(2.2=>Set{String}(),4.1=>Set{String}())
    max_projected_error=0.0
    mktempdir() do directory
        for version in (2.2,4.1),binary in (false,true)
            path=joinpath(directory,"native-projected-$version-$binary.msh")
            write_mixed_msh(
                path,native_projection;version=version,binary=binary)
            reread=read_mixed_msh(path)
            Tessella.Elements.validate(reread).ok || error(
                "Tessella rejected its native periodic MSH $version projection")
            length(reread.periodic_links)==3 || error(
                "Tessella native periodic MSH $version projection lost endpoint or curve links")
            push!(projected_crcs[version],mixed_crc(reread).sha)

            gmsh.clear()
            gmsh.option.setNumber("Mesh.IgnorePeriodicity",0)
            gmsh.open(path)
            projected_master,projected_slaves,projected_masters,projected_affine=
                gmsh.model.mesh.getPeriodicNodes(1,2)
            projected_master==4 || error(
                "Gmsh recovered projected periodic master $projected_master, expected 4")
            length(projected_slaves)==length(projected_masters)==5 || error(
                "Gmsh recovered $(length(projected_slaves)) projected periodic pairs, expected 5")
            projected_affine==expected_affine || error(
                "Gmsh changed Tessella's projected periodic affine transform")
            projected_point_two=gmsh.model.mesh.getPeriodicNodes(0,2)
            projected_point_three=gmsh.model.mesh.getPeriodicNodes(0,3)
            projected_point_two[1]==1 && length(projected_point_two[2])==1 ||
                error("Gmsh lost Tessella's projected Point[2] -> Point[1] link")
            projected_point_three[1]==4 && length(projected_point_three[2])==1 ||
                error("Gmsh lost Tessella's projected Point[3] -> Point[4] link")
            projected_tags,projected_values,_=gmsh.model.mesh.getNodes()
            projected_coordinates=Dict{Int,NTuple{3,Float64}}()
            for (index,tag) in enumerate(projected_tags)
                projected_coordinates[Int(tag)]=(
                    projected_values[3index-2],projected_values[3index-1],
                    projected_values[3index])
            end
            for (slave_tag,master_tag) in
                zip(projected_slaves,projected_masters)
                slave=projected_coordinates[Int(slave_tag)]
                master=projected_coordinates[Int(master_tag)]
                max_projected_error=max(max_projected_error,
                    hypot(slave[1]-master[1]-1.0,
                          slave[2]-master[2],slave[3]-master[3]))
            end
        end
    end
    all(length(crcs)==1 for crcs in values(projected_crcs)) || error(
        "native periodic MSH projection CRC depends on file mode")
    only(projected_crcs[2.2])==
        "506ae0fac8562df49231df71f3b12d7259ba44b3fb5618a064a15f97698951a0" ||
        error("native projected MSH2 CRC changed")
    only(projected_crcs[4.1])==
        "cf03be1a36427f1ef0fbc4e852996bd65d2630b5ac384fa0267dd14e46ea6280" ||
        error("native projected MSH4 CRC changed")
    max_projected_error<=1e-12 || error(
        "Gmsh native periodic MSH projection error is $max_projected_error")

    # The shared-corner lattice fixture requires a deterministic point-link
    # forest: Gmsh's MSH readers permit one master per periodic point even when
    # two curve transformations meet there.
    gmsh.clear()
    gmsh.open(TWO_DIRECTION_GEO)
    gmsh.model.mesh.generate(2)
    gmsh_two_x=gmsh.model.mesh.getPeriodicNodes(1,2)
    gmsh_two_y=gmsh.model.mesh.getPeriodicNodes(1,3)
    gmsh_two_x[1]==4 && length(gmsh_two_x[2])==3 &&
        gmsh_two_x[4]==expected_affine || error(
        "Gmsh two-direction x-periodic relation changed")
    expected_y_affine=[1.0,0.0,0.0,0.0,
                       0.0,1.0,0.0,1.0,
                       0.0,0.0,1.0,0.0,
                       0.0,0.0,0.0,1.0]
    gmsh_two_y[1]==1 && length(gmsh_two_y[2])==3 &&
        gmsh_two_y[4]==expected_y_affine || error(
        "Gmsh two-direction y-periodic relation changed")

    two_execution=execute_geo(TWO_DIRECTION_GEO;mesh_dim=2)
    two_mesh=two_execution.mesh
    two_mesh===nothing && error("Tessella two-direction `.geo` did not mesh")
    two_crc=mesh_crc(two_mesh).sha
    two_crc=="95ef6d0db94505d4f35ff870af09e952d74a32508a338b3994af347b406e9d05" ||
        error("Tessella two-direction periodic mesh CRC changed to $two_crc")
    two_projection=model_to_mixed(two_execution.model,two_mesh,1)
    length(two_projection.periodic_links)==5 || error(
        "Tessella two-direction projection did not emit three point and two curve links")
    mixed_crc(two_projection).sha==
        "d5fd8bd6ef46c78772792f0cee0c7b19cdd747f1c2932c2a8992760f19e69b20" ||
        error("Tessella two-direction projected mixed CRC changed")
    two_projection_crcs=Dict(2.2=>Set{String}(),4.1=>Set{String}())
    mktempdir() do directory
        for version in (2.2,4.1),binary in (false,true)
            path=joinpath(directory,"two-direction-$version-$binary.msh")
            write_mixed_msh(
                path,two_projection;version=version,binary=binary)
            reread=read_mixed_msh(path)
            length(reread.periodic_links)==5 || error(
                "Tessella two-direction MSH $version projection lost a link")
            push!(two_projection_crcs[version],mixed_crc(reread).sha)
            gmsh.clear()
            gmsh.option.setNumber("Mesh.IgnorePeriodicity",0)
            gmsh.open(path)
            projected_x=gmsh.model.mesh.getPeriodicNodes(1,2)
            projected_y=gmsh.model.mesh.getPeriodicNodes(1,3)
            projected_x[1]==4 && length(projected_x[2])==5 &&
                projected_x[4]==expected_affine || error(
                "Gmsh lost Tessella's projected x-periodic curve")
            projected_y[1]==1 && length(projected_y[2])==5 &&
                projected_y[4]==expected_y_affine || error(
                "Gmsh lost Tessella's projected y-periodic curve")
            projected_points=Dict(
                point=>gmsh.model.mesh.getPeriodicNodes(0,point)
                for point in (2,3,4))
            all(projected_points[point][1]==master &&
                length(projected_points[point][2])==1
                for (point,master) in ((2,1),(3,2),(4,1))) || error(
                "Gmsh lost Tessella's projected shared-corner point-link forest")
        end
    end
    only(two_projection_crcs[2.2])==
        "bac00f74b86af8d1a6b70de445cdb17a16a9513f0fc4a542bd995d9120923a58" ||
        error("Tessella two-direction projected MSH2 CRC changed")
    only(two_projection_crcs[4.1])==
        "d5fd8bd6ef46c78772792f0cee0c7b19cdd747f1c2932c2a8992760f19e69b20" ||
        error("Tessella two-direction projected MSH4 CRC changed")

    segments=Matrix{Int32}(undef,2,2(n-1))
    for i in 1:n-1
        segments[:,i].=(Int32(i),Int32(i+1))
        segments[:,n-1+i].=(Int32(n+i),Int32(n+i+1))
    end
    input=Mesh(tessella_coordinates;segs=segments)
    translation_output=periodic_identify(
        input,(1.0,0.0,0.0),collect(1:n),collect(n+1:2n);atol=1e-11)
    validate(translation_output).ok || error("Tessella periodic output is invalid")
    translation_output.segs==input.segs || error(
        "Tessella periodic operation changed connectivity")
    for i in 1:n
        translation_output.coords[:,n+i]==translation_output.coords[:,i].+
                                              [1.0,0.0,0.0] || error(
            "Tessella periodic pair $i was not snapped exactly")
    end

    # A separate pair of transfinite curves exercises Gmsh's documented
    # row-major 4×4 convention with a +90° rotation about the z axis.
    gmsh.model.add("periodic_rotation")
    p1=gmsh.model.geo.addPoint(1.0,0.0,0.0,0.2)
    p2=gmsh.model.geo.addPoint(2.0,0.0,0.0,0.2)
    p3=gmsh.model.geo.addPoint(0.0,1.0,0.0,0.2)
    p4=gmsh.model.geo.addPoint(0.0,2.0,0.0,0.2)
    master_curve=gmsh.model.geo.addLine(p1,p2)
    slave_curve=gmsh.model.geo.addLine(p3,p4)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.setTransfiniteCurve(master_curve,5)
    gmsh.model.mesh.setTransfiniteCurve(slave_curve,5)
    rotation=[0.0,-1.0,0.0,0.0,
              1.0, 0.0,0.0,0.0,
              0.0, 0.0,1.0,0.0,
              0.0, 0.0,0.0,1.0]
    gmsh.model.mesh.setPeriodic(1,[slave_curve],[master_curve],rotation)
    gmsh.model.mesh.generate(1)
    rotation_master,rotation_slave_tags,rotation_master_tags,gmsh_rotation=
        gmsh.model.mesh.getPeriodicNodes(1,slave_curve)
    rotation_master==master_curve || error(
        "Gmsh rotational periodic master curve is $rotation_master, expected $master_curve")
    length(rotation_slave_tags)==length(rotation_master_tags)==5 || error(
        "Gmsh rotational periodic pair count is $(length(rotation_slave_tags)), expected 5")
    length(gmsh_rotation)==16 || error(
        "Gmsh rotational periodic affine transform is not 4×4")
    maximum(abs.(gmsh_rotation.-rotation))<=1e-14 || error(
        "Gmsh rotational periodic affine transform is $gmsh_rotation")

    rotation_node_tags,rotation_node_coordinates,_=gmsh.model.mesh.getNodes()
    rotation_coordinates=Dict{Int,NTuple{3,Float64}}()
    for (i,tag) in enumerate(rotation_node_tags)
        rotation_coordinates[Int(tag)]=(
            rotation_node_coordinates[3i-2],rotation_node_coordinates[3i-1],
            rotation_node_coordinates[3i])
    end
    all(tag->haskey(rotation_coordinates,Int(tag)),rotation_master_tags) || error(
        "Gmsh omitted a rotational periodic master node from getNodes")
    all(tag->haskey(rotation_coordinates,Int(tag)),rotation_slave_tags) || error(
        "Gmsh omitted a rotational periodic slave node from getNodes")
    rotation_order=sortperm(eachindex(rotation_master_tags);
                           by=i->rotation_coordinates[Int(rotation_master_tags[i])][1])
    nr=length(rotation_order)
    tessella_rotation_coordinates=Matrix{Float64}(undef,3,2nr)
    max_rotation_error=0.0
    for (column,pair_index) in enumerate(rotation_order)
        master=rotation_coordinates[Int(rotation_master_tags[pair_index])]
        slave=rotation_coordinates[Int(rotation_slave_tags[pair_index])]
        expected=(-master[2],master[1],master[3])
        max_rotation_error=max(max_rotation_error,
            hypot(slave[1]-expected[1],slave[2]-expected[2],slave[3]-expected[3]))
        tessella_rotation_coordinates[:,column].=master
        tessella_rotation_coordinates[:,nr+column].=slave
    end
    max_rotation_error<=1e-11 || error(
        "Gmsh rotational periodic node correspondence error is $max_rotation_error")
    rotation_segments=Matrix{Int32}(undef,2,2(nr-1))
    for i in 1:nr-1
        rotation_segments[:,i].=(Int32(i),Int32(i+1))
        rotation_segments[:,nr-1+i].=(Int32(nr+i),Int32(nr+i+1))
    end
    rotation_input=Mesh(tessella_rotation_coordinates;segs=rotation_segments)
    rotation_output=periodic_identify_affine(
        rotation_input,gmsh_rotation,collect(1:nr),collect(nr+1:2nr);atol=1e-11)
    validate(rotation_output).ok || error(
        "Tessella affine-periodic output is invalid")
    rotation_output.segs==rotation_input.segs || error(
        "Tessella affine-periodic operation changed connectivity")
    for i in 1:nr
        master=rotation_output.coords[:,i]
        rotation_output.coords[:,nr+i]==[-master[2],master[1],master[3]] || error(
            "Tessella rotational periodic pair $i was not snapped exactly")
    end

    # Exercise the real MSH2/MSH4 $Periodic payloads emitted by Gmsh in both
    # modes, then require Gmsh to recover the relation from Tessella's rewrite.
    io_crcs=Dict(2.2=>Set{String}(),4.1=>Set{String}())
    mktempdir() do directory
        raw_cases=NamedTuple[]
        gmsh.option.setNumber("Mesh.SaveAll",1)
        for version in (2.2,4.1),binary in (false,true)
            gmsh.option.setNumber("Mesh.MshFileVersion",version)
            gmsh.option.setNumber("Mesh.Binary",binary ? 1 : 0)
            raw=joinpath(directory,"gmsh-periodic-$version-$binary.msh")
            gmsh.write(raw)
            push!(raw_cases,(version=version,binary=binary,path=raw))
        end
        for (index,case) in pairs(raw_cases)
            mixed=read_mixed_msh(case.path)
            Tessella.Elements.validate(mixed).ok || error(
                "Tessella rejected Gmsh periodic MSH $(case.version) metadata")
            if case.version==2.2
                mixed.entity_data===nothing || error(
                    "Tessella fabricated MSH4 metadata for MSH2 input")
                mixed.elementary_entities!==nothing || error(
                    "Tessella dropped MSH2 elementary-entity tags")
            else
                mixed.entity_data!==nothing || error(
                    "Tessella dropped MSH4 entity metadata")
                mixed.elementary_entities===nothing || error(
                    "Tessella fabricated MSH2 metadata for MSH4 input")
            end
            length(mixed.periodic_links)==3 || error(
                "Tessella read $(length(mixed.periodic_links)) periodic links, expected 3")
            curve_link=only(filter(link->link.dim==1,mixed.periodic_links))
            curve_link.slave_entity==slave_curve || error(
                "Tessella changed the periodic slave curve tag")
            curve_link.master_entity==master_curve || error(
                "Tessella changed the periodic master curve tag")
            length(curve_link.slave_nodes)==5 || error(
                "Tessella changed the periodic curve pair count")
            curve_link.affine==Tuple(rotation) || error(
                "Tessella changed the periodic affine transform")
            for i in eachindex(curve_link.slave_nodes)
                slave=mixed.coords[:,curve_link.slave_nodes[i]]
                master=mixed.coords[:,curve_link.master_nodes[i]]
                hypot(slave[1]+master[2],slave[2]-master[1],
                      slave[3]-master[3])<=1e-11 || error(
                    "Tessella read an invalid periodic node pair")
            end
            push!(io_crcs[case.version],mixed_crc(mixed).sha)
            rewritten=joinpath(directory,"tessella-periodic-$index.msh")
            write_mixed_msh(
                rewritten,mixed;version=case.version,binary=!case.binary)
            reread=read_mixed_msh(rewritten)
            mixed_crc(reread).sha==mixed_crc(mixed).sha || error(
                "Tessella periodic MSH $(case.version) rewrite changed its CRC")

            gmsh.clear()
            gmsh.option.setNumber("Mesh.IgnorePeriodicity",0)
            gmsh.open(rewritten)
            checked_master,checked_slaves,checked_masters,checked_affine=
                gmsh.model.mesh.getPeriodicNodes(1,slave_curve)
            checked_master==master_curve || error(
                "Gmsh lost Tessella's periodic master curve")
            length(checked_slaves)==length(checked_masters)==5 || error(
                "Gmsh lost Tessella's periodic node pairs")
            checked_affine==rotation || error(
                "Gmsh changed Tessella's periodic affine transform")

            # Gmsh 4.15.2 stores a parsed MSH2 $Periodic section as a raw model
            # attribute as well as live mesh metadata. Remove that duplicate
            # representation before asking Gmsh to write the live relation.
            attribute_names=gmsh.model.getAttributeNames()
            if case.version==2.2
                "Periodic" in attribute_names || error(
                    "Gmsh did not expose its parsed MSH2 Periodic attribute")
                gmsh.model.removeAttribute("Periodic")
            end
            gmsh_roundtrip=joinpath(directory,"gmsh-rewritten-$index.msh")
            gmsh.option.setNumber("Mesh.MshFileVersion",case.version)
            gmsh.option.setNumber("Mesh.Binary",case.binary ? 1 : 0)
            gmsh.write(gmsh_roundtrip)
            gmsh.clear()
            gmsh.option.setNumber("Mesh.IgnorePeriodicity",0)
            gmsh.open(gmsh_roundtrip)
            final_master,final_slaves,final_masters,final_affine=
                gmsh.model.mesh.getPeriodicNodes(1,slave_curve)
            final_master==master_curve || error(
                "Gmsh rewrite lost Tessella's periodic master curve")
            length(final_slaves)==length(final_masters)==5 || error(
                "Gmsh rewrite lost Tessella's periodic node pairs")
            final_affine==rotation || error(
                "Gmsh rewrite changed Tessella's periodic affine transform")
        end
    end
    all(length(crcs)==1 for crcs in values(io_crcs)) || error(
        "Gmsh ASCII/binary periodic inputs produced format-dependent Tessella CRCs")

    println("GMSH_PARITY_PERIODIC_OK gmsh=$(gmsh.GMSH_API_VERSION) "*
            "translation_pairs=$n translation_error=$max_gmsh_error "*
            "translation_sha=$(mesh_crc(translation_output).sha) "*
            "native_pairs=$(length(native_mapping.slave_nodes)) "*
            "native_difference=$max_native_difference native_sha=$native_crc "*
            "projected_error=$max_projected_error "*
            "projected_msh2_crc=$(only(projected_crcs[2.2])) "*
            "projected_msh4_crc=$(only(projected_crcs[4.1])) "*
            "two_direction_sha=$two_crc "*
            "two_direction_msh2_crc=$(only(two_projection_crcs[2.2])) "*
            "two_direction_msh4_crc=$(only(two_projection_crcs[4.1])) "*
            "rotation_pairs=$nr rotation_error=$max_rotation_error "*
            "rotation_sha=$(mesh_crc(rotation_output).sha) "*
            "msh2_links=3 msh2_crc=$(only(io_crcs[2.2])) "*
            "msh4_links=3 msh4_crc=$(only(io_crcs[4.1]))")
finally
    gmsh.finalize()
end
