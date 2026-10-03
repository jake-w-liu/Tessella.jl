#!/usr/bin/env julia
# Native mixed cache reference maps and keys against the pinned optional oracle.
using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.Elements: MixedMesh, ElementBlock, lagrange_nodes, msh_dimension
using Random
using LinearAlgebra

const _API=Tessella.API
const _A=Float64[2.0 0.25 0.5;0.5 3.0 0.25;0.25 0.5 4.0]
const _ORIGIN=Float64[10,-5,3]
const _POINTS=Dict(1=>[-0.2,0.0,0.0],2=>[0.2,0.3,0.0],
    3=>[-0.2,0.3,0.0],4=>[0.1,0.2,0.3],5=>[-0.2,0.3,-0.4],
    6=>[0.2,0.3,-0.4],7=>[-0.2,0.3,0.4],15=>[0.0,0.0,0.0])
for (quadratic,linear) in ((8,1),(9,2),(10,3),(11,4),(12,5),(13,6),(14,7))
    _POINTS[quadratic]=_POINTS[linear]
end

api_path=get(ENV,"GMSH_JULIA_API","")
isfile(api_path) || error("set GMSH_JULIA_API to the pinned Gmsh 4.15.2 gmsh.jl")
include(api_path)
gmsh.GMSH_API_VERSION=="4.15.2" || error("Gmsh 4.15.2 required")

function close_arrays(native,oracle,label)
    length(native)==length(oracle) || error("$label output count differs")
    all(isapprox.(native,oracle;rtol=2e-12,atol=2e-12)) ||
        error("$label output differs: native=$native oracle=$oracle")
end

gmsh.initialize();_API.initialize()
gmsh.option.setNumber("General.Terminal",0)
checks=0
cases=0
pyramid_strict_failures=0
warped_quad_exceptions=0
point_location_exceptions=0
try
    gmsh.option.getString("General.Version")=="4.15.2" || error("loaded Gmsh 4.15.2 library required")
    rng=MersenneTwister(425152)
    for msh in (1,2,3,4,5,6,7,15,8,9,10,11,12,13,14),warped in (false,true),permuted in (false,true)
        case_filter=get(ENV,"API_MIXED_QUERY_TYPE","")
        isempty(case_filter) || msh==parse(Int,case_filter) || continue
        global cases+=1
        gmsh.clear();gmsh.model.add("mixed_reference")
        reference=lagrange_nodes(msh)
        count=size(reference,2)
        coordinates=_A*reference .+ _ORIGIN
        warped && msh in (3,5,6,7,11,12,13,14) && (coordinates[3,end]+=0.7)
        node_tags=permuted ? randperm(rng,count) : collect(1:count)
        stored=similar(coordinates)
        for local_node in 1:count
            stored[:,node_tags[local_node]]=coordinates[:,local_node]
        end
        native=MixedMesh(stored,[ElementBlock(msh,reshape(Int32.(node_tags),count,1))])
        lock(_API.STATE_LOCK) do
            _API._replace_mesh_cache_locked!(_API._copy_mesh(native))
        end
        dimension=msh_dimension(msh)
        gmsh.model.addDiscreteEntity(dimension,1)
        gmsh.model.mesh.addNodes(dimension,1,collect(1:count),vec(stored))
        gmsh.model.mesh.addElementsByType(1,msh,[1],node_tags)
        for primary in (false,true)
            _API.mesh.get_element_edge_nodes(msh,-1,primary)==
                gmsh.model.mesh.getElementEdgeNodes(msh,-1,primary) ||
                error("type$msh primary=$primary edge nodes differ")
            global checks+=1
            for face_type in (3,4)
                _API.mesh.get_element_face_nodes(msh,face_type,-1,primary)==
                    gmsh.model.mesh.getElementFaceNodes(msh,face_type,-1,primary) ||
                    error("type$msh face$face_type primary=$primary nodes differ")
                global checks+=1
            end
            for fast in (false,true)
                close_arrays(_API.mesh.get_barycenters(msh,-1,fast,primary),
                    gmsh.model.mesh.getBarycenters(msh,-1,fast,primary),
                    "type$msh fast=$fast primary=$primary barycenters")
                global checks+=1
            end
        end
        local_coord=_POINTS[msh]
        local_coord=vcat(local_coord,local_coord./2)
        native_jac=_API.mesh.get_jacobian(1,local_coord)
        oracle_jac=gmsh.model.mesh.getJacobian(1,local_coord)
        for component in 1:3
            close_arrays(native_jac[component],oracle_jac[component],"type$msh Jacobian$component")
            global checks+=1
        end
        for physical in (native_jac[3][1:3],native_jac[3][4:6])
            native_local=_API.mesh.get_local_coordinates_in_element(1,physical...)
            oracle_local=gmsh.model.mesh.getLocalCoordinatesInElement(1,physical...)
            if msh==3 && warped
                # Pinned MElement::xyz2uvw accumulates a nonzero normal w for
                # points on warped bilinear quads; MQuadrangle::isInside then
                # rejects its own mapped point. Native queries retain canonical
                # unused w=0 and certify the actual bilinear surface residual.
                close_arrays(collect(native_local)[1:2],collect(oracle_local)[1:2],"warped quad active coordinates")
                native_local[3]==0.0 || error("warped quad unused w must be zero")
                abs(oracle_local[3])>1e-6 || error("pinned warped quad artifact changed")
            else
                close_arrays(collect(native_local),collect(oracle_local),"type$msh local coordinates")
            end
            native_location=_API.mesh.get_element_by_coordinates(physical...,-1,true)
            if (msh==3 && warped) || msh==15
                rejected=false
                try
                    gmsh.model.mesh.getElementByCoordinates(physical...,-1,true)
                catch err
                    startswith(sprint(showerror,err),"No element found at") || rethrow()
                    rejected=true
                end
                rejected || error("pinned warped quad strict rejection changed")
                native_location[2]==msh || error("native located element must retain its type")
                if msh==15
                    global point_location_exceptions+=1
                else
                    global warped_quad_exceptions+=1
                end
            elseif msh==14
                # In this pinned build the strict quadratic-pyramid octree
                # lookup can omit an affine mapped point. Inverse coordinates
                # independently satisfy the strict reference domain; the direct
                # relaxed scan finds the same cell without changing native tolerance.
                u,v,w=oracle_local
                0<=w<=1 && abs(u)<=1-w && abs(v)<=1-w ||
                    error("quadratic pyramid oracle point is not strictly inside")
                oracle_location=try
                    gmsh.model.mesh.getElementByCoordinates(physical...,-1,true)
                catch err
                    startswith(sprint(showerror,err),"No element found at") || rethrow()
                    global pyramid_strict_failures+=1
                    gmsh.model.mesh.getElementByCoordinates(physical...,-1,false)
                end
                native_location[1:3]==oracle_location[1:3] || error("quadratic pyramid record differs")
                close_arrays(collect(native_location[4:6]),collect(oracle_location[4:6]),"quadratic pyramid coordinates")
            else
                oracle_location=try
                    gmsh.model.mesh.getElementByCoordinates(physical...,-1,true)
                catch err
                    error("type$msh warped=$warped permuted=$permuted point=$physical local=$oracle_local: "*sprint(showerror,err))
                end
                native_location[1:3]==oracle_location[1:3] || error("type$msh location record differs")
                close_arrays(collect(native_location[4:6]),collect(oracle_location[4:6]),"type$msh location coordinates")
            end
            global checks+=3
        end
        for space in ("Lagrange","H1Legendre1","H1Legendre2","H1Legendre3",
                      "HcurlLegendre0","HcurlLegendre1","HcurlLegendre2")
            msh in (7,14) && space!="Lagrange" && continue
            msh==15 && startswith(space,"Hcurl") && continue
            native_keys=_API.mesh.get_keys(msh,space)
            oracle_keys=gmsh.model.mesh.getKeys(msh,space)
            native_keys[1:2]==oracle_keys[1:2] || error(
                "type$msh $space keys differ: native=$(native_keys[1:2]) oracle=$(oracle_keys[1:2])")
            close_arrays(native_keys[3],oracle_keys[3],"type$msh $space key coordinates")
            native_orientation=_API.mesh.get_basis_functions_orientation_for_element(1,space)
            oracle_orientation=gmsh.model.mesh.getBasisFunctionsOrientationForElement(1,space)
            native_orientation==oracle_orientation || error("type$msh $space orientation differs")
            global checks+=3
        end
        if msh!=15
            for quality in ("minEdge","maxEdge")
                close_arrays(_API.mesh.get_element_qualities([1],quality),
                    gmsh.model.mesh.getElementQualities([1],quality),"type$msh $quality")
                global checks+=1
            end
        end
        if msh in 8:14
            # Independently exercise the attached-record route, which must
            # honor the same primary flags and complete P2 node closures.
            _API.model.add("record_query")
            _API.model.add_discrete_entity(dimension,1)
            _API.mesh.add_nodes(dimension,1,collect(1:count),vec(stored))
            _API.mesh.add_elements_by_type(1,msh,[1],node_tags)
            for primary in (false,true)
                _API.mesh.get_element_edge_nodes(msh,-1,primary)==
                    gmsh.model.mesh.getElementEdgeNodes(msh,-1,primary) ||
                    error("record type$msh primary=$primary edge nodes differ")
                global checks+=1
                for face_type in (3,4)
                    _API.mesh.get_element_face_nodes(msh,face_type,-1,primary)==
                        gmsh.model.mesh.getElementFaceNodes(msh,face_type,-1,primary) ||
                        error("record type$msh face$face_type primary=$primary nodes differ")
                    global checks+=1
                end
                for fast in (false,true)
                    close_arrays(_API.mesh.get_barycenters(msh,-1,fast,primary),
                        gmsh.model.mesh.getBarycenters(msh,-1,fast,primary),
                        "record type$msh fast=$fast primary=$primary barycenters")
                    global checks+=1
                end
            end
            _API.model.remove()
        end
    end
    cases>0 || error("no cases selected by API_MIXED_QUERY_TYPE")
    println("API_MIXED_QUERIES_DIFFERENTIAL_OK gmsh=4.15.2 cases=$cases checks=$checks warped_quad_upstream_exceptions=$warped_quad_exceptions point_location_upstream_exceptions=$point_location_exceptions quadratic_pyramid_strict_fallbacks=$pyramid_strict_failures")
finally
    _API.finalize();gmsh.finalize()
end
