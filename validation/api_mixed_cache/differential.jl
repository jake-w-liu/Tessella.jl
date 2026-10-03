#!/usr/bin/env julia
# Actual mixed cell families, order conversion, and geometric entity ownership.
using Pkg
Pkg.activate(joinpath(@__DIR__,"..","..");io=devnull)
using Tessella
using Tessella.Elements: lagrange_nodes,msh_dimension,msh_family,msh_type

binding=get(ENV,"GMSH_JULIA_API","")
isfile(binding) || error("set GMSH_JULIA_API to the pinned Gmsh 4.15.2 gmsh.jl")
include(binding)
gmsh.GMSH_API_VERSION=="4.15.2" || error("Gmsh 4.15.2 required")
api=Tessella.API
matrix=Float64[2 .25 .5;.5 3 .25;.25 .5 4]
shift=[10.,-5.,3.]
checks=0

function oracle_coordinates(tags)
    return hcat((gmsh.model.mesh.getNode(tag)[1] for tag in tags)...)
end

function build_native_squares(count;recombine=true)
    api.initialize()
    for patch in 1:count
        point_offset=4(patch-1)
        for (x,y) in ((0,0),(1,0),(1,1),(0,1))
            api.model.add_point(x,y,0)
        end
        for (a,b) in ((1,2),(2,3),(3,4),(4,1))
            api.model.add_line(a+point_offset,b+point_offset)
        end
        api.model.add_curve_loop(collect(point_offset+1:point_offset+4))
        api.model.add_plane_surface([patch])
        for curve in point_offset+1:point_offset+4
            api.mesh.set_transfinite_curve(curve,3)
        end
        api.mesh.set_transfinite_surface(patch)
        recombine && api.mesh.set_recombine(2,patch)
    end
    return api.mesh.generate(2)
end

function build_oracle_squares(count;recombine=true)
    gmsh.clear();gmsh.model.add("native_squares")
    for patch in 1:count
        offset=4(patch-1)
        for (local_point,(x,y)) in enumerate(((0,0),(1,0),(1,1),(0,1)))
            gmsh.model.geo.addPoint(x,y,0,1,offset+local_point)
        end
        for (local_curve,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
            gmsh.model.geo.addLine(a+offset,b+offset,offset+local_curve)
        end
        gmsh.model.geo.addCurveLoop(collect(offset+1:offset+4),patch)
        gmsh.model.geo.addPlaneSurface([patch],patch)
    end
    gmsh.model.geo.synchronize()
    for curve in 1:4count
        gmsh.model.mesh.setTransfiniteCurve(curve,3)
    end
    for patch in 1:count
        gmsh.model.mesh.setTransfiniteSurface(patch)
        recombine && gmsh.model.mesh.setRecombine(2,patch)
    end
    gmsh.model.mesh.generate(2)
end

gmsh.initialize();gmsh.option.setNumber("General.Terminal",0)
try
    for msh in (1,2,3,4,5,6,7,15)
        gmsh.clear();gmsh.model.add("native_elevation")
        reference=lagrange_nodes(msh);dim=msh_dimension(msh)
        coords=matrix*reference.+shift
        nodes=Int32.(1:size(coords,2))
        linear=MixedMesh(coords,[ElementBlock(msh,reshape(nodes,:,1))])
        owner=(dim,Int32(1))
        class=api._mixed_classification(linear,owner,[owner],fill(owner,length(nodes)),
                                       Dict(owner=>Int32[]),Dict(msh=>Int32[1]))
        elevated,elevated_class=api._mixed_quadratic_cache(linear,class,"oracle")
        gmsh.model.addDiscreteEntity(dim,1)
        gmsh.model.mesh.addNodes(dim,1,Int.(nodes),vec(coords))
        gmsh.model.mesh.addElementsByType(1,msh,[1],Int.(nodes))
        gmsh.model.mesh.setOrder(2)
        target=msh==15 ? 15 : msh_type(msh_family(msh),2)
        native_block=elevated.blocks[1]
        _,oracle_nodes=gmsh.model.mesh.getElementsByType(target)
        native_block.msh==target || error("type $msh elevation family differs")
        actual=elevated.coords[:,vec(native_block.nodes)]
        expected=oracle_coordinates(oracle_nodes)
        size(actual)==size(expected) && all(isapprox.(actual,expected;atol=2e-12,rtol=2e-12)) ||
            error("type $msh elevation map or local node order differs")
        lowered,_=api._mixed_linear_cache(elevated,elevated_class)
        gmsh.model.mesh.setOrder(1)
        _,oracle_linear=gmsh.model.mesh.getElementsByType(msh)
        lowered.coords[:,vec(lowered.blocks[1].nodes)]==oracle_coordinates(oracle_linear) ||
            error("type $msh order-one roundtrip differs")
        global checks+=3
    end
    for recombine in (true,false),patches in (1,2)
        build_native_squares(patches;recombine)
        build_oracle_squares(patches;recombine)
        for order in (1,2)
            api.mesh.set_order(order);gmsh.model.mesh.setOrder(order)
            length(api.mesh.get_nodes()[1])==length(gmsh.model.mesh.getNodes()[1]) ||
                error("$patches independent surfaces order$order node count differs")
            global checks+=1
            for dim in 0:2,(entity_dim,entity) in gmsh.model.getEntities(dim)
                native=api.mesh.get_nodes(entity_dim,entity)[1]
                oracle=gmsh.model.mesh.getNodes(entity_dim,entity)[1]
                length(native)==length(oracle) ||
                    error("surface order$order entity($entity_dim,$entity) node ownership differs")
                global checks+=1
            end
            for surface in 1:patches
                native=api.mesh.get_elements(2,surface)
                oracle=gmsh.model.mesh.getElements(2,surface)
                native[1]==oracle[1] && length.(native[2])==length.(oracle[2]) &&
                    length.(native[3])==length.(oracle[3]) ||
                    error("surface$surface order$order type/connectivity width differs")
                global checks+=1
            end
        end
        api.finalize()
    end
finally
    api.finalize();gmsh.finalize()
end
println("API_MIXED_CACHE_DIFFERENTIAL_OK gmsh=4.15.2 checks=$checks")
