# ── Allocation audit: closure boxing and hot-kernel allocation bounds ──────────
#
# Correctness  : a `Core.Box` in lowered code marks a local that a closure
#                captures after reassignment; Julia then types every use as
#                `Any`, which allocates on every call (measured: a boxed scale
#                in `_quality_frame2` made 2-D refinement allocate 2.49 GB for a
#                3,106-triangle square).  Every method of the package is
#                scanned; boxes are tolerated only in the documented one-shot
#                or exact-fallback paths below.
# Robustness   : the hot kernels are measured after a warm-up call at ordinary,
#                exactly degenerate (coplanar/collinear/right-angle), tiny, and
#                huge scales.
# Completeness : a new hot-path box or allocation fails the gate explicitly.

using Test
using Tessella
using Tessella.Predicates: orient2, orient3, diametral_sign
using Tessella.MeshTypes: tet_dihedral_extrema, tet_radius_edge, tet_circumradius,
                          triangle_area, tet_signed_volume
using Tessella.Mesh2D: _quality_frame2, _circumcenter2, _radius_edge2, _encroaches
using Tessella.Mesh3D: _rb_dedup_surface, _rb_build_regions, _rb_side_float,
                       _rb_edge_pierces
using Tessella.Geometry: box_surface

@testset "closure boxing audit" begin
    # Exact-arithmetic fallbacks, one-shot setup, and recursive local helpers on
    # cold paths; each entry names the boxed locals it tolerates.
    allowed=Dict(
        "Predicates._orient_nd_sos_exact"=>["visit"],
        "SizeField._nearest_point_big"=>["best_id","best"],
        "SizeField._build_geo_field"=>["field","build"],
        "SizeField.StructuredField"=>["data","count","d","o"],
        "SizeField._postview_cell_value_big"=>["b","a"],
        "Model.model_is_inside"=>["values"],
        "Model.model_closest_point"=>["exact"],
        "Model._discrete_surface_locate"=>String[],
        "Model._model_entity_state_targets"=>["add_entity!"],
        "Model._transfinite_automatic_surface!"=>["lc"],
        "Model._model_removal_plan"=>["attempt!"],
        "API._remove_elements"=>["element_tags"],
        "API._partition"=>["count"],
        "CLI.main"=>["dim"],
        "IO._weld_triangles"=>["tol","robust_diagonal","hi","lo"],
        "IO.read_geo_params"=>["control_depth","boundary_layer_fan_elements",
            "mesh_size_from_curvature","geometry_tolerance","background","seed",
            "sfactor","smax","smin","caller","_consume_stmt"],
        "BoundaryLayer._float_fill"=>["m"],
        "Transfinite._exact_plane_frame"=>["best_norm","best"],
        "TransfiniteTriangle._exact_plane_frame"=>["best_norm","best"],
    )
    seen=Set{Module}()
    offending=String[]
    function walk(mod)
        mod in seen && return
        push!(seen,mod)
        for name in names(mod;all=true)
            isdefined(mod,name) || continue
            object=getfield(mod,name)
            if object isa Module && object!==mod && parentmodule(object)===mod
                walk(object)
            elseif object isa Function || object isa Type
                methods_list=try methods(object) catch; continue end
                for method in methods_list
                    method.module===mod || continue
                    info=try Base.uncompressed_ast(method) catch; continue end
                    boxed=String[]
                    for statement in info.code
                        statement isa Expr && statement.head==:(=) || continue
                        rhs=statement.args[2]
                        (rhs isa Expr && rhs.head==:call &&
                         rhs.args[1]==GlobalRef(Core,:Box)) || continue
                        push!(boxed,string(info.slotnames[statement.args[1].id]))
                    end
                    isempty(boxed) && continue
                    method_name=replace(string(method.name),r"^#|#\d+$"=>"")
                    key=string(nameof(mod),".",method_name)
                    tolerated=get(allowed,key,nothing)
                    if tolerated===nothing || !issubset(Set(boxed),Set(tolerated))
                        push!(offending,"$key boxed=$(unique(boxed))")
                    end
                end
            end
        end
    end
    walk(Tessella)
    @test isempty(offending)
    isempty(offending) || foreach(println,offending)
end

@testset "hot kernels allocate nothing after warm-up" begin
    # `measure` specializes on the closure type so the harness itself adds no
    # dynamic dispatch or return boxing to the measurement.
    measure(f)=(f(); @allocated f())
    # Scales inside the expansion magnitude guard (2^-230..2^300 for orient3);
    # beyond it the predicates deliberately use the allocating BigInt path.
    for scale in (1.0,1e-60,1e60,1e-20,1e20)
        a=(0.1,0.2,0.3).*scale;b=(0.9,0.2,0.3).*scale
        c=(0.1,0.7,0.3).*scale;d=(0.4,0.4,0.3).*scale  # exactly coplanar
        e=(0.4,0.4,0.9).*scale
        @test measure(()->orient3(a,b,c,d))==0
        @test measure(()->orient3(a,b,c,e))==0
        @test measure(()->tet_dihedral_extrema(a,b,c,e))==0
        @test measure(()->tet_radius_edge(a,b,c,e))==0
        @test measure(()->tet_circumradius(a,b,c,e))==0
        @test measure(()->triangle_area(a,b,c))==0
        @test measure(()->tet_signed_volume(a,b,c,e))==0
        p=(0.1,0.2).*scale;q=(0.9,0.2).*scale;r=(0.5,0.2).*scale;t=(0.5,0.6).*scale
        right=(0.1,0.6).*scale # right angle at p: exact zero diametral product
        @test measure(()->orient2(p,q,r))==0
        @test measure(()->orient2(p,q,t))==0
        @test measure(()->_quality_frame2(p,q,t))==0
        @test measure(()->_circumcenter2(p,q,t))==0
        @test measure(()->_radius_edge2(p,q,t))==0
        @test measure(()->_encroaches(p,q,t))==0
        @test measure(()->_encroaches(p,q,right))==0
    end
    # diametral test at a lattice right angle and at a grazing point
    @test measure(()->diametral_sign((1.0,0.0),(0.0,1.0),(0.0,0.0)))==0
    @test measure(()->diametral_sign((1.0,0.0),(0.0,1.0),(0.0,2^-60)))==0
    # conformity-gate decisions on a faceted box
    surface=box_surface(0,2,0,1,0,1)
    Px,Py,Pz,facets=_rb_dedup_surface(surface)
    regions=_rb_build_regions(Px,Py,Pz,facets)
    region=regions[1];probe=(0.3,0.4,0.6);on_plane=(Px[1],Py[1],Pz[1])
    @test measure(()->_rb_side_float(region.plane,probe))==0
    @test measure(()->_rb_side_float(region.plane,on_plane))==0
    Qx=vcat(Px,[0.5,0.7]);Qy=vcat(Py,[0.5,0.2]);Qz=vcat(Pz,[-0.5,0.4])
    # distinct names: reusing `p`/`q` here would make the loop closures above
    # capture reassigned (boxed) locals and measure the harness, not the kernel
    edge_p=Int32(length(Qx)-1);edge_q=Int32(length(Qx))
    @test measure(()->_rb_edge_pierces(region,Qx,Qy,Qz,edge_p,edge_q))==0
end
