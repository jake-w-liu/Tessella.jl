import PrecompileTools

function _precompile_meshing_workload()
    Base.Experimental.@force_compile
    mktempdir() do directory
        for (name, source) in (
            ("transfinite.geo", """
                Point(1) = {0,0,0};
                Point(2) = {1,0,0};
                Point(3) = {1,1,0.5};
                Point(4) = {0,1,0};
                Line(1) = {1,2}; Line(2) = {2,3};
                Line(3) = {3,4}; Line(4) = {4,1};
                Curve Loop(1) = {1,2,3,4};
                Surface(1) = {1};
                Transfinite Curve {1,2,3,4} = 5;
                Transfinite Surface {1};
                """),
            ("ordinary.geo", """
                Mesh.MeshSizeMin = 0.25;
                Mesh.MeshSizeMax = 0.25;
                Point(1) = {0,0,0,0.25}; Point(2) = {1,0,0,0.25};
                Point(3) = {1,1,0,0.25}; Point(4) = {0,1,0,0.25};
                Line(1) = {1,2}; Line(2) = {2,3};
                Line(3) = {3,4}; Line(4) = {4,1};
                Curve Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
                """))
            path = joinpath(directory, name)
            write(path, source)
            execute_geo(path; mesh_dim=2)
        end
    end
    # Boundary recovery on a small torus whose quadrilaterals alternate their
    # diagonal: about half of its edges have to be recovered.
    let nu = 10, nv = 6
        id(i, j) = mod(i, nu) * nv + mod(j, nv) + 1
        coords = Matrix{Float64}(undef, 3, nu * nv)
        tris = Matrix{Int32}(undef, 3, 2 * nu * nv)
        column = 0
        for i in 0:nu-1, j in 0:nv-1
            u = 2pi * i / nu; v = 2pi * j / nv
            coords[:, id(i, j)] .= ((2 + 0.6cos(v)) * cos(u), (2 + 0.6cos(v)) * sin(u), 0.6sin(v))
            a = id(i, j); b = id(i + 1, j); c = id(i + 1, j + 1); d = id(i, j + 1)
            first, second = iseven(i + j) ? ((a, b, c), (a, c, d)) : ((a, b, d), (b, c, d))
            tris[:, column+=1] .= first
            tris[:, column+=1] .= second
        end
        BoundaryRecovery.recover_boundary_flips(coords, tris)
    end
    saved_start = API._SESSION_START_TIME[]
    saved_epoch = API._ONELAB_SESSION_EPOCH[]
    saved_lc = API.NATIVE_CONTEXT_LC[]
    saved_forced_bbox = API.NATIVE_CONTEXT_FORCED_BBOX[]
    try
        try
            API.initialize()
            API.option("General.Terminal", 0)
            API.option("Mesh.MeshSizeMin", 0.25)
            API.option("Mesh.MeshSizeMax", 0.25)
            for (tag, (x, y)) in enumerate(((0.0, 0.0), (1.0, 0.0),
                    (1.0, 1.0), (0.0, 1.0)))
                API.model.add_point(x, y, 0.0; tag, meshSize=0.25)
            end
            for (tag, (first, last)) in enumerate(((1, 2), (2, 3), (3, 4), (4, 1)))
                API.model.add_line(first, last; tag)
            end
            API.model.add_curve_loop([1, 2, 3, 4]; tag=1)
            API.model.add_plane_surface([1]; tag=1)
            API.mesh.generate(2)
        finally
            API.finalize()
        end
        try
            API.initialize()
            API.model.occ.add_box(0, 0, 0, 1, 1, 1)
            API.model.occ.synchronize()
            API.mesh.generate(3)
        finally
            API.finalize()
        end
        # Curved faces, trimmed faces and a non-convex Region: without them the
        # first sphere or drilled block of a session spends over a minute
        # compiling. A failure here must never prevent the package from
        # loading, so these two are best effort.
        for build in (
                () -> API.model.occ.add_sphere(0, 0, 0, 1),
                () -> begin
                    block = API.model.occ.add_box(0, 0, 0, 2, 2, 2)
                    drill = API.model.occ.add_cylinder(1, 1, -1, 0, 0, 4, 0.5)
                    API.model.occ.cut([(3, block)], [(3, drill)])
                end)
            try
                API.initialize()
                API.option("General.Terminal", 0)
                build()
                API.model.occ.synchronize()
                API.option("Mesh.MeshSizeMax", 0.5)
                API.mesh.generate(3)
            catch
            finally
                API.finalize()
            end
        end
    finally
        API._SESSION_START_TIME[] = saved_start
        API._ONELAB_SESSION_EPOCH[] = saved_epoch
        API.NATIVE_CONTEXT_LC[] = saved_lc
        API.NATIVE_CONTEXT_FORCED_BBOX[] = saved_forced_bbox
    end
    return nothing
end

PrecompileTools.@compile_workload _precompile_meshing_workload()
