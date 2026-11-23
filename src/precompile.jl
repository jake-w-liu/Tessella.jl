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
        # A sphere fragmented with a box: plane x sphere sections through the
        # poles, pole edges and piece bounds. Without it the first such
        # Boolean of a session runs interpreted for about 26 s under
        # `--compile=min` (the suite) where it takes 9 ms compiled. Best effort.
        try
            API.initialize()
            API.option("General.Terminal", 0)
            API.model.occ.add_box(2, 0, 0, 1, 1, 1)
            API.model.occ.add_sphere(1.7, 0, 0.5, 0.35)
            API.model.occ.fragment([(3, 1)], [(3, 2)])
            API.model.occ.synchronize()
            for (_, v) in API.model.get_entities(3)
                API.model.get_bounding_box(3, v)
            end
        catch
        finally
            API.finalize()
        end
        # Staged generation (`generate(1)`, then 2, then 3) runs other paths
        # than one `generate(3)`: the first staged box of a session spent 19 s
        # compiling them. A torus and a cone (2 s each), and the common calls
        # after a first mesh (element qualities 11 s, `set_order(2)` 4 s, an
        # .msh written and opened again) likewise. Best effort.
        for (build, dims) in (
                (() -> API.model.occ.add_box(0, 0, 0, 1, 1, 1), 1:3),
                (() -> API.model.occ.add_torus(0, 0, 0, 2, 0.5), 3:3),
                (() -> API.model.occ.add_cone(0, 0, 0, 0, 0, 1, 1, 0.3), 3:3))
            try
                API.initialize()
                API.option("General.Terminal", 0)
                build()
                API.model.occ.synchronize()
                API.option("Mesh.MeshSizeMin", 0.3)
                API.option("Mesh.MeshSizeMax", 0.3)
                for dim in dims
                    API.mesh.generate(dim)
                end
            catch
            finally
                API.finalize()
            end
        end
        mktempdir() do directory
            path = joinpath(directory, "queries.msh")
            try
                API.initialize()
                API.option("General.Terminal", 0)
                API.model.occ.add_box(0, 0, 0, 1, 1, 1)
                API.model.occ.synchronize()
                API.option("Mesh.MeshSizeMax", 0.5)
                API.mesh.generate(3)
                _, tags, _ = API.mesh.get_elements(3)
                API.mesh.get_element_qualities(tags[1], "gamma")
                API.write(path)
                API.mesh.set_order(2)
            catch
            finally
                API.finalize()
            end
            try
                API.initialize()
                API.option("General.Terminal", 0)
                API.open(path)
                API.mesh.get_nodes()
            catch
            finally
                API.finalize()
            end
        end
        # The command line: `CLI.main` infers the whole execute/mesh/write tree
        # on its first call (about 40 s), so a 2-D run with physical groups,
        # an extruded 3-D run, an OpenCASCADE box and a `-0` unroll are served
        # from the package image. Best effort, like the curved cases above.
        mktempdir() do directory
            plane = joinpath(directory, "plane.geo")
            write(plane, """
                lc = 0.5;
                Point(1) = {0,0,0,lc}; Point(2) = {1,0,0,lc};
                Point(3) = {1,1,0,lc}; Point(4) = {0,1,0,lc};
                Line(1) = {1,2}; Line(2) = {3,2}; Line(3) = {3,4}; Line(4) = {4,1};
                Curve Loop(1) = {4,1,-2,3};
                Plane Surface(1) = {1};
                Physical Curve(5) = {1,2,4};
                Physical Surface("My surface") = {1};
                """)
            extruded = joinpath(directory, "extruded.geo")
            write(extruded, """
                Include "plane.geo";
                Point(5) = {0.5,0.5,0,0.5};
                Point{5} In Surface{1};
                e[] = Extrude {0,0,0.5} { Surface{1}; Layers{2}; };
                q[] = Translate {2,0,0} { Duplicata { Volume{e[1]}; } };
                Physical Volume("v", 7) = {e[1], q[0]};
                """)
            occ = joinpath(directory, "occ.geo")
            write(occ, """
                SetFactory("OpenCASCADE");
                Box(1) = {0,0,0,1,1,1};
                Mesh.MeshSizeMax = 0.5;
                """)
            for args in ((plane, "-2"), (extruded, "-3"), (occ, "-3"),
                         (plane, "-0"))
                output = joinpath(directory,
                    args[2] == "-0" ? "out.geo_unrolled" : "out.msh")
                try
                    redirect_stdout(devnull) do
                        CLI.main(String[args[1], args[2], "-v", "0", "-nt", "1",
                                        "-o", output])
                    end
                catch
                end
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
