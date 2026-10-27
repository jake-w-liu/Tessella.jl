module CLIFrontendsFunctionalChecks
using Test, Tessella
const A=Tessella.API
const E=Tessella.Elements
const LINE="Point(1)={0,0,0,0.25};Point(2)={1,0,0,0.25};Line(1)={1,2};"
const SQUARE="""
Point(1)={0,0,0,0.4};Point(2)={1,0,0,0.4};
Point(3)={1,1,0,0.4};Point(4)={0,1,0,0.4};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Curve Loop(1)={1,2,3,4};Plane Surface(1)={1};
"""
dim_cells(m,d)=sum(size(b.nodes,2) for b in m.blocks if E.msh_dimension(Int(b.msh))==d;init=0)
total_cells(m)=sum(size(b.nodes,2) for b in m.blocks)

@testset "CLI session actions and conversion" begin
    mktempdir() do directory
        line_geo=joinpath(directory,"line.geo");write(line_geo,LINE)
        line_msh=joinpath(directory,"line.msh")
        Tessella.CLI.main([line_geo,"-1","-o",line_msh])
        baseline=E.read_mixed_msh(line_msh)

        # -refine on a loaded .msh routes through API.mesh.refine.
        refined=joinpath(directory,"refined.msh")
        @test Tessella.CLI.main([line_msh,"-refine","-o",refined])==refined
        product=E.read_mixed_msh(refined)
        @test Tessella.validate(product).ok
        @test dim_cells(product,1)==2*dim_cells(baseline,1)

        # .msh output of a session refinement requires -o (no in-place write).
        @test_throws ArgumentError Tessella.CLI.main([line_msh,"-refine"])
        @test_throws ArgumentError Tessella.CLI.main([line_msh,"-refine","-o",line_msh])

        # -save rewrites the loaded mesh through API.write.
        saved=joinpath(directory,"saved.msh")
        @test Tessella.CLI.main([line_msh,"-save","-o",saved])==saved
        @test Tessella.validate(E.read_mixed_msh(saved)).ok

        # -part assigns partition labels like API.mesh.partition (label mode
        # is session state and serializes no partition entities); the request
        # is real — `API.model.get_number_of_partitions()` reports it.
        square_geo=joinpath(directory,"square.geo");write(square_geo,SQUARE)
        partitioned=joinpath(directory,"partitioned.msh")
        @test Tessella.CLI.main([square_geo,"-2","-part","2","-o",partitioned])==
              partitioned
        product=E.read_mixed_msh(partitioned)
        @test Tessella.validate(product).ok

        # -part on a loaded .msh partitions without meshing.
        unpartitioned=joinpath(directory,"unpartitioned.msh")
        Tessella.CLI.main([square_geo,"-2","-o",unpartitioned])
        repartitioned=joinpath(directory,"repartitioned.msh")
        Tessella.CLI.main([unpartitioned,"-part","3","-o",repartitioned])
        @test Tessella.validate(E.read_mixed_msh(repartitioned)).ok

        # Topology/ghost production is bounded to discrete session models —
        # it rejects loudly on BRep imports rather than silently degrading.
        box_geo=joinpath(directory,"box.geo")
        write(box_geo,"""lc=0.5;
Point(1)={0,0,0,lc};Point(2)={1,0,0,lc};Point(3)={1,1,0,lc};Point(4)={0,1,0,lc};
Point(5)={0,0,1,lc};Point(6)={1,0,1,lc};Point(7)={1,1,1,lc};Point(8)={0,1,1,lc};
Line(1)={1,2};Line(2)={2,3};Line(3)={3,4};Line(4)={4,1};
Line(5)={5,6};Line(6)={6,7};Line(7)={7,8};Line(8)={8,5};
Line(9)={1,5};Line(10)={2,6};Line(11)={3,7};Line(12)={4,8};
Curve Loop(1)={1,2,3,4};Curve Loop(2)={5,6,7,8};
Curve Loop(3)={1,10,-5,-9};Curve Loop(4)={2,11,-6,-10};
Curve Loop(5)={3,12,-7,-11};Curve Loop(6)={4,9,-8,-12};
Plane Surface(1)={1};Plane Surface(2)={2};Plane Surface(3)={3};
Plane Surface(4)={4};Plane Surface(5)={5};Plane Surface(6)={6};
Surface Loop(1)={1,2,3,4,5,6};Volume(1)={1};""")
        box_msh=joinpath(directory,"box.msh")
        Tessella.CLI.main([box_geo,"-3","-o",box_msh])
        @test_throws ArgumentError Tessella.CLI.main(
            [box_msh,"-part","2","-part_topo","-o",partitioned])
        @test_throws ArgumentError Tessella.CLI.main(
            [box_geo,"-3","-part","2","-part_topo","-o",partitioned])

        # -merge combines sources before meshing.
        line2_geo=joinpath(directory,"line2.geo")
        write(line2_geo,
            "Point(11)={0,0,0,0.25};Point(12)={0,1,0,0.25};Line(11)={11,12};")
        merged=joinpath(directory,"merged.msh")
        @test Tessella.CLI.main([line_geo,"-merge",line2_geo,"-1","-o",merged])==
              merged
        merged_mesh=E.read_mixed_msh(merged)
        @test Tessella.validate(merged_mesh).ok
        @test dim_cells(merged_mesh,1)==2*dim_cells(baseline,1)

        # -new restarts the model list; -open replaces the current model.
        restarted=joinpath(directory,"restarted.msh")
        Tessella.CLI.main([line_geo,"-new",line2_geo,"-1","-o",restarted])
        @test dim_cells(E.read_mixed_msh(restarted),1)==dim_cells(baseline,1)
        @test_throws ArgumentError Tessella.CLI.main([line_geo,line2_geo])
        # A trailing marker with no following file is a no-op.
        @test Tessella.CLI.main([line_geo,"-merge"])==line_geo

        # A bare .msh input parses through the session open (no output).
        @test Tessella.CLI.main([line_msh])==line_msh
        # -convert rewrites inputs into binary `_new` files.
        new_path=Tessella.CLI.main(["-convert",line_msh])
        @test new_path==line_msh*"_new" && isfile(new_path)
        product=E.read_mixed_msh(new_path)
        @test Tessella.validate(product).ok
        @test dim_cells(product,1)==dim_cells(baseline,1)

        pos_path=joinpath(directory,"view.pos")
        write(pos_path,"""View "v" {
SP(0,0,0){1};
SP(1,0,0){2};
};
""")
        pos_new=Tessella.CLI.main(["-convert",pos_path])
        @test pos_new==pos_path*"_new" && isfile(pos_new)
        @test length(Tessella.Post.read_pos(pos_new))==1

        # .geo -convert writes the unrolled model script (the `-0` artifact).
        geo_new=Tessella.CLI.main(["-convert",line_geo])
        @test geo_new==line_geo*"_new" && isfile(geo_new)
        @test occursin("Point(",read(geo_new,String))
        @test_throws ArgumentError Tessella.CLI.main(["-convert"])
        @test_throws ArgumentError Tessella.CLI.main(
            ["-convert",joinpath(directory,"missing.msh")])

        # Option persistence: sizing/order flags are session options on the
        # session path and forced numeric options on the standalone path.
        sized=joinpath(directory,"sized.msh")
        Tessella.CLI.main([line_msh,"-refine","-save_all","-o",sized])
        @test Tessella.validate(E.read_mixed_msh(sized)).ok

        # -setnumber / -string inject parser input before the .geo source.
        injected_geo=joinpath(directory,"injected.geo")
        write(injected_geo,"Point(1)={0,0,0,0.25};Point(2)={x,0,0,0.25};Line(1)={1,2};")
        injected_out=joinpath(directory,"injected.msh")
        Tessella.CLI.main([injected_geo,"-setnumber","x","2","-1","-o",injected_out])
        injected_mesh=E.read_mixed_msh(injected_out)
        @test maximum(injected_mesh.coords[1,:])==2.0
        @test_throws ArgumentError Tessella.CLI.main(
            [injected_geo,"-setnumber","x","bad","-1"])

        option_geo=joinpath(directory,"opts.geo")
        write(option_geo,"Mesh.MeshSizeMin = 0.5;\n")
        option_out=joinpath(directory,"opts.msh")
        Tessella.CLI.main([injected_geo,"-setnumber","x","1","-option",option_geo,
                           "-1","-o",option_out])
        option_mesh=E.read_mixed_msh(option_out)
        @test Tessella.validate(option_mesh).ok

        string_out=joinpath(directory,"string.txt")
        Tessella.CLI.main([injected_geo,"-setnumber","x","7",
            "-string","Printf(\"%.17g\", x) > \"$(replace(string_out,'\\'=>'/'))\";"])
        @test parse(Float64,read(string_out,String))==7.0

        list_geo=joinpath(directory,"list.geo")
        write(list_geo,"Point(1)={0,0,0,0.25};Point(2)={xs(1),0,0,0.25};Line(1)={1,2};")
        list_out=joinpath(directory,"list.msh")
        Tessella.CLI.main([list_geo,"-setlist","xs","2","0.5","1.5","-1","-o",list_out])
        @test maximum(E.read_mixed_msh(list_out).coords[1,:])==1.5

        # New numeric/algorithm flags have real mesh-time consumers.
        algo_out=joinpath(directory,"algo.msh")
        Tessella.CLI.main([square_geo,"-2","-algo","del2d","-o",algo_out])
        @test Tessella.validate(E.read_mixed_msh(algo_out)).ok
        Tessella.CLI.main([square_geo,"-2","-algo","front2d","-o",algo_out])
        Tessella.CLI.main([line_geo,"-1","-epslc1d","0.01","-edgelmin","0.001",
            "-clcurv","8","-clextend","1","-smooth","2","-tol","1e-6",
            "-optimize","-optimize_netgen","-optimize_threshold","0.5",
            "-degree","2","-v","3","-nt","1","-nopopup","-o",algo_out])
        @test Tessella.validate(E.read_mixed_msh(algo_out)).ok
        @test_throws ArgumentError Tessella.CLI.main([square_geo,"-2","-algo","bamg"])
        @test_throws ArgumentError Tessella.CLI.main([square_geo,"-2","-algo","delquad"])

        # -bin/-vmsh/-f route through the writer.
        bin_out=joinpath(directory,"bin.msh")
        Tessella.CLI.main([line_geo,"-1","-bin","-o",bin_out])
        @test startswith(read(bin_out,String),"\$MeshFormat\n4.1 1 8")
        v2_out=joinpath(directory,"v2.msh")
        Tessella.CLI.main([line_geo,"-1","-vmsh","2.2","-o",v2_out])
        @test occursin("\$MeshFormat\n2.2",read(v2_out,String))
        f_out=joinpath(directory,"f.msh")
        Tessella.CLI.main([line_geo,"-1","-f","msh2","-o",f_out])
        @test occursin("\$MeshFormat\n2.2",read(f_out,String))

        # Informational flags.
        @test occursin("tessella",Tessella.CLI.main(["-info"]))
        @test occursin("Mesh.MeshSizeMin",Tessella.CLI.main(["-help_options"]))
        pid_out=redirect_stdout(Base.devnull) do
            Tessella.CLI.main([line_geo,"-pid"])
        end
        @test pid_out==line_geo

        # -reclassify/-reparam route to API.mesh.classify_surfaces/create_geometry;
        # exporting a reclassified BRep-imported model is a bounded API path and
        # fails loudly rather than silently writing the wrong thing.
        surface_msh=joinpath(directory,"surface.msh")
        Tessella.CLI.main([square_geo,"-2","-o",surface_msh])
        classified=joinpath(directory,"classified.msh")
        @test_throws ArgumentError Tessella.CLI.main(
            [surface_msh,"-reclassify","40","-o",classified])
        reparamed=joinpath(directory,"reparamed.msh")
        @test_throws ArgumentError Tessella.CLI.main(
            [surface_msh,"-reparam","40","-o",reparamed])

        # Loud rejections: batch conflicts, unsupported and obsolete flags.
        @test_throws ArgumentError Tessella.CLI.main([square_geo,"-2","-refine"])
        @test_throws ArgumentError Tessella.CLI.main([line_msh,"-refine","-save","-o",classified])
        for flag in ("-4","-rand","-ho_min",
                     "-camera","-newreg","-part_weight",
                     "-aniso_max","-numsubedges","-nodb",
                     "-post_pro","-min","-listnumber","-view","-clear")
            @test_throws ArgumentError Tessella.CLI.main([line_geo,"-1",flag])
        end
        # `-0`/`-barycentric_refine` conflict with `-1` as batch actions;
        # `-optimize_ho`/`-save_parametric` are implemented session flags.
        @test_throws ArgumentError Tessella.CLI.main([line_geo,"-1","-0"])
        @test_throws ArgumentError Tessella.CLI.main([line_geo,"-1","-barycentric_refine"])
        @test_throws ArgumentError Tessella.CLI.main([line_geo,"-1","-part"])
        @test_throws ArgumentError Tessella.CLI.main([line_geo,"-1","-part","0"])
        @test_throws ArgumentError Tessella.CLI.main([line_geo,"-1","-part_topo","-part_no_topo"])
        @test_throws ArgumentError Tessella.CLI.main([line_geo,"-1","-setstring","x"])
        @test_throws ArgumentError Tessella.CLI.main([line_geo,"-1","-setlist","xs","3","1","2"])
        @test_throws ArgumentError Tessella.CLI.main([line_geo,"-1","-option","missing.geo"])
        @test_throws ArgumentError Tessella.CLI.main([line_geo,"-1","-"])
        @test_throws ArgumentError Tessella.CLI.main([line_msh,"-setstring","s","v"])
    end
end

@testset "CLI log/cpu/statreport/switch-tags flags" begin
    mktempdir() do directory
        square_geo=joinpath(directory,"square.geo")
        write(square_geo,SQUARE*"Printf(\"cli-log-probe\");\n")
        square_msh=joinpath(directory,"square.msh")
        log_path=joinpath(directory,"run.log")
        rep_path=joinpath(directory,"stats.txt")

        # -log mirrors the Msg channels (Direct rows carry no prefix); -cpu
        # appends the upstream resource suffix to Info/Error rows — the
        # statreport tau index is the one Info this pipeline emits.
        @test Tessella.CLI.main([square_geo,"-2","-o",square_msh,
            "-log",log_path,"-cpu","-statreport",rep_path])==square_msh
        log_text=read(log_path,String)
        @test occursin("cli-log-probe",log_text)
        @test occursin(
            r"Info    : Efficiency index for surface mesh tau=0\s*\(From start: Wall [0-9.e+-]+s, CPU [0-9.e+-]+s, Mem [0-9.e+-]+Mb\)",
            log_text)

        # -statreport truncates with the header + one row; -append adds rows.
        rows=split(chomp(read(rep_path,String)),r"\r?\n")
        @test length(rows)==2 && startswith(rows[1],"2D stats\tname")
        nums=[strip(f) for f in split(rows[2],'\t') if !isempty(strip(f))]
        @test nums[1]=="square" && parse(Int,nums[2])==1 &&
              parse(Int,nums[3])==0 && parse(Int,nums[4])>0
        Tessella.CLI.main([square_geo,"-2","-o",square_msh,
                           "-append_statreport",rep_path])
        @test length(split(chomp(read(rep_path,String)),r"\r?\n"))==3
        # A pure .msh session run reports only non-discrete faces.
        session_rep=joinpath(directory,"session.txt")
        Tessella.CLI.main([square_msh,"-2","-o",
            joinpath(directory,"session.msh"),"-statreport",session_rep])
        session_rows=split(chomp(read(session_rep,String)),r"\r?\n")
        @test length(session_rows)==2

        # -preserve_numbering_msh2 emits the stored element numbering on MSH
        # 2.2 output; the same run without it renumbers sequentially.
        sparse_in=joinpath(directory,"sparse.msh")
        write(sparse_in,"""\$MeshFormat
2.2 0 8
\$EndMeshFormat
\$Nodes
4
1 0 0 0
2 1 0 0
3 0 1 0
4 1 1 0
\$EndNodes
\$Elements
2
5 2 2 9 3 1 2 3
17 2 2 9 3 2 4 3
\$EndElements
""")
        function elem_tags(path)
            lines=split(read(path,String),'\n')
            i=findfirst(==("\$Elements"),chomp.(lines))
            n=parse(Int,strip(lines[i+1]))
            return sort!([parse(Int,split(lines[i+1+j])[1]) for j in 1:n])
        end
        out_seq=joinpath(directory,"seq.msh2")
        out_pres=joinpath(directory,"pres.msh2")
        Tessella.CLI.main([sparse_in,"-save","-o",out_seq])
        Tessella.CLI.main([sparse_in,"-save","-o",out_pres,
                           "-preserve_numbering_msh2"])
        @test elem_tags(out_seq)==[1,2]
        @test elem_tags(out_pres)==[5,17]

        # -switch_tags swaps the MSH 2.2 physical/elementary fields on input:
        # the imported cell's owner entity follows the swapped elementary tag.
        function imported_owners(path,switch)
            A.initialize()
            A.OPTIONS["Mesh.SwitchElementTags"]=switch ? 1.0 : 0.0
            A.open(path)
            mesh=lock(A.STATE_LOCK) do
                A.LAST_MESH[]
            end
            A.finalize()
            return sort!(reduce(vcat,copy.(mesh.entity_data.block_entities)))
        end
        @test imported_owners(sparse_in,false)==Int32[3,3]
        @test imported_owners(sparse_in,true)==Int32[9,9]

        # Flags that stay rejected keep an informative message.
        for (flag,fragment) in (("-4","adaptive meshing"),
                                ("-part_weight","weighted graph partitioning"),
                                ("-swapangle","boundary-edge swap"),
                                ("-match","matching"))
            @test_throws ArgumentError Tessella.CLI.main(
                [square_geo,"-2","-o",square_msh,flag])
        end
    end
end
end
