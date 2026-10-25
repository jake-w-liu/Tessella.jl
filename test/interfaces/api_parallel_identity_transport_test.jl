# PAR01 residuals — Region-fill RNG transports and model-identity/merge
# transports through the detach→mesh→merge cycle.
#
# Residual 1 (RNG-consuming Region fills): explicit-PLC volumes take the
# native-initial Delaunay fill, which consumes the shared `Mesh.RandomSeed`
# ledger — a globally ordered stream, matching upstream's single
# `srand(randomSeed)` per `GenerateMesh`.  `_parallel_volume_prepare` marks
# them `may_consume` and predicts each unit's ledger draw (`3*count+1` for
# the native kernel's perturbation+seed draws over its born-UID input rows);
# detached consumers then enter at speculated positions — the live seed for
# a batch's first consumer, an LCG jump-ahead over the intervening predicted
# draws for later ones, or a completed predecessor's reported exit on the
# cross-phase arm.  The merge-time seed-ledger check admits only a receipt
# whose entry equals the live position, so a shortfall in the predicted draw
# (kernel retries) sends the slot back to the identical serial call.  These
# fixtures pin that byte-identity plus actual adoption of consuming
# produces.
#
# Residual 2 (identity/merge transports): detached work runs on transaction
# copies whose writes are certified by the `foreign` transport certificate
# (curve_params diffs inside the footprint; nothing else may move).  The
# end-to-end pins below run a parallel-generated model through
# `model.set_tag` retagging, `api.merge` file merge and post-merge
# generation, asserting the identity stores (`next_tag`, names, physical
# records, `curve_params`) and the produced mesh are identical to the same
# session driven serially.
#
# Run standalone: `julia -t4 --project=. --compile=min
# test/interfaces/api_parallel_identity_transport_test.jl`.  Under -t1 all
# arms take the serial path and the same checks hold; `_PARALLEL_ADOPTED`
# proves the detached path actually engaged under -t4.

using Test
using Tessella
using Tessella.Elements: mixed_crc, write_mixed_msh

A=Tessella.API

# An explicit PLC volume (native-initial fill — shared-RNG consumer)
# followed by three factory Boxes: the batch detaches the consumer first in
# merge order plus the three boxes.
const EXPLICIT_CUBE(points_base,lines_base,loops_base,surfs_base,shell,vol,
        x0,y0,z0)=
    join((join(("Point($(points_base+i))={$(x0+x),$(y0+y),$(z0+z),0.35}"
        for (i,(x,y,z)) in enumerate(((0,0,0),(1,0,0),(1,1,0),(0,1,0),
            (0,0,1),(1,0,1),(1,1,1),(0,1,1)))),";\n"),
        join(("Line($(lines_base+1))={$(points_base+1),$(points_base+2)}",
            "Line($(lines_base+2))={$(points_base+2),$(points_base+3)}",
            "Line($(lines_base+3))={$(points_base+3),$(points_base+4)}",
            "Line($(lines_base+4))={$(points_base+4),$(points_base+1)}",
            "Line($(lines_base+5))={$(points_base+5),$(points_base+6)}",
            "Line($(lines_base+6))={$(points_base+6),$(points_base+7)}",
            "Line($(lines_base+7))={$(points_base+7),$(points_base+8)}",
            "Line($(lines_base+8))={$(points_base+8),$(points_base+5)}",
            "Line($(lines_base+9))={$(points_base+1),$(points_base+5)}",
            "Line($(lines_base+10))={$(points_base+2),$(points_base+6)}",
            "Line($(lines_base+11))={$(points_base+3),$(points_base+7)}",
            "Line($(lines_base+12))={$(points_base+4),$(points_base+8)}"),";\n"),
        join(("Curve Loop($(loops_base+1))={$(lines_base+1),$(lines_base+2),$(lines_base+3),$(lines_base+4)}",
            "Curve Loop($(loops_base+2))={$(lines_base+5),$(lines_base+6),$(lines_base+7),$(lines_base+8)}",
            "Curve Loop($(loops_base+3))={$(lines_base+1),$(lines_base+10),-$(lines_base+5),-$(lines_base+9)}",
            "Curve Loop($(loops_base+4))={$(lines_base+2),$(lines_base+11),-$(lines_base+6),-$(lines_base+10)}",
            "Curve Loop($(loops_base+5))={$(lines_base+3),$(lines_base+12),-$(lines_base+7),-$(lines_base+11)}",
            "Curve Loop($(loops_base+6))={$(lines_base+4),$(lines_base+9),-$(lines_base+8),-$(lines_base+12)}"),";\n"),
        join(("Plane Surface($(surfs_base+i))={$(loops_base+i)}" for i in 1:6),";\n"),
        "Surface Loop($shell)={$(join((string(surfs_base+i) for i in 1:6),","))}",
        "Volume($vol)={$shell}"),";\n")*";\n"

const CONSUMEBOX=EXPLICIT_CUBE(0,0,0,0,1,1,0.0,0.0,0.0)*"""
SetFactory("OpenCASCADE");
Box(3)={3,0,0, 1,1,1};
Box(4)={3,2.5,0, 1,1,1};
Box(5)={0,2.5,0, 1,1,1};
"""

# Two explicit volumes (both predicted consumers — merge-order tags 1 and 2)
# + three boxes: the first consumer enters at the live ledger position and
# the second at the jump-ahead seed over the first's predicted draw — both
# detach, and the merge-time ledger check keeps either honest.
const TWOCONSUMEBOX=EXPLICIT_CUBE(0,0,0,0,1,1,0.0,0.0,0.0)*
    EXPLICIT_CUBE(20,20,20,20,2,2,0.0,4.0,0.0)*"""
SetFactory("OpenCASCADE");
Box(3)={3,0,0, 1,1,1};
Box(4)={3,2.5,0, 1,1,1};
Box(5)={0,2.5,0, 1,1,1};
"""

# Two eligible batches separated by a transfinite Region (merge-order tag
# 3): batch 1 = {Volume 1 consumer, Box 2}, batch 2 = {Volume 4 consumer,
# Box 5} — each batch's first (and only) consumer must adopt, exercising
# the per-batch ledger rule twice in one traversal.
const TWOBATCH=EXPLICIT_CUBE(0,0,0,0,1,1,0.0,0.0,0.0)*
    EXPLICIT_CUBE(30,30,30,30,3,3,0.0,4.0,0.0)*
    EXPLICIT_CUBE(50,50,50,50,4,4,4.0,4.0,0.0)*"""
Transfinite Curve{31,32,33,34,35,36,37,38,39,40,41,42}=3;
Transfinite Surface{31,32,33,34,35,36};
Transfinite Volume{3}={31,32,33,34,35,36,37,38};
SetFactory("OpenCASCADE");
Box(2)={3,0,0, 1,1,1};
Box(5)={4,6.5,0, 1,1,1};
"""

const MERGEGEO="""
SetFactory("OpenCASCADE");
Box(20)={6.5,0,0, 1,1,1};
Physical Volume("merged")={20};
"""

function _identity_generate(src,mesh_dim,workers)
    path,tempio=mktemp()
    write(tempio,src*"\nMesh $mesh_dim;\n");close(tempio)
    A.finalize();A.initialize()
    A.option("General.Terminal",0)
    A.option("General.NumThreads",workers)
    A._PARALLEL_ADOPTED[]=0
    A._PARALLEL_SEED_ADOPTED[]=0
    ex=execute_geo(path;mesh_dim)
    rm(path)
    actual=ex.model.meshing.mesh_authority.mesh
    msh=joinpath(mktempdir(),"out.msh")
    write_mixed_msh(msh,actual;version=4.1)
    return (mixed_crc(actual),read(msh),validate(actual).ok,
        A._PARALLEL_ADOPTED[],A._PARALLEL_SEED_ADOPTED[],ex.model)
end

function _session_bytes(mesh)
    msh=joinpath(mktempdir(),"out.msh")
    write_mixed_msh(msh,mesh;version=4.1)
    return read(msh)
end

@testset "RNG-consuming Region fills detach and stay byte-identical" begin
    cases=(
        ("consumer+three boxes Mesh3",CONSUMEBOX,3),
        ("two consumers+three boxes Mesh3",TWOCONSUMEBOX,3),
        ("two-batch consumers Mesh3",TWOBATCH,3),
    )
    try
        for (name,src,dim) in cases
            scrc,sbytes,sok,sadopt,sseed,_=_identity_generate(src,dim,1)
            tcrc,tbytes,tok,tadopt,tseed,tmodel=_identity_generate(src,dim,4)
            @test sok
            @test tok
            @test scrc==tcrc
            @test sbytes==tbytes
            @test sadopt==0                       # -t1 never detaches
            @test sseed==0
            @test tadopt>=3                       # consumer+boxes adopted
            @test tseed>=1                        # a consuming produce adopted
            scrc==tcrc && sbytes==tbytes ||
                println("PARALLEL-DIFFERS: $name")
            # Eligibility pin: the explicit-PLC Region is predicted a
            # consumer (`may_consume`), the factory Boxes are not; the
            # predicted-draw table assigns the PLC fill its `3*count+1`
            # lower bound and the Boxes a provable zero.
            plan=A._parallel_volume_prepare(tmodel,
                sort!(collect(keys(tmodel.volumes))),
                tmodel.meshing.mesh_authority,"test";
                size_field=nothing,scope=nothing,extruded=nothing,
                hybrid_batch=nothing,staged=false)
            @test plan!==nothing
            prepared,pos_draws=plan
            works=Any[w for w in prepared if w!==nothing]
            @test any(w->w.may_consume,works)
            @test any(w->!w.may_consume,works)
            @test all(i->prepared[i]===nothing ||
                          (prepared[i].may_consume ==
                           (pos_draws[i]===nothing || pos_draws[i]>0)),
                eachindex(pos_draws))
        end
    finally
        A.finalize()
    end
end

@testset "identity and merge transports observe the same post-detach state" begin
    # Nondefault identity state exercised through the detach cycle: a
    # parallel-generated model is retagged, renamed, merged with a file and
    # grown again — every identity store and the produced mesh must match
    # the serially-produced session exactly.  The session model is built
    # through `model.add_*` so `model.set_tag`/`api.merge`/`add_box` all run
    # inside the same live session the detached produce merged into.
    function run(workers)
        A.finalize();A.initialize()
        A.option("General.Terminal",0)
        A.option("General.NumThreads",workers)
        # Explicit PLC volume — the native-initial fill consumes the shared
        # RNG ledger — plus two factory boxes (rng-free).
        p=[A.model.add_point(x,y,z;meshSize=0.35)
            for (x,y,z) in ((0,0,0),(1,0,0),(1,1,0),(0,1,0),
                            (0,0,1),(1,0,1),(1,1,1),(0,1,1))]
        l=[A.model.add_line(p[a],p[b]) for (a,b) in
            ((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),
             (1,5),(2,6),(3,7),(4,8))]
        s=[A.model.add_plane_surface([A.model.add_curve_loop(loop_)])
            for loop_ in ([l[1],l[2],l[3],l[4]],[l[5],l[6],l[7],l[8]],
                [l[1],l[10],-l[5],-l[9]],[l[2],l[11],-l[6],-l[10]],
                [l[3],l[12],-l[7],-l[11]],[l[4],l[9],-l[8],-l[12]])]
        v1=A.model.add_volume([A.model.add_surface_loop(s)])
        A.model.add_box(3,0,0,1,1,1)
        A.model.add_box(0,3,0,1,1,1)
        A._PARALLEL_ADOPTED[]=0
        A.mesh.generate(3)
        mesh=A.CURRENT[].meshing.mesh_authority.mesh
        adopted=A._PARALLEL_ADOPTED[]
        first=(mixed_crc(mesh),_session_bytes(mesh),validate(mesh).ok,adopted)
        A.model.set_tag(3,v1,10)
        A.model.set_entity_name(3,10,"domain")
        A.model.set_entity_name(2,s[5],"shell")
        A.mesh.generate(3)
        mesh=A.CURRENT[].meshing.mesh_authority.mesh
        regen=(mixed_crc(mesh),_session_bytes(mesh),validate(mesh).ok)
        merge_dir=mktempdir()
        merge_path=joinpath(merge_dir,"merge.geo")
        write(merge_path,MERGEGEO)
        A.merge(merge_path)
        rm(merge_dir;recursive=true)
        post=A.CURRENT[]
        # A tag-free Box must allocate the next volume-dimension tag —
        # `next_tag` is per-dimension like upstream's `getMaxTag(dim)`, so
        # the volume allocator ignores higher curve/surface counters.  A
        # detached produce that minted a model-entity tag would leak past
        # `before[4]` and push this allocation higher.
        before=copy(post.next_tag)
        new_tag=A.model.add_box(0,-3,0,1,1,1)
        A.mesh.generate(3)
        final_mesh=A.CURRENT[].meshing.mesh_authority.mesh
        return (first=first,regen=regen,
            final=(mixed_crc(final_mesh),_session_bytes(final_mesh),
                validate(final_mesh).ok),
            next_tag=copy(A.CURRENT[].next_tag),before=before,
            entity_names=copy(post.entity_names),
            physical=copy(post.physical),
            physical_names=copy(post.physical_names),
            curve_params=copy(post.curve_params),
            volumes=sort!(collect(keys(post.volumes))),
            new_tag=new_tag)
    end
    try
        serial=run(1);threaded=run(4)
        @test serial.first[1]==threaded.first[1]
        @test serial.first[2]==threaded.first[2]
        @test serial.first[3] && threaded.first[3]
        @test serial.first[4]==0           # -t1 never detaches
        @test threaded.first[4]>=3         # consumer+2 boxes adopted
        @test serial.regen[1]==threaded.regen[1]
        @test serial.regen[2]==threaded.regen[2]
        @test serial.regen[3] && threaded.regen[3]
        @test serial.next_tag==threaded.next_tag
        @test serial.entity_names==threaded.entity_names
        @test serial.entity_names[(3,10)]=="domain"
        @test serial.entity_names[(2,5)]=="shell"
        @test serial.physical==threaded.physical
        @test serial.physical_names==threaded.physical_names
        @test "merged" in values(serial.physical_names)
        @test serial.curve_params==threaded.curve_params
        @test 20 in serial.volumes
        @test serial.new_tag==threaded.new_tag
        @test serial.new_tag==serial.before[4]+1
        @test serial.final[1]==threaded.final[1]
        @test serial.final[2]==threaded.final[2]
        @test serial.final[3] && threaded.final[3]
    finally
        A.finalize()
    end
end
