module PluginDistanceRangeTests
using Test,Tessella
const A=Tessella.API

function fixture(dir,scale;triangle=false,higher=false,builder=false,target_higher=false)
    points=triangle ? [(0.,0.,0.),(scale,0.,0.),(0.,scale,0.),(.2scale,.2scale,scale)] :
        [(0.,0.,0.),(scale,0.,0.),(.5scale,scale,0.)]
    higher && (points=[(0.,0.,0.),(scale,0.,0.),(0.,scale,0.),
        (.5scale,0.,.2scale),(.5scale,.5scale,.1scale),(0.,.5scale,.1scale)])
    ids=UInt64[10i+1 for i in eachindex(points)]
    source_count=length(points)
    if target_higher
        if triangle
            append!(points,[(.5scale,0.,.2scale),(.5scale,.5scale,.1scale),(0.,.5scale,.1scale)])
        else
            push!(points,(.5scale,0.,.5scale))
        end
        ids=UInt64[10i+1 for i in eachindex(points)]
    end
    target_dim=triangle ? 2 : 1
    target_msh=target_higher ? (triangle ? 9 : 8) : (triangle ? 2 : 1)
    target_owner=triangle ? 7 : 1
    source_dim=triangle ? 3 : 2
    source_msh=triangle ? 4 : higher ? 9 : 2
    source_owner=triangle ? 8 : 7
    target_ids=target_higher ? [ids[1:target_dim+1];ids[source_count+1:end]] : ids[1:target_dim+1]
    source_ids=ids[1:source_count]
    if builder
        A.model.add("builder")
        A.model.add_discrete_entity(target_dim,target_owner)
        A.model.add_discrete_entity(source_dim,source_owner)
        A.mesh.add_nodes(source_dim,source_owner,ids,collect(Iterators.flatten(points)))
        A.mesh.add_elements_by_type(target_owner,target_msh,[201],target_ids)
        A.mesh.add_elements_by_type(source_owner,source_msh,[301],source_ids)
    else
        path=joinpath(dir,"distance.msh")
        open(path,"w") do io
            println(io,raw"$MeshFormat");println(io,"2.2 0 8");println(io,raw"$EndMeshFormat")
            println(io,raw"$Nodes");println(io,length(points))
            for (id,p) in zip(ids,points);println(io,id," ",join(p," "));end
            println(io,raw"$EndNodes");println(io,raw"$Elements");println(io,2)
            println(io,"201 ",target_msh," 2 0 ",target_owner," ",join(target_ids," "))
            println(io,"301 ",source_msh," 2 0 ",source_owner," ",join(source_ids," "))
            println(io,raw"$EndElements")
        end
        A.open(path)
    end
    return ids,points
end

function session(f)
    A.initialize()
    try
        A.option("General.Terminal",0)
        mktempdir(f)
    finally
        A.finalize()
    end
end

function output()
    previous=A.view.get_tags()
    A.plugin.run("Distance")
    return only(setdiff(A.view.get_tags(),previous))
end

@testset "Distance primitives preserve finite Euclidean distances across scales" begin
    for triangle in (false,true),scale in (1.,1e-80,1e-160,1e-200,1e100)
        session() do dir
            fixture(dir,scale;triangle)
            before=(A.mesh.get_nodes(),A.mesh.get_elements())
            types,data=A.view.get_list_data(output())
            count=triangle ? 4 : 3
            @test types==[triangle ? "SS" : "ST"]
            @test length(only(data))==4count
            @test only(data)[3count+1:end]≈[zeros(count-1);scale] rtol=8eps(Float64) atol=0.
            @test (A.mesh.get_nodes(),A.mesh.get_elements())==before
        end
    end
end

@testset "Distance reads builder targets and emits primary higher-order records" begin
    # Native SDK4.15.2 emits one first-order ST, 12values, for this Tri6
    # source with a straight Line2 distance target. The source stays curved.
    for builder in (false,true),higher in (false,true)
        session() do dir
            ids,_=fixture(dir,1.;builder,higher)
            before=(A.mesh.get_nodes(),A.mesh.get_elements())
            source=A.view.add("retained source")
            values=higher ? [0.,1.,1.,.25,.5,.25] : [0.,1.,1.]
            A.view.add_model_data(source,0,A.model.get_current(),"NodeData",ids,[[v] for v in values])
            data_before=A.view.get_model_data(source,0)
            tag=output();types,data=A.view.get_list_data(tag)
            @test types==["ST"]
            coords=higher ? [0.,1.,0.,0.,0.,1.,0.,0.,0.] : [0.,1.,.5,0.,0.,1.,0.,0.,0.]
            @test only(data)≈[coords;[0.,0.,1.]]
            @test A.view.get_number(tag,"NbTimeStep")==1.
            @test (A.mesh.get_nodes(),A.mesh.get_elements())==before
            @test A.view.get_model_data(source,0)==data_before
            higher && (@test only(A.view.probe(source,.2,.2,.16)[1])≈.08)
        end
    end
end

@testset "Higher-order Distance targets keep the native primary geometry contract" begin
    # SDK4.15.2 Line3/Tri6 target controls measure the primary chord/face;
    # the ordinary view still interpolates the complete curved target basis.
    for triangle in (false,true),builder in (false,true)
        session() do dir
            ids,points=fixture(dir,1.;triangle,builder,target_higher=true)
            before=(A.mesh.get_nodes(),A.mesh.get_elements())
            source=A.view.add("curved target source")
            data=[[i>(triangle ? 4 : 3) ? 7. : 0.] for i in eachindex(ids)]
            A.view.add_model_data(source,0,A.model.get_current(),"NodeData",ids,data)
            view_before=A.view.get_model_data(source,0)
            tag=output();types,records=A.view.get_list_data(tag)
            count=triangle ? 4 : 3
            @test types==[triangle ? "SS" : "ST"]
            @test only(records)[3count+1:end]≈[zeros(count-1);1.]
            @test A.view.get_number(tag,"NbTimeStep")==1.
            @test (A.mesh.get_nodes(),A.mesh.get_elements())==before
            @test A.view.get_model_data(source,0)==view_before
            @test A.view.probe(source,points[triangle ? 5 : 4]...)[1]≈[7.]
        end
    end
end

@testset "Distance kernels share robust geometry without mesh-size clipping" begin
    a=(0.,0.,0.);b=(1.,0.,0.);c=(0.,1.,0.)
    @test A._plugin_point_seg_dist(.5,1e100,0.,a,b)==1e100
    @test A._plugin_point_tri_dist(.25,.25,1e100,a,b,c)==1e100
    for offset in (0.,1e16,-1e16)
        triangle=((offset,0.,0.),(offset+4,4.,0.),(offset,0.,4.))
        @test A._plugin_point_tri_dist(offset+2,1.,1.,triangle...)≈sqrt(.5) rtol=8eps(Float64)
    end
    @test A._plugin_point_tri_dist(.5,1.,0.,a,b,(2.,0.,0.))==1.
    @test A._plugin_point_seg_dist(0.,1.,0.,a,a)==1.
    before=(precision(BigFloat),rounding(BigFloat))
    setprecision(BigFloat,80) do
        setrounding(BigFloat,RoundUp) do
            @test A._plugin_point_tri_dist(.2e-200,.2e-200,1e-200,
                a,(1e-200,0.,0.),(0.,1e-200,0.))≈1e-200 rtol=8eps(Float64)
            @test (precision(BigFloat),rounding(BigFloat))==(80,RoundUp)
        end
    end
    @test (precision(BigFloat),rounding(BigFloat))==before
    if Base.JLOptions().compile_enabled!=3
        segment()=A._plugin_point_seg_dist(.5,1.,0.,a,b)
        triangle_query()=A._plugin_point_tri_dist(.25,.25,.5,a,b,c)
        segment();triangle_query()
        @test (@allocated segment())==0
        @test (@allocated triangle_query())==0
    end
end
end
