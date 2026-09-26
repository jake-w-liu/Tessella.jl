using Test
using Tessella
using Tessella.Model: add_discrete_entity!, add_discrete_nodes!,
                      add_discrete_elements!, classify_surfaces!

@testset "discrete Point coordinate updates preserve entity identity" begin
    model=GeoModel()
    tag=add_discrete_entity!(model,0,1)
    add_discrete_nodes!(model,0,tag,[7,8],[1.,2,3,9,8,7])
    @test Tessella.Model.set_point_coordinates!(model,tag,4.,5.,6.)===nothing
    @test Tessella.Model.model_value(model,0,tag,Float64[])==[4.,5.,6.]
    @test !haskey(model.points,tag)
    @test model.discrete[(0,tag)].node_coords[:,2]==[9.,8.,7.]
    before=deepcopy(model.discrete)
    @test_throws ArgumentError Tessella.Model.set_point_coordinates!(model,tag,NaN,0,0)
    @test model.discrete==before
    empty_tag=add_discrete_entity!(model,0,2)
    before=deepcopy(model.discrete)
    @test_throws ArgumentError Tessella.Model.set_point_coordinates!(model,empty_tag,0,0,0)
    @test model.discrete==before
    @test isempty(model.points)
end

@testset "atomic discrete mesh record insertion" begin
    model=GeoModel()
    add_discrete_entity!(model,1,1)
    for (tag,x) in ((101,0.),(102,1.),(103,2.))
        add_point!(model,x,0,0;tag)
    end
    add_line!(model,101,102;tag=2)
    add_line!(model,102,103;tag=3)
    add_discrete_nodes!(model,1,1,[10,20],[0.,0,0,1,0,0],[0.,1.])
    add_discrete_nodes!(model,1,2,[30,40],[2.,0,0,3,0,0])
    add_discrete_elements!(model,1,1,[1],[[50]],[[10,20]])
    add_discrete_elements!(model,1,2,[1],[[60]],[[30,40]])
    state()=deepcopy((model.discrete,model.meshing.attached))

    before=state()
    for entity in (1,2,3), bad in (0,-1,true,big(typemax(Int32))+1,big(2)^100)
        @test_throws ArgumentError add_discrete_nodes!(model,1,entity,
            Any[70,bad],[7.,0,0,8,0,0])
        @test state()==before
    end
    for (entity,conflict) in ((1,30),(2,10),(3,10))
        @test_throws ArgumentError add_discrete_nodes!(model,1,entity,
            [70,conflict],[7.,0,0,8,0,0],[0.,1.])
        @test state()==before
    end
    # An attempted update must also wait until ownership validation completes.
    @test_throws ArgumentError add_discrete_nodes!(model,1,1,[10,30],
        [9.,0,0,8,0,0],[0.5,1.])
    @test state()==before
    @test_throws ArgumentError add_discrete_nodes!(model,1,3,[70],[0.,NaN,0])
    @test state()==before
    @test_throws ArgumentError add_discrete_nodes!(model,1,3,[70],[0.,0,0],[0.,1.])
    @test state()==before

    for (entity,connectivity) in ((1,[10,20]),(2,[30,40]))
        for bad in (0,-1,true,big(typemax(Int32))+1,big(2)^100,50,60)
            @test_throws ArgumentError add_discrete_elements!(model,1,entity,
                [1],[Any[70,bad]],[vcat(connectivity,connectivity)])
            @test state()==before
        end
        @test_throws ArgumentError add_discrete_elements!(model,1,entity,
            [1],[[70,70]],[vcat(connectivity,connectivity)])
        @test state()==before
        @test_throws ArgumentError add_discrete_elements!(model,1,entity,
            [1,1],[[70],[71]],[connectivity,[999,999]])
        @test state()==before
    end
    for (types,tags,nodes) in (([1],[[70]],[[1,2]]),
                               ([999],[[70]],[[1,2]]),
                               ([1],[[70]],[[1]]))
        @test_throws ArgumentError add_discrete_elements!(model,1,3,types,tags,nodes)
        @test state()==before
    end

    # Repeated tags on their current owner retain the documented update behavior.
    add_discrete_nodes!(model,1,1,[10,10,typemax(Int32)],
        [4.,0,0,5,0,0,6,0,0],[0.25,0.5,0.75])
    record=model.discrete[(1,1)]
    @test record.node_tags==Int32[10,20,typemax(Int32)]
    @test record.node_coords==[5. 1 6;0 0 0;0 0 0]
    @test record.node_params==reshape([0.5,1.,0.75],1,3)
    # Element dimension remains independent of its classification entity.
    add_discrete_elements!(model,1,1,[2],[[typemax(Int32)]],
                           [[10,20,typemax(Int32)]])
    @test record.element_tags==Int32[50,typemax(Int32)]
    @test record.element_types==Int32[1,2]
    add_discrete_nodes!(model,1,3,[70],[7.,0,0])
    @test model.meshing.attached[(1,3)].node_tags==Int32[70]
end

@testset "boundary turning angles conserve chain edges" begin
    split=Tessella.Model._split_chains_by_angle
    coords=Dict{Int32,NTuple{3,Float64}}(
        1=>(0.,0.,0.),2=>(1.,0.,0.),3=>(2.,0.,0.),
        4=>(2.,1.,0.),5=>(1.,1.,0.),6=>(0.,1.,0.))
    straight=Tuple{Vector{Int32},Bool}[(Int32[1,2,3],false)]
    @test split(straight,coords,pi/4)==straight
    corner=Tuple{Vector{Int32},Bool}[(Int32[1,2,3,4],false)]
    @test split(corner,coords,pi/4)==[(Int32[1,2,3],false),(Int32[3,4],false)]
    edges(chains)=sort!([minmax(chain[i],chain[i+1])
        for (chain,_) in chains for i in 1:length(chain)-1])
    # Move and reverse the seam through corners and collinear vertices.
    for reversed in (false,true), shift in 0:5
        vertices=circshift(Int32[1,2,3,4,5,6],shift)
        reversed && reverse!(vertices)
        chain=vcat(vertices,vertices[1])
        input=Tuple{Vector{Int32},Bool}[(chain,true)]
        result=split(input,coords,pi/4)
        @test length(result)==4
        @test edges(result)==edges(input)
        @test all(!closed && first(part) in (1,3,4,6) &&
                  last(part) in (1,3,4,6) for (part,closed) in result)
        @test split(input,coords,Float64(pi))==input
    end
end

@testset "surface classification follows adjacent face normals" begin
    for offset in (0,8,16)
        model=GeoModel()
        tag=add_discrete_entity!(model,2,1)
        heights=[0.,0.,tand(10.),tand(10.)+tand(20.)]
        points=[(Float64(x),y,heights[x+1]) for x in 0:3 for y in (0.,1.)]
        add_discrete_nodes!(model,2,tag,collect(1:8).+offset,
                            reduce(vcat,collect.(points)))
        connectivity=Int[]
        for k in 1:3
            a=2k-1+offset;b=a+2;c=a+3;d=a+1
            append!(connectivity,(a,b,c,a,c,d))
        end
        add_discrete_elements!(model,2,tag,[2],[collect(1:6)],[connectivity])
        classify_surfaces!(model;angle=15pi/180,boundary=false)
        surfaces=[record for (key,record) in model.discrete if key[1]==2]
        @test length(surfaces)==1
        @test only(surfaces).element_tags==Int32.(1:6)
    end

    model=GeoModel()
    tag=add_discrete_entity!(model,2,1)
    add_discrete_nodes!(model,2,tag,[1,2,3,4],
                        [0.,0,0,1,0,0,1,1,0,0,1,0])
    add_discrete_elements!(model,2,tag,[2],[[1,2]],[[1,2,3,1,3,4]])
    classify_surfaces!(model;angle=pi/4,curve_angle=pi/4)
    segments=sort!([minmax(nodes...) for (key,record) in model.discrete
                   if key[1]==1 for nodes in record.element_nodes])
    @test segments==[(1,2),(1,4),(2,3),(3,4)]
    @test count(key->key[1]==0,keys(model.discrete))==4
    @test count(key->key[1]==1,keys(model.discrete))==4
    @test count(key->key[1]==2,keys(model.discrete))==1
end
