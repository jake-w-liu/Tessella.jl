using Test,Tessella
const _RECLASS_API=Tessella.API
include(joinpath(@__DIR__,"../artifacts/reclassify_primary_records_v1.jl"))

function _reclass_square()
    for (tag,(x,y)) in enumerate(((0.,0.),(1.,0.),(1.,1.),(0.,1.)))
        _RECLASS_API.model.add_point(x,y,0;tag)
    end
    for (tag,(a,b)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        _RECLASS_API.model.add_line(a,b;tag)
    end
    _RECLASS_API.model.add_curve_loop([1,2,3,4];tag=1)
    _RECLASS_API.model.add_plane_surface([1];tag=1)
    _RECLASS_API.model.add_plane_surface([1];tag=2)
end

function _reclass_fixture(name)
    name=="empty" && return nothing
    if name=="discrete_dimension_priority_and_orphan"
        for dim in 0:3
            _RECLASS_API.model.add_discrete_entity(dim,1)
        end
        _RECLASS_API.mesh.add_nodes(3,1,[11,12,13,14,15],
            [0.,0,0,1,0,0,0,1,0,0,0,1,9,9,9])
        for (msh,labels,cells) in ((4,[21],[11,12,13,14]),(2,[22],[11,12,13]),
                (1,[23],[11,12]),(15,[24],[11]))
            _RECLASS_API.mesh.add_elements_by_type(1,msh,labels,cells)
        end
        return nothing
    end
    _reclass_square()
    if name=="quadratic_actual_support_classification"
        _RECLASS_API.mesh.add_nodes(2,1,collect(101:107),
            [0.,0,0,1,0,0,0,1,0,.5,0,0,.5,.5,0,0,.5,0,9,9,9],
            [.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4])
        for (msh,labels,cells) in ((9,[211],[101,102,103,104,105,106]),
                (8,[212],[101,102,104]),(15,[213],[101]))
            _RECLASS_API.mesh.add_elements_by_type(1,msh,labels,cells)
        end
        return nothing
    end
    params=name=="native_unparameterized" ? Float64[] :
        [.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4]
    _RECLASS_API.mesh.add_nodes(2,1,collect(101:107),
        [0.,0,0,1,0,0,1,1,0,0,1,0,.4,.4,0,.6,.7,0,2,2,0],params)
    _RECLASS_API.mesh.add_elements_by_type(1,2,[111,112],[101,102,105,102,103,105])
    _RECLASS_API.mesh.add_elements_by_type(1,1,[121],[101,102])
    _RECLASS_API.mesh.add_elements_by_type(1,15,[131],[101])
    name=="native_same_dimension_first_entity" &&
        _RECLASS_API.mesh.add_elements_by_type(2,2,[113],[103,104,105])
    return nothing
end

@testset "Reclassification matches independent actual Gmsh records" begin
    for row in _RECLASSIFY_PRIMARY_RECORDS
        @testset "$(row["name"])" begin
            _RECLASS_API.initialize()
            try
                _reclass_fixture(row["name"])
                for repeat in 1:2
                    @test _RECLASS_API.mesh.reclassify_nodes()===nothing
                    for entry in row["after"]
                        dim=entry["dim"];tag=entry["tag"]
                        labels,coords,parameters=_RECLASS_API.mesh.get_nodes(dim,tag)
                        types,elements,connectivity=_RECLASS_API.mesh.get_elements(dim,tag)
                        @test labels==entry["nodes"]
                        @test coords==entry["coordinates"]
                        @test parameters==entry["parameters"]
                        @test types==entry["types"]
                        @test elements==entry["elements"]
                        @test connectivity==entry["connectivity"]
                    end
                    for node in row["nodes"]
                        coords,parameters,dim,owner=_RECLASS_API.mesh.get_node(node["tag"])
                        @test coords==node["coordinates"]
                        @test parameters==node["parameters"]
                        @test (dim,owner)==(node["dim"],node["owner"])
                    end
                end
            finally
                _RECLASS_API.finalize()
            end
        end
    end
end

@testset "Reclassified nodes retain their original parameter type" begin
    _RECLASS_API.initialize()
    try
        _reclass_fixture("native_surface_parameters_to_curve_point")
        _RECLASS_API.mesh.reclassify_nodes()
        _RECLASS_API.mesh.set_node(102,[4.,5,6],[.65])
        @test _RECLASS_API.mesh.get_node(102)==([4.,5,6],[.65,.4],1,1)
        _RECLASS_API.mesh.relocate_nodes(1,1)
        @test _RECLASS_API.mesh.get_node(102)==([.65,0.,0.],[.65,.4],1,1)
        _RECLASS_API.mesh.set_node(101,[4.,5,6],[.75,.85])
        @test _RECLASS_API.mesh.get_node(101)==([4.,5,6],[.75,.85],0,1)
        _RECLASS_API.mesh.relocate_nodes(0,1)
        @test _RECLASS_API.mesh.get_node(101)==([0.,0.,0.],[.75,.85],0,1)
        @test Tessella.MeshTypes.validate(_RECLASS_API.LAST_MESH[]).ok
    finally
        _RECLASS_API.finalize()
    end
end

@testset "One owner can contain different original parameter types" begin
    _RECLASS_API.initialize()
    try
        _reclass_square()
        _RECLASS_API.mesh.add_nodes(2,1,[11],[.2,0.,0.],[.25,.5])
        _RECLASS_API.mesh.add_nodes(1,1,[12],[.6,0.,0.],[.75])
        _RECLASS_API.mesh.add_nodes(1,1,[13],[.8,0.,0.])
        _RECLASS_API.mesh.add_elements_by_type(1,1,[21,22],[11,12,12,13])
        _RECLASS_API.mesh.reclassify_nodes()
        @test _RECLASS_API.mesh.get_nodes(1,1)==(UInt64[11,12,13],
            [.2,0.,0.,.6,0.,0.,.8,0.,0.],[.25,.75])
        @test _RECLASS_API.mesh.get_node(11)==([.2,0.,0.],[.25,.5],1,1)
        @test _RECLASS_API.mesh.get_node(12)==([.6,0.,0.],[.75],1,1)
        @test _RECLASS_API.mesh.get_node(13)==([.8,0.,0.],Float64[],1,1)
        _RECLASS_API.mesh.relocate_nodes(1,1)
        @test _RECLASS_API.mesh.get_node(13)==([.8,0.,0.],Float64[],1,1)
        @test _RECLASS_API.mesh.get_node(11)==([.25,0.,0.],[.25,.5],1,1)
        @test _RECLASS_API.mesh.get_node(12)==([.75,0.,0.],[.75],1,1)
    finally
        _RECLASS_API.finalize()
    end
end

@testset "setNode updates the stored vertex slots independently of owner dimension" begin
    for original in (Float64[],[.25],[.25,.75]),supplied in
            (Float64[],[.6],[.6,.8],[.6,.8,.9])
        _RECLASS_API.initialize()
        try
            dim=length(original)==2 ? 2 : 1
            _RECLASS_API.model.add_discrete_entity(dim,1)
            _RECLASS_API.mesh.add_nodes(dim,1,[111],[1.,2,3],original)
            _RECLASS_API.mesh.set_node(111,[4.,5,6],supplied)
            expected=copy(original)
            for index in 1:min(length(expected),length(supplied))
                expected[index]=supplied[index]
            end
            @test _RECLASS_API.mesh.get_node(111)==([4.,5,6],expected,dim,1)
        finally
            _RECLASS_API.finalize()
        end
    end
end
