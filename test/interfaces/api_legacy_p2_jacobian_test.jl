using Test, Tessella
using Tessella.Elements: MixedMesh, ElementBlock

const _LP2J=Tessella.API

function _lp2j_source(dimension)
    base="""
    Mesh.TransfiniteTri=1;
    Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={0,1,0,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,1};
    Curve Loop(1)={1,2,3};Plane Surface(1)={1};
    Transfinite Curve{:}=2;Transfinite Surface{1};
    """
    dimension==3 && (base*="""
    Extrude{0,0,1}{Surface{1};Layers{3};Recombine;QuadTriNoNewVerts;}
    """)
    return base
end

@testset "Legacy quadratic entity filtering precedes task slicing" begin
    try
        _LP2J.initialize()
        _LP2J.model.add_discrete_entity(2,10)
        _LP2J.model.add_discrete_entity(2,20)
        points=zeros(3,12)
        cells=reshape(Int32.(1:12),3,4)
        for column in 1:4
            offset=3(column-1)
            points[:,offset+1]=[2column,0.,0.]
            points[:,offset+2]=[2column+1,0.,0.]
            points[:,offset+3]=[2column,1.,0.]
        end
        fixture=Mesh(points;tris=cells)
        entities=Tuple{Int,Int32}[(2,10),(2,20)]
        owners=Int32[10,20,10,20]
        class=_LP2J._MeshClassification(fixture,(2,Int32(10)),entities,
            [(2,owner) for owner in repeat(owners;inner=3)],
            Dict{Tuple{Int,Int32},Vector{Int32}}(),Int32[],owners,Int32[])
        lock(_LP2J.STATE_LOCK) do
            _LP2J._replace_mesh_cache_locked!(fixture,class)
        end
        _LP2J.mesh.set_order(2)
        uvw=[.125,.25,0.]
        for (entity,tags) in ((10,(1,3)),(20,(2,4)))
            full=_LP2J.mesh.get_jacobians(9,uvw,entity)
            scalar=map(tag->_LP2J.mesh.get_jacobian(tag,uvw),tags)
            @test full==ntuple(i->vcat((result[i] for result in scalar)...),3)
            for task in 0:1
                @test _LP2J.mesh.get_jacobians(9,uvw,entity,task,2)==scalar[task+1]
            end
        end
    finally
        _LP2J.finalize()
    end
end

# Differentiate the barycentric polynomials directly using exact represented
# coordinates. This reference does not call the production nodal engine.
function _lp2j_reference(vertices,uvw,dimension)
    p=Rational{BigInt}.(vertices)
    u,v,w=Rational{BigInt}.(uvw)
    lambda=dimension==2 ? (1-u-v,u,v) : (1-u-v-w,u,v,w)
    gradients=dimension==2 ? ((-1,-1,0),(1,0,0),(0,1,0)) :
        ((-1,-1,-1),(1,0,0),(0,1,0),(0,0,1))
    edges=dimension==2 ? ((1,2),(2,3),(3,1)) :
        ((1,2),(2,3),(3,1),(1,4),(3,4),(2,4))
    values=[l*(2l-1) for l in lambda]
    derivatives=[ntuple(k->(4lambda[i]-1)*gradients[i][k],3)
                 for i in eachindex(lambda)]
    for (a,b) in edges
        push!(values,4lambda[a]*lambda[b])
        push!(derivatives,ntuple(k->4*(lambda[a]*gradients[b][k]+
                                      lambda[b]*gradients[a][k]),3))
    end
    mapped=[Float64(sum(p[k,i]*values[i] for i in eachindex(values))) for k in 1:3]
    columns=ntuple(axis->ntuple(k->sum(p[k,i]*derivatives[i][axis]
                              for i in eachindex(values)),3),3)
    a,b,c=columns
    cross=(a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
    if dimension==2
        magnitude=sqrt(sum(Float64(x)^2 for x in cross))
        return vcat(Float64.(collect(a)),Float64.(collect(b)),
                    Float64.(collect(cross))./magnitude),magnitude,mapped
    end
    determinant=a[1]*(b[2]*c[3]-b[3]*c[2])-
        a[2]*(b[1]*c[3]-b[3]*c[1])+a[3]*(b[1]*c[2]-b[2]*c[1])
    return vcat((Float64.(collect(column)) for column in columns)...),
           Float64(determinant),mapped
end

function _lp2j_nodes(connectivity)
    return reduce(hcat,(first(_LP2J.mesh.get_node(tag)) for tag in connectivity))
end

@testset "Legacy P2 Jacobians evaluate actual retained nodal maps" begin
    for dimension in (2,3)
        try
            _LP2J.initialize()
            mktempdir() do directory
                path=joinpath(directory,"legacy.geo")
                write(path,_lp2j_source(dimension))
                _LP2J.open_geo!(path;mesh_dim=0)
                primary=_LP2J.mesh.generate(dimension)
                @test primary isa Mesh
                _LP2J.mesh.set_order(2)
                msh=dimension==2 ? 9 : 11
                old_msh=dimension==2 ? 2 : 4
                tags,flat=_LP2J.mesh.get_elements_by_type(msh,1)
                n=dimension==2 ? 6 : 10
                cells=reshape(flat,n,:)
                @test !isempty(tags)
                midpoint=cells[dimension==2 ? 4 : 5,1]
                location,_,_,_=_LP2J.mesh.get_node(midpoint)
                delta=dimension==2 ? [.0,.0,.05] :
                    0.05 .* (_lp2j_nodes(cells[:,1])[:,4]-_lp2j_nodes(cells[:,1])[:,1])
                _LP2J.mesh.set_node(midpoint,location+delta,Float64[])
                uvw=[.125,.25,dimension==2 ? 0. : .125,
                     .25,.125,dimension==2 ? 0. : .25]
                bulk=_LP2J.mesh.get_jacobians(msh,uvw)
                @test length.(bulk)==(18length(tags),2length(tags),6length(tags))
                @test _LP2J.mesh.get_jacobians(msh,uvw,1)==bulk
                for (cell,tag) in enumerate(tags)
                    actual=_LP2J.mesh.get_jacobian(tag,uvw)
                    @test actual==(bulk[1][18cell-17:18cell],
                                   bulk[2][2cell-1:2cell],bulk[3][6cell-5:6cell])
                    vertices=_lp2j_nodes(cells[:,cell])
                    mixed=MixedMesh(vertices,[ElementBlock(msh,
                        reshape(Int32.(1:n),n,1))])
                    @test actual==Tessella.API.mesh_jacobians(mixed,msh,uvw)
                    nodal_keys=_LP2J.mesh.get_keys_for_element(tag,"Lagrange")
                    for space in ("Lagrange1","Lagrange3","GradLagrange1")
                        @test _LP2J.mesh.get_keys_for_element(tag,space)==nodal_keys
                    end
                    for point in 1:2
                        reference=_lp2j_reference(vertices,uvw[3point-2:3point],dimension)
                        @test actual[1][9point-8:9point]≈reference[1] atol=2e-14 rtol=2e-14
                        @test actual[2][point]≈reference[2] atol=2e-14 rtol=2e-14
                        @test actual[3][3point-2:3point]≈reference[3] atol=2e-14 rtol=2e-14
                    end
                end
                tasks=[_LP2J.mesh.get_jacobians(msh,uvw,1,task,4) for task in 0:3]
                @test ntuple(i->vcat((task[i] for task in tasks)...),3)==bulk
                empty=(Float64[],Float64[],Float64[])
                @test _LP2J.mesh.get_jacobians(msh,uvw,-1,4,4)==empty
                @test _LP2J.mesh.get_jacobians(msh,Float64[])==empty
                @test _LP2J.mesh.get_jacobians(old_msh,uvw)==empty
                @test _LP2J.mesh.get_jacobians(old_msh,uvw,1)==empty
                points,weights=_LP2J.mesh.get_integration_points(msh,"Gauss4")
                @test length.(_LP2J.mesh.get_jacobians(msh,points))==
                    (9length(weights)*length(tags),length(weights)*length(tags),
                     3length(weights)*length(tags))
                before=(_LP2J.mesh.get_nodes(),_LP2J.mesh.get_elements())
                cache=_LP2J.LAST_MESH[];overlay=_LP2J.LAST_MESH_HIGH_ORDER[]
                for call in (
                    ()->_LP2J.mesh.get_jacobians(msh,[0,0]),
                    ()->_LP2J.mesh.get_jacobians(old_msh,[NaN,0,0]),
                    ()->_LP2J.mesh.get_jacobians(msh,uvw,999),
                    ()->_LP2J.mesh.get_jacobians(msh,uvw,-1,-1,1),
                    ()->_LP2J.mesh.get_jacobian(first(tags),[Inf,0,0]))
                    @test_throws ArgumentError call()
                    @test (_LP2J.mesh.get_nodes(),_LP2J.mesh.get_elements())==before
                    @test _LP2J.LAST_MESH[]===cache && _LP2J.LAST_MESH_HIGH_ORDER[]===overlay
                end
                bulk[1][1]=99.
                @test first(_LP2J.mesh.get_jacobian(first(tags),uvw))[1]!=99.
                _LP2J.mesh.reverse_elements([first(tags)])
                reversed_cells=_LP2J.mesh.get_element(first(tags))[2]
                reversed_nodes=_lp2j_nodes(reversed_cells)
                for point in 1:2
                    reference=_lp2j_reference(reversed_nodes,uvw[3point-2:3point],dimension)
                    actual=_LP2J.mesh.get_jacobian(first(tags),uvw)
                    @test actual[1][9point-8:9point]≈reference[1] atol=2e-14 rtol=2e-14
                    @test actual[2][point]≈reference[2] atol=2e-14 rtol=2e-14
                end
                _LP2J.mesh.set_order(1)
                @test _LP2J.mesh.get_jacobians(msh,uvw)==empty
                @test !isempty(_LP2J.mesh.get_jacobians(old_msh,uvw)[1])
            end
        finally
            _LP2J.finalize()
        end
    end
end

@testset "Legacy tet overlay keeps published lower linear Jacobians and tag offsets" begin
    try
        _LP2J.initialize()
        fixture=Mesh([0. 1 0 0;0 0 1 0;0 0 0 1];
            segs=reshape(Int32[1,2],2,1),tris=reshape(Int32[1,2,3],3,1),
            tets=reshape(Int32[1,2,3,4],4,1))
        class=_LP2J._MeshClassification(fixture,(3,Int32(7)),[(3,Int32(7))],
            fill((3,Int32(7)),4),Dict{Tuple{Int,Int32},Vector{Int32}}(),
            Int32[9],Int32[8],Int32[7])
        lock(_LP2J.STATE_LOCK) do
            _LP2J._replace_mesh_cache_locked!(fixture,class)
        end
        _LP2J.mesh.set_order(2)
        @test _LP2J.mesh.get_element(1)[1]==1
        @test _LP2J.mesh.get_element(2)[1]==2
        @test _LP2J.mesh.get_element(3)[1]==11
        uvw=[.125,.25,.125]
        for (msh,tag) in ((1,1),(2,2))
            reference=Tessella.MeshReferenceGeometry.mesh_jacobian(fixture,tag,uvw)
            @test _LP2J.mesh.get_jacobian(tag,uvw)==reference
            @test _LP2J.mesh.get_jacobians(msh,uvw)==reference
        end
        @test _LP2J.mesh.get_jacobians(4,uvw)==(Float64[],Float64[],Float64[])
        @test _LP2J.mesh.get_jacobian(3,uvw)==_LP2J.mesh.get_jacobians(11,uvw)
    finally
        _LP2J.finalize()
    end
end

@testset "Node queries select the actual parent family on raw and mixed caches" begin
    for (linear,quadratic) in ((1,8),(2,9),(3,10),(4,11),(5,12),(6,13),(7,14))
        try
            _LP2J.initialize()
            _LP2J.option("Mesh.Renumber",0)
            dimension=Int(Tessella.Elements.msh_dimension(quadratic))
            coordinates=Tessella.Elements.lagrange_nodes(quadratic)
            tags=UInt64.(1001:1000+size(coordinates,2))
            _LP2J.model.add_discrete_entity(dimension,7)
            _LP2J.mesh.add_nodes(dimension,7,tags,vec(coordinates))
            _LP2J.mesh.add_elements_by_type(7,quadratic,[7001],tags)
            expected=(tags,vec(coordinates),Float64[])
            for requested in (linear,quadratic)
                @test _LP2J.mesh.get_nodes_by_element_type(requested,7,false)==expected
            end
            _LP2J.mesh.generate(0)
            @test _LP2J.mesh.get() isa MixedMesh
            for requested in (linear,quadratic)
                @test _LP2J.mesh.get_nodes_by_element_type(requested,7,false)==expected
            end
            @test _LP2J.mesh.get_elements_by_type(linear)==(UInt64[],UInt64[])
            @test _LP2J.mesh.get_elements_by_type(quadratic)==(UInt64[7001],tags)
            nodal_keys=_LP2J.mesh.get_keys(quadratic,"Lagrange",7)
            @test nodal_keys==(zeros(Int32,length(tags)),tags,vec(coordinates))
            for space in ("Lagrange1","Lagrange3","GradLagrange1")
                @test _LP2J.mesh.get_keys(quadratic,space,7)==nodal_keys
                @test _LP2J.mesh.get_keys_for_element(7001,space)==nodal_keys
            end
            @test _LP2J.mesh.get_number_of_keys(quadratic,"Lagrange1")==
                Tessella.Elements.msh_spec(linear).nnodes
            @test _LP2J.mesh.get_number_of_keys(quadratic,"Lagrange3")==
                size(Tessella.Elements.lagrange_nodes(Tessella.Elements.msh_type(
                    Tessella.Elements.msh_spec(linear).family,3)),2)
        finally
            _LP2J.finalize()
        end
    end
end
