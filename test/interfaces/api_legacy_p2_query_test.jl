using Test, Tessella, LinearAlgebra

const _LP2Q=Tessella.API

function _lp2q_source(dimension)
    source="""
    Mesh.TransfiniteTri=1;
    Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={0,1,0,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,1};
    Curve Loop(1)={1,2,3};Plane Surface(1)={1};
    Transfinite Curve{:}=2;Transfinite Surface{1};
    """
    dimension==3 && (source*="""
    Extrude{0,0,1}{Surface{1};Layers{1};Recombine;QuadTriNoNewVerts;}
    """)
    return source
end

_lp2q_points(nodes)=reduce(hcat,(first(_LP2Q.mesh.get_node(node)) for node in nodes))

# Exact barycentric polynomials give an independent forward map. In particular,
# these expectations do not query the implementation's reference-map engine.
function _lp2q_map(points,q,dimension)
    p=Rational{BigInt}.(points)
    u,v,w=Rational{BigInt}.(q)
    lambda=dimension==2 ? (1-u-v,u,v) : (1-u-v-w,u,v,w)
    edges=dimension==2 ? ((1,2),(2,3),(3,1)) :
        ((1,2),(2,3),(3,1),(1,4),(3,4),(2,4))
    values=[l*(2l-1) for l in lambda]
    append!(values,[4lambda[a]*lambda[b] for (a,b) in edges])
    return [Float64(sum(p[axis,node]*values[node] for node in eachindex(values)))
            for axis in 1:3]
end

function _lp2q_orientation(nodes)
    remaining=sort(collect(nodes))
    rank=0
    for (position,node) in enumerate(nodes)
        index=findfirst(==(node),remaining)
        rank+=(index-1)*factorial(length(nodes)-position)
        deleteat!(remaining,index)
    end
    return Int32(rank)
end

function _lp2q_patterns(dimension,primary)
    edges=dimension==2 ? ((1,2,4),(2,3,5),(3,1,6)) :
        ((1,2,5),(2,3,6),(3,1,7),(4,1,8),(4,3,9),(4,2,10))
    faces=dimension==2 ? ((1,2,3,4,5,6),) :
        ((1,3,2,7,6,5),(1,2,4,5,10,8),
         (1,4,3,8,9,7),(4,2,3,10,6,9))
    return primary ? (Tuple(edge[1:2] for edge in edges),
                      Tuple(face[1:3] for face in faces)) : (edges,faces)
end

_lp2q_pattern_nodes(cells,patterns)=UInt64[
    cells[slot,column] for column in axes(cells,2)
    for pattern in patterns for slot in pattern]

function _lp2q_snapshot()
    return (_LP2Q.LAST_MESH[],_LP2Q.LAST_MESH_HIGH_ORDER[],
            _LP2Q.mesh.get_nodes(),_LP2Q.mesh.get_elements(),
            _LP2Q.NODE_TAG_MAX[],_LP2Q.ELEMENT_TAG_MAX[])
end

function _lp2q_unchanged(before)
    @test _LP2Q.LAST_MESH[]===before[1]
    @test _LP2Q.LAST_MESH_HIGH_ORDER[]===before[2]
    @test (_LP2Q.mesh.get_nodes(),_LP2Q.mesh.get_elements(),
           _LP2Q.NODE_TAG_MAX[],_LP2Q.ELEMENT_TAG_MAX[])==before[3:6]
end

@testset "Legacy simplex P2 queries use actual public cells" begin
    for dimension in (2,3),reverse in (false,true)
        try
            _LP2Q.initialize()
            mktempdir() do directory
                path=joinpath(directory,"legacy_query.geo")
                write(path,_lp2q_source(dimension))
                _LP2Q.open_geo!(path;mesh_dim=0)
                @test _LP2Q.mesh.generate(dimension) isa Mesh
                _LP2Q.mesh.set_order(2)
                msh,old_msh,nprimary,nn=dimension==2 ? (9,2,3,6) : (11,4,4,10)
                tags,flat=_LP2Q.mesh.get_elements_by_type(msh,1)
                reverse && _LP2Q.mesh.reverse_elements(tags)
                tags,flat=_LP2Q.mesh.get_elements_by_type(msh,1)
                cells=reshape(flat,nn,:)
                @test !isempty(tags)
                first_cell=copy(cells[:,1])
                primary=_lp2q_points(first_cell[1:nprimary])
                midpoint=first_cell[nprimary+1]
                # Moving edge AB's midpoint toward the final corner yields
                # x=A+sum(q_i E_i)+.2*u*lambda_A*E_last. Its full determinant
                # is D*(1-.2u), strictly oriented throughout the simplex.
                old=first(_LP2Q.mesh.get_node(midpoint))
                _LP2Q.mesh.set_node(midpoint,old+.05*(primary[:,end]-primary[:,1]))
                actual=_lp2q_points(first_cell)
                q=dimension==2 ? [.2,.2,0.] : [.2,.2,.2]
                point=_lp2q_map(actual,q,dimension)
                located=_LP2Q.mesh.get_element_by_coordinates(point...,dimension,true)
                @test located[1:3]==(first(tags),Int32(msh),first_cell)
                @test collect(located[4:6])≈q atol=2e-12 rtol=2e-12
                @test _LP2Q.mesh.get_elements_by_coordinates(point...,dimension,true)==[first(tags)]
                @test collect(_LP2Q.mesh.get_local_coordinates_in_element(first(tags),point...))≈q atol=2e-12 rtol=2e-12
                if dimension==3 && !reverse
                    # This top-face point lies inside the linear skeleton but
                    # outside the actual curved cell. The exact same Tet10
                    # payload is rejected by pinned Gmsh 4.15.2 strict search.
                    outside=_lp2q_map(actual,[.25,.25,-.025],3)
                    @test_throws ArgumentError _LP2Q.mesh.get_element_by_coordinates(outside...,3,true)
                    @test_throws ArgumentError _LP2Q.mesh.get_elements_by_coordinates(outside...,3,true)
                end

                points=reduce(hcat,(_lp2q_points(cells[:,column]) for column in axes(cells,2)))
                @test _LP2Q.mesh.get_nodes_by_element_type(msh,1,false)==
                    (flat,vec(points),Float64[])
                packed=_LP2Q.mesh.get_nodes_by_element_type(msh,1,true)
                expected_parameters=reduce(vcat,(begin
                    _,parameters,owner_dimension,_=_LP2Q.mesh.get_node(node)
                    owner_dimension in (1,2) ? parameters : Float64[]
                end for node in flat);init=Float64[])
                @test packed[1:2]==(flat,vec(points))
                @test packed[3]≈expected_parameters atol=2e-12 rtol=2e-12
                keys=_LP2Q.mesh.get_keys(msh,"Lagrange",1)
                @test keys==(zeros(Int32,length(flat)),flat,vec(points))
                @test _LP2Q.mesh.get_keys(msh,"Lagrange",1,false)==
                    (keys[1],keys[2],Float64[])
                @test _LP2Q.mesh.get_number_of_keys(msh,"Lagrange")==nn
                @test _LP2Q.mesh.get_keys_information(keys[1],keys[2],msh,"Lagrange")==
                    fill((Int32(0),Int32(-1)),length(flat))
                for (column,tag) in enumerate(tags)
                    @test _LP2Q.mesh.get_keys_for_element(tag,"Lagrange")==
                        (zeros(Int32,nn),copy(cells[:,column]),vec(_lp2q_points(cells[:,column])))
                    @test _LP2Q.mesh.get_basis_functions_orientation_for_element(tag,"Lagrange")==0
                    @test _LP2Q.mesh.get_basis_functions_orientation_for_element(tag,"H1Legendre2")==
                        _lp2q_orientation(cells[1:nprimary,column])
                end
                @test _LP2Q.mesh.get_basis_functions_orientation(msh,"Lagrange",1)==zeros(Int32,length(tags))
                @test _LP2Q.mesh.get_basis_functions_orientation(msh,"H1Legendre2",1)==
                    [_lp2q_orientation(cells[1:nprimary,column]) for column in axes(cells,2)]
                for primary_only in (false,true),fast in (false,true)
                    count=primary_only ? nprimary : nn
                    expected=reduce(vcat,(begin
                        xyz=_lp2q_points(cells[1:count,column])
                        vec(sum(xyz;dims=2))./(fast ? 1 : count)
                    end for column in axes(cells,2)))
                    @test _LP2Q.mesh.get_barycenters(msh,1,fast,primary_only)≈expected atol=2e-14 rtol=2e-14
                    tasks=[_LP2Q.mesh.get_barycenters(msh,1,fast,primary_only,task,4) for task in 0:3]
                    @test vcat(tasks...)≈expected atol=2e-14 rtol=2e-14
                end
                for primary_only in (false,true)
                    edges,faces=_lp2q_patterns(dimension,primary_only)
                    expected_edges=_lp2q_pattern_nodes(cells,edges)
                    expected_faces=_lp2q_pattern_nodes(cells,faces)
                    @test _LP2Q.mesh.get_element_edge_nodes(msh,1,primary_only)==expected_edges
                    @test _LP2Q.mesh.get_element_face_nodes(msh,3,1,primary_only)==expected_faces
                    @test _LP2Q.mesh.get_element_face_nodes(msh,4,1,primary_only)==UInt64[]
                    @test vcat((_LP2Q.mesh.get_element_edge_nodes(msh,1,primary_only,task,4) for task in 0:3)...)==expected_edges
                    @test vcat((_LP2Q.mesh.get_element_face_nodes(msh,3,1,primary_only,task,4) for task in 0:3)...)==expected_faces
                end
                tasks=[_LP2Q.mesh.get_basis_functions_orientation(msh,"Lagrange",1,task,4) for task in 0:3]
                @test vcat(tasks...)==zeros(Int32,length(tags))

                # Exact-type queries omit the replaced primary cells. Gmsh's
                # nodes-by-type query instead selects the parent family and
                # visits every actual interpolation node, regardless of order.
                @test _LP2Q.mesh.get_elements_by_type(old_msh)==(UInt64[],UInt64[])
                @test _LP2Q.mesh.get_nodes_by_element_type(old_msh,-1,false)==
                    (flat,vec(points),Float64[])
                @test _LP2Q.mesh.get_keys(old_msh,"Lagrange")== (Int32[],UInt64[],Float64[])
                @test _LP2Q.mesh.get_basis_functions_orientation(old_msh,"Lagrange")==Int32[]
                @test _LP2Q.mesh.get_barycenters(old_msh,-1,false,false)==Float64[]
                @test _LP2Q.mesh.get_element_edge_nodes(old_msh)==UInt64[]
                @test _LP2Q.mesh.get_element_face_nodes(old_msh,3)==UInt64[]

                if dimension==3
                    determinant=det(primary[:,2:4].-primary[:,1])
                    @test abs(determinant)>0
                    expected_min=minimum((determinant,.8determinant))
                    expected_max=maximum((determinant,.8determinant))
                    @test only(_LP2Q.mesh.get_element_qualities([first(tags)],"minDetJac"))≈expected_min atol=2e-13 rtol=2e-13
                    @test only(_LP2Q.mesh.get_element_qualities([first(tags)],"maxDetJac"))≈expected_max atol=2e-13 rtol=2e-13
                    # Gmsh's volume measure intentionally uses primary corners.
                    @test only(_LP2Q.mesh.get_element_qualities([first(tags)],"volume"))≈determinant/6 atol=2e-14 rtol=2e-14
                    @test vcat((_LP2Q.mesh.get_element_qualities(tags,"minDetJac",task,4) for task in 0:3)...)==
                        _LP2Q.mesh.get_element_qualities(tags,"minDetJac")
                else
                    area=norm(cross(primary[:,2]-primary[:,1],primary[:,3]-primary[:,1]))/2
                    # Integrating the exact determinant 1-.2u over the unit
                    # triangle gives 7/15; the area is unchanged by winding.
                    @test only(_LP2Q.mesh.get_element_qualities([first(tags)],"volume"))≈
                        area*14/15 atol=2e-13 rtol=2e-13
                end
                before=_lp2q_snapshot()
                for call in (
                    ()->_LP2Q.mesh.get_nodes_by_element_type(msh,999,false),
                    ()->_LP2Q.mesh.get_keys(msh,"Lagrange",999),
                    ()->_LP2Q.mesh.get_barycenters(msh,1,false,false,-1,1),
                    ()->_LP2Q.mesh.get_element_edge_nodes(msh,1,false,0,0),
                    ()->_LP2Q.mesh.get_local_coordinates_in_element(first(tags),NaN,0,0))
                    @test_throws ArgumentError call()
                    _lp2q_unchanged(before)
                end
                keys[2][1]=typemax(UInt64);keys[3][1]=99.
                @test _LP2Q.mesh.get_keys_for_element(first(tags),"Lagrange")[2]==first_cell
                @test _LP2Q.mesh.get_keys_for_element(first(tags),"Lagrange")[3]==vec(actual)
                _LP2Q.mesh.set_order(1)
                @test _LP2Q.mesh.get_keys(msh,"Lagrange")== (Int32[],UInt64[],Float64[])
                @test _LP2Q.mesh.get_nodes_by_element_type(msh,-1,false)==
                    _LP2Q.mesh.get_nodes_by_element_type(old_msh,-1,false)
                @test length(_LP2Q.mesh.get_nodes_by_element_type(msh,-1,false)[1])==nprimary*length(tags)
                @test !isempty(_LP2Q.mesh.get_keys(old_msh,"Lagrange")[2])
                @test _LP2Q.mesh.get_element(first(tags))[1]==old_msh
            end
        finally
            _LP2Q.finalize()
        end
    end
end

@testset "Legacy Tet10 preserves independently published lower linear cells" begin
    try
        _LP2Q.initialize()
        fixture=Mesh([0. 1 0 0;0 0 1 0;0 0 0 1];
            segs=reshape(Int32[1,2],2,1),tris=reshape(Int32[1,2,3],3,1),
            tets=reshape(Int32[1,2,3,4],4,1))
        lock(_LP2Q.STATE_LOCK) do
            _LP2Q._replace_mesh_cache_locked!(fixture)
        end
        _LP2Q.mesh.set_order(2)
        @test _LP2Q.mesh.get_elements()[1]==[1,2,11]
        for (msh,tag,nodes) in ((1,1,UInt64[1,2]),(2,2,UInt64[1,2,3]))
            @test _LP2Q.mesh.get_elements_by_type(msh)==(UInt64[tag],nodes)
            @test _LP2Q.mesh.get_keys_for_element(tag,"Lagrange")==
                (zeros(Int32,length(nodes)),nodes,vec(fixture.coords[:,Int.(nodes)]))
            @test _LP2Q.mesh.get_keys(msh,"Lagrange")==_LP2Q.mesh.get_keys_for_element(tag,"Lagrange")
            @test _LP2Q.mesh.get_basis_functions_orientation(msh,"Lagrange")==Int32[0]
        end
        @test _LP2Q.mesh.get_elements_by_type(4)==(UInt64[],UInt64[])
        located=_LP2Q.mesh.get_element_by_coordinates(.2,.2,.2,-1,true)
        @test located[1:2]==(UInt64(3),Int32(11))
        @test _LP2Q.mesh.get_element_by_coordinates(.2,.2,0,2,true)[1:2]==(UInt64(2),Int32(2))
    finally
        _LP2Q.finalize()
    end
end

@testset "Actual P2 queries retain supported sparse external labels" begin
    for dimension in (2,3)
        try
            _LP2Q.initialize()
            _LP2Q.option("Mesh.Renumber",0)
            _LP2Q.model.add_discrete_entity(dimension,50)
            msh=dimension==2 ? 9 : 11
            corners=dimension==2 ? [0. 1 0;0 0 1;0 0 0] :
                [0. 1 0 0;0 0 1 0;0 0 0 1]
            edges=dimension==2 ? ((1,2),(2,3),(3,1)) :
                ((1,2),(2,3),(3,1),(1,4),(3,4),(2,4))
            mids=reduce(hcat,((corners[:,a]+corners[:,b])/2 for (a,b) in edges))
            mids[:,1]+=.05*(corners[:,end]-corners[:,1])
            points=hcat(corners,mids)
            labels=UInt64[500,100,900,200,701,302,803,404,605,106][1:size(points,2)]
            _LP2Q.mesh.add_nodes(dimension,50,labels,vec(points))
            _LP2Q.mesh.add_elements_by_type(50,msh,[700],labels)
            _LP2Q.mesh.generate(0)
            q=dimension==2 ? [.2,.2,0.] : [.2,.2,.2]
            point=_lp2q_map(points,q,dimension)
            located=_LP2Q.mesh.get_element_by_coordinates(point...,dimension,true)
            @test located[1:3]==(UInt64(700),Int32(msh),labels)
            @test collect(located[4:6])≈q atol=2e-12 rtol=2e-12
            @test collect(_LP2Q.mesh.get_local_coordinates_in_element(700,point...))≈q atol=2e-12 rtol=2e-12
            @test _LP2Q.mesh.get_nodes_by_element_type(msh,50,false)==(labels,vec(points),Float64[])
            @test _LP2Q.mesh.get_keys(msh,"Lagrange",50)==
                (zeros(Int32,length(labels)),labels,vec(points))
            @test _LP2Q.mesh.get_keys_for_element(700,"Lagrange")==_LP2Q.mesh.get_keys(msh,"Lagrange",50)
            @test _LP2Q.mesh.get_basis_functions_orientation(msh,"H1Legendre2",50)==
                [_lp2q_orientation(labels[1:dimension+1])]
        finally
            _LP2Q.finalize()
        end
    end
end
