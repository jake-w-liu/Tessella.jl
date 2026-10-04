using Test, Tessella, LinearAlgebra

const _LP2QUAL=Tessella.API

function _lp2qual_source(dimension)
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

@testset "Warped Tri6 quality uses each upstream Jacobian convention" begin
    try
        _LP2QUAL.initialize()
        mktempdir() do directory
            path=joinpath(directory,"warped_triangle_quality.geo")
            write(path,_lp2qual_source(2))
            _LP2QUAL.open_geo!(path;mesh_dim=0)
            _LP2QUAL.mesh.generate(2)
            _LP2QUAL.mesh.set_order(2)
            tags,nodes=_LP2QUAL.mesh.get_elements_by_type(9,1)
            @test length(tags)==1
            midpoint=first(_LP2QUAL.mesh.get_node(nodes[4]))
            _LP2QUAL.mesh.set_node(nodes[4],midpoint+[0.,0.,.05])
            quality(name)=only(_LP2QUAL.mesh.get_element_qualities(tags,name))
            # Pinned Gmsh4.15.2 on these exact six nodes. The polynomial is
            # F=(u,v,.2u*(1-u-v)): projected determinant1, Gram norm>1.
            @test quality("volume")≈.5033115062813058 rtol=2e-12
            @test quality("minSICN")≈.8658589084635143 rtol=2e-12
            @test quality("minSIGE")≈.90385683908723 rtol=2e-12
            @test quality("minIsotropy")≈.8490445135141556 rtol=2e-12
            for name in ("minSJ","minDetJac","maxDetJac")
                @test quality(name)≈1. atol=2e-12
            end
            jacobians=_LP2QUAL.mesh.get_jacobian(tags[1],[.2,.2,0.])
            @test only(jacobians[2])≈sqrt(1.008) rtol=2e-12
        end
    finally
        _LP2QUAL.finalize()
    end
end

function _lp2qual_sampled_shape(edges,dimension,name)
    # Direct derivative of F(q)=A+E*q+.2*u*(1-sum(q))*E_last.
    vertices=dimension==2 ? ([0.,0.],[1.,0.],[0.,1.]) :
        ([0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.])
    pairs=dimension==2 ? ((1,2),(2,3),(3,1)) :
        ((1,2),(2,3),(3,1),(1,4),(3,4),(2,4))
    points=vcat(collect(vertices),[(vertices[a]+vertices[b])/2 for (a,b) in pairs])
    minimum(points) do q
        matrix=copy(edges)
        matrix[:,1]+=.2*(1-sum(q)-q[1])*edges[:,end]
        for column in 2:dimension
            matrix[:,column]-=.2q[1]*edges[:,end]
        end
        a,b=matrix[:,1],matrix[:,2]
        if name=="minSICN"
            ideal=dimension==2 ? hcat(a,(2b-a)/sqrt(3)) :
                hcat(a,(2b-a)/sqrt(3),sqrt(1.5)*matrix[:,3]-(a+b)/sqrt(6))
            if dimension==2
                return 2norm(cross(ideal[:,1],ideal[:,2]))/sum(abs2,ideal)
            end
            return 3det(ideal)/(norm(ideal)*norm(det(ideal)*inv(ideal)))
        end
        if dimension==2
            lengths=(norm(a),norm(b),norm(b-a))
            x,y,z=lengths
            return (2/sqrt(3))*norm(cross(a,b))*(1/(x*y)+1/(x*z)+1/(y*z))/3
        end
        c=matrix[:,3]
        lengths=(norm(a),norm(b),norm(c),norm(b-a),norm(c-a),norm(c-b))
        # Twelve products in Gmsh's inverse-gradient-error definition.
        triples=((1,6,2),(1,6,3),(1,6,4),(1,6,5),
                 (2,5,1),(2,5,3),(2,5,4),(2,5,6),
                 (3,4,1),(3,4,2),(3,4,5),(3,4,6))
        return sqrt(2)*det(matrix)*sum(1/prod(lengths[i] for i in triple)
                                         for triple in triples)/12
    end
end

@testset "Legacy P2 quality measures evaluate retained geometry" begin
    for dimension in (2,3),reverse in (false,true)
        try
            _LP2QUAL.initialize()
            mktempdir() do directory
                path=joinpath(directory,"quadratic_quality.geo")
                write(path,_lp2qual_source(dimension))
                _LP2QUAL.open_geo!(path;mesh_dim=0)
                _LP2QUAL.mesh.generate(dimension)
                _LP2QUAL.mesh.set_order(2)
                kind,primary_count,node_count=dimension==2 ? (9,3,6) : (11,4,10)
                tags,nodes=_LP2QUAL.mesh.get_elements_by_type(kind,1)
                reverse && _LP2QUAL.mesh.reverse_elements(tags)
                tags,nodes=_LP2QUAL.mesh.get_elements_by_type(kind,1)
                cells=reshape(nodes,node_count,:)
                cell=cells[:,1];tag=tags[1]
                points=reduce(hcat,(first(_LP2QUAL.mesh.get_node(node)) for node in cell))
                edges=points[:,2:primary_count].-points[:,1]
                primary_det=dimension==2 ? norm(cross(edges[:,1],edges[:,2])) : det(edges)
                _LP2QUAL.mesh.set_node(cell[primary_count+1],
                    points[:,primary_count+1]+.05edges[:,end])
                quality(name)=only(_LP2QUAL.mesh.get_element_qualities([tag],name))
                det_extrema=sort([primary_det,.8primary_det])
                @test quality("minDetJac")≈det_extrema[1] rtol=2e-12
                @test quality("maxDetJac")≈det_extrema[2] rtol=2e-12
                @test quality("minSJ")≈(primary_det>0 ? .8 : -1.0) atol=2e-12
                @test quality("volume")≈(dimension==2 ? 7primary_det/15 : primary_det/6) rtol=2e-12
                for name in ("minSICN","minSIGE")
                    @test quality(name)≈_lp2qual_sampled_shape(edges,dimension,name) rtol=2e-12
                end
                isotropy=quality("minIsotropy")
                @test primary_det<0 ? isotropy==0 : 0<isotropy<=1
                repeated=_LP2QUAL.mesh.get_element_qualities([tag,tag],"volume")
                @test repeated==fill(quality("volume"),2)
                @test _LP2QUAL.mesh.get_element_qualities([tag,tag],"volume",1,2)==[quality("volume")]
                # Midpoint-dependent qualities must return to their original
                # values once the retained geometry is lowered to its corners.
                _LP2QUAL.mesh.set_order(1)
                linear=only(_LP2QUAL.mesh.get_element_qualities([tag],"minSJ"))
                @test linear==(primary_det>0 ? 1.0 : -1.0)
                @test only(_LP2QUAL.mesh.get_element_qualities([tag],"volume"))≈
                    primary_det/(dimension==2 ? 2 : 6) rtol=2e-12
            end
        finally
            _LP2QUAL.finalize()
        end
    end
end
