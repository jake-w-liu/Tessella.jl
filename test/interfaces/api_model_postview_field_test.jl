module ModelPostViewFieldChecks
using Test, Tessella, LinearAlgebra
const A=Tessella.API
const F=Tessella.SizeField
const E=Tessella.Elements

function session(body)
    A.initialize()
    try
        A.option("General.Terminal",0)
        body()
    finally
        A.finalize()
    end
end

function model_cell(name,msh;coords=E.lagrange_nodes(msh))
    A.model.add(name)
    dim=E.msh_dimension(msh)
    A.model.add_discrete_entity(dim,7)
    tags=UInt64[100+100i for i in axes(coords,2)]
    A.mesh.add_nodes(dim,7,tags,vec(coords))
    A.mesh.add_elements_by_type(7,msh,[1301],tags)
    return tags,coords
end

function background(view;crop_negative=true,use_closest=true)
    field=A.mesh.field.add("PostView")
    A.mesh.field.set_number(field,"ViewTag",view)
    A.mesh.field.set_number(field,"CropNegativeValues",crop_negative ? 1 : 0)
    A.mesh.field.set_number(field,"UseClosest",use_closest ? 1 : 0)
    A.mesh.field.set_as_background_mesh(field)
    return lock(A.STATE_LOCK) do
        A._session_size_field_locked(A._model_locked();raw=true)
    end
end

@testset "model PostView preserves source identity, aliases and initial step" begin
    session() do
        tags,coords=model_cell("background",2)
        view=A.view.add("sizes")
        A.view.add_model_data(view,0,"background","NodeData",tags,[[1.],[3.],[5.]])
        A.view.add_model_data(view,2,"background","NodeData",tags,[[9.],[9.],[9.]])
        A.view.set_number(view,"TimeStep",2)
        before=A.view.get_model_data(view,0)
        # Identical UIDs in another model must not rebind the background data.
        model_cell("target",2;coords=coords.+[10.,0.,0.])
        target=A.mesh.get_nodes()
        field=background(view)
        @test F.field_value(field,.25,.25,0.)≈2.5
        @test A.model.get_current()=="target"
        @test A.mesh.get_nodes()==target
        @test A.view.get_model_data(view,0)==before
        alias=A.view.add_alias(view)
        A.view.remove(view)
        @test F.field_value(background(alias),.25,.25,0.)≈2.5
        A.view.add_model_data(alias,0,"background","NodeData",tags,[[2.],[2.],[2.]])
        @test F.field_value(field,.25,.25,0.)≈2.5
        @test F.field_value(background(alias),.25,.25,0.)≈2.
        A.model.set_current("background")
        A.mesh.set_node(tags[1],[4.,4.,0.])
        @test F.field_value(field,.25,.25,0.)≈2.5
        A.model.set_current("target")
        late=A.view.add("only later data")
        A.view.add_model_data(late,2,"background","NodeData",tags,[[.25],[.25],[.25]])
        @test F.field_value(background(late),.25,.25,0.)==F.GMSH_MAX_SIZE
    end
end

@testset "model data kinds, vector and tensor fields" begin
    session() do
        tags,_=model_cell("data",4)
        p=(.1,.2,.3)
        scalar=A.view.add("element constant")
        A.view.add_model_data(scalar,0,"data","ElementData",[1301],[[.375]])
        @test F.field_value(background(scalar),p...)≈.375
        nodal=A.view.add("element nodes")
        A.view.add_model_data(nodal,0,"data","ElementNodeData",[1301],[[1.,3.,4.,5.]])
        @test F.field_value(background(nodal),p...)≈3.
        vector=A.view.add("vector")
        A.view.add_model_data(vector,0,"data","NodeData",tags,[[3.,4.,0.] for _ in tags])
        @test F.field_value(background(vector),p...)≈5.
        tensor=A.view.add("metric")
        metric=[16.,0.,0.,0.,4.,0.,0.,0.,1.]
        A.view.add_model_data(tensor,0,"data","ElementNodeData",[1301],[repeat(metric,4)],0.,9)
        field=background(tensor)
        @test field isa F.PostViewAnisoField
        actual=F.metric_at(field,p...)
        @test [actual.m11,actual.m22,actual.m33]≈[16.,4.,1.]
        @test (actual.m12,actual.m13,actual.m23)==(0.,0.,0.)
    end
end

@testset "curved model cells retain precedence over flat boundaries" begin
    session() do
        tags,_=model_cell("background",9)
        A.model.add_discrete_entity(1,5)
        A.mesh.add_elements_by_type(5,1,[1302],tags[1:2])
        view=A.view.add("dimension precedence")
        A.view.add_model_data(view,0,"background","ElementData",[1301,1302],[[.25],[.5]])
        @test F.field_value(background(view),.25,0.,0.)≈.25
    end
end

@testset "full curved geometry and complete pyramid interpolation" begin
    session() do
        # Curved quadratic line x=u, y=1-u^2. The interior node is required
        # to locate its midpoint; the endpoint chord lies on y=0 instead.
        linecoords=copy(E.lagrange_nodes(8))
        linecoords[2,:].=1 .- linecoords[1,:].^2
        tags,_=model_cell("curve",8;coords=linecoords)
        view=A.view.add("curved sizes")
        values=[[2+linecoords[1,i]^2] for i in axes(linecoords,2)]
        A.view.add_model_data(view,0,"curve","NodeData",tags,values)
        @test F.field_value(background(view;use_closest=false),0.,1.,0.)≈2.

        # This curve's maximum is between all stored/order-2 sample nodes.
        # Its actual point (.5,1,0) lies outside their y-range [-8,0].
        bulge=copy(E.lagrange_nodes(8))
        bulge[2,:].=4 .* bulge[1,:].*(1 .- bulge[1,:])
        btags,_=model_cell("bulging curve",8;coords=bulge)
        bview=A.view.add("bulging sizes")
        A.view.add_model_data(bview,0,"bulging curve","NodeData",btags,
            [[2+bulge[1,i]^2] for i in axes(bulge,2)])
        @test F.field_value(background(bview;use_closest=false),.5,1.,0.)≈2.25
        @test F.field_value(background(bview),.5,1.,0.)≈2.25
        shifted=bulge.+[1e15,1e15,0.]
        stags,_=model_cell("shifted curve",8;coords=shifted)
        shiftedview=A.view.add("large offset cancellation")
        A.view.add_model_data(shiftedview,0,"shifted curve","NodeData",stags,
            [[2+bulge[1,i]^2] for i in axes(bulge,2)])
        @test F.field_value(background(shiftedview;use_closest=false),1e15+.5,1e15+1,0.)≈2.25

        # A non-planar bilinear quad must keep its actual isoparametric map.
        quad=[-1. 1. 1. -1.;-1. -1. 1. 1.;.3 -.3 .3 -.3]
        qtags,_=model_cell("warped",3;coords=quad)
        qview=A.view.add("warped values")
        A.view.add_model_data(qview,0,"warped","NodeData",qtags,
            [[2+quad[1,i]+quad[2,i]+quad[3,i]] for i in 1:4])
        @test F.field_value(background(qview;use_closest=false),.2,.3,.018)≈2.518 atol=1e-10

        # Cubic Bergot pyramid: a quadratic value term distinguishes full
        # interpolation from corner-only sampling. The point is in its domain.
        ptags,pcoords=model_cell("pyramid",118)
        pview=A.view.add("pyramid values")
        polynomial(x,y,z)=2+x*x+2y*y+3z*z
        A.view.add_model_data(pview,0,"pyramid","NodeData",ptags,
            [[polynomial(pcoords[:,i]...)] for i in axes(pcoords,2)])
        p=(.1,-.2,.3)
        @test F.field_value(background(pview;use_closest=false),p...)≈polynomial(p...) atol=1e-10
        constant=A.view.add("curved geometry constant")
        A.view.add_model_data(constant,0,"pyramid","ElementData",[1301],[[.125]])
        @test F.field_value(background(constant;use_closest=false),p...)≈.125
    end
end

@testset "sparse, invalid and cropped model rows" begin
    session() do
        tags,_=model_cell("sparse",2)
        sparse=A.view.add("partial nodes")
        A.view.add_model_data(sparse,0,"sparse","NodeData",tags[1:2],[[1.],[2.]])
        @test F.field_value(background(sparse),.25,.25,0.)==F.GMSH_MAX_SIZE
        bad=A.view.add("unknown UID")
        A.view.add_model_data(bad,0,"sparse","NodeData",[999999],[[1.]])
        before=A.mesh.get_nodes()
        @test_throws ArgumentError background(bad)
        @test A.mesh.get_nodes()==before
        short=A.view.add("short element row")
        A.view.add_model_data(short,0,"sparse","ElementNodeData",[1301],[[1.,2.]])
        @test_throws ArgumentError background(short)
        negative=A.view.add("negative constant")
        A.view.add_model_data(negative,0,"sparse","ElementData",[1301],[[-2.]])
        @test F.field_value(background(negative),.25,.25,0.)==F.GMSH_MAX_SIZE
        @test F.field_value(background(negative;crop_negative=false),.25,.25,0.)≈-2.
    end
end

@testset "model PostView drives the native target mesher" begin
    session() do
        tags,_=model_cell("background",2;coords=[-1. 3. -1.;-1. -1. 3.;0. 0. 0.])
        view=A.view.add("quarter size")
        A.view.add_model_data(view,0,"background","NodeData",tags,[[.25] for _ in tags])
        A.model.add("target line")
        first=A.model.geo.add_point(0.,0.,0.)
        last=A.model.geo.add_point(1.,0.,0.)
        A.model.geo.add_line(first,last)
        A.model.geo.synchronize()
        background(view)
        A.option("Mesh.MeshSizeFromPoints",0)
        A.option("Mesh.MeshSizeExtendFromBoundary",0)
        A.mesh.generate(1)
        _,coordinates,_=A.mesh.get_nodes()
        @test sort(reshape(coordinates,3,:)[1,:])≈[0.,.25,.5,.75,1.]
        @test A.model.get_current()=="target line"
    end
end
end
