module OCCNativeHyperbolaChecks
using Test,Tessella
using LinearAlgebra: norm
const A=Tessella.API
const M=Tessella.Model
const FIXTURES=joinpath(@__DIR__,"..","fixtures")
const LOWER=-.7
const UPPER=.8

point(t)=[2cosh(t),sinh(t),0.]
tangent(t)=[2sinh(t),cosh(t),0.]
second(t)=[2cosh(t),sinh(t),0.]

function line_support(nodes,cells)
    coordinates=reshape(nodes[2],3,:)
    points=Dict(uid=>Tuple(coordinates[:,row]) for (row,uid) in enumerate(nodes[1]))
    return Set(minmax(points[c[i]],points[c[i+1]]) for c in cells[3] for i in 1:2:length(c))
end

function session(f,name="brep_native_hyperbola.brep")
    A.initialize()
    try
        imported=A.model.occ.import_shapes(joinpath(FIXTURES,name))
        @test imported==(name=="brep_native_hyperbola.brep" ? [(1,1)] : [(2,1)])
        A.model.occ.synchronize()
        return f(A.model.occ)
    finally
        A.finalize()
    end
end

@testset "A native BREP hyperbola retains its actual hyperbolic parameter" begin
    session() do occ
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(2,1,0,0)
        @test A.model.get_type(1,1)=="Hyperbola"
        lo,hi=A.model.get_parametrization_bounds(1,1)
        @test collect(lo)==[LOWER] && collect(hi)==[UPPER]
        geometry=A.CURRENT[].curve_geometry[1].curve
        inner=M._ModelOCCNURBSCurve(:bspline,geometry.curve,geometry.nets)
        for t in (LOWER,-.13,0.,.35,UPPER)
            @test A.model.get_value(1,1,[t])≈point(t) atol=3e-14
            @test A.model.get_derivative(1,1,[t])≈tangent(t) atol=3e-14
            @test collect(M._model_occ_curve_d2(geometry,t)[3])≈second(t) atol=3e-14
            a,b=t-.001<=LOWER ? (t,t+.001) :
                t+.001>=UPPER ? (t-.001,t) : (t-.001,t+.001)
            @test A.model.get_second_derivative(1,1,[t])≈
                (tangent(b).-tangent(a))./(b-a) atol=2e-11
            parameter,ds,dds=M._occ_conic_reparam(geometry.reparam,t)
            p,d,dd=M._model_occ_curve_d2(inner,parameter)
            @test collect(p)≈point(t) atol=3e-14
            @test ds.*collect(d)≈tangent(t) atol=3e-14
            @test dds.*collect(d).+(ds^2).*collect(dd)≈second(t) atol=4e-13
        end
        # Independent 60-digit integration of sqrt(4sinh(t)^2+cosh(t)^2).
        @test occ.get_mass(1,1)≈2.0951434453736509710044315769459 rtol=1e-11
        @test collect(occ.get_center_of_mass(1,1))≈
            [2.2394755927147707789248287569951,.082357717432982042789280534307,0.] atol=2e-11
        inertia=[.53663981426191689839,-.06121953968636170811,0.,
                 -.06121953968636170811,.07630662726709414076,0.,
                 0.,0.,.61294644152901103915]
        @test occ.get_matrix_of_inertia(1,1)≈inertia rtol=1e-10 atol=1e-12
        expected_box=[2-1e-7,sinh(LOWER)-1e-7,-1e-7,
            2cosh(UPPER)+1e-7,sinh(UPPER)+1e-7,1e-7]
        @test collect(occ.get_bounding_box(1,1))≈expected_box atol=1e-14
        @test A.model.get_bounding_box(1,1)==occ.get_bounding_box(1,1)
        xyz,parameter=A.model.get_closest_point(1,1,point(.21).+[0.,0.,.4])
        @test xyz≈point(.21) atol=3e-10
        @test collect(parameter)≈[.21] atol=3e-10
        A.option("Mesh.MeshSizeMin",.3);A.option("Mesh.MeshSizeMax",.3)
        A.mesh.generate(1)
        tags,xyz,parameters=A.mesh.get_nodes(1,1)
        @test !isempty(tags) && length(tags)==length(parameters)
        @test A.model.get_value(1,1,parameters)≈xyz atol=2e-12
        families,_,cells=A.mesh.get_elements(1,1)
        @test families==Int32[1] && !isempty(only(cells))
        tags,xyz,_=A.mesh.get_nodes();coordinates=reshape(xyz,3,:)
        rows=Dict(uid=>row for (row,uid) in enumerate(tags))
        @test maximum(norm(coordinates[:,rows[c[i]]].-coordinates[:,rows[c[i+1]]])
            for c in cells for i in 1:2:length(c))<=.3+1e-9
    end
end

@testset "Hyperbola affine images and rational conversion retain held geometry" begin
    session() do occ
        @test occ.copy([(1,1)])==[(1,2)]
        occ.translate([(1,2)],3.,-2.,1.)
        occ.dilate([(1,2)],0.,0.,0.,2.,3.,4.)
        occ.mirror([(1,2)],1.,0.,0.,0.)
        occ.synchronize()
        for t in (-.43,.11,.61)
            p=point(t)
            @test A.model.get_value(1,1,[t])≈p atol=3e-14
            @test A.model.get_value(1,2,[t])≈[-2(p[1]+3),3(p[2]-2),4.] atol=2e-13
            @test A.model.get_derivative(1,2,[t])≈[-2,3,4].*tangent(t) atol=2e-13
        end
        box=occ.get_bounding_box(1,2)
        @test all(box[axis]<=A.model.get_value(1,2,[t])[axis]<=box[axis+3]
            for t in range(LOWER,UPPER;length=101),axis in 1:3)
        held=A.model.get_value(1,2,[.21])
        occ.remove([(1,1)],true);occ.synchronize()
        @test A.model.get_value(1,2,[.21])==held
        occ.convert_to_nurbs([(1,2)]);occ.synchronize()
        @test A.model.get_type(1,2)=="BSpline"
        middle=(LOWER+UPPER)/2;width=UPPER-LOWER;q=tanh(width/4)
        new_parameter=middle+width/2*tanh((.21-middle)/2)/q
        @test A.model.get_value(1,2,[new_parameter])≈held atol=3e-13
        z=(2*(.21-LOWER)/width-1)*q
        source_parameter=middle+2atanh(z)
        derivative=4q/(width*(1-z^2))
        @test A.model.get_derivative(1,2,[.21])≈
            [-2,3,4].*tangent(source_parameter).*derivative atol=3e-13
    end
end

@testset "Explicit hyperbola conversion exposes the actual rational parameter" begin
    session() do occ
        held=A.model.get_value(1,1,[.21])
        occ.convert_to_nurbs([(1,1)]);occ.synchronize()
        @test A.model.get_type(1,1)=="BSpline"
        lo,hi=A.model.get_parametrization_bounds(1,1)
        @test collect(lo)==[LOWER] && collect(hi)==[UPPER]
        @test A.model.get_value(1,1,[.21])≈
            [2.0414317239075026,.20459929336141333,0.] atol=3e-14
        @test A.model.get_derivative(1,1,[.21])≈
            [.3933371162572466,.9811505872352092,0.] atol=3e-14
        geometry=A.CURRENT[].curve_geometry[1].curve
        @test geometry isa M._ModelOCCNURBSCurve
        @test geometry.curve.degree==2 && length(geometry.curve.controls)==3
        closest,parameter=A.model.get_closest_point(1,1,held)
        @test closest≈held atol=3e-11
        @test abs(only(parameter)-.21)>1e-3
    end
end

@testset "A hyperbolic generatrix remains exact on its native extruded face" begin
    session("brep_native_hyperbola_extrusion.brep") do occ
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(4,4,1,0)
        for (u,v) in ((LOWER,0.),(-.13,.3),(.35,.7),(UPPER,1.))
            @test A.model.get_value(2,1,[u,v])≈point(u).+[0.,0.,v] atol=3e-13
            @test A.model.get_derivative(2,1,[u,v])≈vcat(tangent(u),[0.,0.,1.]) atol=3e-13
            @test A.model.get_second_derivative(2,1,[u,v])≈
                vcat(second(u),zeros(6)) atol=6e-12
        end
        A.option("Mesh.MeshSizeMin",.6);A.option("Mesh.MeshSizeMax",.6)
        A.option("Mesh.MeshOnlyEmpty",1);A.option("Mesh.Algorithm",1)
        A.mesh.generate(1)
        prior_tags,prior_xyz,_=A.mesh.get_nodes();prior_coords=reshape(prior_xyz,3,:)
        _,_,groups=A.mesh.get_elements(1)
        lower_edges=Set(minmax(c[i],c[i+1]) for c in groups for i in 1:2:length(c))
        A.mesh.generate(2)
        tags,xyz,_=A.mesh.get_nodes();coordinates=reshape(xyz,3,:)
        rows=Dict(uid=>row for (row,uid) in enumerate(tags))
        @test all(haskey(rows,uid) && coordinates[:,rows[uid]]==prior_coords[:,index]
            for (index,uid) in enumerate(prior_tags))
        _,xyz,parameters=A.mesh.get_nodes(2,1)
        @test !isempty(parameters)
        @test A.model.get_value(2,1,parameters)≈xyz atol=3e-12
        families,_,groups=A.mesh.get_elements(2,1)
        @test families==Int32[2] && !isempty(only(groups))
        incidence=Dict{Tuple{UInt64,UInt64},Int}()
        for c in groups,i in 1:3:length(c),(a,b) in ((1,2),(2,3),(3,1))
            edge=minmax(c[i+a-1],c[i+b-1]);incidence[edge]=get(incidence,edge,0)+1
        end
        @test all(count->count in (1,2),values(incidence))
        @test Set(edge for (edge,count) in incidence if count==1)==lower_edges
    end
end

@testset "Whole-face conversion exposes rational UV and consistent boundary curves" begin
    session("brep_native_hyperbola_extrusion.brep") do occ
        occ.convert_to_nurbs([(2,1)]);occ.synchronize()
        @test A.model.get_type(2,1)=="BSpline surface"
        @test all(A.model.get_type(1,tag)=="BSpline" for (_,tag) in A.model.get_entities(1))
        @test A.model.get_value(2,1,[.21,.3])≈
            [2.0414317239075026,.20459929336141333,.3] atol=4e-14
        @test A.model.get_derivative(2,1,[.21,.3])≈
            [.393337116257247,.981150587235209,0.,0.,0.,1.] atol=5e-14
        geometry=A.CURRENT[].surface_geometry[1]
        @test size(geometry.patch.surface.controls)==(3,2)
        @test geometry.patch.surface.weights[:,1]≈[1.,cosh(.75),1.]
        for curve in keys(geometry.pcurves),t in range(
                A.CURRENT[].curve_geometry[curve].t0,A.CURRENT[].curve_geometry[curve].t1;length=7)
            uv=M._occ_pcurve_d0(geometry.pcurves[curve].fwd,t)
            @test A.model.get_value(2,1,collect(uv))≈A.model.get_value(1,curve,[t]) atol=4e-13
        end
        A.option("Mesh.MeshSizeMin",.6);A.option("Mesh.MeshSizeMax",.6)
        A.option("Mesh.Algorithm",1)
        A.mesh.generate(2)
        _,xyz,uv=A.mesh.get_nodes(2,1)
        @test !isempty(uv)
        @test A.model.get_value(2,1,uv)≈xyz atol=4e-12
    end
end

@testset "Retained owner faces keep their physical pcurves after edge conversion" begin
    session("brep_native_hyperbola_extrusion.brep") do occ
        held=A.model.get_value(2,1,[.21,.3])
        occ.convert_to_nurbs([(1,1)]);occ.synchronize()
        @test A.model.get_value(2,1,[.21,.3])==held
        @test A.model.get_type(1,1)=="BSpline"
        geometry=A.CURRENT[].surface_geometry[1]
        for t in range(LOWER,UPPER;length=9)
            uv=M._occ_pcurve_d0(geometry.pcurves[1].fwd,t)
            @test A.model.get_value(2,1,collect(uv))≈A.model.get_value(1,1,[t]) atol=4e-13
        end
        A.option("Mesh.MeshSizeMin",.6);A.option("Mesh.MeshSizeMax",.6)
        A.option("Mesh.Algorithm",1)
        A.mesh.generate(2)
        _,xyz,uv=A.mesh.get_nodes(2,1)
        @test !isempty(uv)
        @test A.model.get_value(2,1,uv)≈xyz atol=4e-12
    end
end

@testset "Conversion clears its mesh closure and preserves unrelated products" begin
    for dimension in (1,2)
        name=dimension==1 ? "brep_native_hyperbola.brep" : "brep_native_hyperbola_extrusion.brep"
        session(name) do occ
            a=occ.add_point(10.,0.,0.);b=occ.add_point(11.,0.,0.)
            foreign=occ.add_line(a,b);occ.synchronize()
            A.option("Mesh.MeshSizeMin",.6);A.option("Mesh.MeshSizeMax",.6)
            A.option("Mesh.Algorithm",1);A.mesh.generate(dimension)
            foreign_nodes=deepcopy(A.mesh.get_nodes(1,foreign,true))
            foreign_cells=deepcopy(A.mesh.get_elements(1,foreign))
            prior=deepcopy(A.mesh.get_nodes())
            @test !isempty(A.mesh.get_nodes(dimension,1)[1])
            @test_throws ArgumentError occ.convert_to_nurbs([(dimension,999)])
            @test A.mesh.get_nodes()==prior
            occ.convert_to_nurbs([(dimension,1)]);occ.synchronize()
            @test isempty(A.mesh.get_elements(dimension,1)[1])
            @test A.mesh.get_nodes(1,foreign,true)==foreign_nodes
            @test A.mesh.get_elements(1,foreign)==foreign_cells
            @test Set(A.mesh.get_nodes()[1])==Set(foreign_nodes[1])
            A.option("Mesh.MeshOnlyEmpty",1);A.mesh.generate(dimension)
            _,xyz,parameters=A.mesh.get_nodes(dimension,1)
            @test !isempty(parameters)
            @test A.model.get_value(dimension,1,parameters)≈xyz atol=4e-12
            current_nodes=A.mesh.get_nodes(1,foreign,true)
            current_cells=A.mesh.get_elements(1,foreign)
            @test current_nodes[2:3]==foreign_nodes[2:3]
            @test current_cells[1]==foreign_cells[1]
            @test line_support(current_nodes,current_cells)==line_support(foreign_nodes,foreign_cells)
        end
    end
end

@testset "Uncached native node records clear without exposing pending CAD" begin
    session() do occ
        a=occ.add_point(10.,0.,0.);b=occ.add_point(11.,0.,0.)
        foreign=occ.add_line(a,b);occ.synchronize()
        A.mesh.add_nodes(1,1,UInt64[101],point(.21),[.21])
        A.mesh.add_nodes(1,foreign,UInt64[201],[10.5,0.,0.],[.5])
        held=A.mesh.get_nodes(1,foreign)
        pending=occ.add_point(20.,0.,0.)
        @test !((0,pending) in A.model.get_entities(0))
        occ.convert_to_nurbs([(1,1)])
        @test !((0,pending) in A.model.get_entities(0))
        occ.synchronize()
        @test (0,pending) in A.model.get_entities(0)
        @test isempty(A.mesh.get_nodes(1,1)[1])
        @test A.mesh.get_nodes(1,foreign)==held
        @test Set(A.mesh.get_nodes()[1])==Set(UInt64[201])
    end
end

@testset "Partial edge replacement regenerates a face over retained lower products" begin
    session("brep_native_hyperbola_extrusion.brep") do occ
        A.option("Mesh.MeshSizeMin",.6);A.option("Mesh.MeshSizeMax",.6)
        A.option("Mesh.Algorithm",1);A.mesh.generate(2)
        retained_nodes=deepcopy(A.mesh.get_nodes(1,3,true))
        retained_cells=deepcopy(A.mesh.get_elements(1,3))
        occ.convert_to_nurbs([(1,1)]);occ.synchronize()
        @test isempty(A.mesh.get_elements(1,1)[1])
        @test isempty(A.mesh.get_elements(2,1)[1])
        @test A.mesh.get_nodes(1,3,true)==retained_nodes
        @test A.mesh.get_elements(1,3)==retained_cells
        A.option("Mesh.MeshOnlyEmpty",1);A.mesh.generate(2)
        _,xyz,parameters=A.mesh.get_nodes(2,1)
        @test !isempty(parameters)
        @test A.model.get_value(2,1,parameters)≈xyz atol=4e-12
        current_nodes=A.mesh.get_nodes(1,3,true)
        current_cells=A.mesh.get_elements(1,3)
        @test current_nodes[2:3]==retained_nodes[2:3]
        @test line_support(current_nodes,current_cells)==line_support(retained_nodes,retained_cells)
    end
end

@testset "Replacing CAD resets an already-cleared completion status" begin
    session() do occ
        A.option("Mesh.MeshSizeMin",.6);A.option("Mesh.MeshSizeMax",.6)
        A.mesh.generate(1)
        @test !isempty(A.mesh.get_nodes(1,1)[1])
        A.mesh.clear()
        @test isempty(A.mesh.get_nodes()[1])
        occ.convert_to_nurbs([(1,1)]);occ.synchronize()
        A.option("Mesh.MeshOnlyEmpty",1);A.mesh.generate(1)
        _,xyz,parameters=A.mesh.get_nodes(1,1)
        @test !isempty(parameters)
        @test A.model.get_value(1,1,parameters)≈xyz atol=3e-12
    end
end
end
