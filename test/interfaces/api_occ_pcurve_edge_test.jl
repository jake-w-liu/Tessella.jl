module OCCNativePcurveEdgeChecks
using Test,Tessella
using LinearAlgebra: norm,det
const A=Tessella.API
const M=Tessella.Model
const FIXTURES=joinpath(@__DIR__,"..","fixtures")

function session(f,name;all_entities=false)
    A.initialize()
    try
        out=A.model.occ.import_shapes(joinpath(FIXTURES,
            "brep_native_pcurve_edge_$name.brep");highestDimOnly=!all_entities)
        @test Set(out)==(all_entities ? Set([(0,1),(0,2),(1,1)]) : Set([(1,1)]))
        A.model.occ.synchronize()
        return f(A.model.occ)
    finally
        A.finalize()
    end
end

point(kind,t)=kind=="helix" ? [cos(2π*t),sin(2π*t),t] : [t*t,t,0.]
tangent(kind,t)=kind=="helix" ? [-2π*sin(2π*t),2π*cos(2π*t),1.] : [2t,1.,0.]
acceleration(kind,t)=kind=="helix" ? [-4π^2*cos(2π*t),-4π^2*sin(2π*t),0.] : [2.,0.,0.]

@testset "Native pcurve-only edges own exact composed geometry" begin
    for kind in ("plane","helix")
        session(kind;all_entities=kind=="helix") do occ
            @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(2,1,0,0)
            lo,hi=A.model.get_parametrization_bounds(1,1)
            @test collect(lo)==[0.] && collect(hi)==[1.]
            for t in (0.,.13,.37,.79,1.)
                @test A.model.get_value(1,1,[t])≈point(kind,t) atol=2e-14
                @test A.model.get_derivative(1,1,[t])≈tangent(kind,t) atol=2e-13
                @test collect(M._occ_brep_edge_d2(A.CURRENT[].curve_geometry[1],t)[3])≈
                    acceleration(kind,t) atol=1e-12
                # The public OCC second-derivative query retains its
                # documented finite difference of first derivatives.
                a,b=t-.001<=0 ? (t,t+.001) : t+.001>=1 ? (t-.001,t) : (t-.001,t+.001)
                @test A.model.get_second_derivative(1,1,[t])≈
                    (tangent(kind,b).-tangent(kind,a))./(b-a) atol=2e-11
            end
            length_oracle=kind=="helix" ? sqrt(4π^2+1) : sqrt(5)/2+asinh(2)/4
            @test occ.get_mass(1,1)≈length_oracle rtol=1e-10
            bounds=A.model.get_bounding_box(1,1)
            @test bounds==occ.get_bounding_box(1,1)
            @test all(bounds[i]-1e-14<=point(kind,t)[i]<=bounds[i+3]+1e-14
                for t in range(0.,1.;length=101),i in 1:3)
            query=point(kind,.3).+(kind=="helix" ?
                [.1cos(.6π),.1sin(.6π),0.] : [0.,0.,.25])
            closest,parameter=A.model.get_closest_point(1,1,query)
            @test closest≈point(kind,.3) atol=2e-10
            @test collect(parameter)≈[.3] atol=2e-10
            A.option("Mesh.MeshSizeMin",.4)
            A.option("Mesh.MeshSizeMax",.4)
            A.mesh.generate(1)
            node_tags,xyz,parameters=A.mesh.get_nodes(1,1)
            @test !isempty(node_tags) && length(parameters)==length(node_tags)
            @test A.model.get_value(1,1,parameters)≈xyz atol=2e-12
            families,_,cells=A.mesh.get_elements(1,1)
            @test families==Int32[1] && !isempty(only(cells))
            all_tags,all_xyz,_=A.mesh.get_nodes()
            coordinates=reshape(all_xyz,3,:)
            rows=Dict(uid=>row for (row,uid) in enumerate(all_tags))
            @test maximum(norm(coordinates[:,rows[cell[i]]].-
                coordinates[:,rows[cell[i+1]]]) for cell in cells for i in 1:2:length(cell))<=.4+1e-9
        end
    end
end

@testset "Composed edges retain parameters through copy and affine transport" begin
    session("helix") do occ
        @test occ.copy([(1,1)])==[(1,2)]
        occ.translate([(1,2)],3.,-2.,1.)
        occ.dilate([(1,2)],0.,0.,0.,2.,3.,4.)
        occ.mirror([(1,2)],1.,0.,0.,0.)
        occ.synchronize()
        for t in (.11,.41,.83)
            original=point("helix",t)
            expected=[-2(original[1]+3),3(original[2]-2),4(original[3]+1)]
            @test A.model.get_value(1,1,[t])≈original atol=2e-14
            @test A.model.get_value(1,2,[t])≈expected atol=2e-13
            @test A.model.get_derivative(1,2,[t])≈[-2,3,4].*tangent("helix",t) atol=2e-12
        end
        bounds=occ.get_bounding_box(1,2)
        @test all(bounds[i]<=A.model.get_value(1,2,[t])[i]<=bounds[i+3]
            for t in range(0.,1.;length=51),i in 1:3)
        held=A.model.get_value(1,2,[.41])
        occ.remove([(1,1)],true)
        occ.synchronize()
        @test A.model.get_value(1,2,[.41])==held
    end
end

@testset "Finite sequential affine images do not require finite matrix products" begin
    surface=(kind=2,center=(0.,0.,0.),axis=(0.,0.,1.),
        X=(1.,0.,0.),Y=(0.,1.,0.),a=1e-200,b=0.)
    pc=(kind=6,degree=1,periodic=false,poles=[(0.,0.),(2π,1e-200)],
        weights=[1.,1.],knots=[0.,0.,1.,1.])
    curve=M._occ_brep_trace_curve(surface,pc,0.,1.,"finite affine images").curve
    linear=(1e200,0.,0.,0.,1e200,0.,0.,0.,1e200)
    transform=M._AffineTransform(linear,(0.,0.,0.),NTuple{16,Float64}[],
        (linear=linear,loc=(0.,0.,0.)))
    image=M._model_occ_transform_nurbs_curve(curve,transform,"finite affine images")
    image=M._model_occ_transform_nurbs_curve(image,transform,"finite affine images")
    p,d,dd=M._model_occ_curve_d2(image,.25)
    @test all(isfinite,(p...,d...,dd...))
    @test collect(p)./1e200≈[0.,1.,.25] atol=1e-13
    @test collect(d)./1e200≈[-2π,0.,1.] atol=1e-12
    @test collect(dd)./1e200≈[0.,-4π^2,0.] atol=1e-11
end

@testset "A solid's pcurve-only edge retains its actual lower mesh in 3D" begin
    A.initialize()
    try
        occ=A.model.occ
        @test occ.import_shapes(joinpath(FIXTURES,"brep_native_pcurve_edge_box.brep"))==[(3,1)]
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(8,12,6,1)
        @test A.model.get_type(1,1)=="Unknown"
        @test occ.get_mass(3,1)≈1. atol=1e-12
        A.option("Mesh.MeshSizeMin",2.)
        A.option("Mesh.MeshSizeMax",2.)
        A.option("Mesh.MeshOnlyEmpty",1)
        A.option("Mesh.Algorithm",1)
        A.mesh.generate(1)
        prior_tags,prior_xyz,_=A.mesh.get_nodes()
        prior_coordinates=reshape(prior_xyz,3,:)
        _,_,prior_groups=A.mesh.get_elements(1)
        prior_edges=Set(minmax(cells[i],cells[i+1]) for cells in prior_groups for i in 1:2:length(cells))
        A.mesh.generate(3)
        tags,xyz,_=A.mesh.get_nodes()
        coordinates=reshape(xyz,3,:)
        rows=Dict(uid=>row for (row,uid) in enumerate(tags))
        @test all(haskey(rows,uid) && coordinates[:,rows[uid]]==prior_coordinates[:,index]
            for (index,uid) in enumerate(prior_tags))
        _,_,current_groups=A.mesh.get_elements(1)
        @test Set(minmax(cells[i],cells[i+1]) for cells in current_groups
            for i in 1:2:length(cells))==prior_edges
        families,_,groups=A.mesh.get_elements(3)
        @test families==Int32[4]
        incidence=Dict{NTuple{3,UInt64},Int}()
        volumes=Float64[]
        for cells in groups,i in 1:4:length(cells)
            tet=cells[i:i+3];p=coordinates[:,rows[tet[1]]]
            push!(volumes,det(hcat((coordinates[:,rows[tet[j]]].-p for j in 2:4)...))/6)
            for face in ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
                key=Tuple(sort(tet[collect(face)]))
                incidence[key]=get(incidence,key,0)+1
            end
        end
        @test !isempty(volumes) && all(>(0),volumes)
        @test sum(volumes)≈1. rtol=1e-12
        _,_,groups=A.mesh.get_elements(2)
        facets=Set(Tuple(sort(cells[i:i+2])) for cells in groups for i in 1:3:length(cells))
        @test facets==Set(face for (face,count) in incidence if count==1)
    finally
        A.finalize()
    end
end

@testset "A cylinder's composed rim remains shared by its actual Faces" begin
    A.initialize()
    try
        occ=A.model.occ
        @test occ.import_shapes(joinpath(FIXTURES,"brep_native_pcurve_edge_cylinder.brep"))==[(3,1)]
        occ.synchronize()
        @test Tuple(length(A.model.get_entities(d)) for d in 0:3)==(2,3,3,1)
        @test A.model.get_type(1,1)=="Unknown"
        @test occ.get_mass(3,1)≈π/2 rtol=1e-12
        for t in (.2,1.7,5.4)
            @test A.model.get_value(1,1,[t])≈[3+.5cos(t),.5sin(t),2.] atol=2e-14
        end
        A.option("Mesh.MeshSizeMin",.8)
        A.option("Mesh.MeshSizeMax",.8)
        A.option("Mesh.MeshOnlyEmpty",1)
        A.option("Mesh.Algorithm",1)
        A.mesh.generate(1)
        prior_tags,prior_xyz,_=A.mesh.get_nodes()
        prior_coordinates=reshape(prior_xyz,3,:)
        _,_,curve_groups=A.mesh.get_elements(1,1)
        @test length(unique(only(curve_groups)))>=3
        A.mesh.generate(2)
        tags,xyz,_=A.mesh.get_nodes()
        coordinates=reshape(xyz,3,:)
        rows=Dict(uid=>row for (row,uid) in enumerate(tags))
        @test all(haskey(rows,uid) && coordinates[:,rows[uid]]==prior_coordinates[:,index]
            for (index,uid) in enumerate(prior_tags))
        data=A.CURRENT[].meshing.mesh_authority.mesh.entity_data
        @test all(owner===nothing || owner[1]!=1 ||
            parameter!==nothing && length(parameter)==1
            for (owner,parameter) in zip(data.node_entities,data.node_parametric))
        families,_,groups=A.mesh.get_elements(2)
        @test families==Int32[2]
        cells=reduce(vcat,groups;init=UInt64[])
        surface=Tessella.Mesh(coordinates;
            tris=reshape(Int32[rows[uid] for uid in cells],3,:))
        @test Tessella.validate(surface).ok
        @test M._plc_self_intersection3(surface)===nothing
        incidence=Dict{Tuple{UInt64,UInt64},Int}()
        for i in 1:3:length(cells),pair in ((1,2),(2,3),(3,1))
            edge=minmax(cells[i+pair[1]-1],cells[i+pair[2]-1])
            incidence[edge]=get(incidence,edge,0)+1
        end
        @test all(==(2),values(incidence))
        @test length(tags)-length(incidence)+length(cells)÷3==2
        for (_,face) in A.model.get_entities(2)
            _,world,uv=A.mesh.get_nodes(2,face)
            @test A.model.get_value(2,face,uv)≈world atol=2e-12
        end
    finally
        A.finalize()
    end
end

@testset "Writer-rounded cylindrical seams retain their ordinary 3D edges" begin
    A.initialize()
    try
        mktempdir() do directory
            source=read(joinpath(FIXTURES,"brep_native_analytic.brep"),String)
            path=joinpath(directory,"cylinder.brep")
            write(path,replace(source,r"\+1 0\s*$"=>"+31 0\n"))
            @test A.model.occ.import_shapes(path)==[(3,1)]
            A.model.occ.synchronize()
            A.option("Mesh.MeshSizeMin",.8)
            A.option("Mesh.MeshSizeMax",.8)
            A.option("Mesh.MeshOnlyEmpty",1)
            A.option("Mesh.Algorithm",1)
            A.mesh.generate(1)
            A.mesh.generate(2)
            tags,xyz,_=A.mesh.get_nodes()
            rows=Dict(uid=>row for (row,uid) in enumerate(tags))
            _,_,groups=A.mesh.get_elements(2)
            cells=reduce(vcat,groups;init=UInt64[])
            surface=Tessella.Mesh(reshape(xyz,3,:);
                tris=reshape(Int32[rows[uid] for uid in cells],3,:))
            @test Tessella.validate(surface).ok
            @test M._plc_self_intersection3(surface)===nothing
            for (_,face) in A.model.get_entities(2)
                _,world,uv=A.mesh.get_nodes(2,face)
                @test A.model.get_value(2,face,uv)≈world atol=2e-12
            end
        end
    finally
        A.finalize()
    end
end

@testset "Plane witnesses use the actual trimmed copied curve" begin
    A.initialize()
    try
        @test A.model.occ.import_shapes(joinpath(FIXTURES,
            "brep_native_trimmed_planar_curve.brep"))==[(2,1)]
        A.model.occ.synchronize()
        @test A.model.get_value(1,1,[.5])≈[.5,-.1,0.] atol=1e-14
        @test A.model.occ.get_mass(2,1)≈16/15 rtol=1e-12
        @test M._model_surface_plane(A.CURRENT[],1,"trimmed planar boundary")!==nothing
        A.option("Mesh.MeshSizeMin",.5)
        A.option("Mesh.MeshSizeMax",.5)
        A.option("Mesh.Algorithm",1)
        A.mesh.generate(2)
        _,xyz,_=A.mesh.get_nodes()
        @test all(==(0.),@view(reshape(xyz,3,:)[3,:]))
        _,world,uv=A.mesh.get_nodes(2,1)
        @test A.model.get_value(2,1,uv)≈world atol=2e-12
    finally
        A.finalize()
    end
end
end
