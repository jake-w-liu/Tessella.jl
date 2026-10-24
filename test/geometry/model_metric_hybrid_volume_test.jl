module ModelMetricHybridVolumeChecks
# MET02: the native hybrid 3-D inserter's tensor arm — metric-composed
# queue radius, metric circumcentre births, metric in-sphere cavity growth
# and normalized-metric close-point admission, mirroring the 2-D
# `_GmshSurfaceMetricCarrier` contract. Exercises the private helpers
# directly, the scoped `refine` entry on a small cube fill, and the public
# API hybrid Volume production under a `MathEvalAniso` background.
using Test,Tessella,LinearAlgebra,Statistics
const M=Tessella.Model
const I=M.ModelNativeHybridInsertion
const MD=Tessella.Mesh3D
const SF=Tessella.SizeField
const N=Tessella.MeshTypes
const A=Tessella.API

# Kuhn 6-tet decomposition of the unit cube — `positive_constructor` fixes
# local orientations, so the raw column order is immaterial.
function cube6(;side=1.0)
    coords=Float64[0 1 1 0 0 1 1 0;
                   0 0 1 1 0 0 1 1;
                   0 0 0 0 1 1 1 1]*side
    tets=Int32[1 1 1 1 1 1;2 2 3 4 5 5;3 6 4 8 6 8;7 7 7 7 7 7]
    return I.positive_constructor(M.Mesh(coords;tets=tets))
end

function run_refine(mesh;metric=nothing,bg=(x,y,z)->1.0,extend=false,minimize=false)
    nn=M.nnodes(mesh);nt=M.ntets(mesh)
    return I.refine(mesh,UInt64.(1:nn),UInt64.(1:nt),NamedTuple[],bg,"test";
        extend,minimize,metric)
end

function signed_volumes(mesh)
    return [LinearAlgebra.det(hcat(mesh.coords[:,c[2]]-mesh.coords[:,c[1]],
        mesh.coords[:,c[3]]-mesh.coords[:,c[1]],
        mesh.coords[:,c[4]]-mesh.coords[:,c[1]]))/6 for c in eachcol(mesh.tets)]
end

# Mean |dx|,|dy|,|dz| over unique tet edges — under a strong m11 the
# x-extent is metric-expensive, so cells elongate along the weak
# directions and the |dx| mean drops.
function edge_means(mesh)
    seen=Set{Tuple{Int32,Int32}}();dx=Float64[];dy=Float64[];dz=Float64[]
    for cell in eachcol(mesh.tets),(a,b) in ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
        key=minmax(cell[a],cell[b])
        key in seen && continue
        push!(seen,key)
        pa,pb=mesh.coords[:,key[1]],mesh.coords[:,key[2]]
        push!(dx,abs(pa[1]-pb[1]));push!(dy,abs(pa[2]-pb[2]));push!(dz,abs(pa[3]-pb[3]))
    end
    return (mean(dx),mean(dy),mean(dz))
end

@testset "native metric floor composes scalar bound per eigendirection" begin
    m=M.Metric3(64.0,1.0,4.0,0.0,0.0,0.0)
    # h=0.4 -> floor 6.25: both eigenvalues below the bound tighten to it,
    # the strongest direction keeps its anisotropy.
    floored=I.native_metric_floor(m,0.4)
    @test sort(collect(SF.metric_eigenvalues(floored)))≈[6.25,6.25,64.0] atol=1e-12
    # A bound looser than every directional size leaves the metric alone.
    @test I.native_metric_floor(m,2.0)===m
    @test I.native_metric_floor(m,0.0)===m
    @test I.native_metric_floor(m,Inf)===m
    # Isotropic input stays isotropic — the scalar-equivalence property.
    iso=I.native_metric_floor(SF.isotropic_metric(3.0),0.5)
    λ=SF.metric_eigenvalues(iso)
    @test all(v->v≈4.0,λ)
    # Anisotropy survives the floor in the strong directions.
    strong=I.native_metric_floor(m,0.45)
    λs=sort(collect(SF.metric_eigenvalues(strong)))
    @test λs[3]≈64.0 atol=1e-10
    # Rotated metric: floor acts in eigen-coordinates, keeping eigenvectors.
    v1=(1.0,0.0,0.0);v2=(0.0,sqrt(0.5),sqrt(0.5));v3=(0.0,-sqrt(0.5),sqrt(0.5))
    rotated=SF.metric_from_axes(100.0,9.0,1.0,v1,v2,v3)
    fr=I.native_metric_floor(rotated,0.6)
    λr=sort(collect(SF.metric_eigenvalues(fr)))
    @test λr≈[2.7777777777777777,9.0,100.0] atol=1e-9
end

@testset "metric circumcentre is M-equidistant and matches Euclidean at identity" begin
    mesh=cube6()
    verts=Tuple(mesh.tets[:,1])
    points=ntuple(k->Tuple(mesh.coords[:,verts[k]]),4)
    T,_,_=MD._gmsh3_triangulation(mesh)
    for field in (SF.MathEvalAnisoField(m11="64",m22="1",m33="1"),
                  SF.MathEvalAnisoField(m11="4",m22="9",m33="25"))
        m=SF.metric_at(field,0.4,0.4,0.4)
        solved=I.native_metric_center(points...,m)
        @test solved!==nothing
        center,r2=solved
        @test r2>0
        for p in points
            d=I._native_metric_quad(m,center[1]-p[1],center[2]-p[2],center[3]-p[3])
            @test d≈r2 rtol=1e-9
        end
    end
    # Identity metric reproduces the Euclidean circumcentre and R^2.
    mid=SF.isotropic_metric(1.0)
    center,r2=I.native_metric_center(points...,mid)
    ecenter=MD._gmsh3_native_center(points...)
    @test collect(center)≈collect(ecenter) atol=1e-9
    er=sum((points[1][k]-ecenter[k])^2 for k in 1:3)
    @test r2≈er rtol=1e-9
    # Co-spherical classification: own metric circumcentre lies strictly
    # inside (d=0 < r2), a far point lies outside.
    @test I.native_insphere(T,verts,center,SF.MathEvalAnisoField(m11="64",m22="1",m33="1"),nothing)
    far=(10.0,10.0,10.0)
    @test !I.native_insphere(T,verts,far,SF.MathEvalAnisoField(m11="64",m22="1",m33="1"),nothing)
end

@testset "normalized radius is scalar-compatible for isotropic input" begin
    mesh=cube6()
    T,_,_=MD._gmsh3_triangulation(mesh)
    verts=Tuple(mesh.tets[:,1])
    sizes=fill(0.7,N.nnodes(mesh));background=fill(0.7,N.nnodes(mesh))
    iso=SF.MathEvalAnisoField(m11="1/0.49",m22="1/0.49",m33="1/0.49")
    metric_r=I.native_metric_radius(T,verts,sizes,background,false,iso,nothing)
    scalar_r=MD._gmsh3_native_radius(T,verts,sizes,background,false)
    @test metric_r≈scalar_r rtol=1e-9
    # The extend_boundary floor tightens the weak directions: propagated
    # sizes below the field's weak-direction h raise the normalized radius
    # — the metric-space `min(lc_bgm, vsize)` composition.
    aniso=SF.MathEvalAnisoField(m11="64",m22="1",m33="1")
    free=I.native_metric_radius(T,verts,sizes,background,false,aniso,nothing)
    tight=fill(0.3,N.nnodes(mesh))
    floored=I.native_metric_radius(T,verts,tight,background,true,aniso,nothing)
    @test floored>free
    # The floored radius equals the metric circumradius under the
    # intersected metric evaluated at the cell centroid.
    points=ntuple(k->Tuple(mesh.coords[:,verts[k]]),4)
    cx,cy,cz=Tuple(sum(p[i] for p in points)/4 for i in 1:3)
    raw=SF.metric_at(aniso,cx,cy,cz)
    mf=I.native_metric_floor(raw,0.3)
    @test SF.metric_eigenvalues(mf)[1]>SF.metric_eigenvalues(raw)[1]
    @test floored≈sqrt(I.native_metric_center(points...,mf)[2]) rtol=1e-9
end

@testset "tensor-arm refine elongates along the weak direction" begin
    mesh=cube6()
    scalar=run_refine(mesh)
    @test M.validate(scalar.mesh).ok
    @test M.nnodes(scalar.mesh)==8 && M.ntets(scalar.mesh)==6
    xfield=SF.MathEvalAnisoField(m11="64",m22="1",m33="1")
    yfield=SF.MathEvalAnisoField(m11="1",m22="64",m33="1")
    xmesh=run_refine(mesh;metric=xfield,bg=(x,y,z)->0.125)
    ymesh=run_refine(mesh;metric=yfield,bg=(x,y,z)->0.125)
    @test M.validate(xmesh.mesh).ok && M.validate(ymesh.mesh).ok
    xv=signed_volumes(xmesh.mesh);yv=signed_volumes(ymesh.mesh)
    @test all(>(0),xv) && all(>(0),yv)
    @test sum(xv)≈1.0 atol=1e-12
    @test sum(yv)≈1.0 atol=1e-12
    # The field reaches the initial queue: the metric arm refines where the
    # scalar propagated-size radius stays below the bound.
    @test M.ntets(xmesh.mesh)>M.ntets(scalar.mesh)
    # Both fields scalarize to h=0.125; only the transported tensor
    # differs, so connectivity or coordinates must differ.
    @test xmesh.mesh.tets!=ymesh.mesh.tets ||
          !isequal(xmesh.mesh.coords,ymesh.mesh.coords)
    ex=edge_means(xmesh.mesh);ey=edge_means(ymesh.mesh)
    @test ex[1]<ex[2] && ey[2]<ey[1]
    @test ex[2]≈ex[3] rtol=1e-9
    @test ey[1]≈ey[3] rtol=1e-9
    # Operation identity and UID clocks survive the tensor arm.
    @test allunique(xmesh.node_numbers) && allunique(xmesh.element_numbers)
    @test length(xmesh.node_numbers)==M.nnodes(xmesh.mesh)
    @test length(xmesh.element_numbers)==M.ntets(xmesh.mesh)
    # Isotropic tensor input reproduces the scalar arm's product.
    iso=run_refine(mesh;metric=SF.MathEvalAnisoField(m11="1",m22="1",m33="1"))
    @test M.validate(iso.mesh).ok
    @test M.nnodes(iso.mesh)==M.nnodes(scalar.mesh) &&
          M.ntets(iso.mesh)==M.ntets(scalar.mesh)
end

@testset "metric arm keeps boundary and extend certificates" begin
    mesh=cube6()
    xfield=SF.MathEvalAnisoField(m11="64",m22="1",m33="1")
    for (extend,minimize) in ((true,false),(true,true))
        r=run_refine(mesh;metric=xfield,bg=(x,y,z)->0.125,extend,minimize)
        @test M.validate(r.mesh).ok
        @test all(>(0),signed_volumes(r.mesh))
        @test Set(first(N.boundary_faces(r.mesh.tets)))==
              Set(first(N.boundary_faces(mesh.tets)))
    end
end

# ── Public API hybrid Volume production under a tensor background ──────
function api_cube()
    points=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
            (0,0,1),(1,0,1),(1,1,1),(0,1,1))
    for (index,p) in enumerate(points)
        A.model.add_point(p[1],p[2],p[3];tag=index,meshSize=1.0)
    end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
           (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (index,ends) in enumerate(edges)
        A.model.add_line(ends[1],ends[2];tag=index)
        A.mesh.set_transfinite_curve(index,2)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),
           (2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (index,lines) in enumerate(loops)
        A.model.add_curve_loop([sign(line)*abs(line) for line in lines];tag=index)
        A.model.add_plane_surface([index];tag=index)
        A.mesh.set_transfinite_surface(index)
        A.mesh.set_recombine(2,index)
    end
    A.model.add_surface_loop(collect(1:6);tag=1)
    A.model.add_volume([1];tag=1)
    return nothing
end

function api_tensor(entries)
    tag=A.mesh.field.add("MathEvalAniso")
    for (name,expression) in entries
        A.mesh.field.set_string(tag,name,expression)
    end
    A.mesh.field.set_as_background_mesh(tag)
    return tag
end

function api_volume_cells(region)
    families,labels,nodes=A.mesh.get_elements(3,region)
    result=Dict{Int,Matrix{UInt64}}()
    for (kind,ids,vertices) in zip(families,labels,nodes)
        width=kind==4 ? 4 : kind==7 ? 5 : error("unexpected volume family $kind")
        result[Int(kind)]=reshape(vertices,width,length(ids))
    end
    return result
end

function api_nodes()
    tags,xyz,_=A.mesh.get_nodes()
    return tags,xyz
end

function tet_edge_dirs(tets,xyz,tags)
    seen=Set{Tuple{UInt64,UInt64}}();dx=Float64[];dy=Float64[];dz=Float64[]
    row=Dict(tag=>i for (i,tag) in enumerate(tags))
    for cell in eachcol(tets),(a,b) in ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
        key=minmax(cell[a],cell[b])
        key in seen && continue
        push!(seen,key)
        pa=xyz[3row[key[1]]-2:3row[key[1]]];pb=xyz[3row[key[2]]-2:3row[key[2]]]
        push!(dx,abs(pa[1]-pb[1]));push!(dy,abs(pa[2]-pb[2]));push!(dz,abs(pa[3]-pb[3]))
    end
    return (mean(dx),mean(dy),mean(dz))
end

function cell_volumes(cells,xyz,tags)
    row=Dict(tag=>i for (i,tag) in enumerate(tags))
    out=Float64[]
    for cell in eachcol(cells[4])
        a,b,c,d=(xyz[3row[cell[k]]-2:3row[cell[k]]] for k in 1:4)
        push!(out,det(hcat(b-a,c-a,d-a))/6)
    end
    return out
end

function session(body)
    A.initialize()
    try;body();finally;A.finalize();end
end

@testset "API hybrid cube consumes a MathEvalAniso background" begin
    counts=Dict{Symbol,Int}()
    dirs=Dict{Symbol,NTuple{3,Float64}}()
    # Upstream defaults all six MathEvalAniso entries to `F2 + Sin(z)`;
    # off-diagonal zeros must be spelled out for the intended diagonal
    # tensor — the partial `:partial` spec below keeps the poisoned
    # `F2 + Sin(z)` defaults (an unresolved field reference evaluates to
    # MAX_LC upstream, so the sampled matrix is not SPD) and 4.15.2's
    # volume mesher ignores it exactly like it ignores 3-D anisotropy.
    zero=[("M12","0"),("M13","0"),("M23","0")]
    for (name,entries) in ((:x,vcat([("M11","64"),("M22","1"),("M33","1")],zero)),
                           (:y,vcat([("M11","1"),("M22","64"),("M33","1")],zero)),
                           (:partial,[("M11","64"),("M22","1"),("M33","1")]),
                           (:dense,[("F","0.125")]))
        session() do
            A.option("General.Terminal",0)
            api_cube()
            if name===:dense
                tag=A.mesh.field.add("MathEval")
                A.mesh.field.set_string(tag,"F","0.125")
                A.mesh.field.set_as_background_mesh(tag)
            else
                api_tensor(entries)
            end
            A.mesh.generate(3)
            cells=api_volume_cells(1)
            @test Set(keys(cells))⊆Set([4,7])
            @test haskey(cells,4) && size(cells[4],2)>0
            tags,xyz=api_nodes()
            vols=cell_volumes(cells,xyz,tags)
            @test all(>(0),vols)
            counts[name]=length(tags)
            dirs[name]=tet_edge_dirs(cells[4],xyz,tags)
        end
    end
    # The tensor arm engages: M11=64/M22=64 contract the strong direction.
    @test counts[:x]>14 && counts[:y]>14
    @test counts[:x]!=counts[:y] || dirs[:x]!=dirs[:y]
    @test dirs[:x][1]<dirs[:x][2] && dirs[:x][1]<dirs[:x][3]
    @test dirs[:y][2]<dirs[:y][1] && dirs[:y][2]<dirs[:y][3]
    # Gmsh-faithful boundary transport: vertex `lc` is edge-derived —
    # a transfinite boundary prescribes its own spacing and neither a
    # scalar field nor a non-SPD metric sample refines the volume fill.
    @test counts[:partial]==14
    @test counts[:dense]==14
    println("HYBRID_METRIC nodes x=",counts[:x]," y=",counts[:y],
            " partial=",counts[:partial]," dense=",counts[:dense],
            " dirx=",dirs[:x]," diry=",dirs[:y]);flush(stdout)
end

@testset "API free cube scalar background parity" begin
    # A scalar field refines the volume only through the boundary mesh
    # upstream (vSizesBGM inherits boundary edge lengths for existing
    # vertices): on a free boundary 4.15.2 gives 716 nodes.
    session() do
        A.option("General.Terminal",0)
        pts=((0,0,0),(1,0,0),(1,1,0),(0,1,0),
             (0,0,1),(1,0,1),(1,1,1),(0,1,1))
        for (i,p) in enumerate(pts)
            A.model.add_point(p[1],p[2],p[3];tag=i,meshSize=1.0)
        end
        edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),
               (7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
        for (i,e) in enumerate(edges)
            A.model.add_line(e[1],e[2];tag=i)
        end
        loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),(2,11,-6,-10),
               (3,12,-7,-11),(4,9,-8,-12))
        for (i,ls) in enumerate(loops)
            A.model.add_curve_loop([sign(l)*abs(l) for l in ls];tag=i)
            A.model.add_plane_surface([i];tag=i)
        end
        A.model.add_surface_loop(collect(1:6);tag=1)
        A.model.add_volume([1];tag=1)
        tag=A.mesh.field.add("MathEval")
        A.mesh.field.set_string(tag,"F","0.125")
        A.mesh.field.set_as_background_mesh(tag)
        A.mesh.generate(3)
        cells=api_volume_cells(1)
        @test Set(keys(cells))⊆Set([4,7])
        @test haskey(cells,4) && size(cells[4],2)>0
        tags,xyz=api_nodes()
        @test all(>(0),cell_volumes(cells,xyz,tags))
        @test length(tags)>=600
        println("HYBRID_METRIC free-scalar nodes=",length(tags));flush(stdout)
    end
end

end
