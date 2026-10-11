using Test,Tessella,SHA,TOML,LinearAlgebra
module NativeSizeSourceRetention
using Test,Tessella,SHA,TOML,LinearAlgebra
const A=Tessella.API
const ARTIFACTS=joinpath(@__DIR__,"..","artifacts","native_size_source_retention")
const PRIMARY=joinpath(ARTIFACTS,"higher_generated_setter_phase_primary_20261009.json")
@assert bytes2hex(sha256(read(PRIMARY)))=="26bfcad8d0f1f0a1edebd8e7c638ebf6adc51356da099158a40a7175fda32bf6"
const CASES=Dict{String,Any}[]
function snapshot()
    records=Dict{String,Any}()
    for (dim,tag) in A.model.get_entities()
        ids,xyz,uv=A.mesh.get_nodes(dim,tag)
        types,labels,cells=A.mesh.get_elements(dim,tag)
        records["$dim:$tag"]=Dict("nodes"=>collect(ids),"coords"=>collect(xyz),"params"=>collect(uv),"types"=>collect(types),
            "elements"=>[collect(v) for v in labels],"cells"=>[collect(v) for v in cells])
    end
    return Dict("records"=>records,"maxnode"=>A.mesh.get_max_node_tag(),"maxelement"=>A.mesh.get_max_element_tag())
end

function square()
    for (i,(x,y)) in enumerate(((0.,0.),(1.,0.),(1.,1.),(0.,1.)))
        A.model.add_point(x,y,0.;meshSize=.4,tag=i)
    end
    for i in 1:4
        A.model.add_line(i,mod1(i+1,4);tag=i)
        A.mesh.set_transfinite_curve(i,3)
    end
    A.model.add_curve_loop([1,2,3,4];tag=1);A.model.add_plane_surface([1];tag=1)
    A.mesh.set_transfinite_surface(1)
end

function cube()
    points=((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),(0.,0.,1.),(1.,0.,1.),(1.,1.,1.),(0.,1.,1.))
    for (i,xyz) in enumerate(points);A.model.add_point(xyz...;meshSize=.5,tag=i);end
    edges=((1,2),(2,3),(3,4),(4,1),(5,6),(6,7),(7,8),(8,5),(1,5),(2,6),(3,7),(4,8))
    for (i,(a,b)) in enumerate(edges)
        A.model.add_line(a,b;tag=i);A.mesh.set_transfinite_curve(i,3)
    end
    loops=((1,2,3,4),(5,6,7,8),(1,10,-5,-9),(2,11,-6,-10),(3,12,-7,-11),(4,9,-8,-12))
    for (i,edges) in enumerate(loops)
        A.model.add_curve_loop(collect(edges);tag=i);A.model.add_plane_surface([i];tag=i)
        A.mesh.set_transfinite_surface(i);A.mesh.set_recombine(2,i)
    end
    A.model.add_surface_loop([-1,2,3,4,5,6];tag=1);A.model.add_volume([1];tag=1)
    A.mesh.set_transfinite_volume(1)
end

function surface_certificate()
    nodes,flat,_=A.mesh.get_nodes()
    xyz=Dict(tag=>Tuple(flat[3i-2:3i]) for (i,tag) in enumerate(nodes))
    _,cells=A.mesh.get_elements_by_type(2,1)
    incidence=Dict{Tuple{UInt64,UInt64},Int}()
    area=zero(Rational{BigInt})
    for cell in eachcol(reshape(cells,3,:))
        p,q,r=[xyz[tag] for tag in cell]
        exact=map(value->Rational{BigInt}(value),[p[1],p[2],q[1],q[2],r[1],r[2]])
        determinant=(exact[3]-exact[1])*(exact[6]-exact[2])-(exact[4]-exact[2])*(exact[5]-exact[1])
        @test determinant>0
        area+=determinant/2
        for i in 1:3
            key=minmax(cell[i],cell[mod1(i+1,3)])
            incidence[key]=get(incidence,key,0)+1
        end
    end
    @test area==1
    @test all(value->value in (1,2),values(incidence))
    _,lines=A.mesh.get_elements_by_type(1)
    @test Set(key for (key,value) in incidence if value==1)==Set(minmax(cell...) for cell in eachcol(reshape(lines,2,:)))
    @test length(xyz)-length(incidence)+length(cells)÷3==1
end

function volume_certificate()
    nodes,flat,_=A.mesh.get_nodes()
    xyz=Dict(tag=>collect(flat[3i-2:3i]) for (i,tag) in enumerate(nodes))
    _,cells=A.mesh.get_elements_by_type(5,1)
    faces=((1,2,3,4),(5,6,7,8),(1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8))
    signs=((-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),(-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1))
    incidence=Dict{NTuple{4,UInt64},Int}();volume=0.
    gauss=(-inv(sqrt(3.)),inv(sqrt(3.)))
    for cell in eachcol(reshape(cells,8,:))
        for face in faces
            key=Tuple(sort!(collect(cell[collect(face)])))
            incidence[key]=get(incidence,key,0)+1
        end
        for u in gauss,v in gauss,w in gauss
            jacobian=zeros(3,3)
            for (tag,(sx,sy,sz)) in zip(cell,signs)
                gradient=[sx*(1+sy*v)*(1+sz*w),sy*(1+sx*u)*(1+sz*w),sz*(1+sx*u)*(1+sy*v)]/8
                jacobian .+=xyz[tag]*transpose(gradient)
            end
            value=det(jacobian)
            @test value>0
            volume+=value
        end
    end
    @test volume≈1 atol=3e-14 rtol=0
    _,quads=A.mesh.get_elements_by_type(3)
    @test Set(key for (key,value) in incidence if value==1)==Set(Tuple(sort!(collect(cell))) for cell in eachcol(reshape(quads,4,:)))
    @test all(value->value in (1,2),values(incidence))
end

function run_generated(dim,onlyempty)
    A.initialize()
    try
        A.option("Mesh.Renumber",0);A.option("Mesh.MeshOnlyEmpty",onlyempty)
        dim==2 ? square() : cube()
        A.mesh.generate(dim)
        before=snapshot();A.mesh.set_size([(0,1)],.125);after_setter=snapshot()
        entry=Dict{String,Any}("name"=>dim==2 ? "TF_square" : "TF_hex","onlyempty"=>onlyempty,
            "before"=>before,"after_setter"=>after_setter)
        push!(CASES,entry)
        @test after_setter==before
        A.mesh.generate(dim);after=snapshot();entry["after_generate"]=after
        lower=[key for key in keys(before["records"]) if parse(Int,first(split(key,':')))<dim]
        for key in lower,field in ("nodes","coords","params","types","cells")
            @test after["records"][key][field]==before["records"][key][field]
        end
        @test after["maxnode"]==(dim==2 ? (onlyempty==0 ? 10 : 9) : (onlyempty==0 ? 28 : 27))
        @test after["maxelement"]==(dim==2 ? (onlyempty==0 ? 44 : 36) : (onlyempty==0 ? 128 : 120))
        dim==2 ? surface_certificate() : volume_certificate()
        A.mesh.generate(dim);repeat=snapshot();entry["after_repeat"]=repeat
        @test repeat["maxnode"]==(dim==2 ? (onlyempty==0 ? 11 : 9) : (onlyempty==0 ? 29 : 27))
        @test repeat["maxelement"]==(dim==2 ? (onlyempty==0 ? 68 : 52) : (onlyempty==0 ? 192 : 176))
        @test Tessella.validate(A.LAST_MESH[]).ok
        println("GENERATED_REGEN_CASE ",entry["name"]," OE",onlyempty," max=",after["maxnode"],"/",after["maxelement"]," repeat=",repeat["maxnode"],"/",repeat["maxelement"])
    finally
        A.finalize()
    end
end

function raw_sources()
    A.initialize()
    try
        A.option("Mesh.Renumber",0)
        A.model.add_point(0.,0.,0.;tag=1);A.model.add_point(1.,0.,0.;tag=2)
        A.model.add_line(1,2;tag=1);A.mesh.set_transfinite_curve(1,3)
        A.mesh.add_nodes(0,1,[101],[.125,0.,0.]);A.mesh.add_elements_by_type(1,15,[401],[101])
        A.model.add_discrete_entity(0,99);A.mesh.add_nodes(0,99,[999],[3.,3.,0.]);A.mesh.add_elements_by_type(99,15,[899],[999])
        A.model.add_discrete_entity(1,98);A.mesh.add_nodes(1,98,[1001,1002],[2.,0.,0.,3.,0.,0.],[0.,1.]);A.mesh.add_elements_by_type(98,1,[898],[1001,1002])
        A.mesh.generate(1);before=snapshot();A.mesh.set_size([(0,1)],.2);after_setter=snapshot()
        entry=Dict{String,Any}("name"=>"raw_sources","before"=>before,"after_setter"=>after_setter,"events"=>Dict{String,Any}[]);push!(CASES,entry)
        @test before==after_setter
        for dim in (1,2,1,2)
            A.mesh.generate(dim);after=snapshot();push!(entry["events"],Dict("dimension"=>dim,"state"=>after))
            for key in ("0:1","0:99","1:98");@test after["records"][key]==before["records"][key];end
            @test A.mesh.get_node(101)[1]==[.125,0.,0.]
            @test A.mesh.get_node(101)[3:4]==(0,1)
        end
    finally
        A.finalize()
    end
end

function imported_sources(version)
    A.initialize()
    try
        A.option("Mesh.Renumber",0)
        source=joinpath(ARTIFACTS,"higher_generated_import_source_$(replace(string(version),'.'=>'_'))_20261009.msh")
        fixture=joinpath(ARTIFACTS,"higher_generated_import_setter_$(replace(string(version),'.'=>'_'))_20261009.geo")
        @test isfile(source)
        A.open_geo!(fixture)
        A.mesh.generate(1);before=snapshot()
        A.mesh.set_size([(0,1)],.125);after_setter=snapshot()
        entry=Dict{String,Any}("name"=>"MSH_import","version"=>version,"before"=>before,"after_setter"=>after_setter)
        push!(CASES,entry)
        @test before==after_setter
        A.mesh.generate(2);after=snapshot();entry["after_generate"]=after
        for key in ("0:99","1:98","2:97")
            @test after["records"][key]==before["records"][key]
        end
        tags,cells=A.mesh.get_elements_by_type(2,97)
        @test length(tags)==1
        points=[Rational{BigInt}.(A.mesh.get_node(node)[1]) for node in cells]
        p,q,r=points
        @test ((q[1]-p[1])*(r[2]-p[2])-(q[2]-p[2])*(r[1]-p[1]))/2==1//2
        @test Tessella.validate(A.LAST_MESH[]).ok
    finally
        A.finalize()
    end
end

function closed_point_exception(onlyempty)
    A.initialize()
    try
        A.option("Mesh.Renumber",0);A.option("Mesh.MeshOnlyEmpty",onlyempty)
        A.open_geo!(joinpath(ARTIFACTS,"higher_closed_point_mesh_fixture_20261009.geo"))
        A.mesh.add_nodes(0,1,[101],[.25,0.,0.]);A.mesh.add_elements_by_type(1,15,[401],[101])
        A.mesh.generate(1)
        entry=Dict{String,Any}("name"=>"closed_Spline_Point","onlyempty"=>onlyempty,"after"=>snapshot());push!(CASES,entry)
        @test A.mesh.get_nodes(0,1)[1]==UInt64[101]
        @test A.mesh.get_node(101)[1]==[1.,0.,0.]
        @test A.mesh.get_elements(0,1)[2]==[UInt64[401]]
        A.mesh.set_node(101,[.25,0.,0.],Float64[]);A.mesh.generate(1)
        entry["after_repeat"]=snapshot()
        @test A.mesh.get_node(101)[1]==(onlyempty==0 ? [1.,0.,0.] : [.25,0.,0.])
        @test A.mesh.get_elements(0,1)[2]==[UInt64[401]]
        tags,flat,_=A.mesh.get_nodes();xyz=Dict(tag=>Tuple(flat[3i-2:3i]) for (i,tag) in enumerate(tags))
        _,lines=A.mesh.get_elements_by_type(1,1)
        @test length(lines)==16
        @test all(cell->xyz[cell[1]]!=xyz[cell[2]],eachcol(reshape(lines,2,:)))
        @test allunique(tags)
    finally
        A.finalize()
    end
end


const POINT_PRIMARY=[(transfinite=false,extra=false,onlyempty=0,maxnode=111,maxelement=449,repeat_node=114,repeat_element=498,interior=[0.312499999999346, 0.3125000000006088, 0.0, 0.6796874999994118, 0.3203124999996483, 0.0, 0.48934520914417934, 0.6335876702335537, 0.0]),
(transfinite=false,extra=false,onlyempty=1,maxnode=111,maxelement=449,repeat_node=111,repeat_element=469,interior=[0.312499999999346, 0.3125000000006088, 0.0, 0.6796874999994118, 0.3203124999996483, 0.0, 0.48934520914417934, 0.6335876702335537, 0.0]),
(transfinite=false,extra=true,onlyempty=0,maxnode=112,maxelement=450,repeat_node=115,repeat_element=499,interior=[0.312499999999346, 0.3125000000006088, 0.0, 0.6796874999994118, 0.3203124999996483, 0.0, 0.48934520914417934, 0.6335876702335537, 0.0]),
(transfinite=false,extra=true,onlyempty=1,maxnode=112,maxelement=450,repeat_node=112,repeat_element=470,interior=[0.312499999999346, 0.3125000000006088, 0.0, 0.6796874999994118, 0.3203124999996483, 0.0, 0.48934520914417934, 0.6335876702335537, 0.0]),
(transfinite=true,extra=false,onlyempty=0,maxnode=109,maxelement=428,repeat_node=110,repeat_element=452,interior=[0.4669211784807215, 0.500000000000488, 0.0]),
(transfinite=true,extra=false,onlyempty=1,maxnode=109,maxelement=428,repeat_node=109,repeat_element=444,interior=[0.4669211784807215, 0.500000000000488, 0.0]),
(transfinite=true,extra=true,onlyempty=0,maxnode=110,maxelement=429,repeat_node=111,repeat_element=453,interior=[0.4669211784807215, 0.500000000000488, 0.0]),
(transfinite=true,extra=true,onlyempty=1,maxnode=110,maxelement=429,repeat_node=110,repeat_element=445,interior=[0.4669211784807215, 0.500000000000488, 0.0])]
@assert bytes2hex(sha256(read(joinpath(ARTIFACTS,"higher_actual_point_surface_setter_primary_20261009_v2.json"))))=="018ca8eb2f402d793a35dc6f98ec421c7cc7639b1f3f7721fe49320e4af5d9ad"

function point_surface_certificate()
    nodes,flat,_=A.mesh.get_nodes()
    xyz=Dict(tag=>Tuple(flat[3i-2:3i]) for (i,tag) in enumerate(nodes))
    _,cells=A.mesh.get_elements_by_type(2,1)
    incidence=Dict{Tuple{UInt64,UInt64},Int}();used=Set{UInt64}()
    area=zero(Rational{BigInt})
    for cell in eachcol(reshape(cells,3,:))
        union!(used,cell)
        p,q,r=[Rational{BigInt}.(xyz[tag]) for tag in cell]
        determinant=(q[1]-p[1])*(r[2]-p[2])-(q[2]-p[2])*(r[1]-p[1])
        @test determinant>0
        area+=determinant/2
        for i in 1:3
            key=minmax(cell[i],cell[mod1(i+1,3)])
            incidence[key]=get(incidence,key,0)+1
        end
    end
    _,lines=A.mesh.get_elements_by_type(1)
    boundary_area=zero(Rational{BigInt})
    for cell in eachcol(reshape(lines,2,:))
        p,q=[Rational{BigInt}.(xyz[tag]) for tag in cell]
        boundary_area+=(p[1]*q[2]-q[1]*p[2])/2
    end
    @test area==boundary_area
    @test area>0
    @test all(value->value in (1,2),values(incidence))
    @test Set(key for (key,value) in incidence if value==1)==Set(minmax(cell...) for cell in eachcol(reshape(lines,2,:)))
    @test length(used)-length(incidence)+length(cells)÷3==1
end

function point_surface(reference)
    A.initialize()
    try
        A.option("Mesh.Renumber",0);A.option("Mesh.MeshOnlyEmpty",reference.onlyempty)
        square()
        reference.transfinite || A.mesh.remove_constraints([(2,1)])
        tags=reference.extra ? [101,102] : [101]
        coordinates=reference.extra ? [.125,0.,0.,.2,0.,0.] : [.125,0.,0.]
        labels=reference.extra ? [401,402] : [401]
        A.mesh.add_nodes(0,1,tags,coordinates);A.mesh.add_elements_by_type(1,15,labels,tags)
        A.mesh.generate(1);before=snapshot()
        @test_throws ArgumentError Tessella.Model._model_require_point_mesh_geometry(A.CURRENT[],"legacy no-input context")
        A.mesh.set_size([(0,1)],.125);setter=snapshot()
        entry=Dict{String,Any}("name"=>"actual_Point_Surface","transfinite"=>reference.transfinite,
            "extra"=>reference.extra,"onlyempty"=>reference.onlyempty,"before"=>before,"after_setter"=>setter)
        push!(CASES,entry)
        @test setter==before
        A.mesh.generate(2);after=snapshot();entry["after_generate"]=after
        @test after["records"]["0:1"]==before["records"]["0:1"]
        for curve in 1:4,field in ("nodes","coords","params","types","cells")
            @test after["records"]["1:$curve"][field]==before["records"]["1:$curve"][field]
        end
        point_surface_certificate()
        for tag in tags
            @test A.mesh.get_node(tag)[3:4]==(0,1)
        end
        if reference.transfinite
            @test after["maxnode"]==reference.maxnode
            @test after["maxelement"]==reference.maxelement
            @test length(A.mesh.get_elements_by_type(2,1)[1])==8
            @test A.mesh.get_nodes(2,1)[2]≈reference.interior atol=3e-11 rtol=0
        end
        A.mesh.generate(2);repeat=snapshot();entry["after_repeat"]=repeat
        @test repeat["records"]["0:1"]==before["records"]["0:1"]
        point_surface_certificate()
        if reference.transfinite
            @test repeat["maxnode"]==reference.repeat_node
            @test repeat["maxelement"]==reference.repeat_element
        end
        @test Tessella.validate(A.LAST_MESH[]).ok
        println("ACTUAL_POINT_SURFACE_CASE TF=",reference.transfinite," extra=",reference.extra," OE=",reference.onlyempty,
            " max=",after["maxnode"],"/",after["maxelement"])
    finally
        A.finalize()
    end
end


const _TOPOLOGY_API=A
const TOPOLOGY_PRIMARY=Dict{String,Any}("entities"=>[Dict{String,Any}("dimension"=>0,"tag"=>2,"type"=>"Point","bounds_bits"=>["3ff0000000000000","0000000000000000","0000000000000000","3ff0000000000000","0000000000000000","0000000000000000"],"node_tags"=>[1],"xyz_bits"=>["3ff0000000000000","0000000000000000","0000000000000000"],"param_bits"=>[],"types"=>[15],"element_tags"=>[[1]],"connectivity"=>[[1]]),Dict{String,Any}("dimension"=>0,"tag"=>5,"type"=>"Point","bounds_bits"=>["0000000000000000","0000000000000000","3ff0000000000000","0000000000000000","0000000000000000","3ff0000000000000"],"node_tags"=>[2],"xyz_bits"=>["0000000000000000","0000000000000000","3ff0000000000000"],"param_bits"=>[],"types"=>[15],"element_tags"=>[[2]],"connectivity"=>[[2]]),Dict{String,Any}("dimension"=>0,"tag"=>7,"type"=>"Point","bounds_bits"=>["0000000000000000","3ff0000000000000","0000000000000000","0000000000000000","3ff0000000000000","0000000000000000"],"node_tags"=>[3],"xyz_bits"=>["0000000000000000","3ff0000000000000","0000000000000000"],"param_bits"=>[],"types"=>[15],"element_tags"=>[[3]],"connectivity"=>[[3]]),Dict{String,Any}("dimension"=>0,"tag"=>10,"type"=>"Point","bounds_bits"=>["0000000000000000","0000000000000000","0000000000000000","0000000000000000","0000000000000000","0000000000000000"],"node_tags"=>[4],"xyz_bits"=>["0000000000000000","0000000000000000","0000000000000000"],"param_bits"=>[],"types"=>[15],"element_tags"=>[[4]],"connectivity"=>[[4]]),Dict{String,Any}("dimension"=>1,"tag"=>3,"type"=>"Line","bounds_bits"=>["0000000000000000","0000000000000000","0000000000000000","3ff0000000000000","3ff0000000000000","0000000000000000"],"node_tags"=>[5],"xyz_bits"=>["3fe0000000003b41","3fdfffffffff897e","0000000000000000"],"param_bits"=>["3fdfffffffff897e"],"types"=>[1],"element_tags"=>[[5,6]],"connectivity"=>[[1,5,5,3]]),Dict{String,Any}("dimension"=>1,"tag"=>4,"type"=>"Line","bounds_bits"=>["0000000000000000","0000000000000000","0000000000000000","0000000000000000","3ff0000000000000","3ff0000000000000"],"node_tags"=>[6],"xyz_bits"=>["0000000000000000","3fe0000000003b41","3fdfffffffff897e"],"param_bits"=>["3fdfffffffff897e"],"types"=>[1],"element_tags"=>[[7,8]],"connectivity"=>[[3,6,6,2]]),Dict{String,Any}("dimension"=>1,"tag"=>6,"type"=>"Line","bounds_bits"=>["0000000000000000","0000000000000000","0000000000000000","0000000000000000","0000000000000000","3ff0000000000000"],"node_tags"=>[],"xyz_bits"=>[],"param_bits"=>[],"types"=>[1],"element_tags"=>[[9]],"connectivity"=>[[4,2]]),Dict{String,Any}("dimension"=>1,"tag"=>8,"type"=>"Line","bounds_bits"=>["0000000000000000","0000000000000000","0000000000000000","3ff0000000000000","0000000000000000","0000000000000000"],"node_tags"=>[],"xyz_bits"=>[],"param_bits"=>[],"types"=>[1],"element_tags"=>[[10]],"connectivity"=>[[4,1]]),Dict{String,Any}("dimension"=>1,"tag"=>11,"type"=>"Line","bounds_bits"=>["0000000000000000","0000000000000000","0000000000000000","0000000000000000","3ff0000000000000","0000000000000000"],"node_tags"=>[],"xyz_bits"=>[],"param_bits"=>[],"types"=>[1],"element_tags"=>[[11]],"connectivity"=>[[3,4]]),Dict{String,Any}("dimension"=>1,"tag"=>12,"type"=>"Line","bounds_bits"=>["0000000000000000","0000000000000000","0000000000000000","3ff0000000000000","0000000000000000","3ff0000000000000"],"node_tags"=>[7],"xyz_bits"=>["3fe0000000003b41","0000000000000000","3fdfffffffff897e"],"param_bits"=>["3fdfffffffff897e"],"types"=>[1],"element_tags"=>[[12,13]],"connectivity"=>[[1,7,7,2]]),Dict{String,Any}("dimension"=>2,"tag"=>21,"type"=>"Plane","bounds_bits"=>["0000000000000000","0000000000000000","0000000000000000","3ff0000000000000","3ff0000000000000","0000000000000000"],"node_tags"=>[],"xyz_bits"=>[],"param_bits"=>[],"types"=>[2],"element_tags"=>[[14,15]],"connectivity"=>[[5,4,1,5,3,4]]),Dict{String,Any}("dimension"=>2,"tag"=>22,"type"=>"Plane","bounds_bits"=>["0000000000000000","0000000000000000","0000000000000000","3ff0000000000000","0000000000000000","3ff0000000000000"],"node_tags"=>[],"xyz_bits"=>[],"param_bits"=>[],"types"=>[2],"element_tags"=>[[16,17]],"connectivity"=>[[7,4,1,7,2,4]]),Dict{String,Any}("dimension"=>2,"tag"=>23,"type"=>"Plane","bounds_bits"=>["0000000000000000","0000000000000000","0000000000000000","3ff0000000000000","3ff0000000000000","3ff0000000000000"],"node_tags"=>[],"xyz_bits"=>[],"param_bits"=>[],"types"=>[2],"element_tags"=>[[18,19,20,21]],"connectivity"=>[[5,7,1,7,6,2,6,5,3,6,7,5]]),Dict{String,Any}("dimension"=>2,"tag"=>24,"type"=>"Plane","bounds_bits"=>["0000000000000000","0000000000000000","0000000000000000","0000000000000000","3ff0000000000000","3ff0000000000000"],"node_tags"=>[],"xyz_bits"=>[],"param_bits"=>[],"types"=>[2],"element_tags"=>[[22,23]],"connectivity"=>[[2,6,4,4,6,3]]),Dict{String,Any}("dimension"=>3,"tag"=>40,"type"=>"Volume","bounds_bits"=>["0000000000000000","0000000000000000","0000000000000000","3ff0000000000000","3ff0000000000000","3ff0000000000000"],"node_tags"=>[],"xyz_bits"=>[],"param_bits"=>[],"types"=>[4],"element_tags"=>[[24,25,26,27]],"connectivity"=>[[4,6,2,7,3,6,4,5,4,5,7,1,6,4,5,7]])],"global_node_tags"=>[1,2,3,4,5,6,7],"global_xyz_bits"=>["3ff0000000000000","0000000000000000","0000000000000000","0000000000000000","0000000000000000","3ff0000000000000","0000000000000000","3ff0000000000000","0000000000000000","0000000000000000","0000000000000000","0000000000000000","3fe0000000003b41","3fdfffffffff897e","0000000000000000","0000000000000000","3fe0000000003b41","3fdfffffffff897e","3fe0000000003b41","0000000000000000","3fdfffffffff897e"],"global_param_bits"=>[],"global_element_types"=>[1,2,4,15],"global_element_tags"=>[[5,6,7,8,9,10,11,12,13],[14,15,16,17,18,19,20,21,22,23],[24,25,26,27],[1,2,3,4]],"global_connectivity"=>[[1,5,5,3,3,6,6,2,4,2,4,1,3,4,1,7,7,2],[5,4,1,5,3,4,7,4,1,7,2,4,5,7,1,7,6,2,6,5,3,6,7,5,2,6,4,4,6,3],[4,6,2,7,3,6,4,5,4,5,7,1,6,4,5,7],[1,2,3,4]],"scalar_nodes"=>[Dict{String,Any}("tag"=>1,"xyz_bits"=>["3ff0000000000000","0000000000000000","0000000000000000"],"params_bits"=>[],"owner"=>[0,2]),Dict{String,Any}("tag"=>2,"xyz_bits"=>["0000000000000000","0000000000000000","3ff0000000000000"],"params_bits"=>[],"owner"=>[0,5]),Dict{String,Any}("tag"=>3,"xyz_bits"=>["0000000000000000","3ff0000000000000","0000000000000000"],"params_bits"=>[],"owner"=>[0,7]),Dict{String,Any}("tag"=>4,"xyz_bits"=>["0000000000000000","0000000000000000","0000000000000000"],"params_bits"=>[],"owner"=>[0,10]),Dict{String,Any}("tag"=>5,"xyz_bits"=>["3fe0000000003b41","3fdfffffffff897e","0000000000000000"],"params_bits"=>["3fdfffffffff897e"],"owner"=>[1,3]),Dict{String,Any}("tag"=>6,"xyz_bits"=>["0000000000000000","3fe0000000003b41","3fdfffffffff897e"],"params_bits"=>["3fdfffffffff897e"],"owner"=>[1,4]),Dict{String,Any}("tag"=>7,"xyz_bits"=>["3fe0000000003b41","0000000000000000","3fdfffffffff897e"],"params_bits"=>["3fdfffffffff897e"],"owner"=>[1,12])],"max_node_tag"=>15,"max_element_tag"=>82,"volume_boundary"=>[[2,21],[2,-22],[2,23],[2,-24]])
@assert bytes2hex(sha256(read(joinpath(ARTIFACTS,"tessella_api_topology_tetrahedron_primary_20261009.json"))))=="9237b93365b86c030fbb7c97707bb1422fc3b2b9e473e64c97cc540190163723"
function _add_api_topology_tetrahedron()
    for (tag,x,y,z) in ((10,0.0,0.0,0.0),(2,1.0,0.0,0.0),
                        (7,0.0,1.0,0.0),(5,0.0,0.0,1.0))
        _TOPOLOGY_API.model.add_point(x,y,z;tag=tag)
    end
    for (tag,first_point,last_point) in
            ((8,10,2),(3,2,7),(11,7,10),(6,10,5),(12,2,5),(4,7,5))
        _TOPOLOGY_API.model.add_line(first_point,last_point;tag=tag)
    end
    for (tag,curves) in ((21,[8,3,11]),(22,[8,12,-6]),
                         (23,[3,4,-12]),(24,[11,6,-4]))
        _TOPOLOGY_API.model.add_curve_loop(curves;tag=tag)
        _TOPOLOGY_API.model.add_plane_surface([tag];tag=tag)
    end
    _TOPOLOGY_API.model.add_surface_loop([21,-22,23,-24];tag=30)
    _TOPOLOGY_API.model.add_volume([30];tag=40)
    return nothing
end


_topology_bits(values)=[string(reinterpret(UInt64,Float64(value));base=16,pad=16) for value in values]
function topology_snapshot()
    entities=Dict{String,Any}[]
    for (dim,tag) in A.model.get_entities()
        ids,xyz,uv=A.mesh.get_nodes(dim,tag)
        types,labels,cells=A.mesh.get_elements(dim,tag)
        push!(entities,Dict("dimension"=>dim,"tag"=>tag,"type"=>A.model.get_type(dim,tag),
            "bounds_bits"=>_topology_bits(A.model.get_bounding_box(dim,tag)),
            "node_tags"=>collect(ids),"xyz_bits"=>_topology_bits(xyz),"param_bits"=>_topology_bits(uv),
            "types"=>collect(types),"element_tags"=>[collect(v) for v in labels],"connectivity"=>[collect(v) for v in cells]))
    end
    ids,xyz,uv=A.mesh.get_nodes();types,labels,cells=A.mesh.get_elements()
    scalar=Dict{String,Any}[]
    for tag in ids
        p,u,dim,entity=A.mesh.get_node(tag)
        push!(scalar,Dict("tag"=>tag,"xyz_bits"=>_topology_bits(p),"params_bits"=>_topology_bits(u),"owner"=>[dim,entity]))
    end
    Dict("entities"=>entities,"global_node_tags"=>collect(ids),"global_xyz_bits"=>_topology_bits(xyz),
        "global_param_bits"=>_topology_bits(uv),"global_element_types"=>collect(types),
        "global_element_tags"=>[collect(v) for v in labels],"global_connectivity"=>[collect(v) for v in cells],
        "scalar_nodes"=>scalar,"max_node_tag"=>A.mesh.get_max_node_tag(),"max_element_tag"=>A.mesh.get_max_element_tag(),
        "volume_boundary"=>[collect(pair) for pair in A.model.get_boundary([(3,40)],false,true,false)])
end
function topology_control()
    A.initialize()
    try
        _add_api_topology_tetrahedron()
        generated=A.mesh.generate(3);actual=topology_snapshot()
        push!(CASES,Dict("name"=>"original_Topology_CAD_helper","actual"=>actual))
        for field in keys(TOPOLOGY_PRIMARY)
            @test actual[field]==TOPOLOGY_PRIMARY[field]
        end
        @test Tessella.validate(generated).ok
        nodes,flat,_=A.mesh.get_nodes();xyz=Dict(tag=>Rational{BigInt}.(flat[3i-2:3i]) for (i,tag) in enumerate(nodes))
        labels,cells=A.mesh.get_elements_by_type(4,40)
        volume=zero(Rational{BigInt});faces=Dict{NTuple{3,UInt64},Int}()
        patterns=((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        for cell in eachcol(reshape(cells,4,:))
            p,q,r,s=[xyz[tag] for tag in cell];u=q-p;v=r-p;w=s-p
            determinant=u[1]*(v[2]*w[3]-v[3]*w[2])-u[2]*(v[1]*w[3]-v[3]*w[1])+u[3]*(v[1]*w[2]-v[2]*w[1])
            @test determinant>0;volume+=determinant/6
            for pattern in patterns
                key=Tuple(sort!(collect(cell[collect(pattern)])));faces[key]=get(faces,key,0)+1
            end
        end
        @test volume==1//6
        _,triangles=A.mesh.get_elements_by_type(2)
        @test Set(face for (face,count) in faces if count==1)==Set(Tuple(sort!(collect(cell))) for cell in eachcol(reshape(triangles,3,:)))
        @test all(value->value in (1,2),values(faces))
        println("ORIGINAL_TOPOLOGY_CONTROL nodes=",length(nodes)," tets=",length(labels)," exact_volume=",volume,
            " max=",actual["max_node_tag"],"/",actual["max_element_tag"])
    finally
        A.finalize()
    end
end

try
    @testset "Setter retains source products and generation phase state" begin
        @testset "Original CAD Topology helper native state" begin
            topology_control()
        end
        for dim in (2,3),onlyempty in (0,1)
            @testset "Generated dimension $dim OnlyEmpty $onlyempty" begin
                run_generated(dim,onlyempty)
            end
        end
        @testset "Caller supplied native Point and discrete sources" begin
            raw_sources()
        end
        for version in (2.2,4.1)
            @testset "MSH$version source retention" begin
                imported_sources(version)
            end
        end
        for reference in POINT_PRIMARY
            @testset "Actual Point Surface TF$(reference.transfinite) extra$(reference.extra) OE$(reference.onlyempty)" begin
                point_surface(reference)
            end
        end
        for onlyempty in (0,1)
            @testset "Closed single-endpoint Point lifecycle OE$onlyempty" begin
                closed_point_exception(onlyempty)
            end
        end
    end
finally
    empty!(CASES)
    A.finalize()
end
end
