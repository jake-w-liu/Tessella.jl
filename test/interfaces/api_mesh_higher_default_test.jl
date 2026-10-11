using Tessella,Test,SHA,TOML
using Tessella.MeshTypes: Mesh,mesh_crc
const _DSQ=Tessella.API

@testset "Higher transaction preserves callback identity without copying userdata" begin
    m=Tessella.Model.GeoModel()
    userdata=zeros(UInt8,1_000_000)
    callback=(dim,tag,x,y,z,lc)->(userdata[1]+=UInt8(1);lc)
    Tessella.Model.set_size_callback!(m,callback)
    staged=Tessella.Model._model_mesh_transaction_copy(m)
    @test staged!==m && staged.meshing!==m.meshing && staged.points!==m.points
    @test staged.meshing.size_callback===callback
    staged.meshing.size_callback(0,1,0.,0.,0.,1.)
    @test userdata[1]==1
    # Warm both operations; the positive control proves captured storage
    # would otherwise be charged to each detached model transaction.
    deepcopy(m);Tessella.Model._model_mesh_transaction_copy(m)
    snapshot_bytes=minimum(@allocated(deepcopy(m)) for _ in 1:3)
    transaction_bytes=minimum(@allocated(Tessella.Model._model_mesh_transaction_copy(m)) for _ in 1:3)
    @test snapshot_bytes>=length(userdata)
    @test transaction_bytes+65536<snapshot_bytes
end

# Literal independent Gmsh4.15.2 owner payloads; no implementation-derived labels.
const _DSQ_PRIMARY_ROWS = [Dict("only_empty"=>0,"order"=>1,"renumber"=>0,"counts"=>Dict("1"=>4,"2"=>2,"15"=>4),"entities"=>[Dict("owner"=>[0,1],"nodes"=>[1],"xyz"=>[0.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[1]],"connections"=>[[1]]),Dict("owner"=>[0,2],"nodes"=>[2],"xyz"=>[1.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[2]],"connections"=>[[2]]),Dict("owner"=>[0,3],"nodes"=>[3],"xyz"=>[1.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[3]],"connections"=>[[3]]),Dict("owner"=>[0,4],"nodes"=>[4],"xyz"=>[0.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[4]],"connections"=>[[4]]),Dict("owner"=>[1,1],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[5]],"connections"=>[[1,2]]),Dict("owner"=>[1,2],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[6]],"connections"=>[[2,3]]),Dict("owner"=>[1,3],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[7]],"connections"=>[[3,4]]),Dict("owner"=>[1,4],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[8]],"connections"=>[[4,1]]),Dict("owner"=>[2,1],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[2],"elements"=>[[9,10]],"connections"=>[[1,2,4,4,2,3]])],"all_node_tags"=>[1,2,3,4],"node_highwater"=>4,"element_highwater"=>10),Dict("only_empty"=>0,"order"=>1,"renumber"=>1,"counts"=>Dict("1"=>4,"2"=>2,"15"=>4),"entities"=>[Dict("owner"=>[0,1],"nodes"=>[1],"xyz"=>[0.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[1]],"connections"=>[[1]]),Dict("owner"=>[0,2],"nodes"=>[2],"xyz"=>[1.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[2]],"connections"=>[[2]]),Dict("owner"=>[0,3],"nodes"=>[3],"xyz"=>[1.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[3]],"connections"=>[[3]]),Dict("owner"=>[0,4],"nodes"=>[4],"xyz"=>[0.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[4]],"connections"=>[[4]]),Dict("owner"=>[1,1],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[5]],"connections"=>[[1,2]]),Dict("owner"=>[1,2],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[6]],"connections"=>[[2,3]]),Dict("owner"=>[1,3],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[7]],"connections"=>[[3,4]]),Dict("owner"=>[1,4],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[8]],"connections"=>[[4,1]]),Dict("owner"=>[2,1],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[2],"elements"=>[[9,10]],"connections"=>[[1,2,4,4,2,3]])],"all_node_tags"=>[1,2,3,4],"node_highwater"=>4,"element_highwater"=>10),Dict("only_empty"=>0,"order"=>2,"renumber"=>0,"counts"=>Dict("8"=>4,"9"=>2,"15"=>4),"entities"=>[Dict("owner"=>[0,1],"nodes"=>[1],"xyz"=>[0.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[1]],"connections"=>[[1]]),Dict("owner"=>[0,2],"nodes"=>[2],"xyz"=>[1.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[2]],"connections"=>[[2]]),Dict("owner"=>[0,3],"nodes"=>[3],"xyz"=>[1.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[3]],"connections"=>[[3]]),Dict("owner"=>[0,4],"nodes"=>[4],"xyz"=>[0.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[4]],"connections"=>[[4]]),Dict("owner"=>[1,1],"nodes"=>[5],"xyz"=>[0.49999999999867184,0.0,0.0],"parameters"=>[0.49999999999867184],"types"=>[8],"elements"=>[[11]],"connections"=>[[1,2,5]]),Dict("owner"=>[1,2],"nodes"=>[6],"xyz"=>[1.0,0.49999999999867184,0.0],"parameters"=>[0.49999999999867184],"types"=>[8],"elements"=>[[12]],"connections"=>[[2,3,6]]),Dict("owner"=>[1,3],"nodes"=>[7],"xyz"=>[0.5000000000013304,1.0,0.0],"parameters"=>[0.4999999999986696],"types"=>[8],"elements"=>[[13]],"connections"=>[[3,4,7]]),Dict("owner"=>[1,4],"nodes"=>[8],"xyz"=>[0.0,0.5000000000013304,0.0],"parameters"=>[0.4999999999986696],"types"=>[8],"elements"=>[[14]],"connections"=>[[4,1,8]]),Dict("owner"=>[2,1],"nodes"=>[9],"xyz"=>[0.5,0.5,0.0],"parameters"=>[0.5,0.5],"types"=>[9],"elements"=>[[15,16]],"connections"=>[[1,2,4,5,9,8,4,2,3,9,6,7]])],"all_node_tags"=>[1,2,3,4,5,6,7,8,9],"node_highwater"=>9,"element_highwater"=>16),Dict("only_empty"=>0,"order"=>2,"renumber"=>1,"counts"=>Dict("8"=>4,"9"=>2,"15"=>4),"entities"=>[Dict("owner"=>[0,1],"nodes"=>[1],"xyz"=>[0.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[1]],"connections"=>[[1]]),Dict("owner"=>[0,2],"nodes"=>[2],"xyz"=>[1.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[2]],"connections"=>[[2]]),Dict("owner"=>[0,3],"nodes"=>[3],"xyz"=>[1.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[3]],"connections"=>[[3]]),Dict("owner"=>[0,4],"nodes"=>[4],"xyz"=>[0.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[4]],"connections"=>[[4]]),Dict("owner"=>[1,1],"nodes"=>[5],"xyz"=>[0.49999999999867184,0.0,0.0],"parameters"=>[0.49999999999867184],"types"=>[8],"elements"=>[[5]],"connections"=>[[1,2,5]]),Dict("owner"=>[1,2],"nodes"=>[6],"xyz"=>[1.0,0.49999999999867184,0.0],"parameters"=>[0.49999999999867184],"types"=>[8],"elements"=>[[6]],"connections"=>[[2,3,6]]),Dict("owner"=>[1,3],"nodes"=>[7],"xyz"=>[0.5000000000013304,1.0,0.0],"parameters"=>[0.4999999999986696],"types"=>[8],"elements"=>[[7]],"connections"=>[[3,4,7]]),Dict("owner"=>[1,4],"nodes"=>[8],"xyz"=>[0.0,0.5000000000013304,0.0],"parameters"=>[0.4999999999986696],"types"=>[8],"elements"=>[[8]],"connections"=>[[4,1,8]]),Dict("owner"=>[2,1],"nodes"=>[9],"xyz"=>[0.5,0.5,0.0],"parameters"=>[0.5,0.5],"types"=>[9],"elements"=>[[9,10]],"connections"=>[[1,2,4,5,9,8,4,2,3,9,6,7]])],"all_node_tags"=>[1,2,3,4,5,6,7,8,9],"node_highwater"=>9,"element_highwater"=>16),Dict("only_empty"=>1,"order"=>1,"renumber"=>0,"counts"=>Dict("1"=>4,"2"=>2,"15"=>4),"entities"=>[Dict("owner"=>[0,1],"nodes"=>[1],"xyz"=>[0.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[1]],"connections"=>[[1]]),Dict("owner"=>[0,2],"nodes"=>[2],"xyz"=>[1.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[2]],"connections"=>[[2]]),Dict("owner"=>[0,3],"nodes"=>[3],"xyz"=>[1.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[3]],"connections"=>[[3]]),Dict("owner"=>[0,4],"nodes"=>[4],"xyz"=>[0.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[4]],"connections"=>[[4]]),Dict("owner"=>[1,1],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[5]],"connections"=>[[1,2]]),Dict("owner"=>[1,2],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[6]],"connections"=>[[2,3]]),Dict("owner"=>[1,3],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[7]],"connections"=>[[3,4]]),Dict("owner"=>[1,4],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[8]],"connections"=>[[4,1]]),Dict("owner"=>[2,1],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[2],"elements"=>[[9,10]],"connections"=>[[1,2,4,4,2,3]])],"all_node_tags"=>[1,2,3,4],"node_highwater"=>4,"element_highwater"=>10),Dict("only_empty"=>1,"order"=>1,"renumber"=>1,"counts"=>Dict("1"=>4,"2"=>2,"15"=>4),"entities"=>[Dict("owner"=>[0,1],"nodes"=>[1],"xyz"=>[0.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[1]],"connections"=>[[1]]),Dict("owner"=>[0,2],"nodes"=>[2],"xyz"=>[1.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[2]],"connections"=>[[2]]),Dict("owner"=>[0,3],"nodes"=>[3],"xyz"=>[1.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[3]],"connections"=>[[3]]),Dict("owner"=>[0,4],"nodes"=>[4],"xyz"=>[0.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[4]],"connections"=>[[4]]),Dict("owner"=>[1,1],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[5]],"connections"=>[[1,2]]),Dict("owner"=>[1,2],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[6]],"connections"=>[[2,3]]),Dict("owner"=>[1,3],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[7]],"connections"=>[[3,4]]),Dict("owner"=>[1,4],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[1],"elements"=>[[8]],"connections"=>[[4,1]]),Dict("owner"=>[2,1],"nodes"=>[],"xyz"=>[],"parameters"=>[],"types"=>[2],"elements"=>[[9,10]],"connections"=>[[1,2,4,4,2,3]])],"all_node_tags"=>[1,2,3,4],"node_highwater"=>4,"element_highwater"=>10),Dict("only_empty"=>1,"order"=>2,"renumber"=>0,"counts"=>Dict("8"=>4,"9"=>2,"15"=>4),"entities"=>[Dict("owner"=>[0,1],"nodes"=>[1],"xyz"=>[0.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[1]],"connections"=>[[1]]),Dict("owner"=>[0,2],"nodes"=>[2],"xyz"=>[1.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[2]],"connections"=>[[2]]),Dict("owner"=>[0,3],"nodes"=>[3],"xyz"=>[1.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[3]],"connections"=>[[3]]),Dict("owner"=>[0,4],"nodes"=>[4],"xyz"=>[0.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[4]],"connections"=>[[4]]),Dict("owner"=>[1,1],"nodes"=>[5],"xyz"=>[0.49999999999867184,0.0,0.0],"parameters"=>[0.49999999999867184],"types"=>[8],"elements"=>[[11]],"connections"=>[[1,2,5]]),Dict("owner"=>[1,2],"nodes"=>[6],"xyz"=>[1.0,0.49999999999867184,0.0],"parameters"=>[0.49999999999867184],"types"=>[8],"elements"=>[[12]],"connections"=>[[2,3,6]]),Dict("owner"=>[1,3],"nodes"=>[7],"xyz"=>[0.5000000000013304,1.0,0.0],"parameters"=>[0.4999999999986696],"types"=>[8],"elements"=>[[13]],"connections"=>[[3,4,7]]),Dict("owner"=>[1,4],"nodes"=>[8],"xyz"=>[0.0,0.5000000000013304,0.0],"parameters"=>[0.4999999999986696],"types"=>[8],"elements"=>[[14]],"connections"=>[[4,1,8]]),Dict("owner"=>[2,1],"nodes"=>[9],"xyz"=>[0.5,0.5,0.0],"parameters"=>[0.5,0.5],"types"=>[9],"elements"=>[[15,16]],"connections"=>[[1,2,4,5,9,8,4,2,3,9,6,7]])],"all_node_tags"=>[1,2,3,4,5,6,7,8,9],"node_highwater"=>9,"element_highwater"=>16),Dict("only_empty"=>1,"order"=>2,"renumber"=>1,"counts"=>Dict("8"=>4,"9"=>2,"15"=>4),"entities"=>[Dict("owner"=>[0,1],"nodes"=>[1],"xyz"=>[0.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[1]],"connections"=>[[1]]),Dict("owner"=>[0,2],"nodes"=>[2],"xyz"=>[1.0,0.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[2]],"connections"=>[[2]]),Dict("owner"=>[0,3],"nodes"=>[3],"xyz"=>[1.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[3]],"connections"=>[[3]]),Dict("owner"=>[0,4],"nodes"=>[4],"xyz"=>[0.0,1.0,0.0],"parameters"=>[],"types"=>[15],"elements"=>[[4]],"connections"=>[[4]]),Dict("owner"=>[1,1],"nodes"=>[5],"xyz"=>[0.49999999999867184,0.0,0.0],"parameters"=>[0.49999999999867184],"types"=>[8],"elements"=>[[5]],"connections"=>[[1,2,5]]),Dict("owner"=>[1,2],"nodes"=>[6],"xyz"=>[1.0,0.49999999999867184,0.0],"parameters"=>[0.49999999999867184],"types"=>[8],"elements"=>[[6]],"connections"=>[[2,3,6]]),Dict("owner"=>[1,3],"nodes"=>[7],"xyz"=>[0.5000000000013304,1.0,0.0],"parameters"=>[0.4999999999986696],"types"=>[8],"elements"=>[[7]],"connections"=>[[3,4,7]]),Dict("owner"=>[1,4],"nodes"=>[8],"xyz"=>[0.0,0.5000000000013304,0.0],"parameters"=>[0.4999999999986696],"types"=>[8],"elements"=>[[8]],"connections"=>[[4,1,8]]),Dict("owner"=>[2,1],"nodes"=>[9],"xyz"=>[0.5,0.5,0.0],"parameters"=>[0.5,0.5],"types"=>[9],"elements"=>[[9,10]],"connections"=>[[1,2,4,5,9,8,4,2,3,9,6,7]])],"all_node_tags"=>[1,2,3,4,5,6,7,8,9],"node_highwater"=>9,"element_highwater"=>16)]

function _dsq_model(row)
    _DSQ.initialize()
    _DSQ.option("Mesh.Renumber",row["renumber"])
    _DSQ.option("Mesh.MeshOnlyEmpty",row["only_empty"])
    _DSQ.option("Mesh.ElementOrder",row["order"])
    for (tag,(x,y)) in enumerate(((0.,0.),(1.,0.),(1.,1.),(0.,1.)))
        _DSQ.model.add_point(x,y,0.;tag,meshSize=1.)
    end
    for (tag,(first,last)) in enumerate(((1,2),(2,3),(3,4),(4,1)))
        _DSQ.model.add_line(first,last;tag)
        _DSQ.mesh.set_transfinite_curve(tag,2)
    end
    _DSQ.model.add_curve_loop([1,2,3,4];tag=1)
    _DSQ.model.add_plane_surface([1];tag=1);_DSQ.mesh.set_transfinite_surface(1)
end

function _dsq_check_stage(row)
    for entity in row["entities"]
        dim,tag=entity["dim"],entity["tag"]
        nodes,xyz,parameters=_DSQ.mesh.get_nodes(dim,tag,false,true)
        @test nodes==entity["nodes"]
        @test xyz≈Float64.(entity["xyz"]) atol=2e-11 rtol=0
        @test parameters≈Float64.(entity["parameters"]) atol=2e-11 rtol=0
        types,tags,connections=_DSQ.mesh.get_elements(dim,tag)
        @test types==entity["types"]
        @test tags==entity["elements"]
        @test connections==entity["connections"]
    end
    types,tags,connections=_DSQ.mesh.get_elements()
    nodes=_DSQ.mesh.get_nodes()[1]
    @test nodes==row["all_node_tags"]
    @test Dict(string(kind)=>length(labels) for (kind,labels) in zip(types,tags))==row["counts"]
    @test allunique(nodes) && allunique(vcat(tags...))
    @test Set(vcat(connections...)) ⊆ Set(nodes)
    @test _DSQ.mesh.get_max_node_tag()==row["node_highwater"]
    @test _DSQ.mesh.get_max_element_tag()==row["element_highwater"]
    authority=_DSQ.CURRENT[].meshing.mesh_authority
    @test authority.mesh===_DSQ.LAST_MESH[]
    @test _DSQ.LAST_MESH_CLASS[].node_parametric==authority.mesh.entity_data.node_parametric
    geometry=_DSQ.mesh.get()
    @test geometry isa Mesh
    @test authority.mesh.coords[:,authority.carrier_rows]==geometry.coords
end

@testset "Ordinary complete products preserve native carrier and exact lifecycle identities" begin
    path=joinpath(@__DIR__,"../artifacts/onlyempty_higher_lifecycle_primary.toml")
    # The saved payload is fixed; checkout newline conversion is portable.
    @test bytes2hex(sha256(replace(read(path,String),"\r\n"=>"\n")))==
        "c29bebd84e64081326e150429d2a5dcb421a1da2b1559802a7a6f43ea1d991d0"
    primary=TOML.parsefile(path)
    @test primary["gmsh_version"]=="4.15.2"
    try
        _dsq_model(Dict("renumber"=>0,"only_empty"=>0,"order"=>1))
        original=_DSQ.mesh.generate(2)
        _DSQ.mesh.renumber_nodes([1,2,3,4],[41,11,31,21])
        _DSQ.mesh.renumber_elements(collect(1:10),collect(101:110))
        _dsq_check_stage(primary["stages"][1])
        @test mesh_crc(_DSQ.mesh.get())==mesh_crc(original)
        _DSQ.mesh.set_order(2)
        _dsq_check_stage(primary["stages"][2])
        @test mesh_crc(_DSQ.mesh.get())==mesh_crc(original)
        before=_DSQ.LAST_MESH[]
        _DSQ.mesh.set_recombine(2,1,1.)
        _DSQ.mesh.recombine()
        _dsq_check_stage(primary["stages"][3])
        @test _DSQ.LAST_MESH[]===before
        @test mesh_crc(_DSQ.mesh.get())==mesh_crc(original)
        _DSQ.mesh.set_order(1)
        _dsq_check_stage(primary["stages"][4])
        @test mesh_crc(_DSQ.mesh.get())==mesh_crc(original)
        refined=_DSQ.mesh.refine()
        _dsq_check_stage(primary["stages"][5])
        @test refined isa Mesh
        @test mesh_crc(refined)==mesh_crc(Tessella.Refine.refine_uniform(original))
        @test mesh_crc(_DSQ.mesh.get())==mesh_crc(refined)
        surface_tags=_DSQ.mesh.get_elements(2,1)[2][1]
        @test sum(_DSQ.mesh.get_element_qualities(surface_tags,"volume"))≈1 atol=2e-11 rtol=0
    finally
        _DSQ.finalize()
    end
end

@testset "Ordinary actual API recombination preserves lower products and candidate allocation" begin
    for angle in (1.,2.)
        try
            _dsq_model(Dict("renumber"=>0,"only_empty"=>0,"order"=>1))
            _DSQ.option("Mesh.RecombinationAlgorithm",0)
            original=_DSQ.mesh.generate(2)
            lower=[(_DSQ.mesh.get_nodes(dim,tag,false,true),_DSQ.mesh.get_elements(dim,tag))
                   for (dim,tag) in ((0,1),(0,2),(0,3),(0,4),(1,1),(1,2),(1,3),(1,4))]
            _DSQ.mesh.set_recombine(2,1,angle)
            _DSQ.mesh.recombine()
            @test [(_DSQ.mesh.get_nodes(dim,tag,false,true),_DSQ.mesh.get_elements(dim,tag))
                   for (dim,tag) in ((0,1),(0,2),(0,3),(0,4),(1,1),(1,2),(1,3),(1,4))]==lower
            types,tags,connections=_DSQ.mesh.get_elements(2,1)
            @test types==Int32[angle==1 ? 2 : 3]
            @test tags==[angle==1 ? UInt64[9,10] : UInt64[12]]
            @test connections==[angle==1 ? UInt64[1,2,4,4,2,3] : UInt64[2,3,4,1]]
            @test _DSQ.mesh.get_max_element_tag()==(angle==1 ? 11 : 12)
            @test _DSQ.mesh.get_max_node_tag()==4
            if angle==1
                @test _DSQ.mesh.get() isa Mesh
                @test mesh_crc(_DSQ.mesh.get())==mesh_crc(original)
            end
        finally
            _DSQ.finalize()
        end
    end
end

@testset "Ordinary default higher generation publishes complete actual phase products" begin
    for row in _DSQ_PRIMARY_ROWS
        try
            _dsq_model(row);view=_DSQ.mesh.generate(2)
            authority=_DSQ.CURRENT[].meshing.mesh_authority
            @test authority.mesh===_DSQ.LAST_MESH[]
            if row["only_empty"]==0
                @test view isa Mesh
                @test mesh_crc(_DSQ.mesh.get())==mesh_crc(view)
                @test authority.mesh.coords[:,authority.carrier_rows]==view.coords
            end
            for entity in row["entities"]
                dim,tag=entity["owner"]
                nodes,xyz,parameters=_DSQ.mesh.get_nodes(dim,tag,false,true)
                @test nodes==entity["nodes"]
                # Native Gmsh line interpolation differs from the exact
                # affine midpoint by about 1e-12; use the established
                # carrier/UV tolerance while keeping every identity exact.
                @test xyz≈Float64.(entity["xyz"]) atol=2e-11 rtol=0
                @test parameters≈Float64.(entity["parameters"]) atol=2e-11 rtol=0
                types,tags,connections=_DSQ.mesh.get_elements(dim,tag)
                @test types==entity["types"]
                @test tags==entity["elements"]
                @test connections==entity["connections"]
            end
            types,tags,connections=_DSQ.mesh.get_elements()
            @test Dict(string(t)=>length(v) for (t,v) in zip(types,tags))==row["counts"]
            nodes=_DSQ.mesh.get_nodes()[1]
            @test nodes==row["all_node_tags"]
            @test allunique(nodes) && allunique(vcat(tags...))
            @test Set(vcat(connections...)) ⊆ Set(nodes)
            @test _DSQ.mesh.get_max_node_tag()==row["node_highwater"]
            @test _DSQ.mesh.get_max_element_tag()==row["element_highwater"]
        finally
            _DSQ.finalize()
        end
    end
end

@testset "Ordinary complete products retain the fixed native geometry CRC" begin
    try
        _DSQ.initialize()
        _DSQ.model.add_box(0,0,0,1,1,1;tag=1)
        _DSQ.model.add_point(0.25,0.25,0.25;tag=101)
        _DSQ.model.add_point(0.75,0.75,0.75;tag=102)
        view=_DSQ.mesh.generate(3)
        @test view isa Mesh
        @test mesh_crc(view).sha=="ec97caf92e736f8590fe61176430799f103bbe3804756c0577ab56c1d6f4a52e"
        @test mesh_crc(_DSQ.mesh.get())==mesh_crc(view)
        @test _DSQ.mesh.get_elements(0,101)[1]==Int32[15]
        @test _DSQ.mesh.get_elements(0,102)[1]==Int32[15]
        _DSQ.model.add_point(0.5,0.5,0.5;tag=103)
        repeated=_DSQ.mesh.generate(3)
        @test repeated isa Mesh
        @test mesh_crc(repeated)==mesh_crc(view)
        # Independent Gmsh4.15.2 added-Point rows: old status3 skips Mesh0D.
        @test isempty(_DSQ.mesh.get_nodes(0,103)[1])
        @test isempty(_DSQ.mesh.get_elements(0,103)[1])
    finally
        _DSQ.finalize()
    end
end
