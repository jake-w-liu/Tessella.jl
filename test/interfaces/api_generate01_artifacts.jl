module APIGenerate01Artifacts

using Tessella

const api=Tessella.API

function _line(;count=3,law="Progression",coefficient=1.)
    api.model.add_point(0.,0,0;tag=1,meshSize=.2)
    api.model.add_point(1.,0,0;tag=2,meshSize=.2)
    api.model.add_line(1,2;tag=1)
    api.mesh.set_transfinite_curve(1,count,law,coefficient)
end

function _record(name,build;renumber=true,order=1,only_empty=false)
    api.initialize()
    try
        api.option("Mesh.Renumber",renumber ? 1 : 0)
        api.option("Mesh.ElementOrder",order)
        api.option("Mesh.MeshOnlyEmpty",only_empty ? 1 : 0)
        build()
        crc=Tessella.Elements.mixed_crc(api.mesh.get())
        return (name=name,n_nodes=crc.n_nodes,n_blocks=crc.n_blocks,n_cells=crc.n_cells,
                max_node=api.mesh.get_max_node_tag(),
                max_element=api.mesh.get_max_element_tag(),sha=crc.sha)
    finally
        api.finalize()
    end
end

function records()
    result=NamedTuple[]
    push!(result,_record("fresh_generation0",()->begin
        _line()
        api.mesh.generate(0)
    end))
    push!(result,_record("physical_independent_curves",()->begin
        for (tag,x,y) in ((1,0.,0.),(2,1.,0.),(3,1.,0.),(4,2.,0.),(5,4.,1.))
            api.model.add_point(x,y,0;tag)
        end
        for (curve,a,b) in ((1,1,2),(2,3,4))
            api.model.add_line(a,b;tag=curve)
            api.mesh.set_transfinite_curve(curve,3)
        end
        api.model.add_physical_group(1,[1];tag=7,name="selected curve")
        api.model.add_physical_group(0,[5];tag=8,name="orphan point")
        api.mesh.generate(1)
    end))
    push!(result,_record("raw_coincident_point_vertices",()->begin
        api.model.add_discrete_entity(0,1)
        api.mesh.add_nodes(0,1,[901,902],[0.,0,0,0.,0,0])
        api.mesh.generate(1)
    end;renumber=false))
    push!(result,_record("raw_coincident_disconnected_chains",()->begin
        api.model.add_discrete_entity(1,1)
        api.mesh.add_nodes(1,1,[500,100,300,700,200,900],
                           [0.,0,0,.5,0,0,1.,0,0,0.,0,0,.5,0,0,1.,0,0])
        api.mesh.add_elements_by_type(1,1,[400,100,300,200],
                                     [500,100,100,300,700,200,200,900])
        api.mesh.generate(1)
    end;renumber=false))
    push!(result,_record("native_sparse_onlyempty_quadratic",()->begin
        _line()
        api.mesh.add_nodes(0,1,[11],[0.,0,0])
        api.mesh.add_nodes(0,2,[22],[1.,0,0])
        api.mesh.add_nodes(1,1,[33],[.5,0,0])
        api.mesh.add_elements_by_type(1,15,[101],[11])
        api.mesh.add_elements_by_type(2,15,[102],[22])
        api.mesh.add_elements_by_type(1,1,[201,202],[11,33,33,22])
        api.mesh.generate(1)
    end;renumber=false,order=2,only_empty=true))
    for (name,law,count,coefficient) in (("progression","Progression",6,2.),
                                        ("bump","Bump",7,2.),
                                        ("beta","Beta",7,1.2))
        push!(result,_record("graded_"*name,()->begin
            _line(;count,law,coefficient)
            api.mesh.generate(1)
        end))
    end
    push!(result,_record("periodic_full_quadratic",()->begin
        for (tag,x,y) in ((1,0.,0.),(2,1.,0.),(3,0.,2.),(4,1.,2.))
            api.model.add_point(x,y,0;tag)
        end
        api.model.add_line(1,2;tag=1)
        api.model.add_line(3,4;tag=2)
        for curve in 1:2;api.mesh.set_transfinite_curve(curve,3);end
        api.mesh.set_periodic(1,[2],[1],[1.,0,0,0,0,1.,0,2.,0,0,1.,0,0,0,0,1.])
        api.mesh.generate(1)
    end;order=2))
    return result
end

end
