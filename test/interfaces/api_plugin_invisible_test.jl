# `plugin("Invisible")` — Gmsh 4.15.2 `Invisible.cpp`: a non-degenerate
# XMin..ZMax box marks elements — `Inside=0` hides an element only when NO
# vertex sits inside the box, `Inside=1` hides it when ANY vertex does.
# `DeleteElements` removes the whole current invisible set (per entity) and
# `ReverseElements` flips its winding. Marks share the visibility slot that
# `mesh.get_visibility` reports; a degenerate box leaves marks untouched.
module PluginInvisibleChecks
using Test,Tessella
const A=Tessella.API

function _square()
    A.model.add("sq")
    for (x,y) in ((0.,0.),(1.,0.),(1.,1.),(0.,1.))
        A.model.occ.add_point(x,y,0.)
    end
    for i in 1:4; A.model.occ.add_line(i,mod(i,4)+1); end
    A.model.occ.add_curve_loop([1,2,3,4])
    A.model.occ.add_plane_surface([1])
    A.model.occ.synchronize()
    A.mesh.generate(2)
    return nothing
end

# Independent reimplementation of the upstream predicate over ALL entities.
function _expected_hidden(xmin,xmax,ymin,ymax,zmin,zmax,inside)
    ntags,coords,_=A.mesh.get_nodes(-1,-1,false,false)
    pos=Dict(UInt64(t)=>(coords[3(i-1)+1],coords[3(i-1)+2],coords[3(i-1)+3])
             for (i,t) in enumerate(ntags))
    hidden=Set{UInt64}()
    for ent in A.model.get_entities(-1)
        types,etags,ntags_e=A.mesh.get_elements(ent[1],ent[2])
        for (ti,etype) in enumerate(types)
            nvert=A.mesh.get_element_properties(etype)[4]
            for (k,etag) in enumerate(etags[ti])
                inbox=false
                for v in 1:nvert
                    p=pos[UInt64(ntags_e[ti][(k-1)*nvert+v])]
                    if xmin<=p[1]<=xmax && ymin<=p[2]<=ymax && zmin<=p[3]<=zmax
                        inbox=true;break
                    end
                end
                visible=inside==0 ? inbox : !inbox
                visible || push!(hidden,UInt64(etag))
            end
        end
    end
    return hidden
end

function _all_element_tags()
    tags=UInt64[]
    for ent in A.model.get_entities(-1)
        _,etags,_=A.mesh.get_elements(ent[1],ent[2])
        for v in etags; append!(tags,UInt64.(v)); end
    end
    return tags
end

@testset "Invisible marks, deletes and reverses the upstream invisible set" begin
    A.initialize()
    try
        _square()
        alltags=_all_element_tags()
        @test !isempty(alltags)

        # Inside=1 over the top-right quadrant: elements with ANY vertex
        # inside are hidden across every entity (surface tris plus boundary
        # lines and points), exactly like the upstream all-entities scan.
        A.plugin.set_number("Invisible","DeleteElements",0.0)
        A.plugin.set_number("Invisible","XMin",0.7);A.plugin.set_number("Invisible","XMax",2.0)
        A.plugin.set_number("Invisible","YMin",0.7);A.plugin.set_number("Invisible","YMax",2.0)
        A.plugin.set_number("Invisible","ZMin",-1.0);A.plugin.set_number("Invisible","ZMax",1.0)
        A.plugin.set_number("Invisible","Inside",1.0)
        A.plugin.run("Invisible")
        expected=_expected_hidden(0.7,2.0,0.7,2.0,-1.0,1.0,1)
        @test !isempty(expected)
        vis=A.mesh.get_visibility(collect(alltags))
        hidden=Set(t for (t,v) in zip(alltags,vis) if v==0)
        @test hidden==expected

        # A degenerate box leaves the mark set untouched; DeleteElements then
        # removes exactly the still-invisible elements, per entity — upstream
        # `GModel::removeInvisibleElements` semantics.
        before=length(alltags)
        A.plugin.set_number("Invisible","DeleteElements",1.0)
        for o in ("XMin","YMin","ZMin","XMax","YMax","ZMax","Inside")
            A.plugin.set_number("Invisible",o,0.0)
        end
        A.plugin.run("Invisible")
        remaining=_all_element_tags()
        @test length(remaining)==before-length(expected)
        @test isempty(intersect(Set(remaining),expected))
        @test all(==(Int32(1)),A.mesh.get_visibility(collect(remaining)))
    finally
        A.finalize()
    end

    A.initialize()
    try
        _square()
        # Inside=0: only elements with NO vertex inside the box are hidden.
        A.plugin.set_number("Invisible","DeleteElements",0.0)
        A.plugin.set_number("Invisible","XMin",0.9);A.plugin.set_number("Invisible","XMax",2.0)
        A.plugin.set_number("Invisible","YMin",0.9);A.plugin.set_number("Invisible","YMax",2.0)
        A.plugin.set_number("Invisible","ZMin",-1.0);A.plugin.set_number("Invisible","ZMax",1.0)
        A.plugin.set_number("Invisible","Inside",0.0)
        A.plugin.run("Invisible")
        expected=_expected_hidden(0.9,2.0,0.9,2.0,-1.0,1.0,0)
        alltags=_all_element_tags()
        vis=A.mesh.get_visibility(collect(alltags))
        hidden=Set(t for (t,v) in zip(alltags,vis) if v==0)
        @test hidden==expected
        @test !isempty(hidden)

        # ReverseElements flips the invisible set's winding (triangles read
        # back with swapped 2nd/3rd vertices) without removing anything.
        types,etags,ntags_before=A.mesh.get_elements(2,1)
        ti=findfirst(==(2),types)
        pre=Dict(UInt64(t)=>copy(ntags_before[ti][(k-1)*3+1:(k-1)*3+3])
                 for (k,t) in enumerate(etags[ti]))
        hidtris=sort!(collect(intersect(hidden,Set(keys(pre)))))
        @test !isempty(hidtris)
        A.plugin.set_number("Invisible","ReverseElements",1.0)
        A.plugin.run("Invisible")
        _,etags2,ntags_after=A.mesh.get_elements(2,1)
        post=Dict(UInt64(t)=>copy(ntags_after[ti][(k-1)*3+1:(k-1)*3+3])
                  for (k,t) in enumerate(etags2[ti]))
        @test Set(keys(post))==Set(keys(pre))
        for t in hidtris
            a,b=pre[t],post[t]
            @test a[1]==b[1] && a[2]==b[3] && a[3]==b[2]
        end
        for t in setdiff(keys(pre),hidtris)
            @test pre[t]==post[t]
        end
    finally
        A.finalize()
    end
end

end # module
