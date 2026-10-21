# `view.probe` on mixed-element caches — the located quadrangle/hexahedron/
# prism/pyramid interpolates with its own first-order corner basis and the
# view rows are keyed by public (file) tags, not dense cache indices.
module ViewProbeMixedChecks
using Test, Tessella
const A=Tessella.API

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

field(x,y,z)=x+2.0y+3.0z

# Unit-square quadrangle on z = 0 with non-identity public tags.
const QUAD_MSH="""
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Entities
0 0 1 0
7 0 0 0 1 1 0 0 0
\$EndEntities
\$Nodes
1 4 101 104
2 7 0 4
101
102
103
104
0 0 0
1 0 0
1 1 0
0 1 0
\$EndNodes
\$Elements
1 1 51 51
2 7 3 1
51 101 102 103 104
\$EndElements
"""

# Unit-cube hexahedron, a prism offset to x in [3, 4] and a pyramid offset to
# x in [5, 7] on one volume entity — the cells are spatially disjoint so each
# probe point has exactly one owning element.
const VOLUME_MSH="""
\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Entities
0 0 0 1
9 0 -1 0 7 1 1 0 0
\$EndEntities
\$Nodes
1 19 201 225
3 9 0 19
201
202
203
204
205
206
207
208
211
212
213
214
215
216
221
222
223
224
225
0 0 0
1 0 0
1 1 0
0 1 0
0 0 1
1 0 1
1 1 1
0 1 1
3 0 0
4 0 0
3 1 0
3 0 1
4 0 1
3 1 1
5 -1 0
7 -1 0
7 1 0
5 1 0
6 0 1
\$EndNodes
\$Elements
3 3 52 54
3 9 5 1
52 201 202 203 204 205 206 207 208
3 9 6 1
53 211 212 213 214 215 216
3 9 7 1
54 221 222 223 224 225
\$EndElements
"""

const QUAD_NODES=[101,102,103,104]
const QUAD_CORNERS=((0.0,0.0,0.0),(1.0,0.0,0.0),(1.0,1.0,0.0),(0.0,1.0,0.0))
const VOLUME_CELLS=(
    (52,[201,202,203,204,205,206,207,208],
     ((0.,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)),
     ((0.3,0.4,0.6),(0.8,0.2,0.5))),
    (53,[211,212,213,214,215,216],
     ((3.,0,0),(4,0,0),(3,1,0),(3,0,1),(4,0,1),(3,1,1)),
     ((3.2,0.3,0.7),(3.25,0.25,0.5))),
    (54,[221,222,223,224,225],
     ((5.,-1,0),(7,-1,0),(7,1,0),(5,1,0),(6,0,1)),
     ((6.1,-0.2,0.3),(6.0,0.0,1.0))))

function _merge_text(text,name)
    A.model.add(name)
    mktempdir() do directory
        path=joinpath(directory,"input.msh")
        Base.write(path,text)
        A.merge(path)
    end
    return A.model.get_current()
end

@testset "view.probe interpolates quadrangle view data" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        mname=_merge_text(QUAD_MSH,"quad probe")
        values=[field(p...) for p in QUAD_CORNERS]
        points=((0.3,0.4),(0.5,0.5),(0.25,0.75))
        v=A.view.add("nodal")
        A.view.add_model_data(v,0,mname,"NodeData",
                              QUAD_NODES,[[value] for value in values])
        for (x,y) in points
            @test A.view.probe(mname,v,x,y,0.0)[1]≈field(x,y,0.0) atol=1e-9
        end
        v2=A.view.add("element nodal")
        A.view.add_model_data(v2,0,mname,"ElementNodeData",
                              [51],[values])
        for (x,y) in points
            @test A.view.probe(mname,v2,x,y,0.0)[1]≈field(x,y,0.0) atol=1e-9
        end
        v3=A.view.add("per element")
        A.view.add_model_data(v3,0,mname,"ElementData",[51],[[42.5]])
        @test A.view.probe(mname,v3,0.3,0.4,0.0)==[42.5]
        @test occursin("no element contains",
            _errtext(()->A.view.probe(mname,v,5.0,5.0,0.0)))
    finally
        A.finalize()
    end
end

@testset "view.probe interpolates hexahedron, prism and pyramid data" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        mname=_merge_text(VOLUME_MSH,"volume probe")
        for (etag,tags,corners,points) in VOLUME_CELLS
            values=[field(p...) for p in corners]
            v=A.view.add("nodal $etag")
            A.view.add_model_data(v,0,mname,"NodeData",
                                  tags,[[value] for value in values])
            v2=A.view.add("element nodal $etag")
            A.view.add_model_data(v2,0,mname,"ElementNodeData",
                                  [etag],[values])
            for p in points
                expected=field(p...)
                @test A.view.probe(mname,v,p...)[1]≈expected atol=1e-9
                @test A.view.probe(mname,v2,p...)[1]≈expected atol=1e-9
            end
        end
        v3=A.view.add("per element")
        A.view.add_model_data(v3,0,mname,"ElementData",
                              [52,53,54],[[52.0],[53.0],[54.0]])
        @test A.view.probe(mname,v3,0.3,0.4,0.6)==[52.0]
        @test A.view.probe(mname,v3,3.2,0.3,0.7)==[53.0]
        @test A.view.probe(mname,v3,6.1,-0.2,0.3)==[54.0]
    finally
        A.finalize()
    end
end

@testset "view.probe keeps interpolating simplex mesh data" begin
    A.initialize()
    try
        A.option("General.Terminal",0)
        A.option("Mesh.MeshSizeMax",0.5)
        A.model.occ.add_rectangle(0,0,0,1,1)
        A.model.occ.synchronize()
        A.mesh.generate(2)
        mname=A.model.get_current()
        ntags,ncoords,_=A.mesh.get_nodes(-1,-1,true,false)
        values=[ncoords[3i-2]+2.0*ncoords[3i-1] for i in eachindex(ntags)]
        v=A.view.add("nodal")
        A.view.add_model_data(v,0,mname,"NodeData",ntags,values)
        for (x,y) in ((0.3,0.4),(0.5,0.5),(0.1,0.8))
            @test A.view.probe(mname,v,x,y,0.0)[1]≈x+2.0*y atol=1e-10
        end
    finally
        A.finalize()
    end
end
end
