module TransfiniteTriangleOrientationTests

using Test,Tessella
using Tessella.Elements: ElementBlock,MixedMesh
const M=Tessella.Model

function source(;tri=1,sign=1,count=2,pins=(),reverse=false,recombine=false,axis=:xy)
    third=axis===:xy ? "0,1,0" : "0,0,1"
    loop=sign==1 ? "1,2,3" : "-3,-2,-1"
    corners=isempty(pins) ? "" : "={$(join(pins,','))}"
    return """
    Mesh.TransfiniteTri=$tri;
    Point(1)={0,0,0,1};Point(2)={1,0,0,1};Point(3)={$third,1};
    Line(1)={1,2};Line(2)={2,3};Line(3)={3,1};
    Curve Loop(1)={$loop};Plane Surface(1)={1};
    Transfinite Curve{:}=$count;Transfinite Surface{1}$corners;
    $(reverse ? "Reverse Surface{1};" : "")
    $(recombine ? "Recombine Surface{1};" : "")
    """
end

function geometry(text)
    mktempdir() do directory
        path=joinpath(directory,"triangle.geo");write(path,text)
        return execute_geo(path;mesh_dim=0)
    end
end

blocks(mesh::Mesh)=(ElementBlock(2,mesh.tris),)
blocks(mesh::MixedMesh)=mesh.blocks

# Literal signed corner determinant against the analytically known CAD normal.
# Exact rationals cover the represented coordinates; no mesher normal or
# production orientation helper is used to determine the expected winding.
function relative_signs(mesh,normal)
    signs=Int[]
    for block in blocks(mesh),cell in eachcol(block.nodes)
        a,b,c=(Rational{BigInt}.(mesh.coords[:,node]) for node in cell[1:3])
        u=b.-a;v=c.-a
        cross=(u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
        push!(signs,sign(sum(cross[d]*normal[d] for d in 1:3)))
    end
    return signs
end

function cases()
    result=NamedTuple[]
    for tri in (0,1),sign in (-1,1),count in (2,5)
        push!(result,(;tri,sign,count))
    end
    for sign in (-1,1)
        for pins in ((1,2,3),(1,3,2))
            push!(result,(;sign,pins))
        end
        push!(result,(;sign,reverse=true))
        push!(result,(;sign,recombine=true,count=5))
        push!(result,(;sign,axis=:xz))
    end
    return result
end

@testset "Public transfinite triangle winding matches CAD before user Reverse" begin
    # Independently pinned Gmsh4.15.2: meshGFaceTransfinite chooses its lattice
    # through unsigned GEdgeLoop corners; Generator then calls orientMeshGFace
    # to align cells with the geometric normal, followed by reverseMesh.
    fixtures=cases()
    @test length(fixtures)==18
    for fixture in fixtures
        @testset "$(fixture)" begin
            execution=geometry(source(;fixture...))
            mesh=mesh_model_surface(execution.model,1)
            normal=get(fixture,:axis,:xy)===:xy ? (0,0,fixture.sign) : (0,-fixture.sign,0)
            expected=get(fixture,:reverse,false) ? -1 : 1
            @test M.model_normal(execution.model,1,[0.,0.])==collect(Float64,normal)
            @test validate(mesh).ok
            @test all(==(expected),relative_signs(mesh,normal))
            tri=get(fixture,:tri,1);count=get(fixture,:count,2)
            ncells=get(fixture,:recombine,false) ? 10 : count==2 ? 1 : tri==1 ? 16 : 28
            @test sum(block->size(block.nodes,2),blocks(mesh))==ncells
            spec,signed,sides,_,_,_,_,_,_=M._transfinite_surface_sides(execution.model,1,"test")
            tri==1 && isempty(spec.corners) &&
                ((_,sides)=M._transfinite_chained_tri(execution.model,signed,sides,1,"test"))
            kernel=if tri==1
                get(fixture,:recombine,false) ?
                    Tessella.TransfiniteTriangle.mesh_transfinite_triangle_patch(sides...;arrangement=spec.arrangement) :
                    Tessella.TransfiniteTriangle.mesh_transfinite_triangle(sides...;arrangement=spec.arrangement)
            else
                Tessella.TransfiniteTriangle.mesh_transfinite_triangle_collapsed(sides...;
                    arrangement=spec.arrangement,allow_corner_rotation=isempty(spec.corners))
            end
            # Public orientation changes cells only. Internal interpolation and
            # its node ordering remain available verbatim to volume generators.
            @test mesh.coords==kernel.coords
            @test Dict(b.msh=>sort!([sort!(collect(c)) for c in eachcol(b.nodes)]) for b in blocks(mesh))==
                  Dict(b.msh=>sort!([sort!(collect(c)) for c in eachcol(b.nodes)])
                       for b in blocks(kernel) if b.msh in (2,3))
        end
    end
end

end
