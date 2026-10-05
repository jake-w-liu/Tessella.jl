if !isdefined(@__MODULE__,:QuadTriNoNewB4RecStripCertificates)
    include("quadtri_nonew_b4_rec_strip_certificates.jl")
end
module QuadTriNoNewB4FreePlanTests
using Tessella,Test,TOML,SHA
using ..QuadTriNoNewB4RecStripCertificates
const Math=QuadTriNoNewB4RecStripCertificates
const Model=Tessella.Model
const Free=Model._ExtrudeNoNewB4Free
const WIDTH=Dict(4=>4,5=>8,6=>6,7=>5)
const FACES=Dict(4=>((1,3,2),(1,2,4),(2,3,4),(3,1,4)),
    5=>((1,4,3,2),(5,6,7,8),(1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8)),
    6=>((1,3,2),(4,5,6),(1,2,5,4),(2,3,6,5),(3,1,4,6)),
    7=>((1,4,3,2),(1,2,5),(2,3,5),(3,4,5),(4,1,5)))
const OUTER=((1,2,6,5),(2,3,7,6),(3,4,8,7),(4,1,5,8),(1,2,3,4),(5,6,7,8))
const CUBE=((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),
    (0.,0.,1.),(1.,0.,1.),(1.,1.,1.),(0.,1.,1.),(.5,.5,.5))
rotation(face)=minimum(Tuple(face[mod1(i+k,length(face))] for k in 0:length(face)-1)
    for i in 1:length(face))
unoriented(face)=min(rotation(face),rotation(reverse(face)))
function typed_cell(msh,nodes)
    faces=sort!([unoriented(Tuple(nodes[k] for k in face)) for face in FACES[msh]])
    return (msh,Tuple(sort!(collect(nodes))),Tuple(faces))
end
function template_complex(template)
    sort!([typed_cell(Int(cell.msh),Tuple(Int.(cell.nodes[1:WIDTH[Int(cell.msh)]])))
        for cell in template.cells[1:Int(template.ncells)]])
end
function expected_exterior(states)
    exterior=Set{Tuple}()
    for (face,state) in zip(OUTER,states)
        a,b,c,d=face
        if state==0
            push!(exterior,unoriented(face))
        elseif state==1
            push!(exterior,unoriented((a,b,c)));push!(exterior,unoriented((a,c,d)))
        else
            push!(exterior,unoriented((a,b,d)));push!(exterior,unoriented((b,c,d)))
        end
    end
    return exterior
end

@testset "Deferred free B4 literal factory selection" begin
    path=joinpath(@__DIR__,"..","artifacts","quadtri_nonew_b4_free_factory_oracle.toml")
    data=TOML.parsefile(path)
    @test data["supported_masks"]==253 && data["gmsh_version"]=="4.15.2"
    expected=Dict(Tuple(UInt8.(row["states"]))=>sort!([
        typed_cell(cell[1],Tuple(cell[2:end])) for cell in row["cells"]]) for row in data["factories"])
    @test length(expected)==253
    templates=Model._extrude_nonew_templates()
    for raw in Iterators.product(ntuple(_->0:2,6)...)
        states=Tuple(UInt8.(raw))
        if haskey(expected,states)
            index=Free.FinalFactories.factory_index(states)
            @test templates[Int(index)].faces==states
            @test template_complex(templates[Int(index)])==expected[states]
        else
            @test_throws ArgumentError Free.FinalFactories.factory_index(states)
        end
    end
    @test Free.FinalFactories.factory_index((0x01,0x02,0x02,0x01,0x01,0x01))==235
    @test Free.FinalFactories.factory_index((0x02,0x01,0x01,0x02,0x01,0x01))==260
    @test_throws ArgumentError Free.FinalFactories.factory_index((0x03,0x00,0x00,0x00,0x00,0x00))
end

@testset "Actual retained free B4 fan maps and complete local shell" begin
    for raw in Iterators.product(ntuple(_->0:2,6)...)
        states=Tuple(UInt8.(raw));fan=Free.CenterFan.fan(states)
        diagonal=count(!iszero,states)
        @test Int(fan.ncells)==6+diagonal
        faces=Dict{Tuple,Vector{Tuple}}();edges=Set{Tuple}();used=Set{Int}();total=zero(Math.Q)
        for cell in fan.cells[1:Int(fan.ncells)]
            msh=Int(cell.msh);nodes=Tuple(Int.(cell.nodes[1:WIDTH[msh]]))
            @test allunique(nodes) && count(==(9),nodes)==1
            points=Tuple(Math.exact(CUBE[node]) for node in nodes)
            @test Math.map_positive(points,msh)
            value=Math.cell_volume(points,msh)
            @test value>0
            total+=value;union!(used,nodes)
            for pattern in FACES[msh]
                face=Tuple(nodes[k] for k in pattern)
                push!(get!(faces,unoriented(face),Tuple[]),face)
                for k in eachindex(face)
                    push!(edges,minmax(face[k],face[mod1(k+1,length(face))]))
                end
            end
        end
        @test total==1 && used==Set(1:9)
        exterior=Set{Tuple}()
        for (face,uses) in faces
            @test length(uses) in (1,2)
            if length(uses)==1
                push!(exterior,face)
            else
                @test rotation(uses[1])==rotation(reverse(uses[2]))
            end
        end
        @test exterior==expected_exterior(states)
        @test length(used)-length(edges)+length(faces)-Int(fan.ncells)==1
        @test length(edges)==20+diagonal
        @test count(face->length(face)==4,keys(faces))==6-diagonal
        # This is one local fan's P2 support cardinality, never a global count.
        @test length(used)+length(edges)+count(face->length(face)==4,keys(faces))==35
    end
    @test_throws ArgumentError Free.CenterFan.fan((0x00,0x00,0x00,0x00,0x00,0x03))
end

@testset "Actual interval center size admission before allocation" begin
    params=(layers=[1],heights=[1.],scale_last=false,recombine=true,
        quad_to_tri=:no_new_verts,recomb_laterals=false)
    levels,refs=Model._extrude_nonew_levels(params,"free size boundary";
        source_nodes=12,extra_nodes=0,extra_nodes_per_interval=5,cells_per_interval=60)
    @test levels==[0.,1.] && refs==[(Int32(0),Int32(0))]
    for (nodes,fixed,per,cells) in ((0,0,0,1),(12,-1,0,60),(12,0,-1,60),
            (12,0,5,0),(typemax(Int),0,0,60),(12,typemax(Int),0,60),
            (12,0,typemax(Int),60),(12,0,5,typemax(Int)))
        @test_throws ArgumentError Model._extrude_nonew_levels(params,"invalid bounds";
            source_nodes=nodes,extra_nodes=fixed,extra_nodes_per_interval=per,cells_per_interval=cells)
    end
    large=merge(params,(layers=[588235],))
    @test 12*(BigInt(588235)+1)<=10_000_000
    @test 12*(BigInt(588235)+1)+5*BigInt(588235)>10_000_000
    @test_throws r"node limit" Model._extrude_nonew_levels(large,"free center reserve";
        source_nodes=12,extra_nodes=0,extra_nodes_per_interval=5,cells_per_interval=60)
    overflow=merge(params,(layers=[typemax(Int),1],heights=[.5,1.]))
    @test_throws r"layer count overflows Int" Model._extrude_nonew_levels(overflow,"overflow guard")
    @test params.layers==[1] && params.heights==[1.]
end
end
