if !isdefined(@__MODULE__,:QuadTriNoNewFourQuadStripCertificates)
    include("../geometry/quadtri_nonew_four_quad_strip_certificates.jl")
end

module NoNewFourQuadStripBoundaryTests
using Test,Tessella
using ..QuadTriNoNewFourQuadStripCertificates
const API=Tessella.API
const Cert=QuadTriNoNewFourQuadStripCertificates

function run_tests()
@testset "Four-Quad strip public construction and P2 ownership" begin
    f=Cert.fixture("golden";laterals=false,layers=:one,pins=(1,2,3,4))
    for coherence in (0,1)
        API.initialize()
        try
            mktempdir() do directory
                path=joinpath(directory,"strip.geo")
                write(path,replace(f.source,"Geometry.AutoCoherence=0"=>"Geometry.AutoCoherence=$coherence"))
                opened=API.open_geo!(path;mesh_dim=0);volume=Int(opened.lists["sweep"][2])
                if coherence==0
                    @test_throws ArgumentError API.mesh.generate(3)
                    @test isempty(API.mesh.get_elements()[1])
                else
                    API.mesh.generate(3)
                    types,labels,connectivity=API.mesh.get_elements(3,volume)
                    @test all(msh->msh in (4,5,6,7),types)
                    @test length(unique(vcat(connectivity...)))==20
                    @test Tessella.Elements.validate(Cert.public_volume(API,volume)).ok
                    volumes=API.mesh.get_element_qualities(vcat(labels...),"volume")
                    @test all(>(0),volumes)
                    @test sum(volumes)≈1 atol=2e-11 rtol=0
                    API.mesh.set_order(2)
                    types,labels,connectivity=API.mesh.get_elements(3,volume)
                    @test all(msh->API.mesh.get_element_properties(msh)[3]==2,types)
                    used=unique(vcat(connectivity...))
                    points=[API.mesh.get_node(node)[1] for node in used]
                    @test all(point->all(isfinite,point),points) && allunique(Tuple.(points))
                    for node in used
                        point,parameters,dim,entity=API.mesh.get_node(node)
                        @test dim in (1,2) ? length(parameters)==dim : isempty(parameters)
                        if dim in (1,2)
                            @test maximum(abs.(API.model.get_value(dim,entity,parameters).-point))<=2e-11
                        elseif dim==3
                            @test entity==volume
                        end
                    end
                    for msh in types,edge in eachcol(reshape(API.mesh.get_element_edge_nodes(msh,volume,false),3,:))
                        a,b,middle=(API.mesh.get_node(node)[1] for node in edge)
                        @test middle≈(a+b)/2 atol=2e-11 rtol=0
                    end
                    volumes=API.mesh.get_element_qualities(vcat(labels...),"volume")
                    @test all(>(0),volumes)
                    @test sum(volumes)≈1 atol=2e-11 rtol=0
                end
            end
        finally
            API.finalize()
        end
    end
end
end

get(ENV,"TESSELLA_FOUR_STRIP_DEFINE_ONLY","")=="1" || run_tests()
end
