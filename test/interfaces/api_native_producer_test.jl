module NativeProducerTestsV1
using Test,Tessella
mkpath(joinpath(Base.pkgdir(Tessella),"test","tmp"))
include("../../validation/api_producer/suite_v1.jl")
include("../../validation/api_producer/source_births_v1.jl")
end
