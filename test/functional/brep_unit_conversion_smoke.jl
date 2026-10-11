module BRepUnitConversionSmokeChecks
using Test, Tessella
let metre=import_step(joinpath(@__DIR__,"..","fixtures","brep_unit_metre.step");fill=false),
    inch=import_step(joinpath(@__DIR__,"..","fixtures","brep_unit_inch.step");fill=false)
    @test metre isa Tessella.Mesh && inch isa Tessella.Mesh &&
        validate(metre).ok && validate(inch).ok
    @test vec(maximum(metre.coords,dims=2))==[1000.,1000.,1000.] &&
        maximum(inch.coords)≈25.4
    @test_throws ArgumentError import_step(
        joinpath(@__DIR__,"..","fixtures","brep_unit_negative.step");fill=false)
end
end
