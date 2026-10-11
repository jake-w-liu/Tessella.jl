help=redirect_stdout(devnull) do
    Tessella.CLI.main(["-help"])
end
version=redirect_stdout(devnull) do
    Tessella.CLI.main(["--version"])
end
@test occursin("Usage: tessella",help) && occursin("-clscale",help)
@test version==string(Base.pkgversion(Tessella))
