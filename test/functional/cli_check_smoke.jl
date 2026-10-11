mktempdir() do directory
    file=joinpath(directory,"check.msh")
    text="\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n2\n1 0 0 0\n2 1 0 0\n\$EndNodes\n\$Elements\n1\n1 1 2 0 1 1 2\n\$EndElements\n"
    write(file,text)
    @test Tessella.CLI.main([file,"-check"])==file && read(file,String)==text
    write(file,replace(text," 1 2\n"=>" 1 99\n"))
    @test_throws ArgumentError Tessella.CLI.main([file,"-check"])
end
