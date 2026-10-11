module CLICheckFunctionalChecks
using Test,Tessella
const C=Tessella.CLI
const E=Tessella.Elements
mktempdir() do directory
    mesh=E.MixedMesh([0. 1.;0. 0.;0. 0.],[E.ElementBlock(1,Int32[1;2;;],Int32[0])])
    for (version,binary) in ((2.2,false),(2.2,true),(4.1,false),(4.1,true))
        file=joinpath(directory,"line-$version-$binary.msh")
        E.write_mixed_msh(file,mesh;version,binary)
        bytes=read(file);before=Set(readdir(directory))
        @test C.main([file,"-check"])==file
        @test read(file)==bytes && Set(readdir(directory))==before
    end
    valid=joinpath(directory,"valid.msh");E.write_mixed_msh(valid,mesh;version=2.2)
    @test C.main(["-check",valid])==valid
    duplicate=joinpath(directory,"duplicate.msh")
    text=replace(read(valid,String),"\r\n"=>"\n")
    text=replace(text,"\$Elements\n1\n"=>"\$Elements\n2\n")
    rows=split(text,'\n');element=only(filter(row->startswith(row,"1 1 "),rows))
    text=replace(text,"\$EndElements"=>replace(element,r"^1 "=>"2 ")*"\n\$EndElements")
    write(duplicate,text);bytes=read(duplicate)
    @test_throws ArgumentError C.main([duplicate,"-check"])
    @test read(duplicate)==bytes
    broken=joinpath(directory,"broken.msh")
    write(broken,replace(replace(read(valid,String),"\r\n"=>"\n")," 1 2\n"=>" 1 99\n"))
    @test_throws ArgumentError C.main([broken,"-check"])
    for arguments in ([valid,"-check","-check"],[valid,"-check","-1"],
            [valid,"-check","-o",joinpath(directory,"forbidden.msh")],
            [valid,"-check","-format","msh4"],[valid,"-check","-bin"],
            [valid,"-check","-clscale","1"],["scene.geo","-check"])
        @test_throws ArgumentError C.main(arguments)
    end
    @test !ispath(joinpath(directory,"forbidden.msh"))
end
end
