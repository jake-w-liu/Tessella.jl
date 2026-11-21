module ViewPosLineEndsTests
using Test, Tessella
const A = Tessella.API

# A list-based view saved as parsed `.pos` text. Gmsh 4.15.2 opens the file in
# text mode, so on Windows it writes `View "a" {\r\nSP(0,0,0){1};\r\n};\r\n`
# and `\n` line ends elsewhere. `API.view.write` used to convert the writer's
# already platform-terminated lines a second time and wrote `\r\r\n`.
@testset "A list view written as .pos has one line terminator per line" begin
    A.finalize(); A.initialize()
    try
        A.option("General.Terminal", 0)
        view = A.view.add("a")
        A.view.add_list_data(view, "SP", 1, [0.0, 0.0, 0.0, 1.0])
        alias = A.view.add_alias(view, true)
        nl = Sys.iswindows() ? "\r\n" : "\n"
        expected = "View \"a\" {" * nl * "SP(0,0,0){1};" * nl * "};" * nl
        mktempdir() do directory
            for (tag, name) in ((view, "view.pos"), (alias, "alias.pos"))
                path = joinpath(directory, name)
                A.view.write(tag, path)
                @test read(path, String) == expected
            end
            # Appending adds a second record with the same line ends.
            path = joinpath(directory, "view.pos")
            A.view.write(view, path, true)
            @test read(path, String) == expected * expected
            # The other text formats keep one terminator per line as well.
            for name in ("view.txt", "view.csv")
                path = joinpath(directory, name)
                A.view.write(view, path)
                text = read(path, String)
                @test !occursin("\r\r", text) && endswith(text, nl)
            end
        end
    finally
        A.finalize()
    end
end
end
