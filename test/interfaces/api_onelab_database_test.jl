module OnelabDatabaseChecks
using Test,Tessella
const A=Tessella.API

@testset "native ONELAB database records, rollback and byte bounds" begin
    mktempdir() do directory
        runner=A._onelab_geo_runner();store=runner.store
        number=A._OnelabWireParameter(:number,"same");number.values_f=[1.5,2.]
        number.min=0.;number.max=10.;number.step=.5
        number.value_labels=Dict(1.5=>"small",2.0=>"big")
        number.label="Number λ";number.attributes["Custom"]="line1\nline2"
        string_record=A._OnelabWireParameter(:string,"same");string_record.values_s=["λ\nsecond line","tail"]
        string_record.choices_s=["a","b"];string_record.kind_name="file"
        A._onelab_store_receive!(store,number,"Gmsh")
        A._onelab_store_receive!(store,string_record,"Gmsh")
        path=joinpath(directory,"parameters.db")
        @test A._onelab_database_write(store,path)===nothing
        records=A._onelab_database_read(path)
        @test [(p.kind,p.name) for p in records]==[(:number,"same"),(:string,"same")]
        @test records[1].values_f==[1.5,2.] && records[1].value_labels==Dict(1.5=>"small",2.0=>"big")
        @test records[1].label=="Number λ" && records[1].attributes["Custom"]=="line1\nline2"
        @test records[2].values_s==["λ\nsecond line","tail"] && records[2].kind_name=="file"
        target=A._onelab_geo_runner()
        @test A._onelab_database_load!(target.store,path)===nothing
        @test target.store.parameters[(:number,"same")].values_f==[1.5,2.]
        @test target.store.parameters[(:string,"same")].values_s==["λ\nsecond line","tail"]
        # Native files terminate each counted record with NUL, then LF or
        # CRLF. Embedded newlines in values remain part of the counted data.
        bytes=IOBuffer();write(bytes,"ONELAB database control\r\n")
        for record in records
            wire=A._onelab_wire_encode(record)
            write(bytes,string(ncodeunits(wire)+1)," ",wire,UInt8(0),"\r\n")
        end
        write(path,take!(bytes))
        @test [p.values_s for p in A._onelab_database_read(path)]==[p.values_s for p in records]

        before=copy(target.store.parameters)
        write(path,"ONELAB database malformed\n1000000000 short\n")
        @test_throws ArgumentError A._onelab_database_load!(target.store,path)
        @test target.store.parameters==before
        @test_throws ArgumentError A._onelab_database_read(joinpath(directory,"missing.db"))
        for malformed in ("not a database\n","ONELAB database\n-1 x\n",
                "ONELAB database\n0 \n","ONELAB database\n3 xy\0",
                "ONELAB database\n3 xyz\n","ONELAB database\n3 xy\0x")
            write(path,malformed)
            @test_throws ArgumentError A._onelab_database_read(path)
        end
        write(path,"ONELAB database\n1000000000 short\n")
        try A._onelab_database_read(path) catch end
        allocation=@allocated try A._onelab_database_read(path) catch error
            error isa ArgumentError || rethrow()
        end
        @test allocation<65536
        @test_throws ArgumentError A._onelab_database_read(path;max_bytes=8)

        # A late invalid numerical record leaves every prior staged record
        # unpublished, rather than partially adopting the file.
        first_record=A._OnelabWireParameter(:string,"new record")
        first_record.values_s=["unpublished"]
        invalid=replace(A._onelab_wire_encode(records[1]),"1.5"=>"NaN")
        bytes=IOBuffer();write(bytes,"ONELAB database\n")
        for wire in (A._onelab_wire_encode(first_record),invalid)
            write(bytes,string(ncodeunits(wire)+1)," ",wire,UInt8(0),UInt8('\n'))
        end
        write(path,take!(bytes))
        @test_throws ArgumentError A._onelab_database_load!(target.store,path)
        @test target.store.parameters==before

        large=A._OnelabWireParameter(:string,"oversized");large.values_s=[repeat("x",1024*1024)]
        A._onelab_store_receive!(store,large,"Gmsh")
        write(path,"retained output")
        @test_throws ArgumentError A._onelab_database_write(store,path;max_bytes=128)
        allocation=@allocated try A._onelab_database_write(store,path;max_bytes=128) catch error
            error isa ArgumentError || rethrow()
        end
        @test allocation<65536
        @test read(path,String)=="retained output"
        @test sort(readdir(directory))==["parameters.db"]
        @test_throws ArgumentError A._onelab_database_write(store,directory)
        bounded=A._OnelabWireParameter(:number,"bounded")
        bounded.values_f=fill(1.23456789,10000)
        @test_throws ArgumentError A._onelab_wire_encode(bounded;max_bytes=128)
        allocation=@allocated try A._onelab_wire_encode(bounded;max_bytes=128) catch error
            error isa ArgumentError || rethrow()
        end
        @test allocation<65536
    end
end

@testset "CLI ONELAB database load and save options" begin
    mktempdir() do directory
        path=joinpath(directory,"load.geo");database=joinpath(directory,"load.db")
        runner=A._onelab_geo_runner();record=A._OnelabWireParameter(:number,"Input/N")
        record.values_f=[7.,8.];A._onelab_store_receive!(runner.store,record,"Gmsh")
        A._onelab_database_write(runner.store,database)
        write(path,"Solver.AutoLoadDatabase=1; n=DefineNumber(3,Name \"Input/N\");"*
            "Printf(\"%g\",n) >> \"values.txt\";\n")
        @test Tessella.CLI.main([path,"-run"])==path
        @test readlines(joinpath(directory,"values.txt"))==["3","7"]
        @test only(A._onelab_database_read(database)).values_f==[7.,8.]
        rm(joinpath(directory,"values.txt"))
        write(path,"n=DefineNumber(3,Name \"Input/N\"); Printf(\"%g\",n) >> \"values.txt\";\n")
        @test Tessella.CLI.main([path,"-run","-setnumber","Solver.AutoLoadDatabase","1"])==path
        @test readlines(joinpath(directory,"values.txt"))==["3","7"]
        write(database,"retained output")
        write(path,"Solver.AutoSaveDatabase=0; n=DefineNumber(3,Name \"Input/N\");\n")
        @test Tessella.CLI.main([path,"-run"])==path
        @test read(database,String)=="retained output"
        write(path,"n=DefineNumber(3,Name \"Input/N\"); Exit 0;\n")
        @test Tessella.CLI.main([path,"-run"])==path
        @test read(database,String)=="retained output"
    end
end
end
