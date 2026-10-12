module PartitionLowerBRepFunctionalChecks
using Test, LinearAlgebra, Tessella
const A=Tessella.API
const E=Tessella.Elements
const CORNERS=[0. 1. 0. 0. 0.;0. 0. 1. 0. 0.;0. 0. 0. 1. -1.]
const EXTERIOR=UInt64[11 12 13 11 13 12;12 13 11 13 12 11;14 14 14 15 15 15]
const TETS=UInt64[11 11;12 13;13 12;14 15]
const FACES=((1,3,2),(1,2,4),(2,3,4),(3,1,4))
# Curve tag, actual Line2 cells (stored directions), declared signed Point bounds.
const CURVES=[(41,[(11,12)],[1,-2]),(42,[(12,13)],[2,-3]),(43,[(13,11)],[3,-1]),
    (44,[(11,14)],[1,-4]),(45,[(12,14)],[2,-4]),(46,[(13,14)],[3,-4]),
    (47,[(11,15)],[1,-5]),(48,[(12,15)],[2,-5]),(49,[(13,15)],[3,-5])]
# Surface tag, EXTERIOR column, declared signed Curve boundary.
const SURFACES=[(31,1,[41,45,-44]),(32,2,[42,46,-45]),(33,3,[43,44,-46]),
    (34,4,[-43,49,-47]),(35,5,[-42,48,-49]),(36,6,[-41,47,-48])]

function fixture!(;surfaces=SURFACES,curves=CURVES,points=1:5,point_cells=collect(points),
                  region_bounds=[31,32,33,34,35,36],physical=true)
    A.finalize();A.initialize();A.model.add("lower BRep")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    for tag in points
        A.model.add_discrete_entity(0,tag)
    end
    for (tag,_,bounds) in curves
        A.model.add_discrete_entity(1,tag,bounds)
    end
    for (tag,column,bounds) in surfaces
        A.model.add_discrete_entity(2,tag,bounds)
    end
    A.model.add_discrete_entity(3,11,region_bounds)
    if physical
        A.model.add_physical_group(1,[41,44];tag=210,name="weld")
        A.model.add_physical_group(0,[1];tag=110,name="ports")
    end
    for tag in unique(point_cells)
        A.mesh.add_nodes(0,tag,[10+tag],CORNERS[:,tag])
    end
    for (tag,column,_) in surfaces
        A.mesh.add_elements_by_type(tag,2,[630+column],EXTERIOR[:,column])
    end
    uid=701
    for (tag,cells,_) in curves
        for (a,b) in cells
            A.mesh.add_elements_by_type(tag,1,[uid],[a,b]);uid+=1
        end
    end
    for (position,tag) in enumerate(point_cells)
        A.mesh.add_elements_by_type(tag,15,[800+position],[10+tag])
    end
    A.mesh.add_elements_by_type(11,4,[501,502],vec(TETS))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return uid-1
end
# Tri3 parts follow their carrier Tet4; Line2/Point15 parts are ignored.
partition!(last_line=709,last_point=805)=A.mesh.partition(2,
    [631:636...,701:last_line...,801:last_point...,501,502],
    vcat(Int32[1,1,1,2,2,2],fill(Int32(1),last_line-700),
         fill(Int32(1),last_point-800),Int32[1,2]))
snapshot()=(A.model.get_entities(),A.mesh.get_nodes(),A.mesh.get_elements(),
    A.mesh.get_max_node_tag(),A.mesh.get_max_element_tag(),
    A.model.get_number_of_partitions())
function child(dim,parent)
    matches=[tag for (_,tag) in A.model.get_entities(dim)
        if A.model.get_parent(dim,tag)==(dim,parent)]
    return only(matches)
end
function family(dim,owner,kind)
    types,ids,nodes=A.mesh.get_elements(dim,owner);index=findfirst(==(kind),types)
    index===nothing && return UInt64[],zeros(UInt64,E.msh_num_nodes(kind),0)
    return ids[index],reshape(nodes[index],E.msh_num_nodes(kind),length(ids[index]))
end
support(nodes)=Tuple(sort(collect(nodes)))
cycle_equal(a,b)=any(Tuple(a)==ntuple(i->b[mod1(i+shift,3)],3) for shift in 0:2)
function closed_region(tag)
    _,tets=family(3,tag,4)
    hits=Dict{Tuple,Vector{Tuple}}()
    for cell in eachcol(tets),face in FACES
        uids=ntuple(i->cell[face[i]],3);push!(get!(()->Tuple[],hits,support(uids)),uids)
    end
    exterior=Dict(key=>only(faces) for (key,faces) in hits if length(faces)==1)
    actual=Dict{Tuple,Tuple}()
    for (_,signed) in A.model.get_boundary([(3,tag)],false,true,false)
        _,tris=family(2,abs(signed),2)
        for cell in eachcol(tris)
            uids=signed>0 ? Tuple(cell) : (cell[1],cell[3],cell[2]);key=support(uids)
            haskey(actual,key) && return false
            actual[key]=uids
        end
    end
    keys(actual)==keys(exterior) && all(cycle_equal(actual[key],exterior[key]) for key in keys(actual))
end
function records(mesh)
    data=mesh.entity_data
    Dict(uid=>(Int(block.msh),Int(owner),Tuple(data.external_node_tags[block.nodes[:,column]]))
        for (block,ids,owners) in zip(mesh.blocks,data.external_element_tags,data.block_entities)
        for (column,(uid,owner)) in enumerate(zip(ids,owners)))
end
partition_records(mesh)=Dict(key=>(entity.parent,entity.partitions,entity.boundaries,entity.physical_tags)
    for (key,entity) in mesh.partition_data.entities)
function ghost_records(mesh)
    ids=reduce(vcat,mesh.entity_data.external_element_tags;init=UInt64[])
    [(ids[row.element],row.partition,row.ghost_partitions) for row in mesh.partition_data.ghost_elements]
end

@testset "lower BRep partition retains declared Curve and Point topology" begin
    try
        fixture!();original=A.mesh.get_elements()
        A.option("Mesh.PartitionCreateTopology",1);A.option("Mesh.PartitionCreateGhostCells",1)
        partition!()
        @test A.model.get_number_of_partitions()==2
        point_children=Dict(p=>child(0,p) for p in 1:5)
        curve_children=Dict(c=>child(1,c) for (c,_,_) in CURVES)
        surface_children=Dict(s=>child(2,s) for (s,_,_) in SURFACES)
        interface=only(tag for (_,tag) in A.model.get_entities(2)
            if A.model.get_parent(2,tag)==(3,11))
        region1=only(tag for (_,tag) in A.model.get_entities(3)
            if A.model.get_parent(3,tag)==(3,11) && A.model.get_partitions(3,tag)==[1])
        region2=only(tag for (_,tag) in A.model.get_entities(3)
            if A.model.get_parent(3,tag)==(3,11) && A.model.get_partitions(3,tag)==[2])
        @test A.model.get_partitions(0,point_children[1])==[1,2]
        @test A.model.get_partitions(0,point_children[2])==[1,2]
        @test A.model.get_partitions(0,point_children[3])==[1,2]
        @test A.model.get_partitions(0,point_children[4])==[1]
        @test A.model.get_partitions(0,point_children[5])==[2]
        for (curve,parts) in [(41,[1,2]),(42,[1,2]),(43,[1,2]),
                            (44,[1]),(45,[1]),(46,[1]),(47,[2]),(48,[2]),(49,[2])]
            @test A.model.get_partitions(1,curve_children[curve])==parts
        end
        @test Set(A.model.get_boundary([(1,curve_children[41])],false,true,false))==
            Set([(0,point_children[1]),(0,-point_children[2])])
        @test Set(A.model.get_boundary([(1,curve_children[46])],false,true,false))==
            Set([(0,point_children[3]),(0,-point_children[4])])
        @test Set(A.model.get_boundary([(2,surface_children[31])],false,true,false))==
            Set([(1,curve_children[41]),(1,curve_children[45]),(1,-curve_children[44])])
        @test Set(A.model.get_boundary([(2,surface_children[34])],false,true,false))==
            Set([(1,-curve_children[43]),(1,curve_children[49]),(1,-curve_children[47])])
        @test Set(abs.(last.(A.model.get_boundary([(2,interface)],false,true,false))))==
            Set([curve_children[41],curve_children[42],curve_children[43]])
        @test Set(A.model.get_boundary([(3,region1)],false,true,false))==
            Set([(2,surface_children[31]),(2,surface_children[32]),
                 (2,surface_children[33]),(2,interface)])
        @test Set(A.model.get_boundary([(3,region2)],false,true,false))==
            Set([(2,surface_children[34]),(2,surface_children[35]),
                 (2,surface_children[36]),(2,-interface)])
        @test closed_region(region1) && closed_region(region2)
        @test length(A.model.get_entities(0))==10 && length(A.model.get_entities(1))==18
        for (index,(curve,cells,_)) in enumerate(CURVES)
            ids,nodes=family(1,curve_children[curve],1)
            @test ids==UInt64[700+index] && Tuple(vec(nodes))==only(cells)
        end
        for point in 1:5
            ids,nodes=family(0,point_children[point],15)
            @test ids==UInt64[800+point] && vec(nodes)==UInt64[10+point]
            @test first(A.mesh.get_nodes(0,point_children[point]))==UInt64[10+point]
        end
        @test only(family(2,interface,2)[1])==UInt64(806)
        @test Set(A.model.get_entities_for_physical_group(1,210))==
            Set([41,44,curve_children[41],curve_children[44]])
        @test Set(A.model.get_entities_for_physical_group(0,110))==Set([1,point_children[1]])
        @test A.mesh.get_ghost_elements(3,14)==(UInt64[502],Int32[2])
        @test A.mesh.get_ghost_elements(3,15)==(UInt64[501],Int32[1])
        owned=A.mesh.get()
        @test Tessella.validate(owned).ok
        mktempdir() do directory
            for binary in (false,true)
                path=joinpath(directory,"lower-$(binary).msh")
                E.write_mixed_msh(path,owned;version=4.1,binary)
                imported=E.read_mixed_msh(path)
                @test records(imported)==records(owned)
                @test partition_records(imported)==partition_records(owned)
                @test ghost_records(imported)==ghost_records(owned)
            end
            A.option("Mesh.MshFileVersion",4.1);A.option("Mesh.Binary",1)
            A.option("Mesh.SaveAll",1)
            path=joinpath(directory,"public.msh");A.write(path)
            A.model.add("lower roundtrip");A.open(path)
            for point in 1:5
                @test A.model.get_parent(0,point_children[point])==(0,point)
            end
            @test A.model.get_parent(1,curve_children[43])==(1,43)
            @test A.model.get_partitions(1,curve_children[43])==[1,2]
            @test Set(A.model.get_boundary([(1,curve_children[43])],false,true,false))==
                Set([(0,point_children[3]),(0,-point_children[1])])
            @test Set(A.model.get_boundary([(2,surface_children[35])],false,true,false))==
                Set([(1,-curve_children[42]),(1,curve_children[48]),(1,-curve_children[49])])
            @test closed_region(region1) && closed_region(region2)
            @test family(1,curve_children[47],1)==(UInt64[707],reshape(UInt64[11,15],2,1))
            @test family(0,point_children[4],15)==(UInt64[804],reshape(UInt64[14],1,1))
            A.model.set_current("lower BRep")
        end
        A.mesh.unpartition()
        @test A.model.get_entities()==[(0,1),(0,2),(0,3),(0,4),(0,5),
            (1,41),(1,42),(1,43),(1,44),(1,45),(1,46),(1,47),(1,48),(1,49),
            (2,31),(2,32),(2,33),(2,34),(2,35),(2,36),(3,11)]
        @test A.model.get_number_of_partitions()==0
        @test A.mesh.get_elements()==original
        @test Set(A.model.get_boundary([(1,44)],false,true,false))==Set([(0,1),(0,-4)])
        @test Set(A.model.get_boundary([(2,31)],false,true,false))==Set([(1,41),(1,45),(1,-44)])
        @test A.model.get_entities_for_physical_group(1,210)==[41,44]
        @test A.model.get_entities_for_physical_group(0,110)==[1]
        @test closed_region(11)
    finally
        A.finalize()
    end
end

@testset "imported declared BRep retains signed bounds and partitions" begin
    try
        fixture!();original=A.mesh.get_elements()
        mktempdir() do directory
            A.option("Mesh.MshFileVersion",4.1);A.option("Mesh.SaveAll",1)
            A.option("Mesh.Binary",0)
            path=joinpath(directory,"declared.msh");A.write(path)
            A.model.add("imported lower BRep");A.open(path)
        end
        # Signed declared bounds survive the MSH4 transport.
        @test A.model.get_boundary([(1,41)],false,true,false)==[(0,1),(0,-2)]
        @test A.model.get_boundary([(2,31)],false,true,false)==[(1,41),(1,45),(1,-44)]
        @test A.model.get_boundary([(3,11)],false,true,false)==
            [(2,31),(2,32),(2,33),(2,34),(2,35),(2,36)]
        # Unoriented queries still normalize the stored signs.
        @test Set(A.model.get_boundary([(1,41)],false,false,false))==
            Set([(0,1),(0,2)])
        # The imported declared topology partitions like the generated one.
        A.option("Mesh.PartitionCreateTopology",1)
        A.option("Mesh.PartitionCreateGhostCells",1)
        partition!()
        @test A.model.get_number_of_partitions()==2
        point_children=Dict(p=>child(0,p) for p in 1:5)
        curve_children=Dict(c=>child(1,c) for (c,_,_) in CURVES)
        surface_children=Dict(s=>child(2,s) for (s,_,_) in SURFACES)
        interface=only(tag for (_,tag) in A.model.get_entities(2)
            if A.model.get_parent(2,tag)==(3,11))
        region1=only(tag for (_,tag) in A.model.get_entities(3)
            if A.model.get_parent(3,tag)==(3,11) && A.model.get_partitions(3,tag)==[1])
        region2=only(tag for (_,tag) in A.model.get_entities(3)
            if A.model.get_parent(3,tag)==(3,11) && A.model.get_partitions(3,tag)==[2])
        @test Set(A.model.get_boundary([(1,curve_children[41])],false,true,false))==
            Set([(0,point_children[1]),(0,-point_children[2])])
        @test Set(A.model.get_boundary([(2,surface_children[31])],false,true,false))==
            Set([(1,curve_children[41]),(1,curve_children[45]),(1,-curve_children[44])])
        @test Set(abs.(last.(A.model.get_boundary([(2,interface)],false,true,false))))==
            Set([curve_children[41],curve_children[42],curve_children[43]])
        @test closed_region(region1) && closed_region(region2)
        for (index,(curve,cells,_)) in enumerate(CURVES)
            ids,nodes=family(1,curve_children[curve],1)
            @test ids==UInt64[700+index] && Tuple(vec(nodes))==only(cells)
        end
        for point in 1:5
            ids,nodes=family(0,point_children[point],15)
            @test ids==UInt64[800+point] && vec(nodes)==UInt64[10+point]
        end
        A.mesh.unpartition()
        @test A.model.get_entities()==[(0,1),(0,2),(0,3),(0,4),(0,5),
            (1,41),(1,42),(1,43),(1,44),(1,45),(1,46),(1,47),(1,48),(1,49),
            (2,31),(2,32),(2,33),(2,34),(2,35),(2,36),(3,11)]
        @test A.mesh.get_elements()==original
        @test A.model.get_boundary([(1,44)],false,true,false)==[(0,1),(0,-4)]
        @test closed_region(11)
    finally
        A.finalize()
    end
end

@testset "closed declared Curve bounds a partition interface" begin
    try
        A.finalize();A.initialize();A.model.add("closed curve")
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        A.model.add_discrete_entity(1,41)
        A.model.add_discrete_entity(2,30,[41])
        A.model.add_discrete_entity(2,31,[-41])
        A.model.add_discrete_entity(3,11,[30,31])
        A.mesh.add_nodes(1,41,11:13,vec(CORNERS[:,1:3]))
        A.mesh.add_nodes(3,11,14:15,vec(CORNERS[:,4:5]))
        A.mesh.add_elements_by_type(30,2,[631,632,633],vec(EXTERIOR[:,1:3]))
        A.mesh.add_elements_by_type(31,2,[634,635,636],vec(EXTERIOR[:,4:6]))
        A.mesh.add_elements_by_type(41,1,[701,702,703],vec(UInt64[11 12 13;12 13 11]))
        A.mesh.add_elements_by_type(11,4,[501,502],vec(TETS))
        A.mesh.generate(3)
        A.option("Mesh.PartitionCreateTopology",1);A.option("Mesh.PartitionCreateGhostCells",1)
        A.mesh.partition(2,[631:636...,701:703...,501,502],
            Int32[1,1,1,2,2,2,1,1,1,1,2])
        rim=child(1,41);top=child(2,30);bottom=child(2,31)
        interface=only(tag for (_,tag) in A.model.get_entities(2)
            if A.model.get_parent(2,tag)==(3,11))
        part1=only(tag for (_,tag) in A.model.get_entities(3)
            if A.model.get_parent(3,tag)==(3,11) && A.model.get_partitions(3,tag)==[1])
        part2=only(tag for (_,tag) in A.model.get_entities(3)
            if A.model.get_parent(3,tag)==(3,11) && A.model.get_partitions(3,tag)==[2])
        @test A.model.get_partitions(1,rim)==[1,2]
        @test isempty(A.model.get_boundary([(1,rim)],false,true,false))
        @test Set(first(A.mesh.get_nodes(1,rim)))==Set(UInt64[11,12,13])
        @test A.model.get_boundary([(2,top)],false,true,false)==[(1,rim)]
        @test A.model.get_boundary([(2,bottom)],false,true,false)==[(1,-rim)]
        @test A.model.get_boundary([(2,interface)],false,true,false)==[(1,-rim)]
        @test family(1,rim,1)==(UInt64[701,702,703],UInt64[11 12 13;12 13 11])
        @test family(2,top,2)[1]==UInt64[631,632,633]
        @test family(2,bottom,2)[1]==UInt64[634,635,636]
        @test closed_region(part1) && closed_region(part2)
        A.mesh.unpartition()
        @test A.model.get_entities()==[(1,41),(2,30),(2,31),(3,11)]
        @test A.model.get_boundary([(2,30)],false,true,false)==[(1,41)]
        @test A.model.get_boundary([(2,31)],false,true,false)==[(1,-41)]
        @test closed_region(11)
    finally
        A.finalize()
    end
end

@testset "lower BRep topology violations reject atomically" begin
    try
        # Missing declared Curve leaves a Surface border uncovered.
        fixture!(surfaces=[(31,1,[41,45]),SURFACES[2:end]...])
        before=snapshot()
        caught=try partition!();nothing catch error;error end
        @test caught isa ArgumentError &&
            occursin("do not cover its complete border",sprint(showerror,caught))
        @test snapshot()==before
        # A declared Curve cell that is not on the Surface border is foreign.
        fixture!(curves=[CURVES...,(50,[(14,15)],[4,-5])],
            surfaces=[(31,1,[41,45,-44,50]),SURFACES[2:end]...])
        caught=try partition!(710);nothing catch error;error end
        @test caught isa ArgumentError &&
            occursin("must lie on the actual Surface border",sprint(showerror,caught))
        # An undeclared Curve is unbound.
        fixture!(curves=[CURVES...,(50,[(14,15)],[4,-5])])
        caught=try partition!(710);nothing catch error;error end
        @test caught isa ArgumentError &&
            occursin("unbound lower Curves",sprint(showerror,caught))
        # A shared Line2 support is a duplicate role.
        fixture!(curves=[CURVES[1],(50,[(11,12)],[1,-2]),CURVES[2:end]...])
        caught=try partition!(710);nothing catch error;error end
        @test caught isa ArgumentError &&
            occursin("repeats a Curve cell support",sprint(showerror,caught))
        # A Point bound whose cell is not on the Curve end is inconsistent.
        fixture!(curves=[(41,[(11,12)],[3,-2]),CURVES[2:end]...])
        caught=try partition!();nothing catch error;error end
        @test caught isa ArgumentError &&
            occursin("does not lie on the actual Curve end",sprint(showerror,caught))
        # A declared boundary Point without its Point15 cell has no carrier.
        fixture!()
        A.model.add_discrete_entity(0,6)
        A._model_locked().discrete[(1,47)].boundary=NTuple{2,Int}[(0,1),(0,-6)]
        caught=try partition!();nothing catch error;error end
        @test caught isa ArgumentError &&
            occursin("has no actual Point15 cell",sprint(showerror,caught))
        # A Point carrying two cells is ambiguous.
        fixture!(point_cells=[1,2,3,4,5,5])
        caught=try partition!(709,806);nothing catch error;error end
        @test caught isa ArgumentError &&
            occursin("exactly one actual Point15 cell",sprint(showerror,caught))
        # A single Point bound is not a signed start/end pair.
        fixture!(curves=[(41,[(11,12)],[1]),CURVES[2:end]...])
        caught=try partition!();nothing catch error;error end
        @test caught isa ArgumentError &&
            occursin("signed start and end Point boundaries",sprint(showerror,caught))
        # A branched Curve chain cannot bind endpoints.
        fixture!(curves=[(41,[(11,12),(11,13),(12,13)],[1,-3]),CURVES[4:end]...],
            surfaces=[(31,1,[41,45,-44]),(32,2,[46,-45]),(33,3,[44,-46]),
                (34,4,[49,-47]),(35,5,[48,-49]),(36,6,[47,-48])])
        caught=try partition!();nothing catch error;error end
        @test caught isa ArgumentError &&
            occursin("nonbranching directed chain",sprint(showerror,caught))
        # A closed loop cannot carry declared Point bounds.
        A.finalize();A.initialize();A.model.add("closed bound")
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        A.model.add_discrete_entity(0,1)
        A.model.add_discrete_entity(1,41,[1])
        A.model.add_discrete_entity(2,30,[41])
        A.model.add_discrete_entity(2,31,[-41])
        A.model.add_discrete_entity(3,11,[30,31])
        A.mesh.add_nodes(0,1,[11],CORNERS[:,1])
        A.mesh.add_nodes(1,41,12:13,vec(CORNERS[:,2:3]))
        A.mesh.add_nodes(3,11,14:15,vec(CORNERS[:,4:5]))
        A.mesh.add_elements_by_type(30,2,[631,632,633],vec(EXTERIOR[:,1:3]))
        A.mesh.add_elements_by_type(31,2,[634,635,636],vec(EXTERIOR[:,4:6]))
        A.mesh.add_elements_by_type(41,1,[701,702,703],vec(UInt64[11 12 13;12 13 11]))
        A.mesh.add_elements_by_type(1,15,[801],[11])
        A.mesh.add_elements_by_type(11,4,[501,502],vec(TETS))
        A.mesh.generate(3)
        A.option("Mesh.PartitionCreateTopology",1)
        caught=try
            A.mesh.partition(2,[631:636...,701:703...,801,501,502],
                Int32[1,1,1,2,2,2,1,1,1,1,1,2]);nothing
        catch error;error end
        @test caught isa ArgumentError &&
            occursin("cannot have Point boundaries",sprint(showerror,caught))
    finally
        A.finalize()
    end
end

const MTOP=Tessella.MeshEntityTopology
function mixed_facets(msh)
    out=Tuple{Int,Tuple{Vararg{Int}}}[]
    for pattern in MTOP._simplex_face_patterns(msh,3);push!(out,(2,pattern));end
    for pattern in MTOP._simplex_face_patterns(msh,4);push!(out,(3,pattern));end
    return out
end
# A child's closed boundary: every singly-carried volume facet appears exactly
# once, wound outward, among the cells of its signed boundary surfaces.
function closed_mixed_region(tag)
    types,ids,nodes=A.mesh.get_elements(3,tag)
    hits=Dict{Tuple{Vararg{UInt64}},Vector{Tuple{Vararg{UInt64}}}}()
    for (msh,flat) in zip(types,nodes)
        msh in (4,5,6,7) || continue
        cells=reshape(collect(UInt64,flat),E.msh_num_nodes(msh),:)
        for cell in eachcol(cells),(face_msh,pattern) in mixed_facets(Int(msh))
            uids=ntuple(i->cell[pattern[i]],length(pattern))
            push!(get!(()->Tuple{Vararg{UInt64}}[],hits,_partition_support(uids)),uids)
        end
    end
    exterior=Dict(key=>only(faces) for (key,faces) in hits if length(faces)==1)
    actual=Dict{Tuple{Vararg{UInt64}},Tuple{Vararg{UInt64}}}()
    for (_,signed) in A.model.get_boundary([(3,tag)],false,true,false)
        st,si,sn=A.mesh.get_elements(2,abs(signed))
        for (msh,flat) in zip(st,sn)
            cells=reshape(collect(UInt64,flat),E.msh_num_nodes(msh),:)
            for cell in eachcol(cells)
                uids=signed>0 ? Tuple(cell) :
                    ntuple(i->i==1 ? cell[1] : cell[length(cell)+2-i],length(cell))
                key=_partition_support(uids)
                haskey(actual,key) && return false
                actual[key]=uids
            end
        end
    end
    keys(actual)==keys(exterior) &&
        all(cycle_any(actual[key],exterior[key]) for key in keys(actual))
end
_partition_support(uids)=Tuple(sort(collect(uids)))
cycle_any(a,b)=length(a)==length(b) &&
    any(a==ntuple(i->b[mod1(i+shift,length(b))],length(b)) for shift in 0:length(b)-1)

# Unit cube Hex8 with a pyramid cap on its top quad; surfaces 31-35 are the cube
# sides (Quad4), 36-39 the pyramid laterals (Tri3). Curves 41-44 and points
# 1,2,5,6 declare the border of surface 32.
function mixed_fixture!(;region_bounds=[31,32,33,34,35,36,37,38,39],
                         hex_nodes=[1,2,3,4,5,6,7,8],extra_quad=nothing)
    A.finalize();A.initialize();A.model.add("mixed hex pyr")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    for tag in (1,2,5,6);A.model.add_discrete_entity(0,tag);end
    A.model.add_discrete_entity(1,41,[1,-2]);A.model.add_discrete_entity(1,42,[2,-6])
    A.model.add_discrete_entity(1,43,[6,-5]);A.model.add_discrete_entity(1,44,[5,-1])
    A.model.add_discrete_entity(2,31);A.model.add_discrete_entity(2,32,[41,42,43,44])
    for tag in 33:39;A.model.add_discrete_entity(2,tag);end
    A.model.add_discrete_entity(3,11,region_bounds)
    A.model.add_physical_group(2,[32];tag=220,name="plated")
    coords=[0. 1. 1. 0. 0. 1. 1. 0. .5;
            0. 0. 1. 1. 0. 0. 1. 1. .5;
            0. 0. 0. 0. 1. 1. 1. 1. 1.5]
    A.mesh.add_nodes(0,1,[1],coords[:,1]);A.mesh.add_nodes(0,2,[2],coords[:,2])
    A.mesh.add_nodes(0,5,[5],coords[:,5]);A.mesh.add_nodes(0,6,[6],coords[:,6])
    A.mesh.add_nodes(3,11,[3,4,7,8,9],vec(coords[:,[3,4,7,8,9]]))
    A.mesh.add_elements_by_type(31,3,[601],[1,4,3,2])
    A.mesh.add_elements_by_type(32,3,[602],[1,2,6,5])
    A.mesh.add_elements_by_type(33,3,[603],[1,5,8,4])
    A.mesh.add_elements_by_type(34,3,[604],[2,3,7,6])
    A.mesh.add_elements_by_type(35,3,[605],[3,4,8,7])
    A.mesh.add_elements_by_type(36,2,[606],[5,6,9])
    A.mesh.add_elements_by_type(37,2,[607],[6,7,9])
    A.mesh.add_elements_by_type(38,2,[608],[8,5,9])
    A.mesh.add_elements_by_type(39,2,[609],[7,8,9])
    extra_quad===nothing || A.mesh.add_elements_by_type(31,3,[610],extra_quad)
    for (tag,uid,edge) in ((41,701,[1,2]),(42,702,[2,6]),(43,703,[6,5]),(44,704,[5,1]))
        A.mesh.add_elements_by_type(tag,1,[uid],edge)
    end
    A.mesh.add_elements_by_type(1,15,[801],[1]);A.mesh.add_elements_by_type(2,15,[802],[2])
    A.mesh.add_elements_by_type(5,15,[805],[5]);A.mesh.add_elements_by_type(6,15,[806],[6])
    A.mesh.add_elements_by_type(11,5,[501],hex_nodes)
    A.mesh.add_elements_by_type(11,7,[502],[5,6,7,8,9])
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
end
# Tet4 (A,B,C,D) plus a Prism6 (B',D',C',B,D,C) sharing the Tri3 face B,C,D.
function prism_fixture!()
    A.finalize();A.initialize();A.model.add("tet prism")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    for tag in 31:37;A.model.add_discrete_entity(2,tag);end
    A.model.add_discrete_entity(3,11,[31,32,33,34,35,36,37])
    coords=[0. 1. 0. 0. 2. 1. 1.;
            0. 0. 1. 0. 0. 0. 1.;
            0. 0. 0. 1. 0. 1. 0.]
    A.mesh.add_nodes(3,11,1:7,vec(coords))
    A.mesh.add_elements_by_type(31,2,[601],[1,3,2])
    A.mesh.add_elements_by_type(32,2,[602],[1,2,4])
    A.mesh.add_elements_by_type(33,2,[603],[1,4,3])
    A.mesh.add_elements_by_type(34,2,[604],[5,7,6])
    A.mesh.add_elements_by_type(35,3,[605],[5,6,4,2])
    A.mesh.add_elements_by_type(36,3,[606],[5,2,3,7])
    A.mesh.add_elements_by_type(37,3,[607],[6,7,3,4])
    A.mesh.add_elements_by_type(11,4,[501],[1,2,3,4])
    A.mesh.add_elements_by_type(11,6,[502],[5,6,7,2,4,3])
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
end

@testset "mixed linear families partition declared BRep" begin
    try
        # Hex8 + Pyr5 over Quad4 and Tri3 surfaces, with a declared Curve
        # border on one Quad4 side and a Quad4 partition interface.
        mixed_fixture!();original=A.mesh.get_elements()
        A.option("Mesh.PartitionCreateGhostCells",1)
        A.mesh.partition(2,[601:609...,701:704...,801,802,805,806,501,502],
            Int32[1,1,1,1,1,2,2,2,2,1,1,1,1,1,1,1,1,1,2])
        @test A.model.get_number_of_partitions()==2
        surface_children=Dict(s=>only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(2,s)) for s in 31:39)
        curve_children=Dict(c=>child(1,c) for c in 41:44)
        point_children=Dict(p=>child(0,p) for p in (1,2,5,6))
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        region1=only(tag for (dim,tag) in A.model.get_entities()
            if dim==3 && A.model.get_parent(3,tag)==(3,11) && A.model.get_partitions(3,tag)==[1])
        region2=only(tag for (dim,tag) in A.model.get_entities()
            if dim==3 && A.model.get_parent(3,tag)==(3,11) && A.model.get_partitions(3,tag)==[2])
        @test A.model.get_partitions(2,interface)==[1,2]
        @test family(2,interface,3)==(UInt64[807],reshape(UInt64[5,6,7,8],4,1))
        @test isempty(family(2,interface,2)[1])
        # The interface border binds declared Curve 43 (shared by Quad4 and
        # Tri3 carriers) and three generated seams.
        interface_curves=Set(abs.(last.(A.model.get_boundary([(2,interface)],false,true,false))))
        @test curve_children[43] in interface_curves && length(interface_curves)==4
        seams=setdiff(interface_curves,[curve_children[43]])
        @test Set(A.model.get_parent.(1,collect(seams)))==Set([(2,33),(2,34),(2,35)])
        @test all(A.model.get_partitions(1,seam)==[1,2] for seam in seams)
        @test Set(A.model.get_boundary([(3,region1)],false,true,false))==
            Set((2,surface_children[s]) for s in 31:35) ∪ Set([(2,interface)])
        @test Set(A.model.get_boundary([(3,region2)],false,true,false))==
            Set((2,surface_children[s]) for s in 36:39) ∪ Set([(2,-interface)])
        @test Set(A.model.get_boundary([(2,surface_children[32])],false,true,false))==
            Set([(1,curve_children[41]),(1,curve_children[42]),
                 (1,curve_children[43]),(1,curve_children[44])])
        @test A.model.get_boundary([(2,surface_children[31])],false,true,false)==
            [(1,-curve_children[41])]
        @test A.model.get_boundary([(2,surface_children[36])],false,true,false)==
            [(1,-curve_children[43])]
        @test A.model.get_partitions(1,curve_children[43])==[1,2]
        @test A.model.get_partitions(1,curve_children[41])==[1]
        @test A.model.get_partitions(0,point_children[5])==[1,2]
        @test A.model.get_partitions(0,point_children[1])==[1]
        @test Set(A.model.get_boundary([(1,curve_children[43])],false,true,false))==
            Set([(0,point_children[6]),(0,-point_children[5])])
        @test closed_mixed_region(region1) && closed_mixed_region(region2)
        for (tag,uid,edge) in ((41,701,[1,2]),(42,702,[2,6]),(43,703,[6,5]),(44,704,[5,1]))
            ids,nodes=family(1,curve_children[tag],1)
            @test ids==UInt64[uid] && Tuple(vec(nodes))==Tuple(edge)
        end
        for (point,uid,node) in ((1,801,1),(2,802,2),(5,805,5),(6,806,6))
            ids,nodes=family(0,point_children[point],15)
            @test ids==UInt64[uid] && vec(nodes)==UInt64[node]
        end
        @test A.model.get_entities_for_physical_group(2,220)==[32,surface_children[32]]
        ghost1=only(tag for (tag,part) in A.mesh.get().partition_data.ghost_entities if part==1)
        ghost2=only(tag for (tag,part) in A.mesh.get().partition_data.ghost_entities if part==2)
        @test A.mesh.get_ghost_elements(3,ghost1)==(UInt64[502],Int32[2])
        @test A.mesh.get_ghost_elements(3,ghost2)==(UInt64[501],Int32[1])
        owned=A.mesh.get()
        @test Tessella.validate(owned).ok
        mktempdir() do directory
            for binary in (false,true)
                path=joinpath(directory,"mixed-$(binary).msh")
                E.write_mixed_msh(path,owned;version=4.1,binary)
                imported=E.read_mixed_msh(path)
                @test records(imported)==records(owned)
                @test partition_records(imported)==partition_records(owned)
            end
            A.option("Mesh.MshFileVersion",4.1);A.option("Mesh.Binary",0)
            A.option("Mesh.SaveAll",1)
            path=joinpath(directory,"public.msh");A.write(path)
            A.model.add("mixed roundtrip");A.open(path)
            @test family(2,interface,3)==(UInt64[807],reshape(UInt64[5,6,7,8],4,1))
            @test Set(A.model.get_boundary([(3,region2)],false,true,false))==
                Set((2,surface_children[s]) for s in 36:39) ∪ Set([(2,-interface)])
            @test closed_mixed_region(region1) && closed_mixed_region(region2)
            A.model.set_current("mixed hex pyr")
        end
        A.mesh.unpartition()
        @test A.mesh.get_elements()==original
        @test A.model.get_entities()==[(0,1),(0,2),(0,5),(0,6),
            (1,41),(1,42),(1,43),(1,44),(2,31),(2,32),(2,33),(2,34),(2,35),
            (2,36),(2,37),(2,38),(2,39),(3,11)]
        @test A.model.get_boundary([(2,32)],false,true,false)==[(1,41),(1,42),(1,43),(1,44)]
        # Tet4 + Prism6 share a Tri3 face: a Tri3 interface between mixed
        # volume families, with Tri3 and Quad4 outer surfaces.
        prism_fixture!();original=A.mesh.get_elements()
        A.option("Mesh.PartitionCreateGhostCells",1)
        A.mesh.partition(2,[601:607...,501,502],Int32[1,1,1,2,2,2,2,1,2])
        @test A.model.get_number_of_partitions()==2
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        @test A.model.get_partitions(2,interface)==[1,2]
        @test family(2,interface,2)==(UInt64[608],reshape(UInt64[4,2,3],3,1))
        @test isempty(family(2,interface,3)[1])
        interface_curves=Set(abs.(last.(A.model.get_boundary([(2,interface)],false,true,false))))
        @test length(interface_curves)==3
        @test Set(A.model.get_parent.(1,collect(interface_curves)))==
            Set([(2,31),(2,32),(2,33)])
        @test all(A.model.get_partitions(1,seam)==[1,2] for seam in interface_curves)
        region1=only(tag for (dim,tag) in A.model.get_entities()
            if dim==3 && A.model.get_parent(3,tag)==(3,11) && A.model.get_partitions(3,tag)==[1])
        region2=only(tag for (dim,tag) in A.model.get_entities()
            if dim==3 && A.model.get_parent(3,tag)==(3,11) && A.model.get_partitions(3,tag)==[2])
        surface_children=Dict(s=>only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(2,s)) for s in 31:37)
        @test Set(A.model.get_boundary([(3,region1)],false,true,false))==
            Set((2,surface_children[s]) for s in 31:33) ∪ Set([(2,interface)])
        @test Set(A.model.get_boundary([(3,region2)],false,true,false))==
            Set((2,surface_children[s]) for s in 34:37) ∪ Set([(2,-interface)])
        @test closed_mixed_region(region1) && closed_mixed_region(region2)
        owned=A.mesh.get()
        @test Tessella.validate(owned).ok
        A.mesh.unpartition()
        @test A.mesh.get_elements()==original
        @test A.model.get_entities()==[(2,31),(2,32),(2,33),(2,34),(2,35),(2,36),(2,37),(3,11)]
    finally
        A.finalize()
    end
end

@testset "mixed family topology violations reject atomically" begin
    try
        # A reversed Hex8 fails the positive-volume certificate.
        mixed_fixture!(hex_nodes=[1,4,3,2,5,8,7,6])
        before=snapshot()
        caught=try
            A.mesh.partition(2,[601:609...,701:704...,801,802,805,806,501,502],
                Int32[1,1,1,1,1,2,2,2,2,1,1,1,1,1,1,1,1,1,2]);nothing
        catch error;error end
        @test caught isa ArgumentError &&
            occursin("positive oriented volume cells",sprint(showerror,caught))
        @test snapshot()==before
        # A duplicated Quad4 support is a repeated Surface cell role.
        mixed_fixture!(extra_quad=[1,2,3,4])
        before=snapshot()
        caught=try
            A.mesh.partition(2,[601:610...,701:704...,801,802,805,806,501,502],
                Int32[1,1,1,1,1,2,2,2,2,1,1,1,1,1,1,1,1,1,1,2]);nothing
        catch error;error end
        @test caught isa ArgumentError &&
            occursin("repeats a Surface cell support",sprint(showerror,caught))
        @test snapshot()==before
        # A Region boundary omitting a carried Quad4 leaves facets uncovered.
        mixed_fixture!(region_bounds=[31,32,33,34,35,36,37,38])
        before=snapshot()
        caught=try
            A.mesh.partition(2,[601:609...,701:704...,801,802,805,806,501,502],
                Int32[1,1,1,1,1,2,2,2,2,1,1,1,1,1,1,1,1,1,2]);nothing
        catch error;error end
        @test caught isa ArgumentError &&
            occursin("omits Region facets",sprint(showerror,caught))
        @test snapshot()==before
        # A Prism6 wound backwards is equally rejected before any mutation.
        A.finalize();A.initialize();A.model.add("bad prism")
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        for tag in 31:37;A.model.add_discrete_entity(2,tag);end
        A.model.add_discrete_entity(3,11,[31,32,33,34,35,36,37])
        coords=[0. 1. 0. 0. 2. 1. 1.;0. 0. 1. 0. 0. 0. 1.;0. 0. 0. 1. 0. 1. 0.]
        A.mesh.add_nodes(3,11,1:7,vec(coords))
        A.mesh.add_elements_by_type(31,2,[601],[1,3,2])
        A.mesh.add_elements_by_type(32,2,[602],[1,2,4])
        A.mesh.add_elements_by_type(33,2,[603],[1,4,3])
        A.mesh.add_elements_by_type(34,2,[604],[5,7,6])
        A.mesh.add_elements_by_type(35,3,[605],[5,6,4,2])
        A.mesh.add_elements_by_type(36,3,[606],[5,2,3,7])
        A.mesh.add_elements_by_type(37,3,[607],[6,7,3,4])
        A.mesh.add_elements_by_type(11,4,[501],[1,2,3,4])
        A.mesh.add_elements_by_type(11,6,[502],[5,7,6,2,3,4])
        A.mesh.generate(3)
        A.option("Mesh.PartitionCreateTopology",1)
        before=snapshot()
        caught=try
            A.mesh.partition(2,[601:607...,501,502],Int32[1,1,1,2,2,2,2,1,2]);nothing
        catch error;error end
        @test caught isa ArgumentError &&
            occursin("positive oriented volume cells",sprint(showerror,caught))
        @test snapshot()==before
    finally
        A.finalize()
    end
end

# Second-order fixtures build a conforming P2 mesh literally: edge nodes are
# keyed by their corner pair and quadrangle face nodes by the corner set, so
# volume facets, surface cells and curve cells share actual nodes.
struct _P2Nodes
    cv::Dict{Int,NTuple{3,Float64}}
    map::Dict{Any,Int32}
    coords::Vector{NTuple{3,Float64}}
end
_P2Nodes()=_P2Nodes(Dict{Int,NTuple{3,Float64}}(),Dict{Any,Int32}(),NTuple{3,Float64}[])
function _p2_corner!(p,tag,pos)
    p.cv[tag]=pos
    get!(p.map,(:v,tag)) do;push!(p.coords,pos);Int32(length(p.coords));end
end
function _p2_mid(p,a,b)
    key=minmax(a,b)
    get!(p.map,(:e,key)) do
        push!(p.coords,(p.cv[key[1]].+p.cv[key[2]])./2);Int32(length(p.coords))
    end
end
function _p2_center(p,corners...)
    key=Tuple(sort!(collect(corners)))
    get!(p.map,(:f,key)) do
        push!(p.coords,ntuple(d->sum(p.cv[c][d] for c in corners)/length(corners),3))
        Int32(length(p.coords))
    end
end
function _p2_body(p,id,corners)
    get!(p.map,(:b,id)) do
        push!(p.coords,ntuple(d->sum(p.cv[c][d] for c in corners)/length(corners),3))
        Int32(length(p.coords))
    end
end
_p2_row(p,tag)=p.map[(:v,tag)]
const _P2_TET_EDGES=((1,2),(2,3),(3,1),(4,1),(4,3),(4,2))
const _P2_HEX_EDGES=((1,2),(1,4),(1,5),(2,3),(2,6),(3,4),
                     (3,7),(4,8),(5,6),(5,8),(6,7),(7,8))
const _P2_HEX_FACETS=((1,4,3,2),(1,2,6,5),(1,5,8,4),(2,3,7,6),(3,4,8,7),(5,6,7,8))
const _P2_PRI_EDGES=((1,2),(1,3),(1,4),(2,3),(2,5),(3,6),(4,5),(4,6),(5,6))
_p2_tet(p,c)=Int32[[_p2_row(p,t) for t in c]...;
    [_p2_mid(p,c[a],c[b]) for (a,b) in _P2_TET_EDGES]...]
_p2_hex(p,c,id)=Int32[[_p2_row(p,t) for t in c]...;
    [_p2_mid(p,c[a],c[b]) for (a,b) in _P2_HEX_EDGES]...;
    [_p2_center(p,c[f[1]],c[f[2]],c[f[3]],c[f[4]]) for f in _P2_HEX_FACETS]...;
    _p2_body(p,id,c)]
_p2_pri(p,c)=Int32[[_p2_row(p,t) for t in c]...;
    [_p2_mid(p,c[a],c[b]) for (a,b) in _P2_PRI_EDGES]...;
    [_p2_center(p,c[1],c[2],c[5],c[4]),_p2_center(p,c[1],c[4],c[6],c[3]),
     _p2_center(p,c[2],c[3],c[6],c[5])]...]
_p2_tri(p,c)=Int32[_p2_row(p,c[1]),_p2_row(p,c[2]),_p2_row(p,c[3]),
    _p2_mid(p,c[1],c[2]),_p2_mid(p,c[2],c[3]),_p2_mid(p,c[3],c[1])]
_p2_quad(p,c)=Int32[_p2_row(p,c[1]),_p2_row(p,c[2]),_p2_row(p,c[3]),_p2_row(p,c[4]),
    _p2_mid(p,c[1],c[2]),_p2_mid(p,c[2],c[3]),_p2_mid(p,c[3],c[4]),
    _p2_mid(p,c[4],c[1]),_p2_center(p,c...)]
_p2_line(p,a,b)=Int32[_p2_row(p,a),_p2_row(p,b),_p2_mid(p,a,b)]
function _p2_nodes!(p)
    tags=Int32[];xyz=Float64[]
    for (i,c) in enumerate(p.coords)
        push!(tags,Int32(i));append!(xyz,collect(c))
    end
    A.mesh.add_nodes(3,11,tags,xyz)
end
_p2_edge_uid(p,a,b)=UInt64(p.map[(:e,minmax(a,b))])

# Two Tet10 over Tri6 surfaces; a closed declared Line3 Curve runs along the
# inward-wound base Surface 45 so one interface border edge binds the declared
# Curve and the other two become a generated Line3 seam.
function p2_tets_fixture!(;curve_maker=_p2_line,tet_b=_p2_tet)
    A.finalize();A.initialize();A.model.add("p2 tets")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    p=_P2Nodes()
    for (tag,pos) in ((1,(0.,0.,0.)),(2,(1.,0.,0.)),(3,(0.,1.,0.)),
                      (4,(0.,0.,1.)),(7,(1.,1.,1.)))
        _p2_corner!(p,tag,pos)
    end
    t1=_p2_tet(p,(1,2,3,4));t2=tet_b(p,(2,3,4,7))
    base=_p2_tri(p,(1,3,2))
    base=Int32[base[1],base[3],base[2],base[6],base[5],base[4]]
    s41=hcat(_p2_tri(p,(1,2,4)),_p2_tri(p,(1,4,3)))
    c2=(2,3,4,7)
    s42=hcat([_p2_tri(p,(c2[f[1]],c2[f[2]],c2[f[3]]))
              for f in ((1,2,4),(1,4,3),(4,2,3))]...)
    l46=hcat(curve_maker(p,1,3),curve_maker(p,3,2),curve_maker(p,2,1))
    A.model.add_discrete_entity(1,46,Int[])
    A.model.add_discrete_entity(2,41,Int[])
    A.model.add_discrete_entity(2,42,Int[])
    A.model.add_discrete_entity(2,45,[-46])
    A.model.add_discrete_entity(3,11,[41,-45,42])
    A.model.add_physical_group(1,[46];tag=210,name="rim")
    _p2_nodes!(p)
    A.mesh.add_elements_by_type(41,9,[701,702],vec(s41))
    A.mesh.add_elements_by_type(42,9,[703,704,705],vec(s42))
    A.mesh.add_elements_by_type(45,9,[706],vec(base))
    A.mesh.add_elements_by_type(46,8,[801,802,803],vec(l46))
    A.mesh.add_elements_by_type(11,11,[501,502],vcat(t1,t2))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return p
end
# Two Hex27 share the x=1 Quad9; every interface rim edge becomes a Line3 seam.
function p2_hexes_fixture!()
    A.finalize();A.initialize();A.model.add("p2 hexes")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    p=_P2Nodes()
    h1=(1,2,3,4,5,6,7,8)
    for (i,pos) in enumerate(((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),
                             (0.,0.,1.),(1.,0.,1.),(1.,1.,1.),(0.,1.,1.)))
        _p2_corner!(p,i,pos)
    end
    for (i,pos) in ((21,(2.,0.,0.)),(22,(2.,1.,0.)),(23,(2.,0.,1.)),(24,(2.,1.,1.)))
        _p2_corner!(p,i,pos)
    end
    h2=(2,21,22,3,6,23,24,7)
    c1=_p2_hex(p,h1,1);c2=_p2_hex(p,h2,2)
    s41=hcat([_p2_quad(p,(h1[f[1]],h1[f[2]],h1[f[3]],h1[f[4]]))
              for f in _P2_HEX_FACETS[[1,2,3,5,6]]]...)
    s42=hcat([_p2_quad(p,(h2[f[1]],h2[f[2]],h2[f[3]],h2[f[4]]))
              for f in _P2_HEX_FACETS[[1,2,4,5,6]]]...)
    A.model.add_discrete_entity(2,41,Int[])
    A.model.add_discrete_entity(2,42,Int[])
    A.model.add_discrete_entity(3,11,[41,42])
    _p2_nodes!(p)
    A.mesh.add_elements_by_type(41,10,[701:705...],vec(s41))
    A.mesh.add_elements_by_type(42,10,[706:710...],vec(s42))
    A.mesh.add_elements_by_type(11,12,[501,502],vcat(c1,c2))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return p
end
# A Prism18 and a Tet10 share one Tri6; Surface 41 keeps a Tri6 and Quad9 block.
function p2_tet_prism_fixture!()
    A.finalize();A.initialize();A.model.add("p2 tet prism")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    p=_P2Nodes()
    for (tag,pos) in ((1,(0.,0.,0.)),(2,(1.,0.,0.)),(3,(0.,1.,0.)),
                      (4,(0.,0.,1.)),(5,(1.,0.,1.)),(6,(0.,1.,1.)),
                      (7,(1/3,1/3,-1.)))
        _p2_corner!(p,tag,pos)
    end
    pr=_p2_pri(p,(1,2,3,4,5,6));tt=_p2_tet(p,(1,3,2,7))
    s41_tri=_p2_tri(p,(4,5,6))
    s41_quad=hcat(_p2_quad(p,(1,2,5,4)),_p2_quad(p,(1,4,6,3)),_p2_quad(p,(2,3,6,5)))
    tc=(1,3,2,7)
    s42=hcat([_p2_tri(p,(tc[f[1]],tc[f[2]],tc[f[3]]))
              for f in ((1,2,4),(1,4,3),(4,2,3))]...)
    A.model.add_discrete_entity(2,41,Int[])
    A.model.add_discrete_entity(2,42,Int[])
    A.model.add_discrete_entity(3,11,[41,42])
    _p2_nodes!(p)
    A.mesh.add_elements_by_type(41,9,[701],vec(s41_tri))
    A.mesh.add_elements_by_type(41,10,[702,703,704],vec(s41_quad))
    A.mesh.add_elements_by_type(42,9,[705,706,707],vec(s42))
    A.mesh.add_elements_by_type(11,13,[501],vec(pr))
    A.mesh.add_elements_by_type(11,11,[502],vec(tt))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return p
end

@testset "second-order cells partition declared BRep" begin
    try
        # Tet10 pair: a Tri6 interface, a declared closed Line3 Curve on the
        # interface border and one generated Line3 seam for the rest.
        p=p2_tets_fixture!();original=A.mesh.get_elements()
        A.option("Mesh.PartitionCreateGhostCells",1)
        A.mesh.partition(2,[501,502,701:706...,801:803...],
            Int32[1,2,1,1,2,2,2,1,1,1,1])
        @test A.model.get_number_of_partitions()==2
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        @test A.model.get_partitions(2,interface)==[1,2]
        ids,nodes=family(2,interface,9)
        @test length(ids)==1 && size(nodes)==(6,1)
        @test sort(vec(nodes))==sort(UInt64[p.map[(:v,c)] for c in (2,3,4)] ∪
            [_p2_edge_uid(p,2,3),_p2_edge_uid(p,3,4),_p2_edge_uid(p,4,2)])
        curve_child=child(1,46)
        interface_curves=Set(abs.(last.(A.model.get_boundary([(2,interface)],false,true,false))))
        @test curve_child in interface_curves && length(interface_curves)==2
        seam=only(setdiff(interface_curves,[curve_child]))
        @test A.model.get_parent(1,seam)==(2,41) && A.model.get_partitions(1,seam)==[1,2]
        ids,nodes=family(1,seam,8)
        @test length(ids)==2 && size(nodes)==(3,2)
        @test Set(Tuple(sort!(collect(UInt64.(column)))) for column in
            eachcol(nodes[1:2,:]))==
            Set([Tuple(sort!(UInt64[p.map[(:v,4)],p.map[(:v,2)]])),
                 Tuple(sort!(UInt64[p.map[(:v,3)],p.map[(:v,4)]]))])
        ids,nodes=family(1,curve_child,8)
        @test ids==UInt64[801,802,803] && size(nodes)==(3,3)
        @test all(nodes[3,column] in (_p2_edge_uid(p,1,3),_p2_edge_uid(p,2,3),
                                      _p2_edge_uid(p,1,2)) for column in 1:3)
        @test A.model.get_partitions(1,curve_child)==[1,2]
        @test A.model.get_entities_for_physical_group(1,210)==[46,curve_child]
        region1=only(tag for (dim,tag) in A.model.get_entities()
            if dim==3 && A.model.get_parent(3,tag)==(3,11) &&
               A.model.get_partitions(3,tag)==[1])
        region2=only(tag for (dim,tag) in A.model.get_entities()
            if dim==3 && A.model.get_parent(3,tag)==(3,11) &&
               A.model.get_partitions(3,tag)==[2])
        surface_children=Dict(s=>only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(2,s)) for s in (41,42,45))
        @test Set(A.model.get_boundary([(3,region1)],false,true,false))==
            Set([(2,surface_children[41]),(2,-surface_children[45]),(2,interface)])
        @test Set(A.model.get_boundary([(3,region2)],false,true,false))==
            Set([(2,surface_children[42]),(2,-interface)])
        ids,nodes=family(3,region1,11);@test ids==UInt64[501] && size(nodes)==(10,1)
        ids,nodes=family(3,region2,11);@test ids==UInt64[502] && size(nodes)==(10,1)
        ghost1=only(tag for (tag,part) in A.mesh.get().partition_data.ghost_entities if part==1)
        ghost2=only(tag for (tag,part) in A.mesh.get().partition_data.ghost_entities if part==2)
        @test A.mesh.get_ghost_elements(3,ghost1)==(UInt64[502],Int32[2])
        @test A.mesh.get_ghost_elements(3,ghost2)==(UInt64[501],Int32[1])
        owned=A.mesh.get()
        @test Tessella.validate(owned).ok
        mktempdir() do directory
            for binary in (false,true)
                path=joinpath(directory,"p2-$(binary).msh")
                E.write_mixed_msh(path,owned;version=4.1,binary)
                imported=E.read_mixed_msh(path)
                @test records(imported)==records(owned)
                @test partition_records(imported)==partition_records(owned)
            end
            A.option("Mesh.MshFileVersion",4.1);A.option("Mesh.Binary",0)
            A.option("Mesh.SaveAll",1)
            path=joinpath(directory,"public.msh");A.write(path)
            A.model.add("p2 roundtrip");A.open(path)
            @test family(1,curve_child,8)[1]==UInt64[801,802,803]
            @test Set(A.model.get_boundary([(3,region2)],false,true,false))==
                Set([(2,surface_children[42]),(2,-interface)])
            A.model.set_current("p2 tets")
        end
        A.mesh.unpartition()
        @test A.mesh.get_elements()==original
        @test A.model.get_entities()==[(1,46),(2,41),(2,42),(2,45),(3,11)]
        @test A.model.get_boundary([(2,45)],false,true,false)==[(1,-46)]
        # Hex27 pair: one Quad9 interface and a closed four-cell Line3 seam.
        p2_hexes_fixture!();original=A.mesh.get_elements()
        A.option("Mesh.PartitionCreateGhostCells",1)
        A.mesh.partition(2,[501,502,701:710...],Int32[1,2,1,1,1,1,1,2,2,2,2,2])
        @test A.model.get_number_of_partitions()==2
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        ids,nodes=family(2,interface,10)
        @test length(ids)==1 && size(nodes)==(9,1)
        interface_curves=Set(abs.(last.(A.model.get_boundary([(2,interface)],false,true,false))))
        @test length(interface_curves)==1
        seam=only(interface_curves)
        ids,nodes=family(1,seam,8)
        @test length(ids)==4 && size(nodes)==(3,4)
        @test A.model.get_partitions(1,seam)==[1,2]
        owned=A.mesh.get()
        @test Tessella.validate(owned).ok
        A.mesh.unpartition()
        @test A.mesh.get_elements()==original
        # Tet10 + Prism18 share a Tri6; the Prism Surface child keeps Tri6 and
        # Quad9 blocks and the interface is a single Tri6 cell.
        p2_tet_prism_fixture!();original=A.mesh.get_elements()
        A.mesh.partition(2,[501,502,701:707...],Int32[1,2,1,1,1,1,2,2,2])
        @test A.model.get_number_of_partitions()==2
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        ids,nodes=family(2,interface,9)
        @test length(ids)==1 && size(nodes)==(6,1)
        prism_surface=child(2,41)
        types,blocks,_=A.mesh.get_elements(2,prism_surface)
        @test Set(types)==Set(Int32[9,10])
        owned=A.mesh.get()
        @test Tessella.validate(owned).ok
        A.mesh.unpartition()
        @test A.mesh.get_elements()==original
    finally
        A.finalize()
    end
end

@testset "second-order topology violations reject atomically" begin
    try
        # A declared Line3 Curve carrying a foreign edge node is not conforming.
        foreign(p,a,b)=Int32[_p2_row(p,a),_p2_row(p,b),
            get!(p.map,(:x,a,b)) do
                push!(p.coords,(9.,9.,9.));Int32(length(p.coords))
            end]
        p2_tets_fixture!(curve_maker=foreign)
        before=snapshot()
        caught=try
            A.mesh.partition(2,[501,502,701:706...,801:803...],
                Int32[1,2,1,1,2,2,2,1,1,1,1]);nothing
        catch error;error end
        @test caught isa ArgumentError &&
            occursin("edge node disagrees",sprint(showerror,caught))
        @test snapshot()==before
        # An incomplete serendipity family is outside the supported catalog.
        p=_P2Nodes()
        A.finalize();A.initialize();A.model.add("serendipity")
        A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
        for (tag,pos) in ((1,(0.,0.,0.)),(2,(1.,0.,0.)),(3,(0.,1.,0.)),
                          (4,(0.,0.,1.)),(7,(1.,1.,1.)))
            _p2_corner!(p,tag,pos)
        end
        t1=_p2_tet(p,(1,2,3,4));t2=_p2_tet(p,(2,3,4,7))
        base=_p2_tri(p,(1,3,2))
        s41=hcat(_p2_tri(p,(1,2,4)),_p2_tri(p,(1,4,3)))
        c2=(2,3,4,7)
        s42=hcat([_p2_tri(p,(c2[f[1]],c2[f[2]],c2[f[3]]))
                  for f in ((1,2,4),(1,4,3),(4,2,3))]...)
        q8=Int32[_p2_row(p,1),_p2_row(p,2),_p2_row(p,3),_p2_row(p,4),
            _p2_mid(p,1,2),_p2_mid(p,2,3),_p2_mid(p,3,4),_p2_mid(p,4,1)]
        A.model.add_discrete_entity(2,41,Int[])
        A.model.add_discrete_entity(2,42,Int[])
        A.model.add_discrete_entity(2,45,Int[])
        A.model.add_discrete_entity(2,47,Int[])
        A.model.add_discrete_entity(3,11,[41,45,42,47])
        _p2_nodes!(p)
        A.mesh.add_elements_by_type(41,9,[701,702],vec(s41))
        A.mesh.add_elements_by_type(42,9,[703,704,705],vec(s42))
        A.mesh.add_elements_by_type(45,9,[706],vec(base))
        A.mesh.add_elements_by_type(47,16,[707],vec(q8))
        A.mesh.add_elements_by_type(11,11,[501,502],vcat(t1,t2))
        A.mesh.generate(3)
        A.option("Mesh.PartitionCreateTopology",1)
        before=snapshot()
        caught=try
            A.mesh.partition(2,[501,502,701:707...],
                Int32[1,2,1,1,2,2,2,1,1]);nothing
        catch error;error end
        @test caught isa ArgumentError &&
            occursin("complete cell families only",sprint(showerror,caught))
        @test snapshot()==before
        # A reversed Tet10 fails the positive-volume certificate.
        flipped(p,c)=Int32[[_p2_row(p,t) for t in (c[1],c[3],c[2],c[4])]...;
            [_p2_mid(p,c[a],c[b]) for (a,b) in ((1,3),(3,2),(2,1),(4,1),(4,2),(4,3))]...]
        p2_tets_fixture!(tet_b=flipped)
        before=snapshot()
        caught=try
            A.mesh.partition(2,[501,502,701:706...,801:803...],
                Int32[1,2,1,1,2,2,2,1,1,1,1]);nothing
        catch error;error end
        @test caught isa ArgumentError &&
            occursin("positive oriented volume cells",sprint(showerror,caught))
        @test snapshot()==before
    finally
        A.finalize()
    end
end

# Junction fixtures: tet1 shares face {1,2,3} with tet2 and {1,2,4} with tet3 —
# rim edge {1,2} closes a three-volume junction fan. The declared form splits
# the shared triangle {1,2,4} onto surfaces 30/31 and declares its edges as
# Curves 40-42 with Points 7-9 at its vertices.
const JCORNERS=Float64[0 0 1 0 -1 -0.5;0 0 0 1 -1 -1;0 1 0 0 0.2 0.2]
function junction_fixture!(;declared=false,fan=false)
    A.finalize();A.initialize();A.model.add("junction")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    if fan
        A.model.add_discrete_entity(2,30)
        A.model.add_discrete_entity(3,11,[30])
        fc=Float64[0 0 1 0 -1 0;0 0 0 1 0 -1;0 1 0 0 0 0]
        A.mesh.add_nodes(3,11,collect(1:6),vec(fc))
        A.mesh.add_elements_by_type(30,2,UInt64[601:608...],
            UInt64[1,4,3, 4,2,3, 1,5,4, 5,2,4, 1,6,5, 6,2,5, 1,3,6, 3,2,6])
        A.mesh.add_elements_by_type(11,4,UInt64[501:504...],
            UInt64[1,2,3,4, 1,2,4,5, 1,2,5,6, 1,2,6,3])
    else
        if declared
            for tag in (7,8,9);A.model.add_discrete_entity(0,tag);end
            A.model.add_discrete_entity(1,40,[7,-8])
            A.model.add_discrete_entity(1,41,[8,-9])
            A.model.add_discrete_entity(1,42,[9,-7])
            A.model.add_discrete_entity(2,30,[-40,-41,-42])
            A.model.add_discrete_entity(2,31,[40,41,42])
            A.model.add_discrete_entity(3,11,[30,31])
            A.mesh.add_nodes(3,11,collect(1:6),vec(JCORNERS))
            A.mesh.add_elements_by_type(30,2,UInt64[601:605...],
                UInt64[1,4,3, 4,2,3, 1,5,2, 1,3,5, 5,3,2])
            A.mesh.add_elements_by_type(31,2,UInt64[606:608...],
                UInt64[1,6,4, 1,2,6, 6,2,4])
            A.mesh.add_elements_by_type(40,1,UInt64[701],UInt64[1,2])
            A.mesh.add_elements_by_type(41,1,UInt64[702],UInt64[2,4])
            A.mesh.add_elements_by_type(42,1,UInt64[703],UInt64[4,1])
            A.mesh.add_elements_by_type(7,15,UInt64[801],UInt64[1])
            A.mesh.add_elements_by_type(8,15,UInt64[802],UInt64[2])
            A.mesh.add_elements_by_type(9,15,UInt64[803],UInt64[4])
        else
            A.model.add_discrete_entity(2,30)
            A.model.add_discrete_entity(3,11,[30])
            A.mesh.add_nodes(3,11,collect(1:6),vec(JCORNERS))
            A.mesh.add_elements_by_type(30,2,UInt64[601:608...],
                UInt64[1,4,3, 4,2,3, 1,5,2, 1,3,5, 5,3,2, 1,6,4, 1,2,6, 6,2,4])
        end
        A.mesh.add_elements_by_type(11,4,UInt64[501:503...],
            UInt64[1,2,3,4, 1,3,2,5, 1,2,4,6])
    end
    A.mesh.generate(3);A.option("Mesh.PartitionCreateTopology",1)
end
junction_tags(;declared=false)=declared ?
    UInt64[601:608...;701:703...;801:803...;501:503...] :
    UInt64[601:608...;501:503...]
junction_parts(;declared=false)=declared ?
    Int32[1,1,2,2,2,3,3,3, 1,1,1, 1,1,1, 1,2,3] :
    Int32[1,1,2,2,2,3,3,3, 1,2,3]
bound_abs(dim,tag)=[abs(v) for (_,v) in A.model.get_boundary([(dim,tag)],false,true,false)]
# Same junction in complete second order: Tet10 volumes over Tri6 carriers —
# the emitted junction seam is a Line3 carrying the junction edge's mid node.
function p2_junction_fixture!()
    A.finalize();A.initialize();A.model.add("p2 junction")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    p=_P2Nodes()
    for (tag,pos) in ((1,(0.,0.,0.)),(2,(0.,0.,1.)),(3,(1.,0.,0.)),
                      (4,(0.,1.,0.)),(5,(-1.,-1.,0.2)),(6,(-0.5,-1.,0.2)))
        _p2_corner!(p,tag,pos)
    end
    t1=_p2_tet(p,(1,2,3,4));t2=_p2_tet(p,(1,3,2,5));t3=_p2_tet(p,(1,2,4,6))
    s30=hcat(_p2_tri(p,(1,4,3)),_p2_tri(p,(4,2,3)),
             _p2_tri(p,(1,5,2)),_p2_tri(p,(1,3,5)),_p2_tri(p,(5,3,2)),
             _p2_tri(p,(1,6,4)),_p2_tri(p,(1,2,6)),_p2_tri(p,(6,2,4)))
    A.model.add_discrete_entity(2,30)
    A.model.add_discrete_entity(3,11,[30])
    _p2_nodes!(p)
    A.mesh.add_elements_by_type(30,9,[601:608...],vec(s30))
    A.mesh.add_elements_by_type(11,11,[501,502,503],vcat(t1,t2,t3))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return p
end

@testset "partition junctions emit shared seams and junction points" begin
    try
        # Public Gmsh 4.15.2 emits one seam entity shared by every incident
        # interface on a junction edge plus dim-0 entities at seam-graph
        # junction vertices; the rim certificate only requires every incident
        # volume to contribute exactly two facets on the edge.
        junction_fixture!()
        A.mesh.partition(3,junction_tags(),junction_parts())
        jpoints=sort!(collect(tag for (_,tag) in A.model.get_entities(0)
            if A.model.get_parent(0,tag)==(2,30)))
        @test length(jpoints)==2
        @test all(A.model.get_partitions(0,tag)==[1,2,3] for tag in jpoints)
        point_nodes=Set{UInt64}()
        for tag in jpoints
            _,nodes=family(0,tag,15);union!(point_nodes,vec(nodes))
        end
        @test point_nodes==Set(UInt64[1,2])
        jseam=only(tag for (_,tag) in A.model.get_entities(1)
            if length(A.model.get_partitions(1,tag))==3)
        @test A.model.get_parent(1,jseam)==(2,30)
        _,lines=family(1,jseam,1)
        @test Set(support(lines[:,column]) for column in axes(lines,2))==
            Set([Tuple(UInt64[1,2])])
        bounds=A.model.get_boundary([(1,jseam)],false,true,false)
        @test sort([abs(v) for (_,v) in bounds])==sort(jpoints)
        @test count(v->v>0,(v for (_,v) in bounds))==1
        rims=sort!(collect(tag for (_,tag) in A.model.get_entities(1)
            if length(A.model.get_partitions(1,tag))==2))
        @test [A.model.get_partitions(1,tag) for tag in rims]==[[1,2],[1,3]]
        @test all(sort(bound_abs(1,tag))==sort(jpoints) for tag in rims)
        interfaces=sort!(collect(tag for (_,tag) in A.model.get_entities(2)
            if A.model.get_parent(2,tag)==(3,11)))
        @test length(interfaces)==2
        @test [A.model.get_partitions(2,tag) for tag in interfaces]==[[1,2],[1,3]]
        @test all(jseam in bound_abs(2,tag) for tag in interfaces)
        @test all(closed_region(tag) for (_,tag) in A.model.get_entities(3)
            if A.model.get_parent(3,tag)==(3,11))
        owned=A.mesh.get()
        @test Tessella.validate(owned).ok
        mktempdir() do directory
            for binary in (false,true)
                path=joinpath(directory,"junction-$(binary).msh")
                E.write_mixed_msh(path,owned;version=4.1,binary)
                imported=E.read_mixed_msh(path)
                @test records(imported)==records(owned)
                @test partition_records(imported)==partition_records(owned)
            end
        end
        A.mesh.unpartition()
        @test Set(A.model.get_entities())==Set([(2,30),(3,11)])
        # The four-volume fan has no exterior carriers on its junction edge —
        # the seam adopts the Region parent and the full partition membership.
        junction_fixture!(fan=true)
        A.mesh.partition(4,UInt64[601:608...;501:504...],
            Int32[1,1,2,2,3,3,4,4,1,2,3,4])
        fseam=only(tag for (_,tag) in A.model.get_entities(1)
            if length(A.model.get_partitions(1,tag))==4)
        @test A.model.get_parent(1,fseam)==(3,11)
        fpoints=[tag for (_,tag) in A.model.get_entities(0)]
        @test length(fpoints)==2
        @test all(A.model.get_partitions(0,tag)==[1,2,3,4] for tag in fpoints)
        finterfaces=[tag for (_,tag) in A.model.get_entities(2)
            if A.model.get_parent(2,tag)==(3,11)]
        @test length(finterfaces)==4
        @test all(fseam in bound_abs(2,tag) for tag in finterfaces)
        @test Tessella.validate(A.mesh.get()).ok
        # A declared Curve can itself be the junction edge: every incident
        # interface binds it and its declared Point children serve as junctions.
        junction_fixture!(declared=true)
        A.mesh.partition(3,junction_tags(declared=true),junction_parts(declared=true))
        @test Set(tag for (_,tag) in A.model.get_entities(0)
            if A.model.get_parent(0,tag)[1]==0)==Set([child(0,7),child(0,8),child(0,9)])
        @test A.model.get_partitions(0,child(0,7))==[1,2,3]
        @test A.model.get_partitions(0,child(0,9))==[1,3]
        jcurve=child(1,40)
        @test A.model.get_partitions(1,jcurve)==[1,2,3]
        dinterfaces=[tag for (_,tag) in A.model.get_entities(2)
            if A.model.get_parent(2,tag)==(3,11)]
        @test length(dinterfaces)==2
        @test all(jcurve in bound_abs(2,tag) for tag in dinterfaces)
        @test Tessella.validate(A.mesh.get()).ok
        A.mesh.unpartition()
        @test (0,7) in A.model.get_entities(0) && (1,40) in A.model.get_entities(1)
        # A second-order junction emits a Line3 junction seam carrying the
        # junction edge's actual mid node; junction points stay corner-only.
        p=p2_junction_fixture!()
        A.mesh.partition(3,UInt64[501:503...;601:608...],
            Int32[1,2,3,1,1,2,2,2,3,3,3])
        pseam=only(tag for (_,tag) in A.model.get_entities(1)
            if length(A.model.get_partitions(1,tag))==3)
        types,_,pnodes=A.mesh.get_elements(1,pseam)
        @test types==Int32[8]
        seam_nodes=reshape(only(pnodes),3,1)
        @test support(seam_nodes[1:2])==Tuple(UInt64[1,2])
        @test seam_nodes[3]==_p2_edge_uid(p,1,2)
        @test length([tag for (_,tag) in A.model.get_entities(0)])==2
        @test Tessella.validate(A.mesh.get()).ok
        # Junction fixtures still reject atomically: an empty partition and a
        # Surface cell coinciding with a junction interface both hold.
        junction_fixture!()
        before=snapshot()
        caught=try
            A.mesh.partition(4,junction_tags(),junction_parts());nothing
        catch error;error end
        @test caught isa ArgumentError && snapshot()==before
        junction_fixture!(declared=true)
        A.mesh.add_elements_by_type(31,2,UInt64[609],UInt64[1,2,4])
        A.mesh.generate(3)
        before=snapshot()
        coincident_parts=vcat(junction_parts(declared=true)[1:8],Int32(3),
                              junction_parts(declared=true)[9:end])
        caught=try
            A.mesh.partition(3,UInt64[601:609...;701:703...;801:803...;501:503...],
                coincident_parts);nothing
        catch error;error end
        @test caught isa ArgumentError && snapshot()==before
    finally
        A.finalize()
    end
end

# Two Regions sharing a facet carried by declared Surface 30: T1=(1,2,3,4) in
# Region 11 and T2=(1,2,4,5) in Region 12 tile around node 5; surfaces 31-36
# close the exterior. With rim=true, declared Curves 40-42 and Points 1,2,4
# bound the shared surface.
function coincident_fixture!(;rim=false,one_sided=false)
    A.finalize();A.initialize();A.model.add("coincident")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    if rim
        for tag in (1,2,4);A.model.add_discrete_entity(0,tag);end
        A.model.add_discrete_entity(1,40,[1,-2]);A.model.add_discrete_entity(1,41,[2,-4])
        A.model.add_discrete_entity(1,42,[4,-1])
        A.model.add_discrete_entity(2,30,[40,41,42])
    else
        A.model.add_discrete_entity(2,30)
    end
    for tag in 31:36;A.model.add_discrete_entity(2,tag);end
    A.model.add_discrete_entity(3,11,[30,31,32,33])
    A.model.add_discrete_entity(3,12,one_sided ? [34,35,36] : [-30,34,35,36])
    coords=[0. 1. 0. 0. .5;0. 0. 1. 0. -.5;0. 0. 0. 1. .2]
    A.mesh.add_nodes(3,11,1:5,vec(coords))
    A.mesh.add_elements_by_type(30,2,[601],[1,2,4])
    A.mesh.add_elements_by_type(31,2,[602],[1,3,2])
    A.mesh.add_elements_by_type(32,2,[603],[1,4,3])
    A.mesh.add_elements_by_type(33,2,[604],[4,2,3])
    A.mesh.add_elements_by_type(34,2,[605],[1,2,5])
    A.mesh.add_elements_by_type(35,2,[606],[1,5,4])
    A.mesh.add_elements_by_type(36,2,[607],[5,2,4])
    if rim
        for (tag,uid,edge) in ((40,701,[1,2]),(41,702,[2,4]),(42,703,[4,1]))
            A.mesh.add_elements_by_type(tag,1,[uid],edge)
        end
        for (tag,uid,node) in ((1,801,[1]),(2,802,[2]),(4,804,[4]))
            A.mesh.add_elements_by_type(tag,15,[uid],node)
        end
    end
    A.mesh.add_elements_by_type(11,4,[501],[1,2,3,4])
    A.mesh.add_elements_by_type(12,4,[502],[1,2,4,5])
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
end
coincident_tags(;rim=false)=rim ? UInt64[601:607...;701:703...;801;802;804;501;502] :
    UInt64[601:607...;501;502]
coincident_parts(;rim=false)=rim ?
    Int32[1,1,1,1,2,2,2,1,1,1,1,1,1,1,2] : Int32[1,1,1,1,2,2,2,1,2]

@testset "coincident Surface partition keeps one cell with merged membership" begin
    try
        # A declared Surface cell on a cross-part facet already carries the
        # interface role: its partition child keeps the cell, merges the
        # incident partition membership and binds both Region children with
        # opposite signs. No duplicate interface entity or cell is emitted —
        # validate's duplicate-cell certificate stays intact.
        coincident_fixture!();original=A.mesh.get_elements()
        A.mesh.partition(2,coincident_tags(),coincident_parts())
        @test A.model.get_number_of_partitions()==2
        shared=child(2,30)
        @test A.model.get_partitions(2,shared)==[1,2]
        @test family(2,shared,2)==(UInt64[601],reshape(UInt64[1,2,4],3,1))
        region1=only(tag for (_,tag) in A.model.get_entities(3)
            if A.model.get_parent(3,tag)==(3,11))
        region2=only(tag for (_,tag) in A.model.get_entities(3)
            if A.model.get_parent(3,tag)==(3,12))
        bound1=Set(A.model.get_boundary([(3,region1)],false,true,false))
        bound2=Set(A.model.get_boundary([(3,region2)],false,true,false))
        @test (2,shared) in bound1 && (2,-shared) in bound2
        @test closed_region(region1) && closed_region(region2)
        # No generated interface entity: every dim-2 entity is an original
        # surface child; no seam entities exist without real interfaces.
        @test isempty(A.model.get_entities(1))
        @test isempty([tag for (_,tag) in A.model.get_entities(2)
            if A.model.get_parent(2,tag)[1]==3])
        owned=A.mesh.get()
        @test Tessella.validate(owned).ok
        mktempdir() do directory
            for binary in (false,true)
                path=joinpath(directory,"coincident-$(binary).msh")
                E.write_mixed_msh(path,owned;version=4.1,binary)
                imported=E.read_mixed_msh(path)
                @test records(imported)==records(owned)
                @test partition_records(imported)==partition_records(owned)
            end
        end
        A.mesh.unpartition()
        @test A.mesh.get_elements()==original
        @test Set(A.model.get_entities())==Set([(2,s) for s in 30:36] ∪ [(3,11),(3,12)])
        # Declared Curves and Points on the coincident rim: the shared child
        # binds the declared Curve children, which merge the same [1,2]
        # membership their rim carriers imply.
        coincident_fixture!(rim=true)
        A.mesh.partition(2,coincident_tags(rim=true),coincident_parts(rim=true))
        rshared=child(2,30)
        @test A.model.get_partitions(2,rshared)==[1,2]
        rcurves=[child(1,c) for c in 40:42]
        @test Set(A.model.get_boundary([(2,rshared)],false,true,false))==
            Set([(1,c) for c in rcurves])
        @test all(A.model.get_partitions(1,c)==[1,2] for c in rcurves)
        rpoints=[child(0,p) for p in (1,2,4)]
        @test all(length(A.model.get_partitions(0,p))==2 for p in rpoints)
        @test Tessella.validate(A.mesh.get()).ok
        A.mesh.unpartition()
        @test (1,40) in A.model.get_entities(1) && (0,1) in A.model.get_entities(0)
        # One-sided declaration still rejects atomically at coverage: the
        # shared facet must be covered by each Region's declared boundary.
        coincident_fixture!(one_sided=true)
        before=snapshot()
        caught=try
            A.mesh.partition(2,coincident_tags(),coincident_parts());nothing
        catch error;error end
        @test caught isa ArgumentError &&
            occursin("omits Region facets",sprint(showerror,caught))
        @test snapshot()==before
    finally
        A.finalize()
    end
end

# Third-order fixtures use coordinate-keyed nodes: edge nodes sit at i/3
# fractions, triangle interiors at centroids and quadrangle interiors at
# bilinear (i,j)/3 positions, so volume facets, surface cells and curve cells
# share actual nodes.
struct _P3Nodes
    cv::Dict{Int,NTuple{3,Float64}}
    map::Dict{Any,Int32}
    coords::Vector{NTuple{3,Float64}}
end
_P3Nodes()=_P3Nodes(Dict{Int,NTuple{3,Float64}}(),Dict{Any,Int32}(),NTuple{3,Float64}[])
function _p3_corner!(p,tag,pos)
    p.cv[tag]=pos
    get!(p.map,(:v,tag)) do;push!(p.coords,pos);Int32(length(p.coords));end
end
function _p3_at(p,pos)
    key=(:p,round.(pos;digits=9))
    get!(p.map,key) do;push!(p.coords,pos);Int32(length(p.coords));end
end
_p3_edge(p,a,b,i)=_p3_at(p,ntuple(d->p.cv[a][d]+(p.cv[b][d]-p.cv[a][d])*i/3,3))
_p3_center(p,c)=_p3_at(p,ntuple(d->sum(p.cv[x][d] for x in c)/length(c),3))
function _p3_quad_interior(p,c,i,j)
    a,b,d=c[1],c[2],c[4]
    _p3_at(p,ntuple(dim->p.cv[a][dim]+(p.cv[b][dim]-p.cv[a][dim])*i/3+
        (p.cv[d][dim]-p.cv[a][dim])*j/3,3))
end
_p3_row(p,tag)=p.map[(:v,tag)]
_p3_edge_uid(p,a,b,i)=UInt64(p.map[(:p,round.(ntuple(d->p.cv[a][d]+
    (p.cv[b][d]-p.cv[a][d])*i/3,3);digits=9))])
const _P3_TET_FACES=((1,3,2),(1,2,4),(1,4,3),(4,2,3))
const _P3_PRI_TRIS=((1,3,2),(4,5,6))
const _P3_PRI_QUADS=((1,2,5,4),(1,4,6,3),(2,3,6,5))
const _P3_PYR_EDGES=((1,2),(1,4),(1,5),(2,3),(2,5),(3,4),(3,5),(4,5))
const _P3_PYR_TRIS=((1,2,5),(4,1,5),(2,3,5),(3,4,5))
const _P3_QINT=((1,1),(2,1),(2,2),(1,2))
_p3_tri(p,c)=Int32[[_p3_row(p,x) for x in c]...;
    [_p3_edge(p,c[a],c[b],i) for (a,b) in ((1,2),(2,3),(3,1)) for i in (1,2)]...;
    _p3_center(p,c)]
_p3_quad(p,c)=Int32[[_p3_row(p,x) for x in c]...;
    [_p3_edge(p,c[a],c[b],i) for (a,b) in ((1,2),(2,3),(3,4),(4,1)) for i in (1,2)]...;
    [_p3_quad_interior(p,c,i,j) for (i,j) in _P3_QINT]...]
_p3_line(p,a,b)=Int32[_p3_row(p,a),_p3_row(p,b),_p3_edge(p,a,b,1),_p3_edge(p,a,b,2)]
_p3_tet(p,c)=Int32[[_p3_row(p,x) for x in c]...;
    [_p3_edge(p,c[a],c[b],i) for (a,b) in _P2_TET_EDGES for i in (1,2)]...;
    [_p3_center(p,ntuple(i->c[f[i]],3)) for f in _P3_TET_FACES]...]
function _p3_hex(p,c)
    a,b,d,e=c[1],c[2],c[4],c[5]
    Int32[[_p3_row(p,x) for x in c]...;
        [_p3_edge(p,c[x],c[y],i) for (x,y) in _P2_HEX_EDGES for i in (1,2)]...;
        [[_p3_quad_interior(p,ntuple(i->c[f[i]],4),i,j) for (i,j) in _P3_QINT]
         for f in _P2_HEX_FACETS]...;
        [_p3_at(p,ntuple(dim->p.cv[a][dim]+(p.cv[b][dim]-p.cv[a][dim])*i/3+
            (p.cv[d][dim]-p.cv[a][dim])*j/3+(p.cv[e][dim]-p.cv[a][dim])*k/3,3))
         for k in (1,2) for (i,j) in _P3_QINT]...]
end
function _p3_pri(p,c)
    bot=ntuple(d->sum(p.cv[c[i]][d] for i in 1:3)/3,3)
    top=ntuple(d->sum(p.cv[c[i]][d] for i in 4:6)/3,3)
    Int32[[_p3_row(p,x) for x in c]...;
        [_p3_edge(p,c[a],c[b],i) for (a,b) in _P2_PRI_EDGES for i in (1,2)]...;
        [_p3_center(p,ntuple(i->c[f[i]],3)) for f in _P3_PRI_TRIS]...;
        [[_p3_quad_interior(p,ntuple(i->c[f[i]],4),i,j) for (i,j) in _P3_QINT]
         for f in _P3_PRI_QUADS]...;
        [_p3_at(p,ntuple(d->bot[d]+(top[d]-bot[d])*r/3,3)) for r in (1,2)]...]
end
_p3_pyr(p,c)=Int32[[_p3_row(p,x) for x in c]...;
    [_p3_edge(p,c[a],c[b],i) for (a,b) in _P3_PYR_EDGES for i in (1,2)]...;
    [_p3_center(p,ntuple(i->c[f[i]],3)) for f in _P3_PYR_TRIS]...;
    [_p3_quad_interior(p,(c[1],c[4],c[3],c[2]),i,j) for (i,j) in _P3_QINT]...;
    _p3_center(p,c)]
function _p3_nodes!(p)
    tags=Int32[];xyz=Float64[]
    for (i,c) in enumerate(p.coords)
        push!(tags,Int32(i));append!(xyz,collect(c))
    end
    A.mesh.add_nodes(3,11,tags,xyz)
end
# Two Tet20 share the {2,3,4} Tri10; a closed declared Line4 Curve rims the
# interface facet.
function p3_tets_fixture!(;extras=()->nothing)
    A.finalize();A.initialize();A.model.add("p3 tets")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    p=_P3Nodes()
    for (tag,pos) in ((1,(0.,0.,0.)),(2,(1.,0.,0.)),(3,(0.,1.,0.)),
                      (4,(0.,0.,1.)),(7,(1.,1.,1.)))
        _p3_corner!(p,tag,pos)
    end
    t1=_p3_tet(p,(1,2,3,4));t2=_p3_tet(p,(2,3,4,7))
    tc=(2,3,4,7)
    s41=hcat([_p3_tri(p,ntuple(i->(1,2,3,4)[f[i]],3))
              for f in _P3_TET_FACES[1:3]]...)
    s42=hcat([_p3_tri(p,ntuple(i->tc[f[i]],3)) for f in _P3_TET_FACES[2:4]]...)
    l46=hcat(_p3_line(p,2,3),_p3_line(p,3,4),_p3_line(p,4,2))
    A.model.add_discrete_entity(1,46,Int[])
    A.model.add_discrete_entity(2,41,[-46])
    A.model.add_discrete_entity(2,42,[46])
    A.model.add_discrete_entity(3,11,[41,42])
    _p3_nodes!(p)
    A.mesh.add_elements_by_type(46,26,[801,802,803],vec(l46))
    A.mesh.add_elements_by_type(41,21,[701,702,703],vec(s41))
    A.mesh.add_elements_by_type(42,21,[704,705,706],vec(s42))
    A.mesh.add_elements_by_type(11,29,[501,502],vcat(t1,t2))
    extras()
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return p
end
# Two Hex64 share the x=1 Quad16; every interface rim edge becomes a Line4 seam.
function p3_hexes_fixture!()
    A.finalize();A.initialize();A.model.add("p3 hexes")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    p=_P3Nodes()
    h1=(1,2,3,4,5,6,7,8)
    for (i,pos) in enumerate(((0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.,1.,0.),
                             (0.,0.,1.),(1.,0.,1.),(1.,1.,1.),(0.,1.,1.)))
        _p3_corner!(p,i,pos)
    end
    for (i,pos) in ((21,(2.,0.,0.)),(22,(2.,1.,0.)),(23,(2.,0.,1.)),(24,(2.,1.,1.)))
        _p3_corner!(p,i,pos)
    end
    h2=(2,21,22,3,6,23,24,7)
    c1=_p3_hex(p,h1);c2=_p3_hex(p,h2)
    s41=hcat([_p3_quad(p,ntuple(i->h1[f[i]],4)) for f in _P2_HEX_FACETS[[1,2,3,5,6]]]...)
    s42=hcat([_p3_quad(p,ntuple(i->h2[f[i]],4)) for f in _P2_HEX_FACETS[[1,2,4,5,6]]]...)
    A.model.add_discrete_entity(2,41,Int[])
    A.model.add_discrete_entity(2,42,Int[])
    A.model.add_discrete_entity(3,11,[41,42])
    _p3_nodes!(p)
    A.mesh.add_elements_by_type(41,36,[701:705...],vec(s41))
    A.mesh.add_elements_by_type(42,36,[706:710...],vec(s42))
    A.mesh.add_elements_by_type(11,92,[501,502],vcat(c1,c2))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return p
end
# Two Prism40 share the {1,2,4,5} Quad16; exterior surfaces mix Tri10/Quad16.
function p3_prisms_fixture!()
    A.finalize();A.initialize();A.model.add("p3 prisms")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    p=_P3Nodes()
    for (tag,pos) in ((1,(0.,0.,0.)),(2,(1.,0.,0.)),(3,(0.,1.,0.)),
                      (4,(0.,0.,1.)),(5,(1.,0.,1.)),(6,(0.,1.,1.)),
                      (7,(0.,-1.,0.)),(10,(0.,-1.,1.)))
        _p3_corner!(p,tag,pos)
    end
    pr1=_p3_pri(p,(1,2,3,4,5,6));pr2=_p3_pri(p,(2,1,7,5,4,10))
    pc=(2,1,7,5,4,10)
    s41t=hcat(_p3_tri(p,(1,3,2)),_p3_tri(p,(4,5,6)))
    s41q=hcat(_p3_quad(p,(1,4,6,3)),_p3_quad(p,(2,3,6,5)))
    s42t=hcat(_p3_tri(p,(pc[1],pc[3],pc[2])),_p3_tri(p,(pc[4],pc[5],pc[6])))
    s42q=hcat([_p3_quad(p,ntuple(i->pc[f[i]],4)) for f in _P3_PRI_QUADS[2:3]]...)
    A.model.add_discrete_entity(2,41,Int[])
    A.model.add_discrete_entity(2,42,Int[])
    A.model.add_discrete_entity(3,11,[41,42])
    _p3_nodes!(p)
    A.mesh.add_elements_by_type(41,21,[701,702],vec(s41t))
    A.mesh.add_elements_by_type(41,36,[703,704],vec(s41q))
    A.mesh.add_elements_by_type(42,21,[705,706],vec(s42t))
    A.mesh.add_elements_by_type(42,36,[707,708],vec(s42q))
    A.mesh.add_elements_by_type(11,90,[501,502],vcat(pr1,pr2))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return p
end
# Two Pyr30 share the Quad16 base; the four side facets are Tri10.
function p3_pyrs_fixture!()
    A.finalize();A.initialize();A.model.add("p3 pyrs")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    p=_P3Nodes()
    for (tag,pos) in ((1,(0.,0.,0.)),(2,(1.,0.,0.)),(3,(1.,1.,0.)),(4,(0.,1.,0.)),
                      (5,(0.5,0.5,1.)),(6,(0.5,0.5,-1.)))
        _p3_corner!(p,tag,pos)
    end
    py1=_p3_pyr(p,(1,2,3,4,5));py2=_p3_pyr(p,(4,3,2,1,6))
    yc=(4,3,2,1,6)
    s41=hcat([_p3_tri(p,c) for c in _P3_PYR_TRIS]...)
    s42=hcat([_p3_tri(p,ntuple(i->yc[f[i]],3)) for f in _P3_PYR_TRIS]...)
    A.model.add_discrete_entity(2,41,Int[])
    A.model.add_discrete_entity(2,42,Int[])
    A.model.add_discrete_entity(3,11,[41,42])
    _p3_nodes!(p)
    A.mesh.add_elements_by_type(41,21,[701:704...],vec(s41))
    A.mesh.add_elements_by_type(42,21,[705:708...],vec(s42))
    A.mesh.add_elements_by_type(11,118,[501,502],vcat(py1,py2))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return p
end

# Order >=4 fixtures share one coordinate-keyed builder: node positions come
# from the monomial lattice — volume corners verbatim, edge nodes on the
# canonical edge table, face-interior slots through each facet's own frame so
# declared Surface cells share the exact node set, and remaining slots through
# the volume family's own map.
mutable struct _HONodes
    map::Dict{NTuple{3,Int64},Int32}
    coords::Vector{NTuple{3,Float64}}
end
_HONodes()=_HONodes(Dict{NTuple{3,Int64},Int32}(),NTuple{3,Float64}[])
function _ho_at!(p,c)
    key=ntuple(d->round(Int64,c[d]*1e9),3)
    get!(p.map,key) do;push!(p.coords,c);Int32(length(p.coords));end
end
function _ho_pos(family,order,m,c)
    if family==:tri
        x,y=m;w=(order-x-y,x,y)./order
        return ntuple(d->sum(w[k]*c[k][d] for k in 1:3),3)
    elseif family==:qua
        x,y=m;w=((order-x)*(order-y),x*(order-y),x*y,(order-x)*y)./(order*order)
        return ntuple(d->sum(w[k]*c[k][d] for k in 1:4),3)
    elseif family==:tet
        x,y,z=m;w=(order-x-y-z,x,y,z)./order
        return ntuple(d->sum(w[k]*c[k][d] for k in 1:4),3)
    elseif family==:hex
        x,y,z=m
        w=((order-x)*(order-y)*(order-z),x*(order-y)*(order-z),
           x*y*(order-z),(order-x)*y*(order-z),
           (order-x)*(order-y)*z,x*(order-y)*z,x*y*z,(order-x)*y*z)./(order^3)
        return ntuple(d->sum(w[k]*c[k][d] for k in 1:8),3)
    elseif family==:pri
        x,y,z=m;w=((order-x-y)*(order-z),x*(order-z),y*(order-z),
                  (order-x-y)*z,x*z,y*z)./(order*order)
        return ntuple(d->sum(w[k]*c[k][d] for k in 1:6),3)
    elseif family==:pyr
        x,y,z=m;t=z/order
        base=ntuple(d->(((order-x)*(order-y)*c[1][d]+x*(order-y)*c[2][d]+
                         x*y*c[3][d]+(order-x)*y*c[4][d])/(order*order)),3)
        return ntuple(d->c[5][d]+t*(base[d]-c[5][d]),3)
    end
    error("unknown family $family")
end
_ho_mono(fam,p)=fam==:tri ? E._tri_monomials(p) : fam==:qua ? E._qua_monomials(p) :
    fam==:tet ? E._tet_monomials(p) : fam==:hex ? E._hex_monomials(p) :
    fam==:pri ? E._pri_monomials(p) : E._pyr_monomials(p)
const _HO_TOP=Tessella.MeshEntityTopology
function _ho_volume!(p,lin,fam,order,c)
    monos=_ho_mono(fam,order);n=length(monos)
    tags=Vector{Int32}(undef,n)
    nc=length(c)
    for i in 1:nc;tags[i]=_ho_at!(p,c[i]);end
    edges=_HO_TOP._simplex_edge_patterns(lin);e=order-1
    for (ord,(a,b)) in enumerate(edges),i in 1:e
        pos=ntuple(d->(c[a][d]*(order-i)+c[b][d]*i)/order,3)
        tags[nc+(ord-1)*e+i]=_ho_at!(p,pos)
    end
    offset=nc+e*length(edges)
    for (ffam,patterns) in ((:tri,_HO_TOP._simplex_face_patterns(lin,3)),
                          (:qua,_HO_TOP._simplex_face_patterns(lin,4)))
        interior=ffam==:tri ? (order-1)*(order-2)÷2 : (order-1)^2
        nfc=ffam==:tri ? 3 : 4
        isempty(patterns) && continue
        fmonos=_ho_mono(ffam,order)
        for pat in patterns
            fc=ntuple(i->c[pat[i]],nfc)
            for i in 1:interior
                tags[offset+i]=_ho_at!(p,_ho_pos(ffam,order,fmonos[nfc+nfc*e+i],fc))
            end
            offset+=interior
        end
    end
    for i in offset+1:n
        tags[i]=_ho_at!(p,_ho_pos(fam,order,monos[i],c))
    end
    return tags
end
_ho_facet!(p,fam,order,fc)=
    Int32[_ho_at!(p,_ho_pos(fam,order,m,fc)) for m in _ho_mono(fam,order)]
_ho_line!(p,order,ca,cb)=
    Int32[_ho_at!(p,ntuple(d->(ca[d]*(order-m[1])+cb[d]*m[1])/order,3))
          for m in E._line_monomials(order)]
# UID in canonical facet order for a corner-position frame.
_ho_uid(p,pos)=UInt64(p.map[ntuple(d->round(Int64,pos[d]*1e9),3)])
_ho_facet_uids(p,fam,order,fc)=
    UInt64[_ho_uid(p,_ho_pos(fam,order,m,fc)) for m in _ho_mono(fam,order)]
function _ho_nodes!(p)
    tags=Int32[];xyz=Float64[]
    for (i,c) in enumerate(p.coords)
        push!(tags,Int32(i));append!(xyz,collect(c))
    end
    A.mesh.add_nodes(3,11,tags,xyz)
end
# Node positions behind an emitted UID connectivity, via the fixture map.
_ho_frame(p,uids)=ntuple(i->p.coords[Int(uids[i])],length(uids))

# Two Tet35 share the {2,3,4} Tri15; a closed declared Line5 Curve rims it.
function p4_tets_fixture!(order=4)
    A.finalize();A.initialize();A.model.add("p$(order) tets")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    p=_HONodes()
    pos=Dict(1=>(0.,0.,0.),2=>(1.,0.,0.),3=>(0.,1.,0.),4=>(0.,0.,1.),7=>(1.,1.,1.))
    c1=(pos[1],pos[2],pos[3],pos[4]);c2=(pos[2],pos[3],pos[4],pos[7])
    t1=_ho_volume!(p,4,:tet,order,c1);t2=_ho_volume!(p,4,:tet,order,c2)
    faces=_HO_TOP._TETRAHEDRON_FACES
    s41=hcat([_ho_facet!(p,:tri,order,ntuple(i->c1[f[i]],3)) for f in faces[1:3]]...)
    s42=hcat([_ho_facet!(p,:tri,order,ntuple(i->c2[f[i]],3)) for f in faces[2:4]]...)
    l46=hcat(_ho_line!(p,order,pos[2],pos[3]),_ho_line!(p,order,pos[3],pos[4]),
             _ho_line!(p,order,pos[4],pos[2]))
    A.model.add_discrete_entity(1,46,Int[])
    A.model.add_discrete_entity(2,41,[-46])
    A.model.add_discrete_entity(2,42,[46])
    A.model.add_discrete_entity(3,11,[41,42])
    _ho_nodes!(p)
    A.mesh.add_elements_by_type(46,E.msh_type(:lin,order),[801,802,803],vec(l46))
    A.mesh.add_elements_by_type(41,E.msh_type(:tri,order),[701,702,703],vec(s41))
    A.mesh.add_elements_by_type(42,E.msh_type(:tri,order),[704,705,706],vec(s42))
    A.mesh.add_elements_by_type(11,E.msh_type(:tet,order),[501,502],vcat(t1,t2))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return p
end
# Two Hex125 share the x=1 Quad25.
function p4_hexes_fixture!()
    A.finalize();A.initialize();A.model.add("p4 hexes")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    p=_HONodes()
    pos=Dict(1=>(0.,0.,0.),2=>(1.,0.,0.),3=>(1.,1.,0.),4=>(0.,1.,0.),
             5=>(0.,0.,1.),6=>(1.,0.,1.),7=>(1.,1.,1.),8=>(0.,1.,1.),
             21=>(2.,0.,0.),22=>(2.,1.,0.),23=>(2.,0.,1.),24=>(2.,1.,1.))
    h1=(1,2,3,4,5,6,7,8);h2=(2,21,22,3,6,23,24,7)
    c1=ntuple(i->pos[h1[i]],8);c2=ntuple(i->pos[h2[i]],8)
    x1=_ho_volume!(p,5,:hex,4,c1);x2=_ho_volume!(p,5,:hex,4,c2)
    faces=_HO_TOP._HEXAHEDRON_FACES
    s41=hcat([_ho_facet!(p,:qua,4,ntuple(i->c1[f[i]],4)) for f in faces[[1,2,3,5,6]]]...)
    s42=hcat([_ho_facet!(p,:qua,4,ntuple(i->c2[f[i]],4)) for f in faces[[1,2,4,5,6]]]...)
    A.model.add_discrete_entity(2,41,Int[])
    A.model.add_discrete_entity(2,42,Int[])
    A.model.add_discrete_entity(3,11,[41,42])
    _ho_nodes!(p)
    A.mesh.add_elements_by_type(41,37,[701:705...],vec(s41))
    A.mesh.add_elements_by_type(42,37,[706:710...],vec(s42))
    A.mesh.add_elements_by_type(11,93,[501,502],vcat(x1,x2))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return p
end
# Two Prism75 share {1,2,4,5} Quad25; exteriors mix Tri15/Quad25.
function p4_prisms_fixture!()
    A.finalize();A.initialize();A.model.add("p4 prisms")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    p=_HONodes()
    pos=Dict(1=>(0.,0.,0.),2=>(1.,0.,0.),3=>(0.,1.,0.),
             4=>(0.,0.,1.),5=>(1.,0.,1.),6=>(0.,1.,1.),
             7=>(0.,-1.,0.),10=>(0.,-1.,1.))
    pr1=(1,2,3,4,5,6);pr2=(2,1,7,5,4,10)
    c1=ntuple(i->pos[pr1[i]],6);c2=ntuple(i->pos[pr2[i]],6)
    x1=_ho_volume!(p,6,:pri,4,c1);x2=_ho_volume!(p,6,:pri,4,c2)
    tris=_HO_TOP._PRISM_TRIANGLE_FACES;quads=_HO_TOP._PRISM_QUADRANGLE_FACES
    s41t=hcat(_ho_facet!(p,:tri,4,ntuple(i->c1[tris[1][i]],3)),
              _ho_facet!(p,:tri,4,ntuple(i->c1[tris[2][i]],3)))
    s41q=hcat(_ho_facet!(p,:qua,4,ntuple(i->c1[quads[2][i]],4)),
              _ho_facet!(p,:qua,4,ntuple(i->c1[quads[3][i]],4)))
    s42t=hcat(_ho_facet!(p,:tri,4,ntuple(i->c2[tris[1][i]],3)),
              _ho_facet!(p,:tri,4,ntuple(i->c2[tris[2][i]],3)))
    s42q=hcat(_ho_facet!(p,:qua,4,ntuple(i->c2[quads[2][i]],4)),
              _ho_facet!(p,:qua,4,ntuple(i->c2[quads[3][i]],4)))
    A.model.add_discrete_entity(2,41,Int[])
    A.model.add_discrete_entity(2,42,Int[])
    A.model.add_discrete_entity(3,11,[41,42])
    _ho_nodes!(p)
    A.mesh.add_elements_by_type(41,23,[701,702],vec(s41t))
    A.mesh.add_elements_by_type(41,37,[703,704],vec(s41q))
    A.mesh.add_elements_by_type(42,23,[705,706],vec(s42t))
    A.mesh.add_elements_by_type(42,37,[707,708],vec(s42q))
    A.mesh.add_elements_by_type(11,91,[501,502],vcat(x1,x2))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return p
end
# Two Pyr55 share the Quad25 base; the four side facets are Tri15.
function p4_pyrs_fixture!()
    A.finalize();A.initialize();A.model.add("p4 pyrs")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    p=_HONodes()
    pos=Dict(1=>(0.,0.,0.),2=>(1.,0.,0.),3=>(1.,1.,0.),4=>(0.,1.,0.),
             5=>(0.5,0.5,1.),6=>(0.5,0.5,-1.))
    py1=(1,2,3,4,5);py2=(4,3,2,1,6)
    c1=ntuple(i->pos[py1[i]],5);c2=ntuple(i->pos[py2[i]],5)
    x1=_ho_volume!(p,7,:pyr,4,c1);x2=_ho_volume!(p,7,:pyr,4,c2)
    tris=_HO_TOP._PYRAMID_TRIANGLE_FACES
    s41=hcat([_ho_facet!(p,:tri,4,ntuple(i->c1[f[i]],3)) for f in tris]...)
    s42=hcat([_ho_facet!(p,:tri,4,ntuple(i->c2[f[i]],3)) for f in tris]...)
    A.model.add_discrete_entity(2,41,Int[])
    A.model.add_discrete_entity(2,42,Int[])
    A.model.add_discrete_entity(3,11,[41,42])
    _ho_nodes!(p)
    A.mesh.add_elements_by_type(41,23,[701:704...],vec(s41))
    A.mesh.add_elements_by_type(42,23,[705:708...],vec(s42))
    A.mesh.add_elements_by_type(11,119,[501,502],vcat(x1,x2))
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    return p
end

@testset "third-order cells partition declared BRep" begin
    try
        # Tet20 pair: a Tri10 interface whose whole border is the declared
        # closed Line4 Curve; ghosts and MSH4 round-trip keep P3 blocks.
        p=p3_tets_fixture!();original=A.mesh.get_elements()
        A.option("Mesh.PartitionCreateGhostCells",1)
        A.mesh.partition(2,[501,502,701:706...,801:803...],
            Int32[1,2,1,1,1,2,2,2,1,1,1])
        @test A.model.get_number_of_partitions()==2
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        @test A.model.get_partitions(2,interface)==[1,2]
        ids,nodes=family(2,interface,21)
        @test length(ids)==1 && size(nodes)==(10,1)
        @test sort(vec(nodes))==sort(UInt64[p.map[(:v,c)] for c in (2,3,4)] ∪
            [_p3_edge_uid(p,2,3,i) for i in (1,2)] ∪
            [_p3_edge_uid(p,3,4,i) for i in (1,2)] ∪
            [_p3_edge_uid(p,4,2,i) for i in (1,2)] ∪
            [UInt64(p.map[(:p,round.(ntuple(d->(p.cv[2][d]+p.cv[3][d]+p.cv[4][d])/3,3);digits=9))])])
        curve_child=child(1,46)
        @test Set(abs.(last.(A.model.get_boundary([(2,interface)],false,true,false))))==
            Set([curve_child])
        @test A.model.get_partitions(1,curve_child)==[1,2]
        ids,nodes=family(1,curve_child,26)
        @test length(ids)==3 && size(nodes)==(4,3)
        @test all(sort(vec(nodes[3:4,c]))==sort(UInt64[_p3_edge_uid(p,a,b,i) for i in (1,2)])
            for (c,(a,b)) in enumerate(((2,3),(3,4),(4,2))))
        region1=only(tag for (dim,tag) in A.model.get_entities()
            if dim==3 && A.model.get_parent(3,tag)==(3,11) &&
               A.model.get_partitions(3,tag)==[1])
        region2=only(tag for (dim,tag) in A.model.get_entities()
            if dim==3 && A.model.get_parent(3,tag)==(3,11) &&
               A.model.get_partitions(3,tag)==[2])
        @test Set(abs.(last.(A.model.get_boundary([(3,region1)],false,true,false))))==
            Set([child(2,41),interface])
        @test Set(abs.(last.(A.model.get_boundary([(3,region2)],false,true,false))))==
            Set([child(2,42),interface])
        ids,nodes=family(3,region1,29);@test ids==UInt64[501] && size(nodes)==(20,1)
        ids,nodes=family(3,region2,29);@test ids==UInt64[502] && size(nodes)==(20,1)
        ghost1=only(tag for (tag,part) in A.mesh.get().partition_data.ghost_entities if part==1)
        ghost2=only(tag for (tag,part) in A.mesh.get().partition_data.ghost_entities if part==2)
        @test A.mesh.get_ghost_elements(3,ghost1)==(UInt64[502],Int32[2])
        @test A.mesh.get_ghost_elements(3,ghost2)==(UInt64[501],Int32[1])
        owned=A.mesh.get()
        @test Tessella.validate(owned).ok
        mktempdir() do directory
            for binary in (false,true)
                path=joinpath(directory,"p3-$(binary).msh")
                E.write_mixed_msh(path,owned;version=4.1,binary)
                imported=E.read_mixed_msh(path)
                @test records(imported)==records(owned)
                @test partition_records(imported)==partition_records(owned)
            end
            A.option("Mesh.MshFileVersion",4.1);A.option("Mesh.Binary",0)
            A.option("Mesh.SaveAll",1)
            path=joinpath(directory,"public.msh");A.write(path)
            A.model.add("p3 roundtrip");A.open(path)
            @test A.model.get_partitions(1,curve_child)==[1,2]
            A.model.set_current("p3 tets")
        end
        A.mesh.unpartition()
        @test A.mesh.get_elements()==original
        @test Set(A.model.get_entities())==
            Set([(1,46),(2,41),(2,42),(3,11)])
        # Hex64 pair: one Quad16 interface and a closed four-cell Line4 seam.
        p3_hexes_fixture!();original=A.mesh.get_elements()
        A.mesh.partition(2,[501,502,701:710...],Int32[1,2,1,1,1,1,1,2,2,2,2,2])
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        ids,nodes=family(2,interface,36)
        @test length(ids)==1 && size(nodes)==(16,1)
        interface_curves=Set(abs.(last.(A.model.get_boundary([(2,interface)],false,true,false))))
        @test length(interface_curves)==1
        seam=only(interface_curves)
        @test A.model.get_partitions(1,seam)==[1,2]
        ids,nodes=family(1,seam,26)
        @test length(ids)==4 && size(nodes)==(4,4)
        @test Tessella.validate(A.mesh.get()).ok
        A.mesh.unpartition()
        @test A.mesh.get_elements()==original
        # Prism40 pair: mixed Tri10/Quad16 exteriors around a Quad16
        # interface; Pyr30 pair: the Quad16 base interface plus Tri10 sides.
        for (fixture,tags) in ((p3_prisms_fixture!,701:708),
                               (p3_pyrs_fixture!,701:708))
            fixture();original=A.mesh.get_elements()
            A.mesh.partition(2,[501,502,collect(tags)...],
                Int32[1,2,1,1,1,1,2,2,2,2])
            interface=only(tag for (dim,tag) in A.model.get_entities()
                if dim==2 && A.model.get_parent(2,tag)==(3,11))
            ids,nodes=family(2,interface,36)
            @test length(ids)==1 && size(nodes)==(16,1)
            owned=A.mesh.get()
            @test Tessella.validate(owned).ok
            A.mesh.unpartition()
            @test A.mesh.get_elements()==original
        end
    finally
        A.finalize()
    end
end

@testset "third-order topology violations reject atomically" begin
    try
        # Serendipity Quad12, Tri12 and Tet22 blocks reject at the admission
        # gate even when they ride on a valid Tet20 mesh.
        for (badmsh,dim) in ((39,2),(22,2),(32,3))
            ncell=Int(E.msh_num_nodes(badmsh))
            p=p3_tets_fixture!(extras=()->begin
                base=Int32(A.mesh.get_max_node_tag())
                if base<ncell
                    xyz=Float64[]
                    for tag in base+1:Int32(ncell);append!(xyz,(0.,0.,10+tag));end
                    A.mesh.add_nodes(3,11,collect(base+1:Int32(ncell)),xyz)
                end
                A.mesh.add_elements_by_type(dim==3 ? 11 : 41,badmsh,[900],
                    Int32[collect(1:ncell)...])
            end)
            before=snapshot()
            caught=try
                A.mesh.partition(2,[501,502,701:706...,801:803...,900],
                    Int32[1,2,1,1,1,2,2,2,1,1,1,1]);nothing
            catch error;error end
            @test caught isa ArgumentError &&
                occursin("complete cell families only",sprint(showerror,caught))
            @test snapshot()==before
        end
    finally
        A.finalize()
    end
end

@testset "order >=4 cells partition declared BRep" begin
    try
        # Tet35 pair: a Tri15 interface rimmed by the declared closed Line5
        # Curve; the emitted interface connectivity is the facet's canonical
        # lattice order under whichever carrier orientation was chosen.
        p=p4_tets_fixture!();original=A.mesh.get_elements()
        A.option("Mesh.PartitionCreateGhostCells",1)
        A.mesh.partition(2,[501,502,701:706...,801:803...],
            Int32[1,2,1,1,1,2,2,2,1,1,1])
        @test A.model.get_number_of_partitions()==2
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        @test A.model.get_partitions(2,interface)==[1,2]
        ids,iface_nodes=family(2,interface,23)
        @test length(ids)==1 && size(iface_nodes)==(15,1)
        @test vec(iface_nodes[:,1])==
            _ho_facet_uids(p,:tri,4,_ho_frame(p,vec(iface_nodes[1:3,1])))
        curve_child=child(1,46)
        @test Set(abs.(last.(A.model.get_boundary([(2,interface)],false,true,false))))==
            Set([curve_child])
        @test A.model.get_partitions(1,curve_child)==[1,2]
        ids,nodes=family(1,curve_child,27)
        @test length(ids)==3 && size(nodes)==(5,3)
        rim=_ho_facet_uids(p,:tri,4,_ho_frame(p,vec(iface_nodes[1:3,1])))[1:12]
        @test Set(vec(nodes))==Set(rim)
        region1=only(tag for (dim,tag) in A.model.get_entities()
            if dim==3 && A.model.get_parent(3,tag)==(3,11) &&
               A.model.get_partitions(3,tag)==[1])
        region2=only(tag for (dim,tag) in A.model.get_entities()
            if dim==3 && A.model.get_parent(3,tag)==(3,11) &&
               A.model.get_partitions(3,tag)==[2])
        @test Set(abs.(last.(A.model.get_boundary([(3,region1)],false,true,false))))==
            Set([child(2,41),interface])
        @test Set(abs.(last.(A.model.get_boundary([(3,region2)],false,true,false))))==
            Set([child(2,42),interface])
        ids,nodes=family(3,region1,30);@test ids==UInt64[501] && size(nodes)==(35,1)
        ids,nodes=family(3,region2,30);@test ids==UInt64[502] && size(nodes)==(35,1)
        ghost1=only(tag for (tag,part) in A.mesh.get().partition_data.ghost_entities if part==1)
        ghost2=only(tag for (tag,part) in A.mesh.get().partition_data.ghost_entities if part==2)
        @test A.mesh.get_ghost_elements(3,ghost1)==(UInt64[502],Int32[2])
        @test A.mesh.get_ghost_elements(3,ghost2)==(UInt64[501],Int32[1])
        owned=A.mesh.get()
        @test Tessella.validate(owned).ok
        mktempdir() do directory
            for binary in (false,true)
                path=joinpath(directory,"p4-$(binary).msh")
                E.write_mixed_msh(path,owned;version=4.1,binary)
                imported=E.read_mixed_msh(path)
                @test records(imported)==records(owned)
                @test partition_records(imported)==partition_records(owned)
            end
            A.option("Mesh.MshFileVersion",4.1);A.option("Mesh.Binary",0)
            A.option("Mesh.SaveAll",1)
            path=joinpath(directory,"public.msh");A.write(path)
            A.model.add("p4 roundtrip");A.open(path)
            @test A.model.get_partitions(1,curve_child)==[1,2]
            A.model.set_current("p4 tets")
        end
        A.mesh.unpartition()
        @test A.mesh.get_elements()==original
        @test Set(A.model.get_entities())==
            Set([(1,46),(2,41),(2,42),(3,11)])
        # Hex125 pair: one Quad25 interface and a closed four-cell Line5 seam.
        p=p4_hexes_fixture!();original=A.mesh.get_elements()
        A.mesh.partition(2,[501,502,701:710...],
            Int32[1,2,1,1,1,1,1,2,2,2,2,2])
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        ids,nodes=family(2,interface,37)
        @test length(ids)==1 && size(nodes)==(25,1)
        @test vec(nodes[:,1])==_ho_facet_uids(p,:qua,4,_ho_frame(p,vec(nodes[1:4,1])))
        interface_curves=Set(abs.(last.(A.model.get_boundary([(2,interface)],false,true,false))))
        seam=only(interface_curves)
        @test A.model.get_partitions(1,seam)==[1,2]
        ids,nodes=family(1,seam,27)
        @test length(ids)==4 && size(nodes)==(5,4)
        @test Tessella.validate(A.mesh.get()).ok
        A.mesh.unpartition()
        @test A.mesh.get_elements()==original
        # Prism75 pair: mixed Tri15/Quad25 exteriors around a Quad25
        # interface; Pyr55 pair: the Quad25 base interface plus Tri15 sides.
        for (fixture,tags) in ((p4_prisms_fixture!,701:708),
                               (p4_pyrs_fixture!,701:708))
            fixture();original=A.mesh.get_elements()
            A.mesh.partition(2,[501,502,collect(tags)...],
                Int32[1,2,1,1,1,1,2,2,2,2])
            interface=only(tag for (dim,tag) in A.model.get_entities()
                if dim==2 && A.model.get_parent(2,tag)==(3,11))
            ids,nodes=family(2,interface,37)
            @test length(ids)==1 && size(nodes)==(25,1)
            owned=A.mesh.get()
            @test Tessella.validate(owned).ok
            A.mesh.unpartition()
            @test A.mesh.get_elements()==original
        end
        # Order 5 travels the same path: Tet56 pair, Tri21 interface, Line6 seam.
        p=p4_tets_fixture!(5);original=A.mesh.get_elements()
        A.mesh.partition(2,[501,502,701:706...,801:803...],
            Int32[1,2,1,1,1,2,2,2,1,1,1])
        interface=only(tag for (dim,tag) in A.model.get_entities()
            if dim==2 && A.model.get_parent(2,tag)==(3,11))
        ids,nodes=family(2,interface,25)
        @test length(ids)==1 && size(nodes)==(21,1)
        @test vec(nodes[:,1])==_ho_facet_uids(p,:tri,5,_ho_frame(p,vec(nodes[1:3,1])))
        curve_child=child(1,46)
        ids,nodes=family(1,curve_child,28)
        @test length(ids)==3 && size(nodes)==(6,3)
        @test Tessella.validate(A.mesh.get()).ok
        A.mesh.unpartition()
        @test A.mesh.get_elements()==original
    finally
        A.finalize()
    end
end

# Discrete prism: slave Surface 32 is the translated master Surface 31. All six
# boundary vertices carry Surface UV parameters so the fixture exercises both
# periodic-link and parametric-node transport at once.
function periodic_prism_fixture!()
    A.finalize();A.initialize();A.model.add("periodic prism")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    for tag in 31:35
        A.model.add_discrete_entity(2,tag)
    end
    A.model.add_discrete_entity(3,11,[31,32,33,34,35])
    A.mesh.add_nodes(2,31,[11,12,13],[0.,0.,0., 1.,0.,0., 0.,1.,0.],
        [0.,0., 1.,0., 0.,1.])
    A.mesh.add_nodes(2,32,[14,15,16],[0.,0.,1., 1.,0.,1., 0.,1.,1.],
        [0.,0., 1.,0., 0.,1.])
    A.mesh.add_elements_by_type(31,2,[631],[11,13,12])
    A.mesh.add_elements_by_type(32,2,[632],[14,15,16])
    A.mesh.add_elements_by_type(33,2,[633,634],[11,12,15, 11,15,14])
    A.mesh.add_elements_by_type(34,2,[635,636],[13,11,16, 11,14,16])
    A.mesh.add_elements_by_type(35,2,[637,638],[12,13,16, 12,16,15])
    A.mesh.add_elements_by_type(11,4,[501,502,503],
        [11,12,13,16, 11,12,16,15, 11,14,15,16])
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
end

@testset "periodic links and parametric nodes survive lower BRep partition" begin
    try
        periodic_prism_fixture!()
        # MSH round-trip imports an actual `$Periodic` link slave 32 -> master 31.
        path=joinpath(mktempdir(),"prism.msh");A.write(path)
        section=["\$Periodic","1","2 32 31","16",
            "1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0 -1.0 0.0 0.0 0.0 1.0",
            "3","14 11","15 12","16 13","\$EndPeriodic"]
        lines=readlines(path);append!(lines,section);write(path,join(lines,"\n")*"\n")
        A.model.add("imported periodic prism");A.open(path);rm(path)
        cached=A._cached_mesh_locked("test")
        @test length(cached.periodic_links)==1
        link=only(cached.periodic_links)
        @test (link.dim,link.slave_entity,link.master_entity)==(2,32,31)
        slave_nodes=copy(link.slave_nodes);master_nodes=copy(link.master_nodes)
        @test length(slave_nodes)==3 && link.affine!==nothing
        @test count(!isnothing,cached.entity_data.node_parametric)==6
        # Surface cells follow their carrier Tet4; Tet 502 gets part 2.
        A.mesh.partition(2,[631:638...,501,502,503],
            Int32[1,1,2,1,1,1,1,2,1,2,1])
        cached=A._cached_mesh_locked("test")
        # Periodic links stay bound to the original entities (Gmsh convention):
        # the slave/master parents persist alongside their partition children.
        @test length(cached.periodic_links)==1
        link=only(cached.periodic_links)
        @test (link.dim,link.slave_entity,link.master_entity)==(2,32,31)
        @test link.slave_nodes==slave_nodes && link.master_nodes==master_nodes
        @test link.affine!==nothing
        @test cached.partition_data.num_partitions==2
        @test Tessella.validate(cached).ok
        # Parametric rows keep their UVs and reclassify to same-parent children.
        data=cached.entity_data
        @test count(!isnothing,data.node_parametric)==6
        for row in 1:6
            parameters=data.node_parametric[row]
            @test parameters!==nothing && length(parameters)==data.node_entities[row][1]
            entity=data.node_entities[row]
            @test A.model.get_parent(entity[1],entity[2])==(2,row<=3 ? 31 : 32)
        end
        # `$Periodic` round-trips through the partitioned file.
        part=joinpath(mktempdir(),"prism-partitioned.msh");A.write(part)
        @test occursin("\$Periodic",read(part,String))
        A.model.add("partitioned periodic reread");A.open(part);rm(part)
        reread=A._cached_mesh_locked("test")
        @test length(reread.periodic_links)==1 &&
            reread.partition_data.num_partitions==2 &&
            Tessella.validate(reread).ok
        # Unpartition rebuilds the link against surviving node rows.
        A.mesh.unpartition()
        restored=A._cached_mesh_locked("test")
        @test restored.partition_data===nothing
        @test length(restored.periodic_links)==1 &&
            length(only(restored.periodic_links).slave_nodes)==3
        @test count(!isnothing,restored.entity_data.node_parametric)==6
    finally
        A.finalize()
    end
end

# `General.NumThreads` requests parallel workers for the O(cells) extraction
# stages in both partition producers; `Threads.@threads` caps actual
# concurrency at the runtime pool. Every parallel stage writes disjoint
# per-index slots and merges serially in traversal order, so any worker count
# must produce the identical partitioned mesh.
function bare_tet_fixture!()
    A.finalize();A.initialize();A.model.add("bare tets")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    A.model.add_discrete_entity(3,11);A.model.add_discrete_entity(3,22)
    A.mesh.add_nodes(3,11,[11,12,13,14,15],
        [0.,0.,0., 1.,0.,0., 0.,1.,0., 0.,0.,1., 0.,0.,-1.])
    A.mesh.add_elements_by_type(11,4,[501],[11,12,13,14])
    A.mesh.add_elements_by_type(22,4,[502],[11,13,12,15])
    A.mesh.generate(3)
    A.option("Mesh.PartitionCreateTopology",1)
    A.option("Mesh.PartitionCreateGhostCells",1)
end

@testset "parallel partition production stays deterministic" begin
    try
        last_line=fixture!()
        A.option("Mesh.PartitionCreateGhostCells",1)
        A.option("General.NumThreads",1)
        partition!(last_line)
        serial=A._cached_mesh_locked("test")
        @test Tessella.validate(serial).ok
        serial_records=partition_records(serial)
        serial_path=joinpath(mktempdir(),"serial.msh");A.write(serial_path)
        serial_bytes=read(serial_path);rm(serial_path)
        last_line=fixture!()
        A.option("Mesh.PartitionCreateGhostCells",1)
        A.option("General.NumThreads",3)
        partition!(last_line)
        threaded=A._cached_mesh_locked("test")
        @test Tessella.validate(threaded).ok
        @test partition_records(threaded)==serial_records
        threaded_path=joinpath(mktempdir(),"threaded.msh");A.write(threaded_path)
        @test read(threaded_path)==serial_bytes;rm(threaded_path)
        @test A.model.get_number_of_partitions()==2 &&
            !isempty(threaded.partition_data.ghost_elements)
        # The bare Tet4 producer follows the same determinism rule.
        bare_tet_fixture!()
        A.option("General.NumThreads",1)
        A.mesh.partition(2,[501,502],Int32[1,2])
        bare_one=joinpath(mktempdir(),"bare-serial.msh");A.write(bare_one)
        one=A._cached_mesh_locked("test")
        @test Tessella.validate(one).ok
        bare_bytes=read(bare_one);rm(bare_one)
        bare_tet_fixture!()
        A.option("General.NumThreads",4)
        A.mesh.partition(2,[501,502],Int32[1,2])
        bare_two=joinpath(mktempdir(),"bare-threaded.msh");A.write(bare_two)
        two=A._cached_mesh_locked("test")
        @test Tessella.validate(two).ok
        @test read(bare_two)==bare_bytes;rm(bare_two)
        @test two.partition_data.num_partitions==2 &&
            length(two.partition_data.ghost_elements)==2
    finally
        A.option("General.NumThreads",1)
        A.finalize()
    end
end
end
