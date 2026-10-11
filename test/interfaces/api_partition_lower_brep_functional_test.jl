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
            occursin("repeats a Line2 support",sprint(showerror,caught))
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
end
