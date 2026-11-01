module NativeGraphPartitionChecks
# Functional checks for the native multilevel element-dual graph partitioner
# (the upstream METIS-quality tier, implemented without a METIS dependency):
# balanced weighted parts, lower edge cuts than the simple BFS tier,
# connected partitions, deterministic output, the Mesh.Metis* option
# dispatch, the -1 simple-tier extension and configured element weights.
using Test, Tessella
const A=Tessella.API

# Build a structured L×L×L node grid meshed into 6 tets per unit cube on a
# single discrete volume. The L-shaped domain (cubes only where i<=L/2 or
# j<=L/2) keeps the dual connected while defeating geometric slicing — the
# same family of fixture where the upstream METIS tier shows its advantage
# over `SimplePartition`'s chessboard cuts.
function _tet_grid!(L::Int)
    A.model.add_discrete_entity(3,1)
    nside=L+1
    node(i,j,k)=(k-1)*nside*nside+(j-1)*nside+i
    nodes=Int[];coords=Float64[]
    for k in 1:nside,j in 1:nside,i in 1:nside
        push!(nodes,node(i,j,k))
        append!(coords,[(i-1.0)/L,(j-1.0)/L,(k-1.0)/L])
    end
    A.mesh.add_nodes(3,1,nodes,coords)
    # Six tets around the d-e body diagonal of each cube.
    v(i,j,k)=Int32(node(i,j,k))
    conn=Int32[]
    for k in 1:L,j in 1:L,i in 1:L
        (i<=L÷2 || j<=L÷2) || continue      # L-shaped domain
        a=v(i,j,k);b=v(i+1,j,k);c=v(i,j+1,k);d=v(i+1,j+1,k)
        e=v(i,j,k+1);f=v(i+1,j,k+1);g=v(i,j+1,k+1);h=v(i+1,j+1,k+1)
        # Six positive-volume tets around the shared d-e body diagonal.
        append!(conn,(a,b,d,e, d,b,f,e, d,f,h,e, a,d,c,e, c,d,g,e, g,d,h,e))
    end
    tags=collect(Int32,1:Int32(length(conn)÷4))
    A.mesh.add_elements_by_type(1,4,tags,conn)
    A.mesh.generate(3)
    return nothing
end

function _edgecut(graph,assignment)
    cut=0
    for v in 1:length(graph.vwgt)
        for e in graph.xadj[v]:graph.xadj[v+1]-1
            u=graph.adjncy[e]
            u>v && assignment[u]!=assignment[v] && (cut+=1)
        end
    end
    return cut
end

_components(graph,assignment,p)=length(A._gp_components_in(
    graph,Int32[i for i in 1:length(assignment) if assignment[i]==p]))

# SimplePartition-equivalent baseline: contiguous bins by x-centroid.
function _geo_slice(mesh,np)
    cells=[b for b in mesh.blocks if b.msh==4][1].nodes
    cent=[sum(mesh.coords[1,cells[r,c]] for r in 1:4)/4
          for c in 1:size(cells,2)]
    o=sortperm(cent);asg=zeros(Int32,size(cells,2));q=cld(length(o),np)
    for (rank,e) in enumerate(o);asg[e]=min(np,cld(rank,q));end
    return asg
end

@testset "native multilevel partitioner" begin
    A.finalize();A.initialize();A.model.add("tet grid")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    try
        _tet_grid!(6)             # 972 tets on the L-domain
        mesh=A.mesh.get()
        graph,_=A._gp_build_dual(mesh)
        ntet=sum(size(b.nodes,2) for b in mesh.blocks if b.msh==4)
        @test length(graph.vwgt)==ntet
        for np in (2,3,4,6)
            asg=A._partition_elements(mesh,np)
            @test length(asg)==ntet
            @test sort(unique(asg))==collect(Int32,1:np)
            # Upstream balance: recursive bisection targets essentially equal
            # weighted parts (ufactor=1 per-mille); bounded tolerance here.
            counts=[count(==(p),asg) for p in 1:np]
            @test maximum(counts)-minimum(counts)<=max(2,ceil(Int,0.03*ntet))
            # Quality tier: better cut than the simple BFS tier and far
            # better than geometric slicing on the irregular domain.
            simple=A._partition_elements_simple(mesh,np)
            geo=_geo_slice(mesh,np)
            @test _edgecut(graph,asg)<_edgecut(graph,simple)
            @test _edgecut(graph,asg)<_edgecut(graph,geo)
            # Every partition connected on the dual graph.
            @test all(p->_components(graph,asg,p)==1,1:np)
            # Deterministic.
            @test A._partition_elements(mesh,np)==asg
        end
        # MetisAlgorithm=-1 selects the legacy simple tier exactly.
        A.option("Mesh.MetisAlgorithm",-1)
        @test A._partition_elements(mesh,4)==
              A._partition_elements_simple(mesh,4)
        # MetisAlgorithm=2 runs the multilevel k-way path.
        A.option("Mesh.MetisAlgorithm",2)
        kway=A._partition_elements(mesh,4)
        @test sort(unique(kway))==Int32[1,2,3,4]
        @test all(p->_components(graph,kway,p)==1,1:4)
        @test A._partition_elements(mesh,4)==kway
        A.option("Mesh.MetisAlgorithm",1)
        # Degenerate inputs keep the contiguous quota convention.
        @test A._partition_elements(mesh,1)==ones(Int32,ntet)
        small=A._partition_elements(mesh,ntet+2)
        @test maximum(small)<=ntet+2 && length(small)==ntet
    finally
        A.option("Mesh.MetisAlgorithm",1)
        A.finalize()
    end
end

@testset "configured element weights" begin
    A.finalize();A.initialize();A.model.add("weighted")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    try
        _tet_grid!(4)
        mesh=A.mesh.get()
        graph,families=A._gp_build_dual(mesh)
        # A 2:1 tet weight means part element counts skew ~2:1 while the
        # weighted balance stays tight.
        A.option("Mesh.PartitionTetWeight",2)
        A.option("Mesh.PartitionLineWeight",1)
        A.option("Mesh.PartitionTriWeight",1)
        A.option("Mesh.PartitionQuadWeight",1)
        A.option("Mesh.PartitionHexWeight",1)
        A.option("Mesh.PartitionPrismWeight",1)
        A.option("Mesh.PartitionPyramidWeight",1)
        asg=A._partition_elements(mesh,2)
        @test sort(unique(asg))==Int32[1,2]
        @test all(p->_components(graph,asg,p)==1,1:2)
        vwgt=A._gp_element_weights(families,(1,1,2,1,1,1,1,0))
        weights=[sum(Int(vwgt[i]) for i in 1:length(asg) if asg[i]==p)
                 for p in 1:2]
        @test abs(weights[1]-weights[2])<=2
        # Auto weighting (-1) weights only top-dimensional cells.
        A.option("Mesh.PartitionTetWeight",-1)
        auto=A._gp_element_weights(families,(-1,-1,-1,-1,-1,-1,-1,0))
        @test all(==(1),auto)         # pure-tet mesh: all top-dim
    finally
        for k in ("Mesh.PartitionTetWeight","Mesh.PartitionLineWeight",
                  "Mesh.PartitionTriWeight","Mesh.PartitionQuadWeight",
                  "Mesh.PartitionHexWeight","Mesh.PartitionPrismWeight",
                  "Mesh.PartitionPyramidWeight")
            A.option(k,-1)
        end
        A.finalize()
    end
end

@testset "partition API integration" begin
    A.finalize();A.initialize();A.model.add("api part")
    A.option("General.Terminal",0);A.option("Mesh.Renumber",0)
    try
        _tet_grid!(3)               # 90 tets on the L-domain
        source_ntet=sum(size(b.nodes,2) for b in A.mesh.get().blocks if b.msh==4)
        A.option("Mesh.PartitionCreateTopology",1)
        A.option("Mesh.PartitionCreateGhostCells",1)
        A.mesh.partition(4)
        mesh=A.mesh.get()
        @test mesh.partition_data!==nothing
        @test Int(mesh.partition_data.num_partitions)==4
        # Owned tets across partition children sum to the source count.
        owned=0
        for ((dim,tag),entity) in mesh.partition_data.entities
            dim==3 || continue
            _,etags,_=A.mesh.get_elements(3,tag)
            owned+=sum(length.(etags);init=0)
        end
        @test owned==source_ntet
        # Ghost records exist across the partition seams.
        @test !isempty(mesh.partition_data.ghost_elements)
        A.mesh.unpartition()
        # Mesh.NbPartitions drives partition(-1) like upstream.
        A.option("Mesh.NbPartitions",4)
        A.option("Mesh.PartitionCreateTopology",0)
        A.option("Mesh.PartitionCreateGhostCells",0)
        A.mesh.partition(-1)
        @test Int(A.mesh.get().partition_data.num_partitions)==4
    finally
        A.option("Mesh.NbPartitions",0)
        A.option("Mesh.PartitionCreateTopology",0)
        A.option("Mesh.PartitionCreateGhostCells",0)
        A.finalize()
    end
end
end
