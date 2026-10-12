module OccClosestSmokeChecks
using Test, Tessella
let api=Tessella.API
    api.initialize()
    try
        occ=api.model.occ
        # `occ.get_closest_entities` — registry-side nearest-entity query
        # materialized into a scratch model; the interface test exercises
        # ordering/dedup/n semantics, wire/shell dims, trimmed faces,
        # curved records, transform visibility and error paths.
        b=occ.add_box(0.0,0.0,0.0, 2.0,3.0,4.0)
        s=occ.add_sphere(10.0,0.0,0.0, 1.5)
        tags,dist,coord=occ.get_closest_entities(
            4.0,1.0,2.0,[(3,b),(3,s)],2)
        @test tags==[(3,b),(3,s)]
        @test dist≈[2.0,sqrt(41)-1.5]
        @test coord[1:3]≈[2.0,1.0,2.0]
        tags,dist,coord=occ.get_closest_entities(1.0,1.0,1.0,[(3,b)])
        @test dist==[0.0] && coord≈[1.0,1.0,1.0]
        p=occ.add_point(2.0,0.0,0.0)
        tags,dist,_=occ.get_closest_entities(3.0,0.0,1.0,[(0,p),(3,b)],2)
        @test tags==[(3,b),(0,p)] && dist[1]≈1.0 && dist[2]≈sqrt(2)
        # `occ.get_distance` — the upstream (distance, x1…z2) layout.
        q=occ.add_point(6.0,0.0,0.0)
        d=occ.get_distance(0,q,3,b)
        @test d[1]≈4.0 && d[2:4]==(6.0,0.0,0.0) && d[5:7]==(2.0,0.0,0.0)
        @test occ.get_distance(3,b,3,b)[1]==0.0
        d=occ.get_distance(0,q,3,s)
        @test d[1]≈2.5 && collect(d[5:7])≈[8.5,0.0,0.0]
    finally
        api.finalize()
    end
end
end
