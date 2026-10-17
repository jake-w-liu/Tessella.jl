# `model.occ.defeature` — `BRepAlgoAPI_Defeature`: listed faces are
# removed and the kept faces' planes refill the gap (exact when the
# healed body is the convex half-space intersection).
using Test
using Tessella
const A=Tessella.API

function _errtext(f)
    try
        f()
    catch err
        return sprint(showerror,err)
    end
    return ""
end

@testset "occ.defeature" begin
    A.initialize(); occ=A.model.occ
    try
        # boss-on-box: unit box + 0.4³ boss fused on top
        b=occ.add_box(0.,0.,0.,1.,1.,1.)
        boss=occ.add_box(0.3,0.3,1.0,0.4,0.4,0.4)
        fused,_=occ.fuse([(3,b)],[(3,boss)])
        v=only(t for (d,t) in fused if d==3)
        @test occ.get_mass(3,v)≈1.0+0.4^3 atol=1e-6
        # boss faces: centroid z > 1
        boss_faces=Int[]
        # faces of the fused solid — enumerate via kernel-independent
        # route: every (2,f) whose mass lies on z>1 planes
        for (_,f) in occ.get_entities(2)
            _,_,z=occ.get_center_of_mass(2,f)
            z>1.0001 && push!(boss_faces,f)
        end
        @test !isempty(boss_faces)
        nvol=length(occ.get_entities(3))
        occ.defeature([(3,v)],boss_faces)
        @test length(occ.get_entities(3))==nvol       # in-place
        @test occ.get_mass(3,v)≈1.0 atol=1e-6       # boss healed → box

        # errors stay atomic
        m1=occ.get_mass(3,v)
        @test occursin("unknown",_errtext(()->occ.defeature([(3,99999)],[1])))
        @test occursin("not a face",_errtext(()->occ.defeature([(3,v)],[99999])))
        @test occursin("no faces",_errtext(()->occ.defeature([(3,v)],Int[])))
        @test occ.get_mass(3,v)≈m1
    finally
        A.finalize()
    end
end
