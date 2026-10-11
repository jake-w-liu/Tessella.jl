# Native homology coverage: element wraparound, signed cells, wide-node
# numbering and physical-name publication of computed chain cells.
# (Merged from api_native_homology_{element_wrap_v2,signed_cells,wide_nodes}_test.jl
# and api_homology_physical_names_test.jl.)

module APINativeHomologyElementWrapTestV2
using Tessella,Test
include("../../validation/native_geometry/native_homology_element_wrap_controls_v2.jl")
prefix=get(ENV,"TESSELLA_NATIVE_HOMOLOGY_ELEMENT_WRAP_OUTPUT_PREFIX",joinpath(mktempdir(),"native_homology_element_wrap"))
NativeHomologyElementWrapV2.run(prefix)
end

module APINativeHomologySignedCellsTests
using Tessella,Test
include("../../validation/native_geometry/native_homology_signed_cells_controls_v1.jl")
prefix=get(ENV,"TESSELLA_NATIVE_SIGNED_HOMOLOGY_OUTPUT_PREFIX",joinpath(mktempdir(),"signed_homology"))
NativeHomologySignedCellsV1.run(prefix*".functional.toml")
end

module APINativeHomologyWideNodesTestV2
using Tessella,Test
include("../../validation/native_geometry/native_homology_wide_nodes_controls_v2.jl")
prefix=get(ENV,"TESSELLA_NATIVE_WIDE_HOMOLOGY_OUTPUT_PREFIX",joinpath(mktempdir(),"wide_homology"))
NativeHomologyWideNodesV2.run(prefix*".functional.toml")
end

module ChainPhysicalNamesProbe
using Tessella,Test
const API=Tessella.API
@testset "Published chain cells retain their physical names in the actual mesh" begin
    API.initialize()
    try
        API.model.add_discrete_entity(1,1)
        API.mesh.add_nodes(1,1,[11,12,13,14],[0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.])
        API.mesh.add_elements_by_type(1,1,[101,102,103,104],[11,12,12,13,13,14,14,11])
        API.model.add_physical_group(1,[1];tag=10)
        API.mesh.reclassify_nodes()
        API.mesh.add_homology_request("Homology",[10],[],[1])
        @test API.mesh.compute_homology()==[(1,11)]
        @test API.model.get_physical_name(1,11)=="H_1{10}1"
        @test get(API.LAST_MESH[].physical_names,(1,11),nothing)=="H_1{10}1"
    finally
        API.finalize()
    end
end
end
