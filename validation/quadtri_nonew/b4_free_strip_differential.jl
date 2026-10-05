# Included by differential.jl after its existing saved-capture/support helpers.
# The original raw input and actual primary identities are authoritative.
const QTB4F=QuadTriNoNewB4FreeStripCertificates

function b4_free_strip_saved_case(record;oracle_only=false)
    raw=record["payload"];f=QTB4F.saved_fixture(record);tolerance=2e-11
    raw["input_sha256"]==record["input_literal_sha256"]==record["recipe_sha256"] || error("free B4 literal recipe provenance")
    isempty(raw["warnings_or_errors"]) && !any(occursin("Error",line) for line in raw["gmsh_log"]) || error("free B4 captured diagnostic provenance")
    # Existing exact payload/support validator needs its neutral audit field;
    # no capture field, original public tag, or coordinate bit is changed.
    validated=merge(record,Dict("audit"=>record["retained_capture_audit"]))
    saved_nodes=rect_grid_saved_integrity(validated)
    source,_=two_tri_saved_mesh(raw,"p1",2,1)
    linear,index=two_tri_saved_mesh(raw,"p1",3,1)
    saved=QTB4F.certify(linear,f,source;oracle=true)
    isempty(saved.centers) || error("free B4 original primary capture center observation changed")
    saved.families==Dict(parse(Int,k)=>Int(v) for (k,v) in record["retained_capture_audit"]["families"]) || error("free B4 actual saved family evidence differs")
    two_tri_saved_boundary(raw,index)==Set(keys(saved.boundary)) || error("free B4 saved typed exterior")
    quadratic=two_tri_saved_mesh(raw,"p2",3,1)[1]
    qsaved=QTB4F.certify_quadratic(quadratic,f,source;oracle=true)
    qsaved.ncell==saved.ncell || error("free B4 saved actual P2 cell identity")
    empty_parameters=count(node->node["owner"][1]==2 && isempty(node["stored_parameters"]),values(saved_nodes))
    empty_parameters==record["retained_capture_audit"]["empty_stored_uv"] || error("free B4 stored UV provenance")
    surface_queries=rect_grid_saved_surface_queries(validated)
    empty_parameters==sum(Int(row["empty_stored_uv"]) for row in surface_queries) || error("free B4 saved surface query coverage")
    if !oracle_only
        execution=QTB4F.execute(f.source)
        volume=geo_entity_mesh(execution,3,1);actual_source=geo_entity_mesh(execution,2,1)
        quad_patch_source_match(actual_source,source,tolerance)
        actual=QTB4F.certify(volume,f,actual_source)
        QTNN.surface_faces(execution.mesh_parts,volume)==Set(keys(actual.boundary)) || error("free B4 native typed exterior")
        projection=model_to_mixed(execution.model,volume,3,1)
        nodes=Dict(Int(node["tag"])=>node for node in raw["p1"]["nodes"])
        # Twelve source recipes happen to emit C=0 in each native run. Cell
        # families and pointer-dependent diagonals are independently certified
        # and may differ from the captured primary choice.
        isempty(actual.centers) || error("free B4 native recipe acquired an unexpected center")
        for node in 1:nnodes(volume)
            matches=[tag for tag in keys(index) if maximum(abs.(QTB4F.point(volume,node).-Tuple(Float64.(nodes[tag]["coordinates"]))))<=tolerance]
            length(matches)==1 || error("free B4 actual primary column geometry bijection")
            Tuple(Int.(nodes[only(matches)]["owner"]))==projection.entity_data.node_entities[node] || error("free B4 actual primary carrier identity")
        end
        for chain in raw["source_chains"]
            params=execution.model.curve_params[Int(chain["curve"])];expected=Float64.(chain["stored_owned_parameters"])
            length(params)==length(expected)+2 && maximum(abs.(sort(params[2:end-1]).-sort(expected));init=0.)<=tolerance || error("free B4 sampled native Curve parameters")
        end
        mktempdir() do directory
            path=joinpath(directory,"b4_free_strip_saved.geo");write(path,f.source)
            native_api(path,api->begin
                cache=api.mesh.generate(3);first=QTB4F.public_volume(api,1)
                QTNN.typed_signature(first)==QTNN.typed_signature(volume) || error("free B4 actual API/GEO primary cells")
                data=api.mesh.get_elements(3,1);projected=model_to_mixed(api.CURRENT[],cache,3,1)
                isempty(api.mesh.get_elements(2)[1]) || error("free B4 API3 lower-cell contract")
                api.mesh.set_order(2);p2=QTB4F.public_volume(api,1)
                certificate=QTB4F.certify_quadratic(p2,f,actual_source)
                abs(certificate.total-actual.total)<=QTB4F.Q(tolerance) || error("free B4 actual P2 whole-map integral")
                strip_native_supports(api,cache,data,projected,1,tolerance)==nnodes(p2)-nnodes(first) || error("free B4 complete actual P2 supports")
                for node in api.mesh.get_nodes()[1]
                    p,uv,dim,owner=api.mesh.get_node(node)
                    if dim in (1,2)
                        length(uv)==dim && maximum(abs.(api.model.get_value(dim,owner,uv).-p))<=tolerance || error("free B4 computed carrier parameters")
                    end
                end
                for row in surface_queries,include_boundary in (false,true)
                    ids,xyz,uv=api.mesh.get_nodes(2,Int(row["tag"]),include_boundary,true)
                    length(ids)==row[include_boundary ? "closure" : "owned"] && length(uv)==2length(ids) || error("free B4 actual surface owned/closure query")
                    maximum(abs.(api.model.get_value(2,Int(row["tag"]),uv).-xyz))<=tolerance || error("free B4 computed surface UV")
                end
                api.mesh.set_order(1)
                QTNN.typed_signature(QTB4F.public_volume(api,1))==QTNN.typed_signature(first) || error("free B4 P2/P1 actual cell roundtrip")
            end)
        end
    end
    return empty_parameters
end

function b4_free_strip_saved_catalog(;selected="",oracle_only=false)
    catalog=TOML.parsefile(joinpath(@__DIR__,"..","..","test","artifacts","quadtri_nonew_b4_free_strip_oracle.toml"))
    catalog["gmsh_version"]=="4.15.2" && catalog["coordinate_abs_tolerance"]==2e-11 && catalog["fixture_count"]==12 || error("free B4 original primary catalog contract")
    cases=0;empty_parameters=0
    for record in catalog["fixtures"]
        name=record["payload"]["name"]
        isempty(selected) || selected=="b4_free_strip" || occursin(selected,"b4_free_strip_"*name) || continue
        empty_parameters+=b4_free_strip_saved_case(record;oracle_only);cases+=1
        println("QUADTRI_NONEW_B4_FREE_STRIP_SAVED_OK name=$name observed_centers=0 families_independently_certified=true")
    end
    println("QUADTRI_NONEW_B4_FREE_STRIP_DIFFERENTIAL_OK saved_cases=$cases p2_cases=$cases empty_stored_uv=$empty_parameters")
    return (;cases,empty_parameters)
end
