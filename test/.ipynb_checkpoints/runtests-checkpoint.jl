using Test
using HelperDataFunc
using Cosmology
using Unitful, UnitfulAstro


let
	id_list::Vector{Int} = [1, 2, 3, 4, 5]
	z_list::Vector{Float64} = [0.1, 0.2, 0.3, 1., 1.]
	logmstar_list::Vector{Float64} = [8.0, 8.3, 9.5, 10.1, 11.2]
	test_catalog = Catalog(id_list, z_list, logmstar_list)
	z_edges::Vector{Float64} = [0., 1.5]
	logm_edges::Vector{Array{Float64}} = [[7.5, 12.0]]
	comoving_volumes = ustrip(comoving_volume(u"Mpc^3", cosmology(), 1.5))
	num_density = 5 ./ comoving_volumes ./ 4.5

	temp_result = calc_log_gal_density(test_catalog, z_edges, logm_edges, 41253., 4.5)

	@test temp_result[1].logϕ_obs ≈ [log10.(num_density)] atol = 1e-6
	@test temp_result[1].σs ≈ [1. ./ (log.(10.) .* sqrt.(5)) .* log10.(1. ./ 4.5 ./ comoving_volumes)] atol = 1e-6
	
end