module HelperDataFunc
export cosmos_mass_compl
export calc_redshift_volumes
export calc_double_schechter
export simulate_smf_from_model_bestfit
export read_gal_fits
export my_subcatalog_struct
export calc_mass_bins
export calc_log_gal_density

export Catalog, masterdata

using Unitful, UnitfulAstro
using Cosmology
using DataFrames
using FITSIO

begin 

    @doc """Calculate the mass completeness curve by Weaver+2022"""
    
    function cosmos_mass_compl(z::Float64)
        return log10.(-3.23e7 * (1 + z) + 7.83e7 * (1+z)^2.)
    end
    
    function cosmos_mass_compl(z::Vector)
        return log10.(-3.23e7 .* (1 .+ z) + 7.83e7 .* (1 .+ z)^2.)
    end

end

#######################################

begin
	

	function calc_redshift_volumes(field_area::Float64, z_edges::Vector{Float64}, cosmo::Cosmology.FlatLCDM{Float64} = cosmology())::Vector

    """Calculate the comoving volumes between two redshifts given the redshift bin edges and survey coverage area in deg^2	Input: field area (`Float64`), z edges (`Vector{Float64}`), `cosmology()`	Output: volumes (`Vector{Float64}`)"""
    	
    	volumes = Vector{Float64}(undef, 0)
    	for (z_low, z_high) in zip(z_edges[1:length(z_edges) - 1], z_edges[2:length(z_edges)])
    		@assert z_low <= z_high
    		volume =  field_area ./ 41253. .* (ustrip(comoving_volume(u"Mpc^3", cosmo, z_high)) .- ustrip(comoving_volume(u"Mpc^3", cosmo, z_low)))	
    		append!(volumes, volume)
    	end
    	@assert length(volumes) == length(z_edges) - 1
    	return volumes
	end
end

########################################

begin
	@doc """Calculate log galaxy number density from logM, logMc, logϕ1, logϕ2, α1, and α2"""
	function calc_double_schechter end
	calc_double_schechter(logM, p::Vector) = calc_double_schechter(logM, p...)
	function calc_double_schechter(logM, logMc, logϕ1, logϕ2, α1, α2)
    	m_over_mc = 10. .^ (logM .- logMc)
		return log10.(log(10.) .* exp.(.-m_over_mc) .* (10. .^ logϕ1 .* m_over_mc .^ (α1 .+ 1) + 10. ^ logϕ2 .* m_over_mc .^ (α2 .+ 1)))
	end
end

#####################################

begin

    @doc """Calculate the double Schechter function using the best-fit parameters returned by MCMC chains"""
    function simulate_smf_from_model_bestfit(chain, xfit; chain_id::Integer=1)
    
    	@assert 1<=chain_id<=size(chain,3)
    	# Extract parameters from chain
    	logMc = quantile(chain[:logMc][:,chain_id], 0.5)
    	logϕ1 = quantile(chain[:logϕ1][:,chain_id], 0.5)
    	logϕ2 = quantile(chain[:logϕ2][:,chain_id], 0.5)
    	α1 = quantile(chain[:α1][:,chain_id], 0.5)
    	α2 = quantile(chain[:α2][:,chain_id], 0.5)
    	σ_add = quantile(chain[:σ_add][:,chain_id], 0.5)
    	return calc_double_schechter.(xfit, logMc, logϕ1, logϕ2, α1, α2)
    end

end

##########################




begin
    @doc """Custom `Catalog` struct"""
	struct Catalog{T1<:Number, T2<:Number, T3<:Number,
					  V1<:AbstractVector{T1}, V2<:AbstractVector{T2}, V3<:AbstractVector{T3}}
		id::V1
		z::V2
		logmstar::V3
	end

#	function Catalog(id::V1, z::V2, logmstar::V3) where {
#		T1<:Number, T2<:Number, T3<:Number,
#		V1<:AbstractVector{T1}, V2<:AbstractVector{T2}, V3<:AbstractVector{T3} }
#		@assert length(id) == length(z) == length(logmstar)
#		Catalog(id, z, logmstar)
#	end

	function read_gal_fits(path::String)
		@assert isfile(path)
		interm_file = DataFrames.DataFrame(FITSIO.FITS(path, "r")[2])
		return Catalog(interm_file."id", interm_file."redshift", log10.(interm_file."bayes.stellar.m_star"))
	end

end	

begin
	@doc """My function to read substruct of my `Catalog` struct"""
	function my_subcatalog_struct(catalog::Catalog, idx)
		return Catalog(collect(view(catalog.id, idx)), collect(view(catalog.z, idx)), collect(view(catalog.logmstar, idx)))
	end
end

begin
    @doc """Construct `masterdata` struct with pre-allocated workspace."""
	struct masterdata{T1<:Number, T2<:Number, T3<:Number, T4<:Number, T5<:Number, T6<:Number,
		V3<:AbstractVector{T3}, V4<:AbstractVector{T4}, V5<:AbstractVector{T5}, V6<:AbstractVector{T6}}
		z_low::T1
		z_high::T2
		logmstar::V3
		logϕ_obs::V4
		σs::V5
		predict::V6
	
        function masterdata(z_low::T1, z_high::T2, logmstar::V3, logϕ_obs::V4, σs::V5) where {
    		T1<:Number, T2<:Number, T3<:Number, T4<:Number, T5<:Number,
    		V3<:AbstractVector{T3},V4<:AbstractVector{T4},V5<:AbstractVector{T5} }
    		@assert z_high > z_low >= 0.
    		@assert length(logϕ_obs) == length(logmstar) == length(σs)
    		WorkspaceT = promote_type(T1, T2, T3, T4, T5)
    		ws1 = Vector{WorkspaceT}(undef, length(logmstar))
    		new{Float64, Float64, Float64, Float64, Float64, Float64, Vector{Float64},Vector{Float64},Vector{Float64}, Vector{Float64}}(z_low, z_high, logmstar, logϕ_obs, σs, ws1)
    	end
    end
end

#######################

begin
    @doc """Calculate the logmass bin edges in each redshift bins, by considering the mass completeness curve in COSMOS by Weaver+2022. Note: The correct output: length(logmass) == length(z_center) == length(z_edges) - 1
	
	Input: redshift bin edges (`Vector{Float64}`).
	Output: logmstar edges """
	function calc_mass_bins(z_edges::Vector{Float64}, logmass_grid::Vector{Float64}; mass_compl_curve = cosmos_mass_compl)
		logmstar_edges = Vector{Array{Float64}}(undef, length(z_edges) - 1)
		n = length(z_edges)
		z_centers = (z_edges[1:n-1] .+ z_edges[2:n]) ./ 2
        for i in 1:length(z_centers)
    		z = z_centers[i]
			logmstar_edges_in_zbin = Vector{Float64}
			for i in 1:length(logmass_grid) - 1
				if logmass_grid[1] >= mass_compl_curve(z)
					logmstar_edges_in_zbin = logmass_grid
				elseif logmass_grid[i] - mass_compl_curve(z) < 0. && logmass_grid[i+1] - mass_compl_curve(z) >= 0.
					logmstar_edges_in_zbin = logmass_grid[i:length(logmass_grid)]
				end
			end
			logmstar_edges[i] = logmstar_edges_in_zbin
		end	
		@assert size(logmstar_edges,1) == length(z_edges) - 1
		return logmstar_edges	
	end
end

##########################

begin
	@doc """Calculate the log galaxy number density in each redshift bins. The function will invoke `calc_redshift_volumes()`.
	
	Input: the catalog (after applying the mass completeness curve), z edges (`Vector{Float64}`), logmstar_edges (`Vector{Any}`), and field area (`Float64`)
	Output: `masterdata` struct"""
	function calc_log_gal_density(catalog::Catalog, z_edges::Vector{Float64},  logmstar_edges::Vector{Array{Float64}}, field_area::Float64, logmass_step::Float64)
		@assert size(logmstar_edges, 1) == length(z_edges) - 1
		z_vols = calc_redshift_volumes(field_area, z_edges)
		log_density_gal = Array{masterdata}(undef, length(z_edges) - 1)
		for i in 1:length(z_vols)
			number_gal_in_logmbins_vec = Vector{Float64}(undef, 0)
			len_logmstar_edges = length(logmstar_edges[i])
			for j in 1: len_logmstar_edges - 1
				z_bin_catalogdata = my_subcatalog_struct(catalog,
											(catalog.z .> z_edges[i]) .& (catalog.z .<= z_edges[i+1]))
				number_gal_in_logmbin = sum((z_bin_catalogdata.logmstar .> logmstar_edges[i][j]) .& (z_bin_catalogdata.logmstar .<= logmstar_edges[i][j+1]))
				
                
				append!(number_gal_in_logmbins_vec, number_gal_in_logmbin)
			end
			log_number_density = log10.(number_gal_in_logmbins_vec ./ logmass_step ./ z_vols[i])
			σN = 1. ./ (log.(10.) .* sqrt.(number_gal_in_logmbins_vec)) .* log10.(1. ./ logmass_step ./ z_vols[i])
			#σSED = 0.04

            ### Removing logmass bins that has no galaxies
            noinf_ind = findall(!isinf, log_number_density)

			logmstar_center = (logmstar_edges[i][1:len_logmstar_edges - 1] .+ logmstar_edges[i][2:len_logmstar_edges]) ./ 2

            logmstar_center = logmstar_center[noinf_ind]
            log_number_density = log_number_density[noinf_ind]
            σN = σN[noinf_ind]    

            @assert length(logmstar_center) == length(log_number_density) == length(σN)
			
			log_density_gal[i] =  masterdata(z_edges[i], z_edges[i+1], logmstar_center, log_number_density, σN)
		end
		@assert size(log_density_gal, 1) == size(z_edges, 1) - 1
		return log_density_gal
	end
end

################################

begin
    """Custom `paramSMF` struct to read the final chains"""
	struct paramSMF{T1<:Number, T2<:Number, T3<:Number, T4<:Number, T5<:Number, T6<:Number, T7<:Number, T8<:Number, T9<:Number, T10<:Number, V1<:AbstractVector{T1}, V2<:AbstractVector{T2}, V3<:AbstractVector{T3}, V4<:AbstractVector{T4}, V5<:AbstractVector{T5}, V6<:AbstractVector{T6}, V7<:AbstractVector{T7}, V8<:AbstractVector{T8}, V9<:AbstractVector{T9}, V10<:AbstractVector{T10}}
		logMc::V1
		logMc_err::V2
		logϕ1::V3
		logϕ1_err::V4
		logϕ2::V5
		logϕ2_err::V6
		α1::V7
		α1_err::V8
		α2::V9
		α2_err::V10
	end


    @doc """read the chains and extract the fitting results"""
	function read_param_from_chains(chains)
		params = Array{Vector{Float64}}(undef,length(chains))

		for i in 1:length(chains)
			logMc = quantile(chains[i][:logMc][:,1], 0.5)
			logϕ1 = quantile(chains[i][:logϕ1][:,1], 0.5)
			logϕ2 = quantile(chains[i][:logϕ2][:,1], 0.5)
			α1 = quantile(chains[i][:α1][:,1], 0.5)
			α2 = quantile(chains[i][:α2][:,1], 0.5)
			
			logMc_err = (quantile(chains[i][:logMc][:,1], 0.84) .- quantile(chains[i][:logMc][:,1], 0.16)) ./ 2
			logϕ1_err = (quantile(chains[i][:logϕ1][:,1], 0.84) .- quantile(chains[i][:logϕ1][:,1], 0.16)) ./ 2
			logϕ2_err = (quantile(chains[i][:logϕ2][:,1], 0.84) .- quantile(chains[i][:logϕ2][:,1], 0.16)) ./ 2
			α1_err = (quantile(chains[i][:α1][:,1], 0.84) .- quantile(chains[i][:α1][:,1], 0.16)) ./ 2
			α2_err = (quantile(chains[i][:α2][:,1], 0.84) .- quantile(chains[i][:α2][:,1], 0.16)) ./ 2
			params[i] =  [logMc, logMc_err, logϕ1, logϕ1_err, logϕ2, logϕ2_err, α1, α1_err, α2, α2_err]
		end
		new_params::Matrix{Float64} = hcat(params...)
		return paramSMF(new_params[1,:], new_params[2,:], new_params[3,:], new_params[4,:], new_params[5,:], new_params[6,:], new_params[7,:], new_params[8,:], new_params[9,:], new_params[10,:])
		
	end
end


end