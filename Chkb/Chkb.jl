#################################################################
# - Research Project:
#   Interconnectedness, Performance, Risk, and Response 
#   to Shocks: Evidence from the Banking Network in Chile.
# - Paper:
#   Chanci, Kumbhakar and Bobadilla - 
#   Interconnectedness and Banking Performance
# - File:
#   Chkb.jl
#   This is the module containing the functions for the estimation. 
# - By: 
#   Luis Chanci, 
#   www.luischanci.com
#   2025
#################################################################

#################################################################
#  Notes
#  - Algorithm for estimation. Returns SEs and p-values for rho, 
#    sv, and delta (incl. intercept).
#  - Optimized: 
#    Using parallel computing and I also pre-compute the spatial 
#    matrices to speed up bootstrap.
#################################################################

module Chkb #Chanci-Kumbhakar-Bobadilla

# Initial Settings

# --- i) Load Dependencies --- (double checking they exist)
using Distributed, DataFrames, LinearAlgebra, Distributions
using Optim, Statistics, StatsBase, GLM, Random, RData
using CategoricalArrays, Dates, StatsModels

# --- ii) Export Public Functions ---
export run_analysis, create_config
export prepare_translog_data, estimate_step1_ols,
       prepare_step2_inputs, estimate_step2_gmm,
       calculate_post_estimation, calculate_post_estimation_determinants,
       run_wild_bootstrap, precompute_gmm_matrices
export optimize_rho_value, run_bootstrap_replication 

#################################################################
#  MODULE 0: SOME UTILITY HELPER FUNCTIONS
#################################################################

log1p_safe(x) = ifelse.(x .== 0, log1p.(x), log.(x))   # In case of values = 0

frobenius_norm(M::AbstractMatrix) = sqrt(sum(abs2, M)) # Frobenius norm

compute_inv_W(W::AbstractMatrix, rho::Float64) = inv(I - rho * W) # (I - rho*W)^(-1)

function build_moment_inputs(dates::Vector,            # Build S_t and optionally D_matrix
                             gdf_step2::GroupedDataFrame,
                             include_determinants::Bool,
                             residual_col::Symbol,
                             level_residual_col::Symbol,
                             D_names::Vector{String})
    S_t      = Dict{Any, Matrix{Float64}}()
    D_matrix = include_determinants ? Dict{Any, Matrix{Float64}}() : nothing
    e_list   = Dict{Any, Vector{Float64}}()

    for tt in dates
        df_t       = gdf_step2[(tt,)]
        bar_e_t    = df_t[!, residual_col]
        S_t[tt]    = bar_e_t * bar_e_t'
        e_list[tt] = df_t[!, level_residual_col]
        if include_determinants
            d_cols = Symbol.(D_names)
            if isempty(d_cols)
                D_matrix[tt] = hcat(ones(nrow(df_t)))
            else
                D_matrix[tt] = hcat(ones(nrow(df_t)), Matrix(df_t[!, d_cols]))
            end
        end
    end
    return S_t, D_matrix, e_list
end

function estimate_step1_ols(data_levels::DataFrame,    # For Step 1 - OLS (using "within data")
                            ols_regressors::Vector{Symbol}, 
                            y_var::Symbol)
    rhs_terms_vector = map(term, ols_regressors)
    rhs_expression   = sum(rhs_terms_vector)
    ols_formula      = term(y_var) ~ term(1) + rhs_expression
    ols_model        = lm(ols_formula, data_levels)
    Beta_hat         = Dict(zip(Symbol.(coefnames(ols_model)), coef(ols_model)))
    return (Beta_hat = Beta_hat, ols_formula = ols_formula, lm_obj = ols_model)
end

function precompute_gmm_matrices(W_list, rhos_grid, Dates) # Pre-compute Winv*Winv' for all rho and t
    WW_dict  = Dict{Float64, Dict{Any, Matrix{Float64}}}()
    for rho in rhos_grid
        WW_t = Dict{Any, Matrix{Float64}}()
        for tt in Dates
            Winv_tt  = compute_inv_W(W_list[tt], rho)
            WW_t[tt] = Winv_tt * Winv_tt'
        end
        WW_dict[rho] = WW_t
    end
    return WW_dict
end

# GMM objective for the basic model
function gmm_objective_basic(theta, rho_val, S_t, Dates, WW_dict_for_rho)
    log_sv, log_su = theta
    scalar_term    = exp(2 * log_sv) + exp(2 * log_su) * (1.0 - 2.0 / pi)
    total_norm     = 0.0
    for tt in Dates
        V_theoretical = scalar_term .* WW_dict_for_rho[tt]
        total_norm += frobenius_norm(S_t[tt] .- V_theoretical)
    end
    return total_norm
end

# GMM objective for the model with determinants
function gmm_objective_determinants(theta, rho_val, S_t, D_matrix, Dates, WW_dict_for_rho, W_list)
    log_sv = theta[1]
    delta  = theta[2:end]
    sv2    = exp(2 * log_sv)
    total_norm = 0.0
    for tt in Dates
        d_tt  = D_matrix[tt]
        dd_tt = exp.(d_tt * delta) 
        Winv  = compute_inv_W(W_list[tt], rho_val)
        WW    = WW_dict_for_rho[tt]
        WdW   = Winv * Diagonal(dd_tt .^ 2) * Winv'
        V_th  = (sv2 .* WW) .+ ((1.0 - 2.0 / pi) .* WdW)
        total_norm += frobenius_norm(S_t[tt] .- V_th)
    end
    return total_norm
end

# Optimizer for a single rho
function optimize_rho_value(rho_val::Float64, gmm_data, model_type::Symbol, WW_dict_for_rho::Dict)
    (; S_t, Dates, D_matrix, W_list) = gmm_data
    local obj_fn, initial_theta
    if model_type == :basic
        obj_fn = (theta) -> gmm_objective_basic(theta, rho_val, S_t, Dates, WW_dict_for_rho)
        initial_theta = [1e-3, 1e-3]
    else
        k_d = size(D_matrix[Dates[1]], 2)  # intercept + D_names
        obj_fn = (theta) -> gmm_objective_determinants(theta, rho_val, S_t, D_matrix, Dates, WW_dict_for_rho, W_list)
        initial_theta = fill(log(1e-3), 1 + k_d)
    end
    opt = optimize(obj_fn, initial_theta, NelderMead())
    return (rho=rho_val, obj_fn_value=Optim.minimum(opt), minimizer=Optim.minimizer(opt))
end

# For Step 2 - GMM driver
function estimate_step2_gmm(gmm_data, rhos_grid, model_type::Symbol, use_parallel::Bool; precomputed_WW=nothing)
    if use_parallel
        println("--- 4. Running Step 2 (GMM) Grid Search (Parallel) ---")
    end

    if isnothing(precomputed_WW)
        println("Pre-computing GMM matrices for $(length(rhos_grid)) rho values ...")
        WW_dict = precompute_gmm_matrices(gmm_data.W_list, rhos_grid, gmm_data.Dates)
        println("Pre-computation complete.")
    else
        WW_dict = precomputed_WW
    end

    opt_store_list = use_parallel ?
        pmap(rho_val -> Chkb.optimize_rho_value(rho_val, gmm_data, model_type, WW_dict[rho_val]), rhos_grid) :
        map(rho_val -> Chkb.optimize_rho_value(rho_val, gmm_data, model_type, WW_dict[rho_val]), rhos_grid)

    opt_store_df = DataFrame(opt_store_list)
    min_opt      = opt_store_df[argmin(opt_store_df.obj_fn_value), :]

    local gmm_results
    if model_type == :basic
        log_sv, log_su = min_opt.minimizer
        gmm_results    = (rho=min_opt.rho, sv=exp(log_sv), su=exp(log_su), delta=Float64[], obj_fn_value=min_opt.obj_fn_value)
    else
        log_sv = min_opt.minimizer[1]
        delta  = min_opt.minimizer[2:end]
        gmm_results = (rho=min_opt.rho, sv=exp(log_sv), su=NaN, delta=delta, obj_fn_value=min_opt.obj_fn_value)
    end

    if use_parallel
        println("GMM estimation complete. Optimal rho: $(Base.round(gmm_results.rho, digits=4))")
    end
    return (gmm_results = gmm_results, opt_store_df = opt_store_df)
end

# Build W_list and ID_list by period
function build_W_list(all_W_df::DataFrame, dates::Vector, gdf_step2::GroupedDataFrame)
    W_list     = Dict{Any, Matrix{Float64}}()
    ID_list    = Dict{Any, Vector{String}}()
    for tt in dates
        df_t        = gdf_step2[(tt,)]
        IDs_t       = string.(df_t.ID)
        W_df_t      = filter(row -> row.PERIOD == tt && row.ID in IDs_t, all_W_df)
        W_df_t      = W_df_t[indexin(IDs_t, W_df_t.ID), :]
        Wt_mat      = Matrix{Float64}(W_df_t[!, IDs_t])
        row_sums    = sum(Wt_mat, dims=2)
        Wt_norm     = Wt_mat ./ row_sums
        Wt_norm[isnan.(Wt_norm)] .= 0.0
        W_list[tt]  = Wt_norm
        ID_list[tt] = IDs_t
    end
    return W_list, ID_list
end


#################################################################
#  MODULE 1: DATA PREPARATION
#################################################################

function prepare_translog_data(df::DataFrame, Y_names, X_names, Z_names, R_names, p1_name::String)
    println("--- 1. Preparing Translog Data for the Cost Function---")
    reg_df = DataFrame(ID = df.INS_COD, PERIOD = df.PERIODO)

    for var in Y_names
        reg_df[!, Symbol("lY_", var)] = log1p_safe(df[!, Symbol("Y_", var)])
    end
    p1_vec = df[!, Symbol(p1_name)]
    for var in X_names
        reg_df[!, Symbol("lp_", var)] = log1p_safe(df[!, Symbol("p_", var)] ./ p1_vec)
    end
    for var in Z_names
        reg_df[!, Symbol("lz_", var)] = log1p_safe(df[!, Symbol("X_", var)])
    end
    if !isempty(R_names)
        for var in R_names
            reg_df[!, Symbol("lr_", var)] = log1p_safe(df[!, Symbol(var)])
        end
    end

    translog_vars = names(reg_df, r"^(lY_|lp_|lz_|lr_)")

    for var in translog_vars
        reg_df[!, Symbol(var, "_sq")] = reg_df[!, var] .^ 2
    end
    for i in 1:length(translog_vars), j in (i+1):length(translog_vars)
        reg_df[!, Symbol(translog_vars[i], "_x_", translog_vars[j])] =
            0.5 .* reg_df[!, translog_vars[i]] .* reg_df[!, translog_vars[j]]
    end

    reg_df.lcostw   = log1p_safe(abs.(df.Cost ./ p1_vec))
    regressor_names = names(reg_df, r"^(lY_|lp_|lz_|lr_|lcostw)")
    gdf             = groupby(reg_df, :PERIOD)
    reg_df_with     = transform(gdf, regressor_names .=> (x -> x .- mean(x)) .=> Symbol.(regressor_names .* "_within"))

    data_levels     = hcat(reg_df, select(reg_df_with, Not(names(reg_df))), df[!, r"^d_"], makeunique=true)
    sort!(data_levels, [:PERIOD, :ID])

    ols_regressors  = Symbol.(filter!(x -> x != "lcostw_within", names(data_levels, r"_within$")))
    y_var           = :lcostw_within

    return (data_levels = data_levels, ols_regressors = ols_regressors, y_var = y_var, translog_names = Symbol.(translog_vars))
end

#################################################################
#  MODULE 2: STEP 2 (GMM) PREPARATION (los datos)
#################################################################

function prepare_step2_inputs(data_levels::DataFrame, Beta_hat::Dict, ols_regressors::Vector{Symbol},
                              wd::String, W_to_use::String, model_type::Symbol, D_names::Vector{String})
    println("--- 3. Preparing Data for Step 2 (GMM) Inputs ---")

    beta_names_level = [Symbol(replace(string(s), "_within" => "")) for s in ols_regressors]
    beta_values      = [Beta_hat[s] for s in ols_regressors]
    Xb_level         = Matrix(data_levels[!, beta_names_level]) * beta_values

    df_step2     = copy(data_levels)
    df_step2.eit = data_levels.lcostw .- Xb_level
    df_step2.Xb  = Xb_level

    tokeep = ["OR","PP","PS","PU","PW","QN","RO","RR","SN","SP","SR","SV","TP",
              "TR","TT","TV","UP","UQ","SOU","ST"] # Check again only use main banks
    filter!(row -> row.ID in tokeep, df_step2)
    df_step2.ID = categorical(df_step2.ID, levels=tokeep)
    sort!(df_step2, [:PERIOD, :ID])

    gdf   = groupby(df_step2, :PERIOD)
    transform!(gdf, :eit => mean => :alpha_star, :eit => (x -> x .- mean(x)) => :bar_e)
    Dates = unique(df_step2.PERIOD)

    println("Loading W matrix: $W_to_use")
    file_path = joinpath(wd, "Data/Processed/", W_to_use * ".Rda") # Ojo si cambio la ruta
    all_W_df  = DataFrame(RData.load(file_path)[W_to_use])

    W_list, ID_list       = build_W_list(all_W_df, Dates, gdf)
    S_t, D_matrix, e_list = build_moment_inputs(Dates, gdf, model_type != :basic, :bar_e, :eit, D_names)

    gmm_data = (S_t=S_t, D_matrix=D_matrix, e_list=e_list, W_list=W_list, ID_list=ID_list, Dates=Dates)
    return (gmm_data = gmm_data, df_step2 = df_step2)
end

#################################################################
#  MODULE 3: STEP 2 (GMM) ESTIMATION (la function principal)
#################################################################
# estimate_step2_gmm defined above

#################################################################
#  MODULE 4: POST-ESTIMATION
#################################################################

# Basic model post-estimation
function calculate_post_estimation(df_step2::DataFrame, gmm_data, gmm_results)
    println("--- 5. Running Post-Estimation (Basic Model) ---")
    (; W_list, Dates, ID_list) = gmm_data
    (; rho, sv, su) = gmm_results
    if isnan(su)
        println("Error: `su` is NaN. Cannot run basic post-estimation.")
        return df_step2
    end

    sv2, su2 = sv^2, su^2
    s_star = sqrt((su2 * sv2) / (sv2 + su2))

    Winv_opt = Dict(tt => compute_inv_W(W_list[tt], rho) for tt in Dates)
    Wrho_opt = Dict(tt => (I - rho * W_list[tt]) for tt in Dates)

    results_list = Vector{DataFrame}()
    gdf = groupby(df_step2, :PERIOD)

    for tt in Dates
        df_t = gdf[(tt,)]
        IDs_tt = ID_list[tt]
        df_t = df_t[indexin(IDs_tt, df_t.ID), :]
        Winv_tt, Wrho_tt = Winv_opt[tt], Wrho_opt[tt]

        Eu_t = (sqrt(2.0 / pi) * su * (sum(Winv_tt, dims=2)))[1]
        beta_0_t = df_t.alpha_star .- Eu_t
        vu_t = df_t.bar_e .+ Eu_t
        vu_dot_t = Wrho_tt * vu_t

        mu_star = vu_dot_t .* (su2 / (sv2 + su2))
        a_it = -mu_star ./ s_star

        if s_star < 1e-10
            inef_dot_t = max.(mu_star, 0.0)
        else
            phi = pdf.(Normal(), a_it)
            Phi_r = 1.0 .- cdf.(Normal(), a_it)
            jlms_term = phi ./ Phi_r
            jlms_term[Phi_r .< 1e-10] .= 0.0
            inef_dot_t = mu_star .+ s_star .* jlms_term
            inef_dot_t[isinf.(inef_dot_t)] .= mu_star[isinf.(inef_dot_t)]
            inef_dot_t = max.(inef_dot_t, 0.0)
        end

        Ineff_t = Winv_tt * inef_dot_t
        D_eff_t = diag(Winv_tt) .* inef_dot_t
        push!(results_list, DataFrame(PERIOD=tt, ID=IDs_tt, beta_0=beta_0_t, vu_dot=vu_dot_t, Ineff=Ineff_t, D_eff=D_eff_t))
    end

    df_out = leftjoin(df_step2, vcat(results_list...), on = [:PERIOD, :ID])
    println("Post-estimation complete.")
    return df_out
end

# Determinants model post-estimation:
function calculate_post_estimation_determinants(df_step2::DataFrame, gmm_data, gmm_results)
    println("--- 5. Running Post-Estimation (Determinants Model) ---")
    (; W_list, Dates, ID_list, D_matrix) = gmm_data
    (; rho, sv, delta) = gmm_results
    sv2 = sv^2

    Winv_opt = Dict(tt => compute_inv_W(W_list[tt], rho) for tt in Dates)
    Wrho_opt = Dict(tt => (I - rho * W_list[tt]) for tt in Dates)

    results_list = Vector{DataFrame}()
    gdf          = groupby(df_step2, :PERIOD)

    for tt in Dates
        df_t   = gdf[(tt,)]
        IDs_tt = ID_list[tt]
        df_t   = df_t[indexin(IDs_tt, df_t.ID), :]
        Winv_tt, Wrho_tt = Winv_opt[tt], Wrho_opt[tt]
        d_tt   = D_matrix[tt]

        su_vec_t  = exp.(d_tt * delta)
        su2_vec_t = su_vec_t .^ 2

        Eu_vec_t = sqrt(2.0 / pi) .* (Winv_tt * su_vec_t)
        beta_0_t = df_t.alpha_star .- Eu_vec_t
        vu_t     = df_t.bar_e .+ Eu_vec_t
        vu_dot_t = Wrho_tt * vu_t

        mu_star_t = vu_dot_t .* su2_vec_t ./ (sv2 .+ su2_vec_t)
        s_star_t  = sqrt.(su2_vec_t .* sv2 ./ (sv2 .+ su2_vec_t))
        a_it      = -mu_star_t ./ s_star_t

        inef_dot_t = zeros(length(a_it))
        for i in eachindex(a_it)
            if s_star_t[i] < 1e-10
                inef_dot_t[i] = max(mu_star_t[i], 0.0)
            else
                phi   = pdf(Normal(), a_it[i])
                Phi_r = 1.0 - cdf(Normal(), a_it[i])
                jlms  = (Phi_r < 1e-10) ? 0.0 : (phi / Phi_r)
                inef_dot_t[i] = max(mu_star_t[i] + s_star_t[i] * jlms, 0.0)
            end
        end

        Ineff_t = Winv_tt * inef_dot_t
        D_eff_t = diag(Winv_tt) .* inef_dot_t
        push!(results_list, DataFrame(PERIOD=tt, ID=IDs_tt, beta_0=beta_0_t, vu_dot=vu_dot_t, Ineff=Ineff_t, D_eff=D_eff_t))
    end

    df_out = leftjoin(df_step2, vcat(results_list...), on = [:PERIOD, :ID])
    println("Post-estimation (determinants) complete.")
    return df_out
end

#################################################################
#  MODULE 5: WILD BOOTSTRAP
#################################################################

# First (one) bootstrap replication
function run_bootstrap_replication(
    bb::Int,
    df_b_data,
    gmm_data_const,
    rhos_grid::AbstractRange,
    model_type::Symbol,
    config,
    WW_dict_precomputed
)
    Random.seed!(123 + bb)

    # i) Mammen two-point
    prob_p = (sqrt(5.0) + 1.0) / (2.0 * sqrt(5.0))
    twop   = [-(sqrt(5.0) - 1.0) / 2.0, (sqrt(5.0) + 1.0) / 2.0]
    wts    = Weights([1.0 - prob_p, prob_p])

    df_boot = copy(df_b_data)
    eta_b = sample(twop, wts, nrow(df_boot))
    df_boot.vu_dot_b = df_boot.vu_dot .* eta_b

    # ii) Synthetic dependent variable by period
    gdf           = groupby(df_boot, :PERIOD)
    lcostw_b_list = Vector{DataFrame}()
    for tt in gmm_data_const.Dates
        df_t      = gdf[(tt,)]
        IDs_tt    = gmm_data_const.ID_list[tt]
        df_t      = df_t[indexin(IDs_tt, df_t.ID), :]
        Winv_tt   = compute_inv_W(gmm_data_const.W_list[tt], gmm_data_const.original_rho)
        lcostw_bt = df_t.alpha_star .+ df_t.Xb .+ (Winv_tt * df_t.vu_dot_b)
        push!(lcostw_b_list, DataFrame(PERIOD=tt, ID=IDs_tt, lcostw_b=lcostw_bt))
    end
    df_boot = leftjoin(df_boot, vcat(lcostw_b_list...), on=[:PERIOD, :ID])

    # iii) Re-estimate Step 1 (OLS within)
    transform!(groupby(df_boot, :PERIOD), :lcostw_b => (x -> x .- mean(x)) => :lcostw_b_within)
    ols_regressors_b = config.ols_regressors
    step1_results_b = estimate_step1_ols(df_boot, ols_regressors_b, :lcostw_b_within)

    # iv) Rebuild moments for Step 2
    beta_names_level_b = [Symbol(replace(string(s), "_within" => "")) for s in ols_regressors_b]
    beta_values_b      = [step1_results_b.Beta_hat[s] for s in ols_regressors_b]
    Xb_level_b         = Matrix(df_boot[!, beta_names_level_b]) * beta_values_b
    df_boot.eit_b      = df_boot.lcostw_b .- Xb_level_b
    gdf_b              = groupby(df_boot, :PERIOD)
    transform!(gdf_b, :eit_b => (x -> x .- mean(x)) => :bar_e_b)

    S_t_b, D_matrix_b, _ = build_moment_inputs(gmm_data_const.Dates, gdf_b, model_type != :basic, :bar_e_b, :eit_b, config.D_names)

    gmm_data_b = (
        S_t     = S_t_b,
        D_matrix= D_matrix_b,
        e_list  = nothing,
        W_list  = gmm_data_const.W_list,
        ID_list = gmm_data_const.ID_list,
        Dates   = gmm_data_const.Dates
    )

    # v) Step 2 (serial within worker) - nota: Pass precomputed matrices
    gmm_results_tuple_b = estimate_step2_gmm(gmm_data_b, rhos_grid, model_type, false; precomputed_WW=WW_dict_precomputed)
    gmm_results_b       = gmm_results_tuple_b.gmm_results

    if model_type == :basic
        return [gmm_results_b.rho, gmm_results_b.sv, gmm_results_b.su]
    else
        return [gmm_results_b.rho, gmm_results_b.sv, gmm_results_b.delta...]
    end
end

_delta_labels(D_names::Vector{String}) = isempty(D_names) ? ["delta_const"] :
                                         vcat(["delta_const"], ["delta_" * n for n in D_names])

# Driver for wild bootstrap; returns DataFrame of SEs and p-values
function run_wild_bootstrap(df_out, gmm_data, gmm_results, config, ols_regressors, WW_dict_precomputed)
    println("--- 6. Running Wild Bootstrap (B = $(config.B_reps)) ---")

    gmm_data_const = merge(gmm_data, (original_rho = gmm_results.rho,))
    config_const   = merge(config, (ols_regressors = ols_regressors,))

    df_b_data = select(df_out, :PERIOD, :ID, :vu_dot, :Xb, :alpha_star,
                       ols_regressors...,
                       Not([:eit, :bar_e, :beta_0, :Ineff, :D_eff]))

    Boot_list = pmap(
        bb -> Chkb.run_bootstrap_replication(
              bb, df_b_data, gmm_data_const, config.rhos_grid, config.model_type, config_const, WW_dict_precomputed
        ),
        1:config.B_reps
    )
    Boot_out_matrix = hcat(Boot_list...)'

    # Parameter names and originals
    local param_names::Vector{String}
    local original_params::Vector{Float64}
    if config.model_type == :basic
        param_names = ["rho", "sv", "su"]
        original_params = [gmm_results.rho, gmm_results.sv, gmm_results.su]
    else
        dlabels = _delta_labels(config.D_names)
        param_names = vcat(["rho", "sv"], dlabels)
        original_params = vcat([gmm_results.rho, gmm_results.sv], gmm_results.delta)
    end

    Boot_out_df = DataFrame(Boot_out_matrix, param_names)
    ses = std.(eachcol(Boot_out_df))

    # z-scores with zero-SE guard
    zscores = map((m, s) -> s > 0 ? m / s : NaN, original_params, ses)
    pvals   = map(z -> isnan(z) ? NaN : 2 * (1.0 - cdf(Normal(), abs(z))), zscores)

    se_results = DataFrame(Parameter = param_names, Estimate = original_params, StdError = ses, PValue = pvals)

    println("Bootstrap complete.")
    return (bootstrap_df = Boot_out_df, se_results = se_results)
end

#################################################################
#  MODULE 6: MAIN CONTROLLER FUNCTION (deberia ser el modulo 0)
#################################################################

function create_config() #Uso esta parte para definir el modelo
    return (
        model_name = "",
        wd         = @__DIR__,
        df_path    = joinpath(@__DIR__, "Data/Processed/",  "df.Rda"),
        Y_names    = ["cloan", "rsloan", "hhloan", "secur"],
        X_names    = ["capital", "deposit"],
        Z_names    = ["equity"],
        model_type = :basic,       # or :determinants_v1
        R_names    = String[],     # set ["risk"],
        D_names    = String[],     # set e.g. ["d_size","d_liq",...]
        p1_name    = "p_labor",
        W_to_use   = "all_W_df_Tot",
        rhos_grid  = -0.600:0.001:0.030,
        do_post_estimation = true,
        do_bootstrap = true,
        B_reps     = 50,
        use_parallel_estimation = true
    )
end

function run_analysis(config)
    println("Starting SFA Analysis: $(config.W_to_use)")
    println("Timestamp: $(now())")
    println("-------------------------------------------")

    # i) Load and prepare data
    df_raw = RData.load(config.df_path) |> x -> DataFrame(x["df"])
    prep_results = prepare_translog_data(
        df_raw, config.Y_names, config.X_names, config.Z_names, config.R_names, config.p1_name
    )

    # ii) Step 1 OLS
    step1_results = estimate_step1_ols(prep_results.data_levels, prep_results.ols_regressors, prep_results.y_var)

    # iii) Step 2 GMM prep
    gmm_prep_results = prepare_step2_inputs(
        prep_results.data_levels, step1_results.Beta_hat, prep_results.ols_regressors,
        config.wd, config.W_to_use, config.model_type, config.D_names
    )

    WW_dict_global = precompute_gmm_matrices(gmm_prep_results.gmm_data.W_list, config.rhos_grid, gmm_prep_results.gmm_data.Dates)

    # iv) Step 2 GMM estimation (passing precomputed matrices)
    gmm_results_tuple = estimate_step2_gmm(
        gmm_prep_results.gmm_data, config.rhos_grid, config.model_type, config.use_parallel_estimation;
        precomputed_WW = WW_dict_global
    )
    gmm_results = gmm_results_tuple.gmm_results

    # Guardo los estimates para luego cargarlos en LaTeX
    output_dir = joinpath(config.wd, "outcomes/includes") #Ojo con la ruta en caso de cambio
    mkpath(output_dir)
    model_prefix = isempty(config.model_name) ? "" : config.model_name * "_"

    # Create a temporary DataFrame of results to iterate over
    local main_params_df
    if config.model_type == :basic
        main_params_df = DataFrame(Parameter=["rho", "sv", "su"], Estimate=[gmm_results.rho, gmm_results.sv, gmm_results.su])
    else
        dlabels = _delta_labels(config.D_names)
        main_params_df = DataFrame(Parameter=vcat(["rho", "sv"], dlabels), Estimate=vcat([gmm_results.rho, gmm_results.sv], gmm_results.delta))
    end

    for row in eachrow(main_params_df)
        open(joinpath(output_dir, "$(model_prefix)$(row.Parameter).txt"), "w") do io
            println(io, Base.round(row.Estimate, digits=3))
        end
    end
    println("Main parameter estimates saved to text files.")

    # v) Post-estimation
    df_out = gmm_prep_results.df_step2
    if config.do_post_estimation
        if config.model_type == :basic
            df_out = calculate_post_estimation(df_out, gmm_prep_results.gmm_data, gmm_results)
        elseif config.model_type == :determinants_v1
            df_out = calculate_post_estimation_determinants(df_out, gmm_prep_results.gmm_data, gmm_results)
        else
            println("No post-estimation function defined for model_type: $(config.model_type). Skipping.")
        end
    end

    # vi) Bootstrap
    local bootstrap_results = nothing
    local se_results = nothing
    if config.do_bootstrap
        if !config.do_post_estimation
            println("Warning: Bootstrap requires post-estimation ... computing now.")
            if config.model_type == :basic
                df_out = calculate_post_estimation(df_out, gmm_prep_results.gmm_data, gmm_results)
            elseif config.model_type == :determinants_v1
                df_out = calculate_post_estimation_determinants(df_out, gmm_prep_results.gmm_data, gmm_results)
            end
        end

        boot_tuple = run_wild_bootstrap( # Pass WW_dict_global to bootstrap to avoid repeated inversions
            df_out, gmm_prep_results.gmm_data, gmm_results, config, prep_results.ols_regressors, WW_dict_global
        )
        bootstrap_results = boot_tuple.bootstrap_df
        se_results = boot_tuple.se_results

        if !isnothing(se_results)         # Guardo los SEs to (for LaTeX)
            for row in eachrow(se_results)

                open(joinpath(output_dir, "$(model_prefix)$(row.Parameter)_se.txt"), "w") do io
                    println(io, Base.round(row.StdError, digits=3))
                end

                p_val = row.PValue
                asterisks = ""
                if !isnan(p_val) #Armé esta vaina pa mejor quedarme con los *
                    if p_val <= 0.01
                        asterisks = "***"
                    elseif p_val <= 0.051
                        asterisks = "**"
                    elseif p_val <= 0.10
                        asterisks = "*"
                    end
                end
                open(joinpath(output_dir, "$(model_prefix)$(row.Parameter)_pval.txt"), "w") do io
                    print(io, asterisks) # Use print to avoid extra newline
                end
            end
            println("Standard errors and p-value significance saved to text files.")
        end
    end

    println("-------------------------------------------")
    println("Analysis Complete.")

    return (
        config      = config,
        step1_model = step1_results.lm_obj,
        gmm_data    = gmm_prep_results.gmm_data,
        gmm_results = gmm_results,
        opt_store_df= gmm_results_tuple.opt_store_df,
        df_out      = df_out,
        bootstrap_results = bootstrap_results,
        se_results  = se_results
    )
end

#################################################################
end # End of module Chkb; Chanci
#################################################################
