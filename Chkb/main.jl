#################################################################
# - Research Project:
#   Interconnectedness, Performance, Risk, and Response 
#   to Shocks: Evidence from the Banking Network in Chile.
# - Paper:
#   Chanci, Kumbhakar and Bobadilla - 
#   Interconnectedness and Banking Performance
# - File:
#   main.jl
#   This is the main file for processing the data and conducting 
#   the estimation. It executes the R file dataprocessing.R and the 
#   module Chkb.jl
# - By: 
#   Luis Chanci, 
#   www.luischanci.com
#   2025
#################################################################

# 1.   Initial Software Settings

# 1.1. Where the folder Chkb is set 
#      (the directory that contains the folder named Chkb)

wd = @__DIR__ # Use the directory of the current Chkb folder as the working directory
cd(wd)

# 1.2. Setting the required Julia packages (ej. for parallel computing)
import Pkg

required_pkg = [
    "Distributed", "DataFrames", "LinearAlgebra", "Distributions", "Optim", "Statistics", "Printf",
    "GLM", "Random", "RData", "CategoricalArrays", "Dates", "StatsModels", "RCall", "StatsBase"
]

for Paquete in required_pkg
    is_installed = any(x -> x.name == Paquete, values(Pkg.dependencies()))
    if !is_installed
        println(" Installing '$Paquete' ...")
        try
            Pkg.add(Paquete)
        catch e
            println("Error installing '$Paquete': $e") # Un error puede surgir por firewalls
        end
    end
end

for Paquete in required_pkg
    eval(Meta.parse("using $Paquete"))
end


# 1.3. Generating the data that will be used in the estimation 
#      (Processing the CSV data using the initial R code)
@rput wd

R"""
setwd(wd)
if (file.exists("dataprocessing.R")) {
    source("dataprocessing.R")
    cat("Successfully executed dataprocessing.R\n")
} else {
    stop("Error: dataprocessing.R not found in current directory.")
}
"""

# 1.4. Setting Julia for the estimation
#      (Setting up Julia for parallel computing and loading the module Chkb.jl)
if nprocs() == 1
    addprocs(Sys.CPU_THREADS - 1)
end
@everywhere using DataFrames, LinearAlgebra, Distributions, Optim, Statistics, StatsBase, GLM, Random, RData, CategoricalArrays, Dates

module_path = joinpath(wd, "Chkb.jl")
@everywhere include($module_path) 
@everywhere using .Chkb
using .Chkb

config_defaults = Chkb.create_config();

# define outcomes directory
outcomes_dir = joinpath(wd, "outcomes", "includes")
mkpath(outcomes_dir) # Checking that the directory exists

#################################################################
# Model 1: Basic model without risk or the determinants of the 
#          mean inefficiency
#################################################################

M1_config = (
    model_name = "M1",
    model_type = :basic
)
config = merge(config_defaults, M1_config)

println("\n--- STARTING FULL ANALYSIS (Model 1) ---")
M1 = Chkb.run_analysis(config);

# To display results in Julia REPL
fmt(x) = @sprintf("%.3f", x)
println("\n--- Results for Model 1 ---")
show(transform(M1.se_results, names(M1.se_results, Real) .=> ByRow(fmt); renamecols=false), allrows=true)
println("\n")

# Save M1 results dataframe
M1_df_out = M1.df_out
save_path_m1 = joinpath(outcomes_dir, "M1_df_out.rds")
@rput M1_df_out save_path_m1
R"""
saveRDS(M1_df_out, file = save_path_m1)
"""
println("M1 dataframe saved to: $save_path_m1")


#################################################################
# Model 2: Basic model with risk and without the determinants of  
#          the mean inefficiency
#################################################################

M2_config = (
    model_name = "M2",
    model_type = :basic,
    R_names    = ["risk"]
)

config = merge(config_defaults, M2_config)

println("\n--- STARTING FULL ANALYSIS (Model 2) ---")
M2 = Chkb.run_analysis(config);

println("\n--- Results for Model 2 ---")
show(transform(M2.se_results, names(M2.se_results, Real) .=> ByRow(fmt); renamecols=false), allrows=true)
println("\n")

# Save M2 results dataframe
M2_df_out = M2.df_out
save_path_m2 = joinpath(outcomes_dir, "M2_df_out.rds")
@rput M2_df_out save_path_m2
R"""
saveRDS(M2_df_out, file = save_path_m2)
"""
println("M2 dataframe saved to: $save_path_m2")


#################################################################
# Model 3: Model with risk and the determinants of  
#          the mean inefficiency
#################################################################
M3_config = (
    model_name = "M3",
    model_type = :determinants_v1,
    R_names    = ["risk"],
    D_names    = ["d_ownership", "d_segment"]
)
config = merge(config_defaults, M3_config)

println("\n--- STARTING FULL ANALYSIS (Model 3) ---")
M3 = Chkb.run_analysis(config);

println("\n--- Results for Model 3 ---")
show(transform(M3.se_results, names(M3.se_results, Real) .=> ByRow(fmt); renamecols=false), allrows=true)
println("\n")

# Save M3 results dataframe
M3_df_out = M3.df_out
save_path_m3 = joinpath(outcomes_dir, "M3_df_out.rds")
@rput M3_df_out save_path_m3
R"""
saveRDS(M3_df_out, file = save_path_m3)
"""
println("M3 dataframe saved to: $save_path_m3")

#################################################################
# End