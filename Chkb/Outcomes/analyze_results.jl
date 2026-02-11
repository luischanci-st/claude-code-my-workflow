#################################################################
# - File: analyze_results.jl
# - Description: Generates econometric comparison tables (M1-M3)
#   with significance stars and standard errors.
# - By: Luis Chanci, 2026
#################################################################

using DataFrames, Printf, Dates

wd = @__DIR__
cd(wd)

params_list = [
    ("Rho",         "rho"),
    ("Sigma V",     "sv"),
    ("Sigma U",     "su"),
    ("Delta (cons.)","delta_const"),
    ("Ownership",    "delta_d_ownership"),
    ("Segment",      "delta_d_segment")
]

model_def = Dict(
    "M1" => ["rho", "sv", "su"],
    "M2" => ["rho", "sv", "su"],
    "M3" => ["rho", "sv", "delta_const", "delta_d_ownership", "delta_d_segment"]
)

function get_est(folder, model, param)
    base_path = joinpath(folder, "includes", "$(model)_$(param)")
    coef = parse(Float64, strip(readline(base_path * ".txt")))
    se   = parse(Float64, strip(readline(base_path * "_se.txt")))
    stars = readline(base_path * "_pval.txt")
    return @sprintf("%.3f%s\n(%.3f)", coef, stars, se)
end

pattern = r"Includes.*\((\d+)\s*a\s*(\d+)\)"
folders = filter(x -> isdir(x) && occursin(pattern, x), readdir())

for folder in folders
        m = match(pattern, folder)
        period = m !== nothing ? m.match : folder
        
        println("\n" * "="^60)
        println(" RESULTS FOR SAMPLE: $period")
        println("="^60)
            
        df_table = DataFrame(Parameter = [p[1] for p in params_list])

        for model in ["M1", "M2", "M3"]
            valid_params = model_def[model]
            
            col_data = map(params_list) do p
                suffix = p[2]
                if suffix in valid_params
                    get_est(folder, model, suffix)
                else
                    ""
                end
            end
            df_table[!, model] = col_data
        end

        show(df_table, allrows=true, truncate=100)
        println("\n")
end