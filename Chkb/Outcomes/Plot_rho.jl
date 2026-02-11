using Plots, Printf

wd = @__DIR__
cd(wd)

plotly() # To save at the end as savefig("rho_stability_plot.eps")

#ENV["GKSwstype"] = "100"
#gr() # to save at the end as savefig("rho_stability_plot.eps")


pattern = r"Includes.*\((\d+)\s*a\s*(\d+)\)"
folders = filter(x -> isdir(x) && occursin(pattern, x), readdir())
data    = sort(
            map(folders) do folder
                date     = match(pattern, folder).captures[2]
                rho_path = joinpath(folder, "includes", "M3_rho.txt")
                rho      = parse(Float64, strip(readline(rho_path)))
                se_path  = joinpath(folder, "includes", "M3_rho_se.txt") 
                se       = parse(Float64, strip(readline(se_path)))
                (date = date, rho = rho, ci_margin = 1.96 * se)
            end, 
          by = x -> x.date)

x_vals = [d.date for d in data]
y_vals = [d.rho  for d in data]
y_errs = [d.ci_margin for d in data]

plt = plot(
    x_vals, y_vals;
    seriestype = :scatter,
    yerror     = y_errs,
    marker     = (:circle, 8, :blue),
    line       = (:dash, :gray),
    ylims      = (-0.5, 0.5),
    legend     = false,
    xrotation  = 45,
    grid       = true,
    size       = (1000, 600),
    title      = "Stability of Estimated Rho across Subsamples (with 95% C.I.)",
    xlabel     = "Start Date (End Date fixed at 2018m1)",
    ylabel     = "Estimated Rho"
    )

display(plt)
savefig("rho_stability_plot.eps")

