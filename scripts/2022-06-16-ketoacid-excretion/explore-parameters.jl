using DrWatson
@quickactivate "CoarseGrainedPombe.jl"

using CSV
using Dates
using Distributions
using Random

modelname = "ketoacid-excretion"

savefile = datadir(modelname, string(today()), "results.csv")

# experiment = "sweep_kC"
# experiment = "sweep_kN"
# experiment = "sweep_kC_nores"
# experiment = "sweep_kC_noferm"
experiment = "random_kCkEfkN"



include(srcdir(modelname * ".jl"))
include(srcdir("find_steady_state.jl"))
include(srcdir("optimise_allocation.jl"))
include(srcdir("parameter_sweep.jl"))

p0 = (
  k_Kre =    10.0,    # Maximal keto-acid recycling enzyme rate (h^-1)
  k_Kex =    20.0,    # Maximal keto-acid excretion enzyme rate (h^-1)
  k_C =      10.0,    # Carbon transporter rate (h^-1)
  k_Ef =     15.0,    # Maximal fermentative enzyme rate (h^-1)
  k_Er =      7.5,    # Maximal respiratory enzyme rate (h^-1)
  k_N =      20.0,    # Nitrogen transporter rate (h^-1)
  k_R =       6.46,   # Maximal ribosomal synthesis rate (h^-1)
  k_sat =     0.0167, # Keto-acid Michaelis constant for keto-acid recycling enzyme
  c_sat =     0.0167, # Carbon Michaelis constant for fermentative/respirotory enzymes
  a_sat =     0.0167, # Precursor Michaelis constant for ribosome
  n_sat =     0.0167, # Nitrogen Michaelis constant for fermentative/respiratory enzymes
  Δ_CN  =     0.0,    # Nonlinearity in fermentative/respiratory enzymes
  f_Kre =     0.0,    # Allocation towards keto-acid recycling enzymes
  f_Kex =     0.0,    # Allocation towards keto-acid excretion enzymes
  f_C   =     0.2,    # Allocation towards carbon transporters
  f_Ef  =     0.15,  # sweep_kC, sweep_kN, random_kCkEfkN
  f_Er  =     0.05,  # sweep_kC, sweep_kN, random_kCkEfkN 
  # f_Ef  =     0.2,   # sweep_kC_nores
  # f_Er  =     0.0,   # sweep_kC_nores
  # f_Ef  =     0.0,    # sweep_kC_noferm 
  # f_Er  =     0.2,    # sweep_kC_noferm  
  f_N   =     0.1,   # Allocation towards nitrogen transporters
  f_R   =     0.3,    # Allocation towards ribosomes
  f_Z   =     0.2,    # Allocation towards non-metabolic proteins (Q & X)
  γ_K   =     0.0,    # Proportion of keto-acid produced by nitrogen transporters
  α_Cr  =   24//31,  # Proportion of carbon consumed by fermentative enzymes
  α_Cf  =   48//55,   # Proportion of carbon consumed by respiratory enzymes
)

function force_zeros(parameters, symbols_tozero)
  if symbols_tozero == Symbol[]
      return(parameters)
  end
  f_R = parameters.f_R + sum(parameters[i] for i in symbols_tozero)
  parameters = merge(
      parameters,
      NamedTuple(i => 0 for i in symbols_tozero)
      )
  parameters = merge(parameters, (f_R = f_R,))
  
  return(parameters)    
end

# At least one of E_C or E_Kr must be expressed for usable carbon to be present in the cell
try_zeros = [Symbol[], [:f_Ef], [:f_Er]]  # sweep_kC, sweep_kN, random_kCkEfkN
# try_zeros = [Symbol[],]  # sweep_kC_nores, sweep_kC_noferm
base_allocs = [:f_C, :f_Ef, :f_Er, :f_N, :f_R]  # sweep_kC, sweep_kN, random_kCkEfkN
# base_allocs = [:f_C, :f_Ef, :f_N, :f_R]  # sweep_kC_nores
# base_allocs = [:f_C, :f_Er, :f_N, :f_R]  # sweep_kC_noferm
try_allocs = [setdiff(base_allocs, try_zero) for try_zero in try_zeros]

randomisation_pars = (a = 0.0, b = 20.0, n = 100)  # random_kCkEfkN
dists = [Uniform(randomisation_pars.a, randomisation_pars.b) for i in 1:3]  # random_kCkEfkN
Random.seed!(20220906)  # random_kCkEfkN

# rates = logrange(1e-2, 1e3, 101)  # sweep_kC, sweep_kN, sweep_kC_nores, sweep_kC_noferm
rates = [rand.(dists) for i in 1:100]  # random_kCkEfkN
nt = []
for (iter, rate) in enumerate(rates)
  progress_msg = string(iter)*"/"*string(length(rates))
  @info progress_msg rate

  # p = merge(p0, (k_C = rate,))  # sweep_kC, sweep_kC_nores, sweep_kC_noferm
  # p = merge(p0, (k_N = rate,))  # sweep_kN
  p = merge(p0, (k_C = rate[1], k_Ef = rate[2], k_N = rate[3]))  # random_kCkEfkN
  
  params = [force_zeros(p, try_zero) for try_zero in try_zeros]
  all_results = []    
  for i in eachindex(params)
    p, c = optimise_allocation(
      evolve_concentrations!,
      params[i];
      optimised_allocations = try_allocs[i],
      # min_alloc = 1e-4,
      algorithm = DynamicSS(Rodas5()),
      init_guesser = concentration_guess,
      optim_alg = NelderMead(
        initial_simplex = Optim.AffineSimplexer(a = 0.0, b = -0.1)
      ),
      optim_opts = Optim.Options(
        time_limit = 5.0,
        g_tol = 1e-10,
      ),
      abstol = 1e-10,
      reltol = 1e-10,
    )
    push!(all_results,
      [growth_rate(c, p), p, c]
    )
  end
  
  mu_opt, p_opt, c_opt = argmax(x->x[1], all_results)
      
  push!(nt, (params = p_opt, state = c_opt))
end

df = nt_to_df(nt) |>
  @mutate(
    C_per_N = 14 * _.γ_K / ( 12 * (1 - _.γ_K)),
    experiment = experiment
  ) |> DataFrame


# wsave(
#     savefile,
#     df
# )

CSV.write(
  savefile,
  df,
  append = true
)