include(srcdir("find_steady_state.jl"))

using Optim
using OrdinaryDiffEq
using SteadyStateDiffEq

function optimise_allocation(
  model_fn,
  model_parameters,
  allocation_guess;
  algorithm = DynamicSS(Rodas5()),
  init_guesser = nothing,
  na_val = 1/eps(),
  optim_alg = NelderMead(initial_simplex = Optim.AffineSimplexer(a = 0.0)),
  optim_opts = Optim.Options(),
  restrictor = nothing,
  min_alloc = 0.0,
  kwargs...
  )
@assert !isnothing(init_guesser)

AllocParameterType = typeof(allocation_guess)
n = length(allocation_guess)

active_fraction = sum(model_parameters[i] for i in keys(allocation_guess))
@assert active_fraction ≈ sum(allocation_guess) "Active fractions must agree"
function objective(x)
  if sum(x) > active_fraction || any(x .< min_alloc)
    return na_val
  end    
  modified_pars = AllocParameterType((x..., active_fraction - sum(x)))        
  pars = merge(model_parameters, modified_pars)
  c0 = init_guesser(pars)
  try  # slow when caught!
    ss, rc = find_steady_state(model_fn, c0, pars, algorithm; kwargs...)
    rc == :Success && all(ss .>= 0) || (@warn(rc, pars, ss); return na_val)
  
          
    mu = growth_rate(ss, pars)      
    if !isnothing(restrictor)
      mu *= restrictor(model_parameters)
    end
          
    return min(log(2) / mu, na_val)
      catch error
          rethrow(error)
    return na_val
  end
end

optim = Optim.optimize(
  objective,
  collect(allocation_guess)[1:end-1],
  optim_alg,
  optim_opts
)

optimizer = Optim.minimizer(optim)
modified_parameters = AllocParameterType((optimizer..., active_fraction - sum(optimizer)))
optimised_parameters = merge(model_parameters, modified_parameters)
init_guess = init_guesser(optimised_parameters)
optimised_steady_state, _ = find_steady_state(
  model_fn,
  init_guess,
  optimised_parameters,
  algorithm;
  kwargs...
  )
return optimised_parameters, optimised_steady_state
end


function optimise_allocation(
    model_fn,
    model_parameters;
    optimised_allocations = [:f_R],
    kwargs...
    )
  # Use model parameters to define initial allocation guess
  allocation_guess = NamedTuple(i => model_parameters[i]
                                for i in optimised_allocations
                                )

  return optimise_allocation(
    model_fn,
    model_parameters,
    allocation_guess;
    kwargs...
  )
end