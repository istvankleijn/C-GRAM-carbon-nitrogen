include(srcdir("utils.jl"))
include(srcdir("nt_to_df.jl"))

function parameter_sweep(
    optimiser,
    sweep_key::Symbol,
    model_function,
    optimisation_parameters,
    initial_parameters,
    initial_state
    ;
    sweep_values = logrange(1e-2, 1e2, 101),
    naive_weight = nothing,
    init_guesser = nothing,
    kwargs...
  )
    
  namedtuples = []
  params = initial_parameters
  state = initial_state
  for value in sweep_values
    params = merge(params, (; sweep_key => value))
      
    allocation_guess = (; (k => params[k] for k in optimisation_parameters)...)
      
    if !isnothing(naive_weight)
    naive_guess = sum(allocation_guess)/length(allocation_guess)
    allocation_guess = (; (k => (1 - naive_weight)*params[k] + naive_weight * naive_guess
                            for k in optimisation_parameters)...
                        )
    end
        
    if isnothing(init_guesser)
        init_guesser = (x -> state)
    end
    params, state = optimiser(
        model_function,
        params,
        allocation_guess;
        init_guesser = init_guesser,
        kwargs...
    )
    push!(namedtuples, (params = params, state = state))
  end

  return(nt_to_df(namedtuples))
end