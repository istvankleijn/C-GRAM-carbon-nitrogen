using SteadyStateDiffEq

function find_steady_state(
    model_ode_function,
    initial_condition,
    model_parameters,
    algorithm::SteadyStateDiffEq.SteadyStateDiffEqAlgorithm;
    kwargs...
    )
  problem = SteadyStateProblem(
    model_ode_function,
    initial_condition,
    model_parameters
  )
  solution = solve(
    problem, algorithm;
    isoutofdomain = (u,p,t) -> any(x -> x < 0, u),  # Refuse negative concentrations/amounts
    kwargs...
  )
  return typeof(initial_condition)(solution.u), solution.retcode
end
