using DataFrames
using Query

function nt_to_df(namedtuples)
  floats = @from row in namedtuples begin
    @select merge(convert(NamedTuple, row.state), row.params)
    @collect DataFrame
  end

  df = hcat(floats, DataFrame(namedtuples)) |>
    @mutate(
        mu = growth_rate(_.state, _.params)
    ) |>
    @select(-:params, -:state) |>
    DataFrame

  return(df)
end