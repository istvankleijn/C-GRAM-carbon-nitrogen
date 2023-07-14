using LabelledArrays
using OrdinaryDiffEq

function evolve_concentrations!(du, u, p, t)
  # Auxiliary parameters
  # Stoichiometries
  α_Nf = 1 - p.α_Cf
  α_Nr = 1 - p.α_Cr
  γ_N  = 1 - p.γ_K  

  # Fluxes
  j_Kre  = p.k_Kre * u.e_Kre  * u.k       / ( u.k + p.k_sat )
  j_Kex  = p.k_Kex * u.e_Kex  * u.k       / ( u.k + p.k_sat )
  j_C    = p.k_C   * u.e_C
  j_Af   = p.k_Ef  * u.e_Af   * u.c * u.n / ( (u.c + p.c_sat) * (u.n + p.n_sat) + p.Δ_CN )
  j_Ar   = p.k_Er  * u.e_Ar   * u.c * u.n / ( (u.c + p.c_sat) * (u.n + p.n_sat) + p.Δ_CN )
  j_N    = p.k_N   * u.e_N
  j_R    = p.k_R   * u.r      * u.a       / ( u.a + p.a_sat )

  # Growth rate of total copy number (mass)
  μ = j_C + j_N - j_Kex
  
  # Differentials
  # Keto-acid
  du.k    = - j_Kre - j_Kex +                                   p.γ_K * j_N       - μ * u.k
  # Carbon precursor
  du.c    =   j_Kre         + j_C - p.α_Cf * j_Af - p.α_Cr * j_Ar                 - μ * u.c
  # Metabolite
  du.a    =                                  j_Af          + j_Ar           - j_R - μ * u.a
  # Nitrogen precursor
  du.n    =                         - α_Nf * j_Af   - α_Nr * j_Ar + γ_N * j_N     - μ * u.n
  # Keto-acid recycling enzyme
  du.e_Kre  =                                                       p.f_Kre * j_R - μ * u.e_Kre
  # Keto-acid excretion enzyme
  du.e_Kex  =                                                       p.f_Kex * j_R - μ * u.e_Kex
  # Carbon enzyme
  du.e_C  =                                                           p.f_C * j_R - μ * u.e_C
  # Fermentative enzyme
  du.e_Af =                                                          p.f_Ef * j_R - μ * u.e_Af
  # Respiratory enzyme
  du.e_Ar =                                                          p.f_Er * j_R - μ * u.e_Ar
  # Nitrogen enzyme
  du.e_N  =                                                           p.f_N * j_R - μ * u.e_N
  # Active ribosome
  du.r    =                                                           p.f_R * j_R - μ * u.r
  # Housekeeping protein (cytoplasm)
  du.z    =                                                           p.f_Z * j_R - μ * u.z
  
  # Avoid negative derivatives for negative concentration vectors
  # (might occur due to numerical inaccuracies)
  # Required for PositiveDomain() callback (not currently working)
  # du = (u .>= 0) .* (du .< 0) .* du
  
  # Avoid allocating the final result
  nothing
end
function growth_rate(u, p)
  j_C    = p.k_C   * u.e_C    
  j_N    = p.k_N   * u.e_N
  j_Kex  = p.k_Kex * u.e_Kex  * u.k       / ( u.k + p.k_sat )    
    
  μ = j_C + j_N - j_Kex
  return μ
end

p0 = (
  k_Kre =    10.0,    # Maximal keto-acid recycling enzyme rate (h^-1)
  k_Kex =    10.0,    # Maximal keto-acid excretion enzyme rate (h^-1)
  k_C =      20.0,    # Carbon transporter rate (h^-1)
  k_Ef =     10.0,    # Maximal fermentative enzyme rate (h^-1)
  k_Er =      5.0,    # Maximal respiratory enzyme rate (h^-1)
  k_N =       5.0,    # Nitrogen transporter rate (h^-1)
  k_R =       6.46,   # Maximal ribosomal synthesis rate (h^-1)
  k_sat =     0.0167, # Keto-acid Michaelis constant for keto-acid recycling enzyme
  c_sat =     0.0167, # Carbon Michaelis constant for fermentative/respirotory enzymes
  a_sat =     0.0167, # Precursor Michaelis constant for ribosome
  n_sat =     0.0167, # Nitrogen Michaelis constant for fermentative/respiratory enzymes
  Δ_CN  =     0.0,    # Nonlinearity in fermentative/respiratory enzymes
  f_Kre =     0.05,    # Allocation towards keto-acid recycling enzymes
  f_Kex =     0.05,    # Allocation towards keto-acid excretion enzymes
  f_C   =     0.2,    # Allocation towards carbon transporters
  f_Ef  =     0.1,   # Allocation towards fermentative enzymes
  f_Er  =     0.05,    # Allocation towards respiratory enzymes
  f_N   =     0.1,   # Allocation towards nitrogen transporters
  f_R   =     0.3,    # Allocation towards ribosomes
  f_Z   =     0.2,    # Allocation towards non-metabolic proteins (Q & X)
  γ_K   =     0.0,    # Proportion of keto-acid produced by nitrogen transporters
  α_Cr  =   24//31,  # Proportion of carbon consumed by fermentative enzymes
  α_Cf  =   48//55,   # Proportion of carbon consumed by respiratory enzymes
)

function concentration_guess(parameters; k = 0.0, c = 0.05, a = 0.05, n = 0.05)
  total_protein = 1.0 - k - c - a - n
  LVector(
    k    = k,   # Keto-acid
    c    = c,   # Carbon precursor
    a    = a,   # Metabolite
    n    = n,   # Nitrogen precursor
    e_Kre  = total_protein * parameters.f_Kre,   # Keto-acid recycling enzyme
    e_Kex  = total_protein * parameters.f_Kex,   # Keto-acid excretion enzyme
    e_C  = total_protein * parameters.f_C,   # Carbon transporter
    e_Af = total_protein * parameters.f_Ef,  # Fermentative enzyme
    e_Ar = total_protein * parameters.f_Er,  # Respiratory enzyme
    e_N  = total_protein * parameters.f_N,   # Nitrogen transporter
    r    = total_protein * parameters.f_R,   # Ribosome
    z    = total_protein * (parameters.f_Z),  # Housekeeping (cytoplasm)
  )
end
