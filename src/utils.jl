logrange(min, max, length) = exp.(range(log(min), stop=log(max), length=length))

"""
    circumference(S, V)

Calculate the circumference of a spherocylinder with given surface and volume
by solving the characteristic cubic equation.
For a spherocylinder, the equation has two positive solutions and one negative,
the physical solution is the smallest positive solution.
"""
function circumference(S, V)
  Cpolynomial = Polynomial([12*pi^2*V, -3*pi*S, 0, 1])
  Croots = roots(Cpolynomial)
  C = real(Croots[2])
end
