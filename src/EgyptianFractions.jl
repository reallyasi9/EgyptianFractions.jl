__precompile__()

module EgyptianFractions

export efgreedy, efoddgreedy, efharmonic, engelexpand, efengel

function _prep(r::Rational)
  @assert(den(r) != 0, "denominator must be > 0")
  rem = abs(big(r))
  ef = BigInt[]
  if rem ≥ 1
    f = floor(rem)
    rem -= f
    append!(ef, collect(repeated(1, Int(f))))
  end
  return (rem, ef)
end

function _greedyloop!(ef::Vector{BigInt}, r::Rational{BigInt})
  c = ceil(1//r)
  push!(ef, c)
  return r - 1//c
end

"""
Performs a greedy (Fibonacci–Sylvester) expansion of the given rational number into a sum of fractions of the form `1/a_1 + 1/a_2 + ...`

This function returns only the denominators of the expansion (i.e., only `[a_1, a_2, ...]`).  If the given rational number `r` satisfies `r ≥ 1`, the first `n` elements of the expansion will be 1, where `n = floor(r)`.  If the given rational number `r` satisfies `r < 0`, all returned denominators will be less than 0, so that the relationship

    r == sum(1 .// efgreedy(r))

is always true.
"""
function efgreedy(r::Rational; nmax::Int = typemax(Int))
  (rem, ef) = _prep(r)
  i = nmax
  while rem != 0 && i > 0
    rem = _greedyloop!(ef, rem)
    i -= 1
  end
  if r < 0
    ef .*= -1
  end
  return ef
end

efgreedy(r::Real; nmax::Int = typemax(Int)) = efgreedy(Rational(big(r)), nmax=nmax)

function _oddgreedyloop!(ef::Vector{BigInt}, r::Rational{BigInt})
  c = ceil(1//r)
  if iseven(num(c))
    c += 1
  end
  push!(ef, c)
  return r - 1//c
end

"""
Performs a greedy (Fibonacci–Sylvester) expansion of the given rational number into a sum of fractions of the form `1/a_1 + 1/a_2 + ...`, but only using odd denominators.

It can be shown that, for any rational number `x/y` where `y` is odd, you can write this as a finite sum of fractions where all denominators are odd.  This method sometimes produces an expansion with fewer elements than the traditional greedy algorithm.

This function will throw an exception if the denominator of the rational number supplied is even.

This function returns only the denominators of the expansion (i.e., only `[a_1, a_2, ...]`).  If the given rational number `r` satisfies `r ≥ 1`, the first `n` elements of the expansion will be 1, where `n = floor(r)`.  If the given rational number `r` satisfies `r < 0`, all returned denominators will be less than 0, so that the relationship

    r == sum(1 .// efoddgreedy(r))

is always true.
"""
function efoddgreedy(r::Rational; nmax::Int = typemax(Int))
  @assert isodd(den(r)) "denominator of rational ($(den(r))) must be odd"
  (rem, ef) = _prep(r)
  i = nmax
  while rem != 0 && i > 0
    rem = _oddgreedyloop!(ef, rem)
    i -= 1
  end
  if r < 0
    ef .*= -1
  end
  return ef
end

efoddgreedy(r::Real; nmax::Int = typemax(Int)) = efoddgreedy(Rational(big(r)), nmax=nmax)

"""
Performs a harmonic expansion of the given rational number into a sum of fractions of the form `1/2 + 1/3 + ...`, and concludes the remainder of the expansion using `efgreedy` if the given rational number cannot be represented as a sum in the harmonic sequence.

If the second argument is specified, the expansion will begin using that value as the first denominator.  In other words, while `efharmonic(r)` will return an array beginning `[2, 3, ...]`, `efharmonic(r, 3)` will return an array beginning `[3, 4, ...]`, and so on.

This function will throw an exception if the second argument `f ≤ 1`.

This function returns only the denominators of the expansion (i.e., only `[a_1, a_2, ...]`).  If the given rational number `r` satisfies `r ≥ 1`, the first `n` elements of the expansion will be 1, where `n = floor(r)`.  If the given rational number `r` satisfies `r < 0`, all returned denominators will be less than 0, so that the relationship

    r == sum(1 .// efharmonic(r))

is always true.
"""
function efharmonic(r::Rational, first::Int = 2; nmax::Int = typemax(Int))
  @assert(first ≥ 2, "harmonic series must start at 2 or greater ($first given)")
  (rem, ef) = _prep(r)
  s = Rational{BigInt}(0)
  i = big(first)
  j = nmax
  while s ≤ rem && j > 0
    hh = 1//i
    if s + hh > rem
      break
    end
    push!(ef, i)
    s += hh
    i += 1
    j -= 1
  end
  rem -= s
  while rem != 0 && j > 0
    rem = _greedyloop!(ef, rem)
    j -= 1
  end
  if r < 0
    ef .*= -1
  end
  return ef
end

efharmonic(r::Real, first::Int = 2; nmax::Int = typemax(Int)) = efharmonic(Rational(big(r)), first, nmax=nmax)

function _engelloop!(ef::Vector{BigInt}, r::Rational{BigInt})
  c = ceil(1//r)
  push!(ef, c)
  return r * c - 1
end

"""
Performs an Engel expansion of the given rational number into a sum of fractions of the form `1/a + 1/(a*b) + 1/(a*b*c)...`.

This function returns only the unique denominators of the expansion (i.e., only `[a, b, ...]`).  If the given rational number `r` satisfies `r ≥ 1`, the first `n` elements of the expansion will be 1, where `n = floor(r)`.
"""
function engelexpand(r::Rational; nmax::Int = typemax(Int))
  (rem, ef) = _prep(r)
  i = nmax
  while rem != 0 && i > 0
    rem = _engelloop!(ef, rem)
    i -= 1
  end
  # Convention:  the Engle expansion of a negative number will have ef[1] < 0.
  # This way, all elements of cumprod(ef) will be negative.
  if !isempty(ef) && r < 0
    ef[1] *= -1
  end
  return ef
end

engelexpand(r::Real; nmax::Int = typemax(Int)) = engelexpand(Rational(big(r)), nmax=nmax)

"""
Performs an Engel expansion of the given rational number into a sum of fractions of the form `1/a + 1/(a*b) + 1/(a*b*c)...`.

This function returns only the denominators of the expansion (i.e., only `[a, a*b, ...]`).  If the given rational number `r` satisfies `r ≥ 1`, the first `n` elements of the expansion will be 1, where `n = floor(r)`.  If the given rational number `r` satisfies `r < 0`, all returned denominators will be less than 0, so that the relationship

    r == sum(1 .// efengel(r))

is always true.
"""
function efengel(r::Rational; nmax::Int = typemax(Int))
  ef = engelexpand(r, nmax=nmax)
  cumprod!(ef, ef)
  return ef
end

efengel(r::Real; nmax::Int = typemax(Int)) = efengel(Rational(big(r)), nmax=nmax)

end # module
