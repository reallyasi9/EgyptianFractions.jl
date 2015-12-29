__precompile__()

module EgyptianFractions

export efgreedy, efoddgreedy, efharmonic, engelexpand, efengel

function _prep{T<:Integer}(r::Rational{T})
  # This doesn't matter--just switch to negatives if this condition is not met.
  #@assert(s ≤ r, "sum of given factors ($s) > target ($r)")
  rem = abs(BigInt(num(r))//BigInt(den(r)))
  ef = BigInt[]
  # Later...
  #append!(ef, given)
  if abs(rem) ≥ 1
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
function efgreedy{T<:Integer}(r::Rational{T}, given::Array{T,1} = T[])
  (rem, ef) = _prep(r - sum(1 .// given))
  while rem != 0
    rem = _greedyloop!(ef, rem)
  end
  if r < 0
    ef .*= -1
  end
  if !isempty(given)
    unshift!(ef, given...)
  end
  return ef
end

efgreedy(r::Real) = efgreedy(Rational{BigInt}(r))

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
function efoddgreedy{T<:Integer}(r::Rational{T}, given::Array{T,1} = T[])
  (rem, ef) = _prep(r - sum(1 .// given))
  @assert isodd(den(rem))
  while rem != 0
    rem = _oddgreedyloop!(ef, rem)
  end
  if r < 0
    ef .*= -1
  end
  if !isempty(given)
    unshift!(ef, given...)
  end
  return ef
end

efoddgreedy(r::Real) = efoddgreedy(Rational{BigInt}(r))

"""
Performs a harmonic expansion of the given rational number into a sum of fractions of the form `1/2 + 1/3 + ...`, and concludes the remainder of the expansion using `efgreedy` if the given rational number cannot be represented as a sum in the harmonic sequence.

If the second argument is specified, the expansion will begin using that value as the first denominator.  In other words, while `efharmonic(r)` will return an array beginning `[2, 3, ...]`, `efharmonic(r, 3)` will return an array beginning `[3, 4, ...]`, and so on.

This function will throw an exception if the second argument `f ≤ 1`.

This function returns only the denominators of the expansion (i.e., only `[a_1, a_2, ...]`).  If the given rational number `r` satisfies `r ≥ 1`, the first `n` elements of the expansion will be 1, where `n = floor(r)`.  If the given rational number `r` satisfies `r < 0`, all returned denominators will be less than 0, so that the relationship

    r == sum(1 .// efharmonic(r))

is always true.
"""
function efharmonic{T<:Integer}(r::Rational{T}, first::Integer = 2, given::Array{T,1} = T[])
  @assert first ≥ 2
  (rem, ef) = _prep(r - sum(1 .// given))
  s = Rational{BigInt}(0)
  i = BigInt(first)
  while s ≤ rem
    hh = 1//i
    if s + hh > rem
      break
    end
    push!(ef, i)
    s += hh
    i += 1
  end
  rem -= s
  while rem != 0
    rem = _greedyloop!(ef, rem)
  end
  if r < 0
    ef .*= -1
  end
  if !isempty(given)
    unshift!(ef, given...)
  end
  return ef
end

efharmonic(r::Real, first::Integer = 2) = efharmonic(Rational{BigInt}(r), first, BigInt[])

function _engelloop!(ef::Vector{BigInt}, r::Rational{BigInt})
  c = ceil(1//r)
  push!(ef, c)
  return r * c - 1
end

"""
Performs an Engel expansion of the given rational number into a sum of fractions of the form `1/a + 1/(a*b) + 1/(a*b*c)...`.

This function returns only the unique denominators of the expansion (i.e., only `[a, b, ...]`).  If the given rational number `r` satisfies `r ≥ 1`, the first `n` elements of the expansion will be 1, where `n = floor(r)`.
"""
function engelexpand{T<:Integer}(r::Rational{T})
  (rem, ef) = _prep(r)
  while rem != 0
    rem = _engelloop!(ef, rem)
  end
  if r < 0
    ef .*= -1
  end
  return ef
end

engelexpand(r::Real) = engelexpand(Rational{BigInt}(r))

"""
Performs an Engel expansion of the given rational number into a sum of fractions of the form `1/a + 1/(a*b) + 1/(a*b*c)...`.

This function returns only the denominators of the expansion (i.e., only `[a, a*b, ...]`).  If the given rational number `r` satisfies `r ≥ 1`, the first `n` elements of the expansion will be 1, where `n = floor(r)`.  If the given rational number `r` satisfies `r < 0`, all returned denominators will be less than 0, so that the relationship

    r == sum(1 .// efengel(r))

is always true.
"""
function efengel{T<:Integer}(r::Rational{T}, given::Array{T, 1} = T[])
  ef = engelexpand(r - sum(1 .// given))
  cumprod!(ef, ef)
  if r < 0
    ef .*= -1
  end
  if !isempty(given)
    unshift!(ef, given...)
  end
  return ef
end

efengel(r::Real) = efengel(Rational{BigInt}(r))

end # module
