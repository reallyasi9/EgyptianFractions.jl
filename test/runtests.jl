using EgyptianFractions
using Base.Test

# Random examples are my favorite!
const seed = 7182818284 # first 10 digits of e after the decimal
const rng = MersenneTwister(seed)

const r32 = rand(rng, Int32, 100, 2)

for i in 1:size(r32, 1)
  n = r32[i, 1]
  d = r32[i, 2]
  efg = efgreedy(n//d)
  @test sum(1.//efg) == n//d
end

for i in 1:size(r32, 1)
  n = r32[i, 1]
  d = r32[i, 2]
  if !isodd(d)
    d += 1
  end
  efog = efoddgreedy(n//d)
  @test sum(1.//efog) == n//d
end

for i in 1:size(r32, 1)
  n = r32[i, 1]
  d = r32[i, 2]
  efh = efharmonic(n//d)
  @test sum(1.//efh) == n//d
end

# Engel expansions are problematic, because they get really large really fast.

# Examples from https://en.wikipedia.org/wiki/Egyptian_fraction

# Comparing the size of some fractions
@test efgreedy(4//5) == [2, 4, 20]
@test efgreedy(3//4) == [2, 4]
# Thus
@test 4//5 - 3//4 == 1//20

@test efgreedy(3//11) == [4, 44]
@test efgreedy(2//7) == [4, 28]
# Thus, because
@test 1//44 < 1//28
# therefore
@test 3//11 < 2//7

# Equally distributing objects
@test efgreedy(5//8) == [2, 8]
@test efgreedy(13//12) == [1, 12] # by definition, however
@test efgreedy(13//12, [2, 3]) == [2, 3, 4]

# Later usage
@test efgreedy(5//121) == [25, 757, 763309, 873960180913, 1527612795642093418846225]
@test efgreedy(5//121, [33, 121]) == [33, 121, 363]

# Examples from https://en.wikipedia.org/wiki/Odd_greedy_expansion

@test efoddgreedy(4//23) == [7, 33, 1329, 2353659]

# Fractions with long expansions
@test efgreedy(8//77) == [10, 257, 197890]
@test efoddgreedy(8//77) == [11, 77]
@test sum(1 .// efgreedy(8//77)) == sum(1 .// efoddgreedy(8//77)) == 8//77

let og = efoddgreedy(3//179)
  @test size(og, 1) == 19
  # Dag, yo!
  #@test_approx_eq_eps(og[end], 1.415e439491, 1e-2)
end

# Engel expansions are a little different.
# From https://en.wikipedia.org/wiki/Engel_expansion
@test engelexpand(1175//1000) == [1, 6, 20]
@test efengel(1175//1000) == cumprod([1, 6, 20])

# Engel expansions for some well-known constants
let enpi = engelexpand(π)
  @test enpi[1:9] == [1, 1, 1, 8, 8, 17, 19, 300, 1991]
end

let ensqrt2 = engelexpand(sqrt(2))
  @test ensqrt2[1:9] == [1, 3, 5, 5, 16, 18, 78, 102, 120]
end

let ene = engelexpand(exp(1))
  @test ene[1:9] == [1, 1, 2, 3, 4, 5, 6, 7, 8]
end
