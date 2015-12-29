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

# Greedy expansions from OEIS
let e = efgreedy(π-3)
  # A001466
  # Numbers get too big, errors compound
  @test e[1:3] == [8, 61, 5020, #=128541455=#]
  @test_approx_eq_eps(sum(1 .// e[4:end]), 1//128541455, 1e-12)
end

let e = efgreedy(sqrt(2))
  # A006487
  # Numbers get too big, errors compound
  @test e[1:5] == [1, 3, 13, 253, 218201, #=61323543802=#]
  @test_approx_eq_eps(sum(1 .// e[6:end]), 1//61323543802, 1e-12)
end

let e = efgreedy(1/π)
  # A006524
  # Numbers get too big, errors compound
  @test e[1:4] == [4, 15, 609, 845029, #=1010073215739=#]
  @test_approx_eq_eps(sum(1 .// e[5:end]), 1//1010073215739, 1e-12)
end

let e = efgreedy(eu - 2)
  # A006525
  # Numbers get too big, errors compound
  @test e[1:4] == [2, 5, 55, 9999, #=3620211523=#]
  @test_approx_eq_eps(sum(1 .// e[5:end]), 1//3620211523, 1e-12)
end

let e = efgreedy(exp(-1))
  # A006526
  # Numbers get too big, errors compound
  @test e[1:3] == [3, 29, 15786, #=513429610, 339840390654894740=#]
  @test_approx_eq_eps(sum(1 .// e[4:end]), 1//513429610, 1e-12)
end

# Engel expansions are a little different.
# From https://en.wikipedia.org/wiki/Engel_expansion
@test engelexpand(1175//1000) == [1, 6, 20]
@test efengel(1175//1000) == cumprod([1, 6, 20])

# Engel expansions for some well-known constants
let e = engelexpand(π)
  # A006784
  @test e[1:9] == [1, 1, 1, 8, 8, 17, 19, 300, 1991]
end

let e = engelexpand(1/π)
  # A014012
  # Numbers get too big, errors compound
  @test e[1:6] == [4, 4, 11, 45, 70, 1111, #=4423, 5478, 49340=#]
end

let e = engelexpand(sqrt(2))
  # A028254
  @test e[1:9] == [1, 3, 5, 5, 16, 18, 78, 102, 120]
end

let e = engelexpand(sqrt(3))
  # A028257
  @test e[1:9] == [1, 2, 3, 3, 6, 17, 23, 25, 27]
end

let e = engelexpand(sqrt(5))
  # A059176
  @test e[1:9] == [1, 1, 5, 6, 13, 16, 16, 38, 48]
end

let e = engelexpand(sqrt(10))
  # A059177
  @test e[1:9] == [1, 1, 1, 7, 8, 12, 20, 86, 94]
end

let e = engelexpand(golden)
  # A028259
  @test e[1:9] == [1, 2, 5, 6, 13, 16, 16, 38, 48]
end

let e = engelexpand(eulergamma)
  # A059177
  # Numbers get too big, errors compound
  @test e[1:7] == [2, 7, 13, 19, 85, 2601, 9602, #=46268, 4812284=#]
end

let e = engelexpand(2^(1/3))
  # A059178
  # Numbers get too big, errors compound
  @test e[1:7] == [1, 4, 26, 32, 58, 1361, 4767, #=22303, 134563=#]
end

let e = engelexpand(3^(1/3))
  # A059179
  @test e[1:9] == [1, 3, 4, 4, 5, 8, 9, 14, 63]
end

let e = engelexpand(log(2))
  # A059180
  # Numbers get too big, errors compound
  @test e[1:7] == [2, 3, 7, 9, 104, 510, 1413, #=2386, 40447=#]
end

let e = engelexpand(log(3))
  # A059181
  # Numbers get too big, errors compound
  @test e[1:7] == [1, 11, 12, 60, 108, 139, 176, #=1228, 1356=#]
end

let e = engelexpand(log(10))
  # A059182
  @test e[1:9] == [1, 1, 4, 5, 20, 30, 48, 74, 265]
end

let e = engelexpand(1/log(2))
  # A059183
  @test e[1:9] == [1, 3, 4, 4, 5, 5, 5, 6, 47]
end

let e = engelexpand(1/log(10))
  # A059184
  # Numbers get too big, errors compound
  @test e[1:8] == [3, 4, 5, 18, 27, 37, 415, 445, #=1812=#]
end

let e = engelexpand(π^2)
  # A059185
  # Expansion gets too small, errors compound
  @test e[1:20] == [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 9, 28, 45, 72, 111, #=329, 415, 846, 1488=#]
end

let e = engelexpand(zeta(2))
  # A059186
  # Numbers get too big, errors compound
  @test e[1:8] == [1, 2, 4, 7, 9, 22, 35, 79, #=2992=#]
end

let e = engelexpand(exp(π))
  # A059196
  # Expansion gets too small, errors compound
  @test e[1:28] == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 232, 238, 428, #=1103=#]
end

let e = engelexpand(π^(eu))
  # A059197
  # Expansion gets too small, errors compound
  @test e[1:29] == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 8, 17, 111, 236, 419, #=2475=#]
end

let e = engelexpand(exp(eulergamma))
  # A059199
  @test e[1:9] == [1, 2, 2, 9, 9, 15, 84, 256, 278]
end

let e = engelexpand(-log(log(2)))
  # A059200
  # Numbers get too big, errors compound
  @test e[1:8] == [3, 11, 11, 23, 62, 66, 466, 1450, #=7617=#]
end

let e = engelexpand(catalan)
  # A054543
  @test e[1:9] == [2, 2, 2, 4, 4, 5, 5, 12, 13]
end

const khintchine = BigFloat(2.685452001065306445309714835481795693820382293994462953051152345557218)
let e = engelexpand(khintchine)
  # A054544
  @test e[1:9] == [1, 1, 2, 3, 9, 70, 117, 503, 648]
end

let e = engelexpand(sqrt(π))
  # A059187
  # Numbers get too big, errors compound
  @test e[1:8] == [1, 2, 2, 12, 13, 90, 121, 3457, #=7372=#]
end

let e = engelexpand(zeta(3))
  # A053980
  # Numbers get too big, errors compound
  @test e[1:5] == [1, 5, 98, 127, 923, #=5474, 16490, 25355, 37910=#]
end

let e = engelexpand(gamma(1/3))
  # A059188
  # Numbers get too big, errors compound
  @test e[1:8] == [1, 1, 2, 3, 14, 33, 57, 236, #=6280=#]
end

let e = engelexpand(gamma(2/3))
  # A059189
  # Numbers get too big, errors compound
  @test e[1:7] == [1, 3, 17, 17, 50, 79, 796, #=3687, 7074=#]
end

let e = engelexpand(eulergamma^2)
  # A059190
  @test e[1:9] == [4, 4, 4, 4, 4, 6, 23, 26, 126]
end

let e = engelexpand(1/eulergamma)
  # A059191
  @test e[1:9] == [1, 2, 3, 3, 6, 10, 20, 46, 226]
end

let e = engelexpand(log(1/eulergamma))
  # A059192
  # Numbers get too big, errors compound
  @test e[1:7] == [2, 11, 12, 13, 53, 348, 5263, #=9960, 17040=#]
end

let e = engelexpand(exp(-1))
  # A059193
  # Expansion gets too small, errors compound
  @test e[1:8] == [3, 10, 28, 54, 88, 130, 180, 238, #=304=#]
end

let e = engelexpand(exp(-2))
  # A059194
  # Expansion gets too small, errors compound
  @test e[1:7] == [8, 13, 14, 21, 87, 92, 119, #=444, 472=#]
end

let e = engelexpand(log(π))
  # A059195
  # Numbers get too big, errors compound
  @test e[1:6] == [1, 7, 77, 107, 150, 167, #=7091, 27852, 31790=#]
end

let e = engelexpand(eu)
  # A028310
  @test e[1:9] == [1, 1, 2, 3, 4, 5, 6, 7, 8]
end

let e = engelexpand(exp(1/2))
  # A004277
  @test e[1:9] == [1, 2, 4, 6, 8, 10, 12, 14, 16]
end
