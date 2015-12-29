# EgyptianFractions

Linux/Mac: [![Build Status](https://travis-ci.org/reallyasi9/EgyptianFractions.jl.svg?branch=master)](https://travis-ci.org/reallyasi9/EgyptianFractions.jl)

Windows: [![Build status](https://ci.appveyor.com/api/projects/status/0mi0m282d5rbu2p0?svg=true)](https://ci.appveyor.com/project/reallyasi9/egyptianfractions-jl)

Everything else: ![Build status](https://img.shields.io/badge/test%20it-in%20production-brightgreen.svg)

This package will decompose rational numbers into [Egyptian fractions](https://en.wikipedia.org/wiki/Egyptian_fraction).

## What the heck are Egyptian fractions?

Egyptian fractions are representations of rational numbers as sums of fractions where the numerator is always 1.  Here is a quick example (using Julia's rational number representation):

    3//4 == 1//2 + 1//4

That's pretty easy to check.  A more complex example is something like this:

    5//121 == 1//33 + 1//121 + 1//363

That result is not particularly obvious, but if you do the math, it turns out to be correct:

    (5*3)//(121*3) == 15//363 == (1*11)//(33*11) + (1*3)//(121*3) + 1//363

Egyptian fractions are called that because that is how the ancient Egyptians wrote down fractions that were more complex than the simple fractions supported by their hieroglyph script.  Ancient Egyptian hieroglyphs had a way to denote common fractions like 1/2, 2/3, and 3/4, but relied on sums of integer recipricals for uncommon fractions.

Egyptian fractions follow these simple rules:

1. Each fraction in the sum must be the reciprical of a positive integer (*i.e.*, *1/n*, where *n* is an integer > 0).
2. No denominator can appear more than once.

Because of this, one can specify an Egyptian fraction as a list of denominators *n_1*, *n_2*, *...*, where each *n* is an integer, and every *n* will be unique.

## Extended Egyptian fractions

Proper Egyptian fractions (those that follow the two rules above) can only really represent rational numbers greater than 0.  If you want to extend the concept of Egyptian fractions to represent all rational numbers, you need to change the rules a little:

1. Each fraction in the sum must be the reciprical of ~~a positive~~ **any** integer (*i.e.*, *1/n*, where *n* is an integer ~~> 0~~).
2. No denominator **except 1** can appear more than once.
3. The Egyptian fraction representation of 0 is the empty set.

The modification of the first rule means you can represent rational numbers < 0 by simply using all negative denominators.  The sum of these (negative) fractions will equal the original (negative) rational number.  The modification of the second rule means you can represent rational numbers > 1 or < -1 by repeating *1/1* for each integral unit of the original rational number, then representing the remainder as a proper Egyptian fraction.  For example, if

    3//4 == 1//2 + 1//4

then a demonstration of the first modified rule would be

    -3//4 == (-1//2) + (-1//4)

and a demonstration of the second modified rule would be

    11//4 == 1//1 + 1//1 + 1//2 + 1//4

From these extended rules, any rational number can be represented as an Egyptian-style fraction.

It should be noted here that, because the harmonic series *1 + 1/2 + 1/3 + ...* does not converge, it is possible to represent any rational number > 0 using the definition of the proper Egyptian fractions; however, the harmonic series dinverges extremely slowly, so the modified rule 2 is adopted in this package for the sake of making runtimes more reasonable.

## Usage

To use this package, simply call from Julia:

```julia
using EgyptianFractions
```

There are several methods for generating Egyptian fractions from a rational number.  The simplest is the "[greedy](https://en.wikipedia.org/wiki/Greedy_algorithm_for_Egyptian_fractions)" algorithm:

```julia
efgreedy(3//4)  # [2, 4]
```

All of the exported functions act on `Rational`s and return a `Vector{BigInt}`.  These integers represent the denominators of the Egyptian fraction expansion of the given rational, as demonstrated:

```julia
sum(1 .// efgreedy(3//4)) == 3//4  # true
```

You can also pass any `Real`, and the function will simply convert the number into a Rational as best it can:
```julia
efgreedy(0.75)  # [2, 4]
```

**Be careful using the `Real` versions of these functions**:  Julia's conversion from `Real` to `Rational` may not be completely accurate due to how floating point numbers are represented by your CPU archetecture, and you may end up with results that are not quite sensible, as demonstrated:

```julia
efgreedy(2//3)  # [2, 6]
efgreedy(2/3)   # A 16-element vector, the last entry approximately 4.47e+14155
```

**You have been warned.**

Expansions of negative numbers work:

```julia
efgreedy(-3//4)  # [-2, -4]
```

So do expansions of numbers > 1 or < -1:

```julia
efgreedy(11//4)  # [1, 1, 2, 4]
```

The expansions are done in such a way that the following identity should always hold:

```julia
if isa(x, Rational)
  sum(1 .// efgreedy(x)) == x  # true
end
```

## Types of expansions included

As mentioned before, there are several ways to generate Egyptian fractions.  Besides the "greedy" algorithm, this package also includes an "[odd greedy](https://en.wikipedia.org/wiki/Odd_greedy_expansion)" algorithm, a "harmonic/greedy" algorithm, and one that performs an [Engel expansion](https://en.wikipedia.org/wiki/Engel_expansion) (also known as an "Egyptian product").

### Greedy expansion

This algorithm chooses greedily the largest possible unit fraction that can be used in any representation of the remaining fraction, then recurses on the remaining fraction.

### Odd greedy expansion

### Harmonic/greedy expansion

### Engel expansion

## Todo
