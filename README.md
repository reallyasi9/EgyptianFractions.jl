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

Egyptian fractions are called that because that is how the ancient Egyptians wrote down fractions that were more complex than the simple fractions supported by their hieroglyph script.  Ancient Egyptian hieroglyphs had a way to denote common fractions like 1/2, 2/3, and 3/4, but relied on sums of integer reciprocals for uncommon fractions.

Egyptian fractions follow these simple rules:

1. Each fraction in the sum must be the reciprocal of a positive integer (*i.e.*, *1/n*, where *n* is an integer > 0).
2. No denominator can appear more than once.

Because of this, one can specify an Egyptian fraction as a list of denominators *n_1*, *n_2*, *...*, where each *n* is an integer, and every *n* will be unique.

## Extended Egyptian fractions

Proper Egyptian fractions (those that follow the two rules above) can only represent rational numbers strictly greater than 0.  If you want to extend the concept of Egyptian fractions to represent all rational numbers, you need to change the rules a little:

1. Each fraction in the sum must be the reciprocal of ~~a positive~~ **any** integer (*i.e.*, *1/n*, where *n* is an integer ~~> 0~~).
2. No denominator **except 1** can appear more than once.
3. The Egyptian fraction representation of 0 is the empty set.

The modification of the first rule means you can represent rational numbers < 0 by simply using all negative denominators.  The sum of these (negative) fractions will equal the original (negative) rational number.  The modification of the second rule means you can represent rational numbers > 1 or < -1 by repeating *1/1* for each integral unit of the original rational number, then representing the remainder as a proper Egyptian fraction.  For example, if

    3//4 == 1//2 + 1//4

then a demonstration of the first modified rule would be

    -3//4 == (-1//2) + (-1//4)

and a demonstration of the second modified rule would be

    11//4 == 1//1 + 1//1 + 1//2 + 1//4

From these extended rules, any rational number can be represented as an Egyptian-style fraction.

It should be noted here that, because the harmonic series *1 + 1/2 + 1/3 + ...* does not converge, it is possible to represent any rational number > 0 using the definition of the proper Egyptian fractions; however, the harmonic series diverges extremely slowly, so the modified rule 2 is adopted in this package for the sake of making runtimes more reasonable.

## Usage

To use this package, simply call from Julia:

```julia
using EgyptianFractions
```

There are several methods for generating Egyptian fractions from a rational number.  The simplest is the "[greedy](https://en.wikipedia.org/wiki/Greedy_algorithm_for_Egyptian_fractions)" algorithm:

```julia
efgreedy(3//4)  # [2, 4]
```

All of the exported functions act on `Rational`s and return a `Vector{BigInt}`.  The returned integers represent the denominators of the Egyptian fraction expansion of the given rational, as demonstrated:

```julia
sum(1 .// efgreedy(3//4)) == 3//4  # true
```

You can also pass any `Real`, and the function will simply convert the number into a `Rational` as best it can:
```julia
efgreedy(0.75)  # [2, 4]
```

**Be careful using the `Real` versions of these functions**:  Julia's conversion from `Real` to `Rational` may not be completely accurate due to how floating point numbers are represented by your CPU architecture, and you may end up with results that are not quite sensible, as demonstrated:

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

This algorithm chooses greedily the largest possible unit fraction that can be used in any representation of the remaining fraction, then recurses on the remaining fraction.  This algorithm tends to produce expansions that have many terms with very large denominators.  For instance:

```julia
efgreedy(5//121)  # [25, 757, 763309, 873960180913, 1527612795642093418846225]
```

When, as seen before, `[33, 121, 163]` would be a much more sensible expansion.

### Odd greedy expansion

This algorithm acts identically to the greedy algorithm, but it greedily chooses the largest possible odd unit fraction.  It only works if the (reduced) denominator of the input rational is odd.  If it were applied to a rational number with an even denominator, the resulting expansion would never converge.  This algorithm may produce expansions that are shorter than those produced using the generic greedy expansion.  For example:

```julia
efgreedy(8//77)     # [10, 257, 197890]
efoddgreedy(8//77)  # [11, 77]
```

However, odd greedy expansions are typically longer than generic greedy expansions, so are not typically useful except in rare situations.

### Harmonic/greedy expansion

Because the harmonic sum *1 + 1/2 + 1/3 + ...* diverges, it is always possible to find an integer *N* for any rational number *r* such that the harmonic sum up to *1/N* is less than or equal to *r*, and the harmonic sum up to *1/(N+1)* is greater than *r*.  In fact, one can start the harmonic series from any point, and that partial harmonic sum to infinity will also diverge.

Given these properties of the harmonic sum, the harmonic/greedy expansion will expand a rational number in a harmonic sum starting from some integer denominator (typically > 2) until the next fraction in the harmonic series would make the partial sum greater than the given rational number, then expands the remainder using generic greedy expansion.  An example may be helpful:

```julia
efgreedy(18//23)       # [2, 4, 31, 2852]
efharmonic(18//23)     # [2, 4, 31, 2852]
efharmonic(18//23, 5)  # [5, 6, 7, 8, 9, 28, 794, 23010120]
#                  ^ Start from the 5th term in the harmonic series
```

The second argument given to the function `efharmonic` tells the algorithm where to begin the harmonic series, which defaults to 2.  Note that:

```julia
if isa(r, Rational)
  efgreedy(r) == efharmonic(r) == efharmonic(r, 2)  # true
end
```

This is because *1/2* is the largest fraction that can ever be chosen for the greedy algorithm, followed by *1/3*, then *1/4*, and so on.

### Engel expansion

The Engel expansion of a rational number is an Egyptian fraction, but with the fractions in the form *1/n_1 + 1/(n_1 * n_2) + 1/(n_1 * n_2 * n_3) + ...*.  Because the denominators are the cumulative products of distinct integers, this type of expansion is sometimes called an "Egyptian product".

This package supplies two functions for representing Engel expansions.  The first is the more typical way of representing Engel expansions in the literature:  instead of returning the full denominators `[n_1, n_1 * n_2, n_1 * n_2 * n_3, ...]`, the function `engelexpand` simply returns the unique *n_i* values.  For instance:

```julia
engelexpand(7//40)  # [6, 20]
```

This means that the full series expansion is `1//6 + 1//(6 * 20)`.  The full series expansion is returned using the `efengel` function, as shown:

```julia
efengel(7//40)  # [6, 120]
```

It should be noted that the denominators in an Engel expansion can quickly become very large, so it is wise to use `engelexpand`.  Note also the relationship:

```julia
if isa(r, Rational)
  efengel(r) == cumprod(engelexpand(r))  # true
end
```

## Terminating expansions

All of the exported functions in this package take an optional named argument `nmax::Int` which, if specified, will terminate the expansion after `nmax` terms.  This is useful if you want an approximation of a `Real` argument that you know cannot be accurately represented using a floating point number.  For instance:

```julia
efgreedy(π-3)          # A length-14 vector with the last element approximately 2.24e+9107
efgreedy(π-3, nmax=4)  # [8, 61, 5020, 128541457]
sum(1 .// efgreedy(π-3)) - sum(1 .// efgreedy(π-3, nmax=4))
                       # Approximately 4.72e-18.  Not bad!
```

## To do

In no particular order, these are the enhancements I would like to make to this package:

- Add an algorithm that attempts to follow the deduced historical rules of Egyptian fractions (as determined from various historical texts).
- Make the modified version of rule 2 optional (add a flag or another set of named functions).
  - This could lead to extremely large expansions for rational numbers > 1, due to how slowly the harmonic series diverges.
- Add an `nmin` argument, using the various term-generation rules to construct larger expansions when needed.
- Fix `nmax` so that it correctly counts the `1//1` terms that are added for rational numbers > 1.

Feel free to suggest other ideas using the [issues](https://github.com/reallyasi9/EgyptianFractions.jl/issues) link or by creating a [pull request](https://github.com/reallyasi9/EgyptianFractions.jl/pulls).
