# Could anyone have been harmed?

A program evaluation reports that a job-training program raised average
earnings by \$500. A city council member asks a different question: did
the program make anyone worse off? The average cannot answer that
question. We would see an average gain of \$500 if every participant
gained \$500. We would see the same average if two thirds of the
participants gained \$1500 and one third lost \$1500.

With data from a randomized experiment,
[`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
gives the council member two numbers, both holding together with 95
percent confidence: at least this many participants were made worse off,
and at most this many could have been. This vignette shows how to
compute both and how to read them.

## What a randomized experiment can tell us about individual effects

Call each participant a unit. Each unit i has two potential outcomes:
$`Y_i(1)`$, the outcome it would have if treated, and $`Y_i(0)`$, the
outcome it would have if not treated. Its individual effect is
$`\tau_i = Y_i(1) - Y_i(0)`$. The unit was harmed if $`\tau_i < 0`$. We
see only one of the two potential outcomes for each unit. So we never
see any $`\tau_i`$. We cannot say which units were harmed.

We can still learn how many units were harmed. Sort the n individual
effects from smallest to largest,
$`\tau_{(1)} \le \tau_{(2)} \le \dots \le
\tau_{(n)}`$. Because treatment was assigned at random, we can test the
hypothesis that the k-th smallest effect is at most some number c,
written $`\tau_{(k)} \le c`$. We reject that hypothesis when the
observed data would be unlikely, under random assignment, if
$`\tau_{(k)}`$ were at most c. The largest c that we reject at the 5
percent level is a 95 percent lower confidence bound for $`\tau_{(k)}`$.

[`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
returns a lower bound for every k from 1 to n. The bounds are built so
that all n of them hold together with 95 percent confidence (see
[`?com_conf_quant_larger_cre`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_conf_quant_larger_cre.md)).
We may therefore look through the whole list and pick out the bounds
above 0. Suppose the lower bound for $`\tau_{(k)}`$ is above 0. The
effects $`\tau_{(k+1)}, \dots, \tau_{(n)}`$ are at least as large as
$`\tau_{(k)}`$. So, with 95 percent confidence, all $`n - k + 1`$
effects from $`\tau_{(k)}`$ to $`\tau_{(n)}`$ are above 0. At least
$`n - k + 1`$ units have positive effects.
[`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
reports this count for the smallest such k.

## Counting harmed units with a minus sign

Lower bounds count units with positive effects. To count harmed units we
need to show that some effects are below 0. A lower bound cannot show
that. To count them, we multiply every outcome by -1 and analyze the new
outcome, $`-Y`$, in place of $`Y`$. Each unit’s effect on $`-Y`$ is
$`(-Y_i(1)) - (-Y_i(0)) = -\tau_i`$. A harmed unit, with $`\tau_i < 0`$,
therefore has a positive effect on $`-Y`$. The count of units with
positive effects on $`-Y`$ is a count of harmed units. In
[`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
the change is a minus sign in the formula: `cmrss(-y ~ z, data = d)`.

When larger outcomes are worse, as with blood pressure or days absent,
harm means $`\tau_i > 0`$. The analysis of y itself then counts harmed
units.

## A simulated program that harms some participants

We simulate an experiment in which we set every unit’s effect, so we can
compare the answers with the truth. Of 200 units, 100 gain 2 points, 40
lose 3 points, and 60 are unaffected. Half the units are assigned to
treatment at random.

``` r

set.seed(20261004)
n <- 200
y0 <- rnorm(n)
tau <- rep(0, n)
who <- sample(n)
tau[who[1:40]] <- -3     # harmed
tau[who[41:140]] <- 2    # helped
z <- sample(rep(c(1, 0), each = n / 2))
sim <- data.frame(z = z, y = ifelse(z == 1, y0 + tau, y0))
diff_means <- mean(sim$y[sim$z == 1]) - mean(sim$y[sim$z == 0])
```

The average effect is $`(100 \times 2 - 40 \times 3) / 200 = 0.4`$
points. In this simulated experiment the treated group’s mean outcome
exceeds the control group’s by 0.51 points. A researcher who reported
only that difference would tell the council member the program helped.
The difference says nothing about the 40 units the program harmed.

We ask the data two questions: how many units were helped, and how many
were harmed. We want both answers to hold together with 95 percent
confidence. The chance that at least one of two statements is wrong is
at most the sum of the two chances. So if we compute each at 97.5
percent confidence, both hold together with at least 95 percent
confidence. In
[`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
that is `alpha = 0.025`, where `alpha` is one minus the confidence
level. To keep this page quick to build, each call below uses 1000
simulated random assignments (`nperm`) and finds each bound to within
0.05 (`tol`). A real analysis should use the defaults of 10,000
assignments and 0.01.

``` r

helped <- cmrss(y ~ z, data = sim, set = "all", alpha = 0.025,
                nperm = nperm, tol = tol)
harmed <- cmrss(-y ~ z, data = sim, set = "all", alpha = 0.025,
                nperm = nperm, tol = tol)
helped
#> cmrss: completely randomized experiment, 200 units, 100 treated
#> Polynomial scores with parameters 2, 11, 66; 1000 simulated assignments
#> With 97.5 percent confidence, at least 18 of 200 units have effects on y above 0.
harmed
#> cmrss: completely randomized experiment, 200 units, 100 treated
#> Polynomial scores with parameters 2, 11, 66; 1000 simulated assignments
#> With 97.5 percent confidence, at least 7 of 200 units have effects on -y above 0.
n_harmed <- sum(harmed$bounds$lower > 0)
n_helped <- sum(helped$bounds$lower > 0)
```

With `set = "all"`,
[`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
bounds the effects of all 200 units, treated and control. It does so by
bounding the treated units and the control units separately, each at
half of `alpha`. The other choices, `set = "treat"` and
`set = "control"`, bound one group at the full `alpha`. Analyzing `-y`,
we find 7 units with positive effects on $`-y`$, so at least 7 units
were harmed. Analyzing `y`, we find at least 18 helped units. A helped
unit was not harmed, so at most $`200 - 18 = 182`$ units could have been
harmed. These 182 units are every unit not shown to be helped. They
include the harmed units, the unaffected units, and helped units that
the data could not detect.

With 95 percent confidence, then, the number of harmed units is between
7 and 182. The true number, 40, lies in that range. The range spans most
of the 200 units, because we never see a unit’s treated and control
outcomes together. The lower end still answers the council member. At
least 7 units were made worse off, although the average effect was
positive.

Each row of `harmed$bounds` is one effect on $`-y`$. Its column `k`
gives the effect’s position in the sorted list. Its column `lower` gives
the effect’s lower bound. The number of harmed units comes from the
first row whose lower bound is above 0:

``` r

b <- harmed$bounds
first <- b[b$lower > 0, ][1, ]
first
#>       k proportion lower
#> 194 194       0.97   0.1
```

There k = 194, so the count is $`200 - 194 + 1 = 7`$ harmed units.

## The teacher professional-development experiment

The `electric_teachers` data come from an experiment with 233 elementary
school teachers in 7 sites. Within each site, the researchers assigned
teachers at random to one of three versions of a professional
development program in science or to a control group. The variable
`TxAny` is 1 for the 164 teachers assigned to any version of the program
and 0 for the 69 in the control group. The outcome `gain` is a teacher’s
score on a test of knowledge about electric circuits after the program
minus the score before it, in percentage points.

Because the researchers randomized within sites, our analysis compares
teachers only with other teachers in their own site. We tell
[`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
the sites by writing `| Site` in the formula. Computing bounds within
sites requires an optimization solver, such as the one in the `highs`
package. As in the simulation, we compute the helped and harmed counts
at `alpha = 0.025` each.

``` r

data(electric_teachers)
t_helped <- cmrss(gain ~ TxAny | Site, data = electric_teachers,
                  set = "all", alpha = 0.025, nperm = nperm, tol = tol,
                  opt.method = "ILP_highs")
#> Warning in cmrss(gain ~ TxAny | Site, data = electric_teachers, set = "all", :
#> 29 of 233 units (12 percent) share the outcome value 10. Tied outcomes are
#> ranked by their order in the data, so reordering the rows can change p-values
#> and confidence bounds (see https://github.com/davidk91919/CMRSS/issues/5).
t_harmed <- cmrss(-gain ~ TxAny | Site, data = electric_teachers,
                  set = "all", alpha = 0.025, nperm = nperm, tol = tol,
                  opt.method = "ILP_highs")
#> Warning in cmrss(-gain ~ TxAny | Site, data = electric_teachers, set = "all", :
#> 29 of 233 units (12 percent) share the outcome value -10. Tied outcomes are
#> ranked by their order in the data, so reordering the rows can change p-values
#> and confidence bounds (see https://github.com/davidk91919/CMRSS/issues/5).
t_helped
#> cmrss: block-randomized experiment, 233 units, 164 treated
#> Polynomial scores with parameters 2, 11, 59; 1000 simulated assignments
#> With 97.5 percent confidence, at least 63 of 233 units have effects on gain above 0.
t_harmed
#> cmrss: block-randomized experiment, 233 units, 164 treated
#> Polynomial scores with parameters 2, 11, 59; 1000 simulated assignments
#> With 97.5 percent confidence, at least 0 of 233 units have effects on -gain above 0.
t_n_helped <- sum(t_helped$bounds$lower > 0)
t_n_harmed <- sum(t_harmed$bounds$lower > 0)
```

Both calls warn that 29 of the 233 teachers share a single value of
`gain`.
[`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
ranks tied outcomes by their order in the data frame. With this many
ties, a different order of the rows could give somewhat different
bounds.

With 95 percent confidence, at least 63 teachers gained from the
program. Analyzing `-gain`, we cannot show that any teacher was harmed.
A failure to show harm is not evidence that no teacher was harmed. From
the analysis of `gain`, as many as 233 - 63 = 170 teachers could have
been harmed.

## Choices that change the answer

With the argument `set`, we choose whose effects are bounded: `"treat"`,
`"control"`, or `"all"`.

With `alpha`, we choose the confidence level, 1 - `alpha`. When we
report the helped and harmed counts together, we compute each at
`alpha / 2`.

With `s`, we choose the rank statistics that
[`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
combines. A rank statistic replaces each outcome by its rank and adds up
a score for the rank of every treated unit. By default
[`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
uses polynomial scores: the unit at rank r among n, counting from the
smallest outcome, gets the score $`(r/(n + 1))^{s - 1}`$. With s = 2 the
score is proportional to the rank, so every rank counts. With a large s
nearly all the score sits on the few largest outcomes. With a small s
the combined test has more power against effects shared by most units.
With a large s it has more power against large effects in a few units.
The default combines three values of s, from 2 up to a maximum set by
the numbers of treated and control units. In the simulation above they
were 2, 11, 66. The section “Choosing s” in
[`?cmrss`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
gives the rule, the simulations behind it, and the alternative
Stephenson scores (`scores = "stephenson"`).

The method is built for outcomes that take many different values. With a
binary outcome, every observed outcome is 0 or 1, so each unit is tied
with every other unit that has the same value. The order of the rows
then decides most of the ranks, and
[`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
warns. Every individual effect is also $`-1`$, $`0`$, or $`1`$, so a
bound can say little more than whether an effect is at least 0 or at
least 1. With a binary outcome there are four kinds of unit:
$`Y_i(1) = Y_i(0) = 1`$, $`Y_i(1) = Y_i(0) = 0`$, $`Y_i(1) = 1`$ and
$`Y_i(0) = 0`$ (helped), and $`Y_i(1) = 0`$ and $`Y_i(0) = 1`$ (harmed).
Methods built on the counts of these four kinds, such as those of Rigdon
and Hudgens (2015, *Statistics in Medicine* 34, <doi:10.1002/sim.6384>)
and Li and Ding (2016, *Statistics in Medicine* 35,
<doi:10.1002/sim.6924>), give exact confidence intervals for the average
effect on a binary outcome.
