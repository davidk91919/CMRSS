# Add a thread cap to a Gurobi parameter list

Gurobi's default is to use every core on the machine for one solve. That
is the right default for a caller solving one problem at a time, and the
wrong one for a caller that is already running in parallel: a `foreach`
loop over 14 workers, each solving with 14 threads, puts roughly 150
runnable threads on 14 cores and the workers queue behind each other
instead of solving.

## Usage

``` r
add_thread_param(params)
```

## Arguments

- params:

  A list of Gurobi parameters.

## Value

`params`, with `Threads` added when the option is set.

## Details

Only the caller knows whether it is already parallel, so the thread
count is the caller's to set, through
`options(CMRSS.gurobi.threads = 1)`. Forked workers inherit R options,
so setting it once before a parallel loop reaches every worker. Unset,
which is the default, this returns `params` untouched and Gurobi chooses
as it always has.

Threads is a resource setting rather than a modelling one, so capping it
does not change the solution; `tests/testthat/test-gurobi-threads.R`
checks that on a four-block problem.
