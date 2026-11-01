# SNCSTM simulation

[Pipeline.Rmd](Pipeline.Rmd) demonstrates simulation and estimation under
sequential ignorability (SI) and instrumental-variable (IV) settings with two
competing causes. It has two parts: a single simulation for each method, followed
by an optional Monte Carlo simulation that produces the SI/IV comparison table.

## Getting started

Open the project with `SNCSTM-simulation/` as the working directory. Install the
R packages used by the pipeline if needed:

```r
install.packages(c("survival", "knitr", "rmarkdown"))
```

Run the Rmd chunks in order, or render the single-simulation example as HTML:

```r
rmarkdown::render("Pipeline.Rmd")
```

HTML rendering also requires Pandoc. The Monte Carlo chunks have `eval=FALSE`,
so rendering the Rmd runs only the single-simulation section.

## 1. Run one simulation

The `setup` chunk sets:

```r
n <- 10000
seed <- 1
event_rate <- 0.30
```

- `n`: number of individuals in each simulated dataset.
- `seed`: random seed for the single-simulation example.
- `event_rate`: event-rate setting; supported values are 0.07, 0.15, 0.30, and 0.60.

For each method, the pipeline generates data, fits nuisance models once, and
estimates the effects for causes 1 and 2 using `n_grid = 10`.

| Step | SI | IV |
|---|---|---|
| Generate data | `simulation_si()` | `simulation_iv()` |
| Fit nuisance models | `fit_si_nuisance()` | `fit_iv_nuisance()` |
| Estimate both causes | `estimation_si()` | `estimation_iv()` |

Each cause has four targets: cumulative effects B(0.5), B(1), B(1.5), and pooled
beta. The cumulative effects are evaluated at the left limits of these times.
The Rmd displays `si_results` and `iv_results`, each with eight rows:

| Column | Meaning |
|---|---|
| `cause` | Cause 1 or 2 |
| `parameter` | One of the four targets |
| `estimate` | Estimated effect |
| `std` | Analytical standard error |
| `lower`, `upper` | Pointwise 95% confidence interval limits |
| `truth` | True effect under the data-generating model |

SI results are taken from `fit$B_summary` and `fit$pooled`; IV results are taken
from `fit$summary`. The Rmd arranges these existing outputs for display.

## 2. Run Monte Carlo simulations

After running `setup`, go to the Monte Carlo section. In `mc-settings`, the
current settings are:

```r
B <- 5
# B <- 1000
sample_sizes <- c(1000, 5000, 10000)
```

`B` controls the number of replications **per method and sample size**. Each
combination uses seeds `1:B`, matching the original seed sequence. The `seed`
used in the single-simulation example does not change this sequence. The event
rate is taken from `setup` (currently 30%).

With `B = 5`, the pipeline runs 2 methods × 3 sample sizes × 5 seeds = **30
simulated datasets**. To test only the sample size from `setup`, replace the
sample-size line with `sample_sizes <- n`.

To use the original 1000 replications, uncomment `B <- 1000`; that assignment
then overrides `B <- 5`. With all three sample sizes and event rate 30%, this
uses the methods, targets, seeds, and settings of the manuscript simulation
table, except that all analytical variances now use the uncentered IF formula
below. SEs and coverage can therefore differ from the historical table. The existing scripts support integer values of `B` from 1 to 1000.

Select and run the code in these four chunks manually, in order:

| Chunk | What it does |
|---|---|
| `mc-settings` | Sets the replication count, seeds, sample sizes, and output paths |
| `mc-run` | Runs each seed sequentially and saves its estimates |
| `mc-merge` | Combines saved replications and calculates summary statistics |
| `mc-table` | Displays success counts and the SI/IV comparison table |

The run and merge chunks call the existing scripts in `scripts/`. Each run
starts a fresh R process. Rerunning `mc-run` recomputes matching saved runs;
files with different code or settings cause an error. The merge step requires
saved files for all seeds `1:B` and does not run missing simulations.

## 3. Read the comparison table

`mc_table` places SI and IV side by side. For each sample size and cause, it
shows four metrics for the three cumulative effects and pooled beta:

| Metric | Definition |
|---|---|
| Perc. Bias | `100 * abs(mean(estimate - truth) / truth)` |
| MCSD | Standard deviation of the estimates across replications |
| Ave. SE | Mean analytical standard error across replications |
| 95% CP | Percentage of intervals containing the true effect |

The saved summaries call signed mean bias `Bias`, MCSD `SEE`, average SE `SD`,
and coverage percentage `CP`. The Rmd converts signed bias to absolute relative
percent bias for display. Percent bias uses two decimal places, MCSD and average
SE use three, and coverage uses one.

Summaries use available finite results. Check `successful_estimates` and
`successful_intervals`, displayed before the table. Five replications are a
workflow check; use the full replication count for the simulation study.

## Saved results and files

Monte Carlo results are written under `results/`. For example, at 30% and
`n = 10000`:

```text
results/
├── results_si_cumulative_30_n10000/
│   ├── simulation_0001.rds
│   ├── ...
│   └── si_cumulative_30_n10000.rds
└── results_iv_table1_30_n10000/
    ├── simulation_0001.rds
    ├── ...
    └── iv_table1_30_n10000.rds
```

Each individual file records the seed, settings, code fingerprints, and eight
target results. Each merged archive contains `settings`, `runs`, `summary`, and
`table1`. The formatted SI/IV comparison is displayed by the Rmd as `mc_table`;
it is not automatically exported as a separate table file.

| File | Purpose |
|---|---|
| [Pipeline.Rmd](Pipeline.Rmd) | Single-simulation examples and optional Monte Carlo workflow |
| [R/SI.R](R/SI.R), [R/IV.R](R/IV.R) | Data generators, estimators, and their required helpers |
| [scripts/run_si_cumulative.R](scripts/run_si_cumulative.R), [scripts/run_iv_table1.R](scripts/run_iv_table1.R) | Run and save one seed |
| [scripts/merge_si_cumulative.R](scripts/merge_si_cumulative.R), [scripts/merge_iv_table1.R](scripts/merge_iv_table1.R) | Summarize saved seeds |
| [Analytical_SE.md](Analytical_SE.md) | Separate technical note on analytical standard errors |

Generated results, logs, R session files, and HTML output are excluded by
`.gitignore`.

## Analytical standard errors

Both [SI.R](R/SI.R) and [IV.R](R/IV.R) recursively compute subject-level
influence functions (IFs), including nuisance estimation uncertainty. All four
targets use `variance = sum(IF^2) / n^2` and `SE = sqrt(variance)`. Contributions
are combined within subjects before squaring; the IFs are not sample-centered.
The same convention applies to internal stage and fixed-nuisance diagnostics.
See [Analytical_SE.md](Analytical_SE.md) for the recursion.

This variance convention changes some historical SEs, confidence intervals, and
coverage, but not data generation, seeds, point estimates, or truth. Use a new
output directory for new runs; do not combine them with historical results.
