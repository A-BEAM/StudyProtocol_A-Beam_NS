# A-BEAM leapfrog trial simulation

Simulation code for the sample-size planning of the A-BEAM study, **App-based Inhibition
Training in Alcohol Use Disorder: Protocol for an Adaptive Multicenter Leapfrog Trial
including Neurophysiological Markers**

* **Code author and contact:** Luzius Brogli (Institute of Psychology, University of Bern), luzius.brogli@unibe.ch
* **Study protocol:** Dr. Natascha Stoffel (corresponding author), natascha.stoffel@unibe.ch
* **Version:** 1.0.0 (24/11/2026)
* **License:** [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
* **How to cite:** Stoffel N, et al. App-based Inhibition Training in Alcohol Use Disorder:
  Protocol for an Adaptive Multicenter Leapfrog Trial including Neurophysiological Markers.

`leapfrog_two_phase_simulation.R` simulates a six-arm, two-phase Bayesian leapfrog trial for
one set of design parameters and estimates false positive and false negative rates, expected
sample sizes and the expected effect of the finally selected arm. It is adapted from a script
by Simon Blackwell (https://osf.io/8mxda/, CC BY 4.0).

## Requirements

R (used: 4.4.2) with `BayesFactor` (used: 0.9.12-4.7) and `reshape2`; `nlme` and
`parallel` ship with R.

```r
install.packages(c("BayesFactor", "reshape2"))
```

## Usage

```bash
Rscript leapfrog_two_phase_simulation.R NMIN_INITIAL NMAX_INITIAL BFFAIL_INITIAL BFSUCCESS_INITIAL \
        NMIN_UPDATED NMAX_UPDATED BFFAIL_UPDATED BFSUCCESS_UPDATED NULL_SCENARIO [N_SIMS] [SEED] [OUTPUT_DIR]
```

| Argument | Meaning |
|---|---|
| `NMIN_*` | minimum n per arm before an arm is analysed |
| `NMAX_*` | maximum n per arm; the arm is dropped when it is exceeded |
| `BFFAIL_*` | futility threshold: drop the arm if BF < value |
| `BFSUCCESS_*` | promotion threshold: the arm becomes the new control if BF > value |
| `NULL_SCENARIO` | `0` = assumed effects (FNR, sample size, improvement); `1` = all true effects 0 (FPR) |
| `N_SIMS`, `SEED`, `OUTPUT_DIR` | optional; defaults 5000, 42, `output` |

Nmin, Nmax and both BF thresholds are set separately for the two parts of the trial:
`*_INITIAL` values apply until the first arm leaves the trial, `*_UPDATED` values to all
arms from then on. Thresholds may be given as fractions (`1/8`) or decimals. Any values can
be passed. Without arguments (e.g. sourced in RStudio) the script runs `DEFAULT_ARGS`,
the design selected for the trial (initial Nmin 70, Nmax 150, BFfail 1/8,
BFsuccess 6; updated Nmin 50, Nmax 150, BFfail 1/10, BFsuccess 3):

```bash
Rscript leapfrog_two_phase_simulation.R 70 150 1/8 6 50 150 1/10 3 0   # assumed effects
Rscript leapfrog_two_phase_simulation.R 70 150 1/8 6 50 150 1/10 3 1   # null scenario
```

The assumed effects, missing-data rate, prior and other fixed settings are set and commented
in section 2 of the script. One simulated trial takes a few seconds; every trial has its own
random-number stream derived from `SEED`, so results are reproducible.

## Output

`<OUTPUT_DIR>/TwoTrial_N<NminInitial>-<NmaxInitial>_NUp<NminUpdated>-<NmaxUpdated>_BF<fail>-<success>_BFUp<fail>-<success>_Final[_FP].RData`
(`_FP` = null scenario; intermediate `_iter<k>` files every 100 trials) with:

* `settings`: all parameters of the run, including seed and package versions
* `summaryStats`: FPR, FNR, mean N (phase 1, phase 2, total) and improvement (true effect
  of the winning arm vs. the original control)
* `allData`: the full course of every simulated trial (arms, sample sizes, Bayes factors,
  decisions); fields are described in `simulate_trial()`
