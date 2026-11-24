################################################################################
#
#  Simulation of a two-phase Bayesian leapfrog trial (A-BEAM)
#
################################################################################
#
#  Study:          A-BEAM: App-based Inhibition Training in Alcohol Use Disorder: Protocol for
#                  an Adaptive Multicenter Leapfrog Trial including Neurophysiological Markers
#  Registration:   German Clinical Trials Register DRKS00040447
#                  (https://drks.de/search/en/trial/DRKS00040447/details); SNCTP HumRes68293
#  Code author:    Luzius Brogli (Institute of Psychology, University of Bern)
#  Code contact:   luzius.brogli@unibe.ch
#  Study protocol: Dr. Natascha Stoffel (corresponding author), natascha.stoffel@unibe.ch
#  Version:        1.0.0  24/11/2026
#  License:        CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).
#  R version:      4.4.2  (BayesFactor 0.9.12-4.7, nlme, reshape2)
#
#  Purpose:
#    Simulation-based optimisation of the design parameters (Nmin, Nmax,
#    BF thresholds) of a six-arm, two-phase leapfrog trial of alcohol-specific
#    inhibition training in alcohol use disorder. The script estimates false
#    positive / false negative rates, expected sample sizes and the expected
#    effect size of the finally selected arm for one parameter set.
#
#  Based on:
#    This script is a modified version of the following script by Simon Blackwell (https://osf.io/8mxda/):
#
#    | Simulation of an entire Leapfrog trial with multiple arms
#    | Script used to inform section 7.1.2 (Global trial simulations) for the study protocol for the study:
#    | "Feasibility randomized controlled trial using a 'leapfrog' design: optimizing an
#    |  internet-delivered cognitive training intervention to reduce anhedonia"
#    | See https://osf.io/8mxda/ and https://clinicaltrials.gov/ct2/show/NCT04791137
#    | Version: 3 22/03/21 (update from previous version: addition of further comments and explanation)
#    | Author: Simon Blackwell (simon.blackwell@rub.de)
#
#  Usage:
#    Rscript leapfrog_two_phase_simulation.R \
#            NMIN_INITIAL NMAX_INITIAL BFFAIL_INITIAL BFSUCCESS_INITIAL \
#            NMIN_UPDATED NMAX_UPDATED BFFAIL_UPDATED BFSUCCESS_UPDATED \
#            NULL_SCENARIO [N_SIMS] [SEED] [OUTPUT_DIR]
#
#    e.g. the parameter set selected for the trial:
#    Rscript leapfrog_two_phase_simulation.R 70 150 1/8 6 50 150 1/10 3 0 5000 42
#
#    Without arguments the script runs DEFAULT_ARGS (section 2a; handy in RStudio).
#    See README.md for the meaning of the arguments and the output file.
#
################################################################################


#### 1. Packages ###############################################################

# Uncomment if packages are installed in a project-local library (e.g. on the cluster)
# .libPaths(new = "R_libs")
suppressPackageStartupMessages({
  library(BayesFactor)  # ttest.tstat(): Bayes factor from a t statistic
  library(nlme)         # gls(): constrained longitudinal data analysis (cLDA)
  library(reshape2)     # melt(): wide -> long format
})


#### 2. Settings ###############################################################

## 2a. Command-line arguments (parameters varied across simulation runs) -------
# Initial phase: Nmin, Nmax and BF thresholds used until the first arm leaves the trial.
# Update phase:  Nmin, Nmax and BF thresholds used for all arms from then on.

# Defaults: the parameter set selected for the trial
# DEFAULT_ARGS: used when the script is run without arguments, as
#   c(NMIN_INITIAL, NMAX_INITIAL, BFFAIL_INITIAL, BFSUCCESS_INITIAL,
#     NMIN_UPDATED, NMAX_UPDATED, BFFAIL_UPDATED, BFSUCCESS_UPDATED, NULL_SCENARIO).
#   Set NULL_SCENARIO to "1" for the null scenario.
# DEFAULT_N_SIMS / DEFAULT_SEED / DEFAULT_OUTPUT_DIR: used when the optional arguments are omitted.
DEFAULT_ARGS       <- c("70", "150", "1/8", "6", "50", "150", "1/10", "3", "0")
DEFAULT_N_SIMS     <- 5000L
DEFAULT_SEED       <- 42L
DEFAULT_OUTPUT_DIR <- "output"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 && !is.null(DEFAULT_ARGS)) {
  message("No arguments given - using DEFAULT_ARGS: ", paste(DEFAULT_ARGS, collapse = " "))
  args <- DEFAULT_ARGS
}
if (length(args) < 9) {
  stop("Usage: Rscript leapfrog_two_phase_simulation.R ",
       "NMIN_INITIAL NMAX_INITIAL BFFAIL_INITIAL BFSUCCESS_INITIAL ",
       "NMIN_UPDATED NMAX_UPDATED BFFAIL_UPDATED BFSUCCESS_UPDATED ",
       "NULL_SCENARIO [N_SIMS] [SEED] [OUTPUT_DIR]")
}

# accepts plain numbers ("0.125") and fractions ("1/8")
parse_number <- function(x) {
  if (grepl("^\\s*[0-9.]+\\s*/\\s*[0-9.]+\\s*$", x)) {
    parts <- as.numeric(strsplit(x, "/")[[1]])
    return(parts[1] / parts[2])
  }
  as.numeric(x)
}

N_MIN_INITIAL      <- parse_number(args[1])  # min. n per arm before the first analysis (initial phase)
N_MAX_INITIAL      <- parse_number(args[2])  # max. n per arm before stopping for futility (initial phase)
BF_FAIL_INITIAL    <- parse_number(args[3])  # futility threshold, initial phase (drop arm if BF < this)
BF_SUCCESS_INITIAL <- parse_number(args[4])  # promotion threshold, initial phase (promote if BF > this)
N_MIN_UPDATED      <- parse_number(args[5])  # min. n per arm before the first analysis (update phase)
N_MAX_UPDATED      <- parse_number(args[6])  # max. n per arm before stopping for futility (update phase)
BF_FAIL_UPDATED    <- parse_number(args[7])  # futility threshold, update phase
BF_SUCCESS_UPDATED <- parse_number(args[8])  # promotion threshold, update phase
NULL_SCENARIO      <- as.logical(as.numeric(args[9]))  # TRUE: all true effects = 0 (to estimate the FPR)
N_SIMS             <- if (length(args) >= 10) as.integer(args[10]) else DEFAULT_N_SIMS  # number of simulated trials
SEED               <- if (length(args) >= 11) as.integer(args[11]) else DEFAULT_SEED
OUTPUT_DIR         <- if (length(args) >= 12) args[12] else DEFAULT_OUTPUT_DIR

stopifnot(!anyNA(c(N_MIN_INITIAL, N_MAX_INITIAL, BF_FAIL_INITIAL, BF_SUCCESS_INITIAL,
                   N_MIN_UPDATED, N_MAX_UPDATED, BF_FAIL_UPDATED, BF_SUCCESS_UPDATED,
                   NULL_SCENARIO, N_SIMS, SEED)),
          N_MIN_INITIAL <= N_MAX_INITIAL, N_MIN_UPDATED <= N_MAX_UPDATED,
          BF_FAIL_INITIAL < 1, BF_FAIL_UPDATED < 1,
          BF_SUCCESS_INITIAL > 1, BF_SUCCESS_UPDATED > 1)

## 2b. Fixed design parameters ---------------------------------------------------
N_START_ARMS <- 3    # arms at the start: control, Alc-IT-6, Alc-IT-12
N_MAX_ARMS   <- 6    # total number of arms over the whole trial

# Analysis settings
PRIOR_RSCALE      <- 0.707  # scale of the Cauchy prior on the effect size (one-sided, delta > 0)
USE_CLDA          <- TRUE   # TRUE: cLDA via nlme::gls; FALSE: Welch t-test on change scores
MIN_RELEVANT_DIFF <- 0      # effect-size difference counted as a "real" improvement (FN / improvement)

# Guard against a numerical failure of BayesFactor::ttest.tstat() for one-sided tests:
# for strongly negative t (treatment clearly worse than control, t < about -6.5) the
# numerical integration breaks down and returns huge BFs in favour of the treatment.
# TRUE: evaluate the BF at t = T_CLAMP instead (an upper bound of the true BF, see
# compute_bayes_factor()). FALSE: use the t statistic unchanged.
GUARD_NEGATIVE_T <- TRUE
T_CLAMP          <- -5

# Data generation. Only the standardised effect sizes matter for the BF; the raw
# scale (mean/SD) is arbitrary but kept so that 'realistic' data could be inspected.
START_MEAN       <- 25    # baseline mean
START_SD         <- 5     # baseline SD (= SD of the change scores)
BASE_CHANGE_D    <- 0.2   # change from baseline to follow-up common to all arms (in SD units)
PCT_MISSING      <- 28.5  # % of follow-up values missing completely at random
# Participants simulated per arm. The control arm recruits during the whole trial and
# its participants are indexed by calendar time, so it needs far more than Nmax values.
DATA_BUFFER      <- 3000
N_DATA_PER_ARM   <- max(N_MAX_INITIAL, N_MAX_UPDATED) + DATA_BUFFER

CHECKPOINT_EVERY <- 100   # save intermediate results / print summary every k simulated trials

## 2c. True effect sizes (Cohen's d relative to the original control arm) -------
# Arm numbers and labels used in the output (events / trialOrder):
#   arm 1  control                                                                    d = 0
#   arm 2  Alc-IT-6                                                                   d = 0.36
#   arm 3  Alc-IT-12                                                                  d = 0.36 + 0.24
#   arm 4  Alc-IT-boost  ("Alc-IT-boost") add-on to the control arm current at entry  + 0.30
#   arm 5  Alc-IT-HIIT   ("Alc-IT-HIIT")  HIIT add-on to the phase-1 winner           + 0.30
#   arm 6  HIIT-control  ("HIIT-control") HIIT only, standalone                       d = 0.10
if (!NULL_SCENARIO) {
  EFFECT_IT6           <- 0.36
  EFFECT_IT12          <- EFFECT_IT6 + 0.24
  INCREMENT_BOOST      <- 0.3
  INCREMENT_HIIT_ADDON <- 0.3
  EFFECT_HIIT_ONLY     <- 0.1
  OUTPUT_SUFFIX        <- ""
} else {
  EFFECT_IT6           <- 0
  EFFECT_IT12          <- 0
  INCREMENT_BOOST      <- 0
  INCREMENT_HIIT_ADDON <- 0
  EFFECT_HIIT_ONLY     <- 0
  OUTPUT_SUFFIX        <- "_FP"
}
INITIAL_EFFECTS <- c(0, EFFECT_IT6, EFFECT_IT12)  # arms 1-3


#### 3. Helper functions #######################################################

# Simulate baseline and change scores for all potential participants of one arm.
generate_arm_data <- function(armEffect) {
  baselines <- rnorm(N_DATA_PER_ARM, START_MEAN, START_SD)
  changes   <- rnorm(N_DATA_PER_ARM, (BASE_CHANGE_D + armEffect) * START_SD, START_SD)
  if (PCT_MISSING > 0) {
    missing <- sample(1:N_DATA_PER_ARM, round(PCT_MISSING * N_DATA_PER_ARM / 100))
    changes[missing] <- NA
  }
  list(baselines = baselines, changes = changes)
}

# Assemble the data for one pairwise comparison (treatment vs. current control).
# Only concurrently randomised control participants are used. Participants are
# indexed by their position in the arm's data vector, i.e. by time since the arm
# entered the trial (one participant per arm and time point).
get_comparison_data <- function(treatArm, ctrlArm, nPerArm, entryTime, allBaselines, allChanges) {
  nTreat <- nPerArm[treatArm]
  if (entryTime[treatArm] >= entryTime[ctrlArm]) {
    # Treatment arm entered at the same time as / after the control arm:
    # control participants recruited since the treatment arm entered
    offset   <- entryTime[treatArm] - entryTime[ctrlArm]
    idxCtrl  <- (offset + 1):(offset + nTreat)
    idxTreat <- 1:nTreat
  } else {
    # Treatment arm entered before the (later promoted) control arm:
    # treatment participants recruited since the control arm entered
    offset   <- entryTime[ctrlArm] - entryTime[treatArm]
    if (nTreat - offset < 1) warning("No overlapping recruitment period for arm ", treatArm)
    idxCtrl  <- 1:(nTreat - offset)
    idxTreat <- (offset + 1):nTreat
  }
  if (max(idxCtrl) > N_DATA_PER_ARM) {
    stop("Simulated data of the control arm exhausted - increase DATA_BUFFER.")
  }

  trialdata <- data.frame(
    Baseline = c(allBaselines[[ctrlArm]][idxCtrl], allBaselines[[treatArm]][idxTreat]),
    C2min1   = c(allChanges[[ctrlArm]][idxCtrl],   allChanges[[treatArm]][idxTreat]),  # change score
    Group    = c(rep("C", length(idxCtrl)), rep("T", length(idxTreat)))
  )
  trialdata$Posttx <- trialdata$Baseline + trialdata$C2min1
  trialdata$id     <- seq_len(nrow(trialdata))
  trialdata
}

# t statistic of the treatment effect from a cLDA model: baseline and post-treatment
# scores as repeated measures, baseline means constrained to be equal across groups,
# separate variances per time point and unstructured within-participant correlation.
clda_tstat <- function(trialdata) {
  trialdata$id <- as.factor(trialdata$id)
  long <- reshape2::melt(trialdata, id.vars = c("id", "Group"))
  long <- long[long$variable %in% c("Baseline", "Posttx"), ]
  X <- model.matrix(~ variable * Group, data = long)
  long$variablePosttx        <- X[, "variablePosttx"]         # time effect
  long$variablePosttx.GroupT <- X[, "variablePosttx:GroupT"]  # treatment effect at follow-up

  mod <- nlme::gls(value ~ variablePosttx + variablePosttx.GroupT,
                   data        = long,
                   weights     = nlme::varIdent(form = ~ 1 | variable),
                   correlation = nlme::corSymm(form = ~ 1 | id),
                   method      = "ML",
                   na.action   = na.exclude)
  list(t  = summary(mod)$tTable["variablePosttx.GroupT", "t-value"],
       df = mod$dims$N - mod$dims$p)
}

# t statistic from a Welch t-test on the change scores (positive = treatment better).
ttest_tstat <- function(trialdata) {
  tt <- t.test(trialdata$C2min1[trialdata$Group == "C"], trialdata$C2min1[trialdata$Group == "T"])
  list(t = -1 * unname(tt$statistic), df = unname(tt$parameter))
}

# One-sided (delta > 0) Bayes factor BF+0 for one comparison.
# The t statistic is converted with the equal-n two-sample formula, choosing n1 = n2
# such that n1 + n2 - 2 equals the residual df of the model (as in the original script).
# Returns list(bf, t); bf = NA if the model could not be fitted.
compute_bayes_factor <- function(trialdata) {
  tryCatch({
    fit <- if (USE_CLDA) clda_tstat(trialdata) else ttest_tstat(trialdata)
    # BF+0 decreases monotonically for negative t, so the BF at T_CLAMP is an upper bound
    # (<= 0.052 for n >= 30) - i.e. still below every futility threshold explored.
    tForBF <- if (GUARD_NEGATIVE_T) max(fit$t, T_CLAMP) else fit$t
    nEquiv <- (fit$df + 2) / 2
    suppressMessages(capture.output(
      bf <- BayesFactor::ttest.tstat(tForBF, n1 = nEquiv, n2 = nEquiv,
                                     nullInterval = c(0, Inf), rscale = PRIOR_RSCALE)
    ))
    list(bf = exp(bf$bf), t = fit$t)
  }, error = function(e) list(bf = NA, t = NA))
}

# Aggregate the outcome measures over all simulated trials.
summarise_simulations <- function(allData) {
  get <- function(field) sapply(allData, `[[`, field)
  nDecisions <- sum(get("nDecisions"))
  improvementSizes <- get("improvementsize")
  list(
    nSims              = length(allData),
    pctFalseNegative   = sum(get("nFalseNegatives")) * 100 / nDecisions,  # FNR (main scenario)
    pctFalsePositive   = sum(get("nFalsePositives")) * 100 / nDecisions,  # FPR (null scenario)
    pctTrialsWithFP    = mean(get("nFalsePositives") > 0) * 100,
    meanImprovement    = mean(improvementSizes),     # d of winner vs. original control
    meanDiscrepancy    = mean(get("discrep")),       # winner vs. best arm present in that trial
    nBestWon           = sum(get("bestwon")),
    nImproved          = sum(get("improvement")),
    meanNTotal         = mean(get("totalTested")),
    meanNPhase1        = mean(get("totalTestPart1")),
    meanNPhase2        = mean(get("totalTestPart2")),
    improvementTable   = table(round(improvementSizes, 2)),
    nBFErrors          = sum(get("nBFErrors")),
    nStrongNegativeT   = sum(get("nStrongNegativeT"))
  )
}

print_summary <- function(s) {
  cat("\n--- Results after", s$nSims, "simulated trials ---\n")
  cat("False negative rate (% of decisions):          ", s$pctFalseNegative, "\n")
  cat("False positive rate (% of decisions):          ", s$pctFalsePositive, "\n")
  cat("Trials with >= 1 false positive (%):           ", s$pctTrialsWithFP, "\n")
  cat("Mean improvement (d of winner vs. control):    ", s$meanImprovement, "\n")
  cat("Mean discrepancy (winner vs. best arm, d):     ", s$meanDiscrepancy, "\n")
  cat("Best arm won (n trials):                       ", s$nBestWon, "\n")
  cat("Improvement over control (n trials):           ", s$nImproved, "\n")
  cat("Mean sample size - total / phase 1 / phase 2:  ", s$meanNTotal, "/", s$meanNPhase1, "/", s$meanNPhase2, "\n")
  cat("BF computation failures / comparisons with t <", T_CLAMP, ":", s$nBFErrors, "/", s$nStrongNegativeT, "\n")
  cat("Frequency of improvement sizes:\n")
  print(s$improvementTable)
}


#### 4. Simulation of one complete trial #######################################
#
# Time advances in discrete steps ("time points"). At each time point every active
# arm recruits one participant; every treatment arm with n >= nMin is analysed
# against the current control arm and may be
#   - dropped  (n > nMax, or BF < BFfail),
#   - promoted (BF > BFsuccess, and the largest such BF at this time point):
#     it becomes the new control, the old control leaves,
#   - continued (otherwise).
# Whenever an arm leaves, the thresholds switch to the update-phase values and a new
# arm may enter: Alc-IT-boost after the first departure (phase 1); Alc-IT-HIIT
# and HIIT-control once only the phase-1 winner is left (phase 2).
# The trial ends when only one arm (the winner) remains.

simulate_trial <- function() {

  ## Data of the initial arms
  allBaselines <- list()
  allChanges   <- list()
  for (arm in seq_along(INITIAL_EFFECTS)) {
    armData <- generate_arm_data(INITIAL_EFFECTS[arm])
    allBaselines[[arm]] <- armData$baselines
    allChanges[[arm]]   <- armData$changes
  }

  ## Trial state. Lists/vectors indexed by time point tp: activeArms, controlArm, BFs, events
  trueEffects   <- INITIAL_EFFECTS          # true effect size of each arm (filled as arms enter)
  activeArms    <- list(seq(1:N_START_ARMS))# arms in the trial at each time point
  controlArm    <- c(1)                     # control arm at each time point
  BFs           <- list()                   # BF of each arm at each time point (-1 = not computed)
  events        <- list()                   # what happened at each time point (text)
  trialOrder    <- c()                      # flat log: (arm, time point, decision) triples
  nPerArm       <- rep(1, N_MAX_ARMS)       # participants per arm in the current phase
  nPerArmPhase1 <- rep(1, 4)                # snapshot of nPerArm at the end of phase 1
  entryTime     <- rep(1, N_MAX_ARMS)       # time point at which each arm entered
  nextArm       <- N_START_ARMS + 1         # number of the next arm to enter
  boostAdded <- FALSE; hiitAddonAdded <- FALSE; hiitOnlyAdded <- FALSE

  nComparisons <- 0; nDecisions <- 0; nFalsePositives <- 0; nFalseNegatives <- 0
  nBFErrors <- 0; nStrongNegativeT <- 0

  # thresholds start with the initial-phase values
  nMin <- N_MIN_INITIAL; nMax <- N_MAX_INITIAL
  BFfail <- BF_FAIL_INITIAL; BFsuccess <- BF_SUCCESS_INITIAL

  # Let a new arm enter at the current time point (modifies the trial state above).
  enter_arm <- function(effect, label) {
    activeArms[[tp + 1]] <<- c(activeArms[[tp + 1]], nextArm)
    entryTime[nextArm]   <<- tp
    events[[tp]]         <<- c(events[[tp]], paste(nextArm, " enters the trial (", label, ")", sep = ""))
    nPerArm[nextArm]     <<- 1
    trueEffects[nextArm] <<- effect
    armData <- generate_arm_data(effect)
    allBaselines[[nextArm]] <<- armData$baselines
    allChanges[[nextArm]]   <<- armData$changes
    trialOrder <<- c(trialOrder, c(nextArm, tp, label))
    nextArm    <<- nextArm + 1
  }

  tp <- 1
  ended <- FALSE
  while (!ended) {
    BFs[[tp]]    <- rep(-1, N_MAX_ARMS)
    events[[tp]] <- list()
    leavingArms  <- vector()   # arms leaving at this time point (incl. a replaced control arm)

    treatmentArms <- activeArms[[tp]][activeArms[[tp]] != controlArm[tp]]

    ## --- Step 1: Bayes factors of all arms due for an analysis ---------------
    # (nMin <= n <= nMax). All arms are compared against the same, current control.
    for (thisArm in treatmentArms) {
      if (nPerArm[thisArm] <= nMax && nPerArm[thisArm] >= nMin) {
        nComparisons <- nComparisons + 1
        trialdata <- get_comparison_data(thisArm, controlArm[tp], nPerArm, entryTime,
                                         allBaselines, allChanges)
        result <- compute_bayes_factor(trialdata)
        BFs[[tp]][thisArm] <- result$bf
        if (is.na(result$bf)) nBFErrors <- nBFErrors + 1
        if (!is.na(result$t) && result$t < T_CLAMP) nStrongNegativeT <- nStrongNegativeT + 1
      }
    }

    # At most one arm is promoted per time point: the one with the largest BF above
    # BFsuccess. Other arms above BFsuccess simply continue and are compared against the
    # new control from the next time point on.
    passed <- treatmentArms[!is.na(BFs[[tp]][treatmentArms]) & BFs[[tp]][treatmentArms] > BFsuccess]
    promotedArm <- if (length(passed) > 0) passed[which.max(BFs[[tp]][passed])] else NA

    ## --- Step 2: decisions ----------------------------------------------------
    for (thisArm in treatmentArms) {
      isTrulyBetter <- trueEffects[thisArm] > trueEffects[controlArm[tp]] + MIN_RELEVANT_DIFF

      if (nPerArm[thisArm] > nMax) {
        # Rule 1: Nmax reached without decisive evidence -> drop (futility)
        events[[tp]] <- c(events[[tp]], paste(thisArm, " hit Nmax", sep = ""))
        leavingArms  <- c(leavingArms, thisArm)
        nDecisions   <- nDecisions + 1
        if (isTrulyBetter) nFalseNegatives <- nFalseNegatives + 1
        trialOrder   <- c(trialOrder, c(thisArm, tp, "dropMax"))

      } else if (nPerArm[thisArm] >= nMin) {
        # Rule 2: enough participants -> decide based on the Bayes factor from step 1
        BF <- BFs[[tp]][thisArm]

        if (is.na(BF)) {
          # BF could not be computed -> keep recruiting
          nPerArm[thisArm] <- nPerArm[thisArm] + 1

        } else if (BF < BFfail) {
          # Rule 2a: evidence for the null -> drop (futility)
          events[[tp]] <- c(events[[tp]], paste(thisArm, " hit BFfail", sep = ""))
          leavingArms  <- c(leavingArms, thisArm)
          nDecisions   <- nDecisions + 1
          if (isTrulyBetter) nFalseNegatives <- nFalseNegatives + 1
          trialOrder   <- c(trialOrder, c(thisArm, tp, "drop"))

        } else if (!is.na(promotedArm) && thisArm == promotedArm) {
          # Rule 2b: evidence for superiority (largest BF > BFsuccess at this time point)
          # -> promote to control; the old control leaves
          events[[tp]] <- c(events[[tp]], paste(thisArm, " hit BFsuccess", sep = ""))
          leavingArms  <- c(leavingArms, controlArm[tp])
          nDecisions   <- nDecisions + 1
          nPerArm[thisArm] <- nPerArm[thisArm] + 1
          controlArm[tp + 1] <- thisArm
          if (!(trueEffects[thisArm] > trueEffects[controlArm[tp]])) nFalsePositives <- nFalsePositives + 1
          trialOrder   <- c(trialOrder, c(thisArm, tp, "promote"))

        } else {
          # Rule 2c: inconclusive, or above BFsuccess but another arm had a larger BF
          # -> keep recruiting
          nPerArm[thisArm] <- nPerArm[thisArm] + 1
        }

      } else {
        # n < nMin: keep recruiting without analysis
        nPerArm[thisArm] <- nPerArm[thisArm] + 1
      }
    }

    ## --- Update the trial for the next time point ----------------------------
    if (length(leavingArms) > 0) {
      # After the first departure, switch to the update-phase thresholds (once)
      if (nMin != N_MIN_UPDATED || nMax != N_MAX_UPDATED ||
          BFfail != BF_FAIL_UPDATED || BFsuccess != BF_SUCCESS_UPDATED) {
        BFfail    <- BF_FAIL_UPDATED
        BFsuccess <- BF_SUCCESS_UPDATED
        events[[tp]] <- c(events[[tp]], paste("BF thresholds updated: BFfail =", BFfail, ", BFsuccess =", BFsuccess))
        nMin <- N_MIN_UPDATED
        nMax <- N_MAX_UPDATED
        events[[tp]] <- c(events[[tp]], paste("Minimum sample size for analysis updated to", nMin))
      }

      activeArms[[tp + 1]] <- activeArms[[tp]][!(activeArms[[tp]] %in% leavingArms)]

      # Effect sizes of add-on arms build on the control arm of the next time point
      newControlEffect <- if (length(controlArm) == tp + 1) trueEffects[controlArm[tp + 1]] else trueEffects[controlArm[tp]]

      for (k in seq_along(leavingArms)) {
        if (!boostAdded) {
          # Phase 1: Alc-IT-boost replaces the first arm that leaves
          boostAdded <- TRUE
          enter_arm(newControlEffect + INCREMENT_BOOST, "Alc-IT-boost")

        } else if (!hiitAddonAdded && length(activeArms[[tp + 1]]) < 2) {
          # End of phase 1 (only the winner is left) -> start phase 2.
          # Sample sizes are counted separately for each phase.
          hiitAddonAdded <- TRUE
          nPerArmPhase1  <- nPerArm
          nPerArm        <- rep(1, N_MAX_ARMS)
          enter_arm(newControlEffect + INCREMENT_HIIT_ADDON, "Alc-IT-HIIT")
        }
        if (!hiitOnlyAdded && hiitAddonAdded) {
          # Phase 2: HIIT-control enters together with Alc-IT-HIIT
          hiitOnlyAdded <- TRUE
          enter_arm(EFFECT_HIIT_ONLY, "HIIT-control")
        }
      }
    } else {
      activeArms[[tp + 1]] <- activeArms[[tp]]
      controlArm[tp + 1]   <- controlArm[tp]
    }

    # If the control arm stays the same, it recruits one participant at this time point
    # (a newly promoted control arm has already been incremented above)
    if (is.na(controlArm[tp + 1]) | controlArm[tp + 1] == controlArm[tp]) {
      controlArm[tp + 1] <- controlArm[tp]
      nPerArm[controlArm[tp]] <- nPerArm[controlArm[tp]] + 1
    }

    ## --- End of trial: one arm left -----------------------------------------
    if (length(activeArms[[tp + 1]]) < 2) {
      winningArm      <- activeArms[[tp + 1]][1]
      improvementsize <- trueEffects[winningArm] - trueEffects[1]
      improvement     <- trueEffects[winningArm] > trueEffects[1] + MIN_RELEVANT_DIFF
      bestwon         <- trueEffects[winningArm] == max(trueEffects)
      discrep         <- if (bestwon) 0 else trueEffects[winningArm] - max(trueEffects)
      totalTestPart1  <- sum(nPerArmPhase1)
      totalTestPart2  <- sum(nPerArm)
      totalTest       <- totalTestPart1 + totalTestPart2
      ended <- TRUE
    }
    tp <- tp + 1
  }

  ## Results of this trial (field names kept from the original script for compatibility)
  list(
    BFs                = list(BFs),            # BF per arm and time point
    cArm               = list(controlArm),     # control arm per time point
    whereAt            = list(nPerArm),        # n per arm in phase 2
    whereAtTrial1      = list(nPerArmPhase1),  # n per arm in phase 1
    whereAtTrial2      = list(nPerArm),        # n per arm in phase 2
    nComparisons       = nComparisons,         # number of interim analyses
    nFalseNegatives    = nFalseNegatives,
    nFalsePositives    = nFalsePositives,
    totalTested        = totalTest,            # total sample size
    totalTestPart1     = totalTestPart1,       # sample size phase 1
    totalTestPart2     = totalTestPart2,       # sample size phase 2
    events             = list(events),
    improvement        = improvement,          # winner better than original control?
    bestwon            = bestwon,              # winner = best arm present in this trial?
    trueds             = list(trueEffects),
    armsIn             = list(activeArms),
    whenEntered        = list(entryTime),
    nDecisions         = nDecisions,           # number of drop/promote decisions (FPR/FNR denominator)
    discrep            = discrep,
    improvementsize    = improvementsize,      # d of winner - d of original control
    maxpossimprovement = max(trueEffects),
    trialOrder         = trialOrder,
    nBFErrors          = nBFErrors,            # BF computations that failed (BF = NA)
    nStrongNegativeT   = nStrongNegativeT      # comparisons with t < T_CLAMP
  )
}


#### 5. Run the simulation #####################################################

settings <- list(
  nMinInitial = N_MIN_INITIAL, nMaxInitial = N_MAX_INITIAL,
  nMinUpdated = N_MIN_UPDATED, nMaxUpdated = N_MAX_UPDATED,
  BFfailInitial = BF_FAIL_INITIAL, BFsuccessInitial = BF_SUCCESS_INITIAL,
  BFfailUpdated = BF_FAIL_UPDATED, BFsuccessUpdated = BF_SUCCESS_UPDATED,
  nullScenario = NULL_SCENARIO, nSims = N_SIMS, seed = SEED,
  priorRscale = PRIOR_RSCALE, useCLDA = USE_CLDA, pctMissing = PCT_MISSING,
  guardNegativeT = GUARD_NEGATIVE_T, tClamp = T_CLAMP,
  rng = "L'Ecuyer-CMRG, one stream per simulated trial",
  initialEffects = INITIAL_EFFECTS, incrementBoost = INCREMENT_BOOST,
  incrementHiitAddon = INCREMENT_HIIT_ADDON, effectHiitOnly = EFFECT_HIIT_ONLY,
  R = R.version.string, BayesFactor = as.character(packageVersion("BayesFactor"))
)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
outputName <- function(tag) {
  fmt <- function(x) format(signif(x, 4))   # e.g. 1/15 -> "0.06667"
  file.path(OUTPUT_DIR, paste0("TwoTrial_N", N_MIN_INITIAL, "-", N_MAX_INITIAL,
                               "_NUp", N_MIN_UPDATED, "-", N_MAX_UPDATED,
                               "_BF", fmt(BF_FAIL_INITIAL), "-", fmt(BF_SUCCESS_INITIAL),
                               "_BFUp", fmt(BF_FAIL_UPDATED), "-", fmt(BF_SUCCESS_UPDATED),
                               "_", tag, OUTPUT_SUFFIX, ".RData"))
}

# Random numbers: every simulated trial k gets its own L'Ecuyer-CMRG stream derived from
# SEED (stream k = parallel::nextRNGStream applied k - 1 times), so each trial can be
# re-simulated on its own and the trials could be run in any order.
# Set SKIP_MAIN_RUN <- TRUE before sourcing this file to load sections 1-4 only.
if (!exists("SKIP_MAIN_RUN") || !isTRUE(SKIP_MAIN_RUN)) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(SEED)
  stream <- .Random.seed
  allData <- list()
  for (sim in 1:N_SIMS) {
    assign(".Random.seed", stream, envir = globalenv())
    allData[[sim]] <- simulate_trial()
    stream <- parallel::nextRNGStream(stream)

    if (sim %% CHECKPOINT_EVERY == 0 && sim < N_SIMS) {
      summaryStats <- summarise_simulations(allData)
      print_summary(summaryStats)
      save(settings, allData, summaryStats, file = outputName(paste0("iter", sim)))
    }
  }

  summaryStats <- summarise_simulations(allData)
  print_summary(summaryStats)
  save(settings, allData, summaryStats, file = outputName("Final"))
}
