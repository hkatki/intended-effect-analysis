##
# Bootstrap for IE analyses that relax assumptions - BMC MRM paper
# HKatki 5 May 2026
##

# setwd("/Users/katkih/hkCurrent/Dropbox/Liquid-Biopsy/2021-MCED-trial-proposal/Conceal-Reveal/R")
setwd("/Users/hormuzd.katki/ACS Guidelines Scr Dropbox/Hormuzd Katki/Conceal-Reveal/R")

# packages <- c("xlsx","scales","lmtest","dplyr","ggplot2","survival","gmodels","coxph.risk","geepack",
#               "MESS","psych","Hmisc","glmnet","boot","zoo", "scales")
# lapply(packages, require, c = T)



##
# Fast simulation + nonparametric bootstrap by G
# Now includes loss of signal through Mtilde.
#
# The original/gold-standard quantities use M:
#   RR_all, RR_pos, RR_neg
#
# The loss-of-signal quantities use Mtilde:
#   RR_pos.loss, RR_neg.loss
#
# Mtilde is generated as follows:
#   - If M = 0, then Mtilde = 0.
#   - If M = 1 and D = 1, then Mtilde ~ Bernoulli(p_Mtilde_pos_given_Mpos_Dpos).
#   - If M = 1 and D = 0, then Mtilde ~ Bernoulli(p_Mtilde_pos_given_Mpos_Dneg).
#   - These retention probabilities do not depend on G.
#
# NEW: unknown ever-positivity
#
# U is generated as follows:
#   - U = 1 means unknown ever-positivity M.
#   - U = 0 means known/observed ever-positivity M.
#   - U ~ Bernoulli(p_U_unknown_given_DG), with 4 probabilities depending on D and G.
#
# The unknown-ever-positivity quantities are:
#   RR_pos.unk, RR_neg.unk
#     calculated only among U = 0, i.e. those with observed M.
#
#   RR_pos.unk.corr, RR_neg.unk.corr
#     calculated among U = 0 but corrected using the U = 1 unknown-positive table.
##


##
# Bootstrap one dataset, preserving the original number of G=0 and G=1
##

bootstrap_superpop_by_G <- function(superpop) {
  
  idx_g0 <- which(superpop[, "G"] == 0)
  idx_g1 <- which(superpop[, "G"] == 1)
  
  boot_idx_g0 <- sample(idx_g0, size = length(idx_g0), replace = TRUE)
  boot_idx_g1 <- sample(idx_g1, size = length(idx_g1), replace = TRUE)
  
  boot_idx <- c(boot_idx_g0, boot_idx_g1)
  
  superpop.boot <- superpop[boot_idx, , drop = FALSE]
  
  return(superpop.boot)
}


##
# NEW: unknown ever-positivity helper
#
# Make a D by G table:
#
#              screen control
#   D+             .       .
#   D-             .       .
#
# Assumes:
#   G = 1 is screen
#   G = 0 is control
#   D = 1 is D+
#   D = 0 is D-
##

make_D_by_G_table_for_unknown <- function(dat) {
  
  if (nrow(dat) == 0) {
    out <- matrix(
      0,
      nrow = 2,
      ncol = 2,
      dimnames = list(
        D = c("D+", "D-"),
        G = c("screen", "control")
      )
    )
    return(out)
  }
  
  D_fac <- factor(
    ifelse(dat[, "D"] == 1, "D+", "D-"),
    levels = c("D+", "D-")
  )
  
  G_fac <- factor(
    ifelse(dat[, "G"] == 1, "screen", "control"),
    levels = c("screen", "control")
  )
  
  out <- table(D_fac, G_fac)
  out <- as.matrix(out)
  storage.mode(out) <- "numeric"
  
  dimnames(out) <- list(
    D = c("D+", "D-"),
    G = c("screen", "control")
  )
  
  return(out)
}


##
# NEW: unknown ever-positivity helper
#
# Calculate RR from a D by G table:
#
#              screen control
#   D+             a       b
#   D-             c       d
#
# RR = P(D+ | screen) / P(D+ | control)
##

RR_from_D_by_G_table_for_unknown <- function(tab) {
  
  tab <- as.matrix(tab)
  
  risk_screen <- tab["D+", "screen"] / sum(tab[, "screen"])
  risk_control <- tab["D+", "control"] / sum(tab[, "control"])
  
  RR <- risk_screen / risk_control
  
  return(RR)
}


##
# NEW: unknown ever-positivity helper
#
# Correct observed never-positive and ever-positive tables using the
# unknown-positive table.
#
# Inputs:
#
#   obs_never:
#     D by G table among U = 0 and M = 0
#
#   obs_ever:
#     D by G table among U = 0 and M = 1
#
#   unknown_pos:
#     D by G table among U = 1
#
# For each D row:
#
#   standard_table = obs_never + obs_ever + unknown_pos
#
#   unknown_frac_screen =
#     unknown_pos_screen / standard_screen
#
#   unknown_frac_control =
#     unknown_pos_control / standard_control
#
#   correction_factor =
#     (1 - unknown_frac_screen) / (1 - unknown_frac_control)
#
# The screen columns are unchanged.
# The control columns are multiplied by the row-specific correction factor.
##

correct_positive_tables_simple <- function(obs_never,
                                           obs_ever,
                                           unknown_pos,
                                           screen_col = "screen",
                                           control_col = "control",
                                           round_output = FALSE) {
  
  obs_never   <- as.matrix(obs_never)
  obs_ever    <- as.matrix(obs_ever)
  unknown_pos <- as.matrix(unknown_pos)
  
  if (!identical(dim(obs_never), dim(obs_ever))) {
    stop("obs_never and obs_ever must have the same dimensions.")
  }
  
  if (!identical(dim(obs_never), dim(unknown_pos))) {
    stop("obs_never and unknown_pos must have the same dimensions.")
  }
  
  if (is.null(rownames(obs_never)) ||
      is.null(rownames(obs_ever)) ||
      is.null(rownames(unknown_pos))) {
    stop("All three input tables must have row names.")
  }
  
  if (is.null(colnames(obs_never)) ||
      is.null(colnames(obs_ever)) ||
      is.null(colnames(unknown_pos))) {
    stop("All three input tables must have column names.")
  }
  
  required_cols <- c(screen_col, control_col)
  
  if (!all(required_cols %in% colnames(obs_never))) {
    stop("obs_never must contain screen_col and control_col.")
  }
  
  if (!all(required_cols %in% colnames(obs_ever))) {
    stop("obs_ever must contain screen_col and control_col.")
  }
  
  if (!all(required_cols %in% colnames(unknown_pos))) {
    stop("unknown_pos must contain screen_col and control_col.")
  }
  
  safe_divide_local <- function(num, den) {
    out <- num / den
    out[den == 0] <- NA_real_
    out
  }
  
  standard_table <- obs_never + obs_ever + unknown_pos
  
  unknown_frac_screen <- safe_divide_local(
    unknown_pos[, screen_col],
    standard_table[, screen_col]
  )
  
  unknown_frac_control <- safe_divide_local(
    unknown_pos[, control_col],
    standard_table[, control_col]
  )
  
  observed_frac_screen <- 1 - unknown_frac_screen
  observed_frac_control <- 1 - unknown_frac_control
  
  correction_factor <- safe_divide_local(
    observed_frac_screen,
    observed_frac_control
  )
  
  corrected_never <- obs_never
  corrected_ever  <- obs_ever
  
  ## Screen columns are unchanged.
  corrected_never[, screen_col] <- obs_never[, screen_col]
  corrected_ever[,  screen_col] <- obs_ever[,  screen_col]
  
  ## Control columns are corrected row by row.
  corrected_never[, control_col] <-
    obs_never[, control_col] * correction_factor
  
  corrected_ever[, control_col] <-
    obs_ever[, control_col] * correction_factor
  
  if (round_output) {
    corrected_never <- round(corrected_never)
    corrected_ever  <- round(corrected_ever)
  }
  
  out <- list(
    standard_trial_analysis_table = standard_table,
    unknown_frac_screen = unknown_frac_screen,
    unknown_frac_control = unknown_frac_control,
    observed_frac_screen = observed_frac_screen,
    observed_frac_control = observed_frac_control,
    correction_factor = correction_factor,
    corrected_never_positive_table = corrected_never,
    corrected_ever_positive_table = corrected_ever
  )
  
  return(out)
}


##
# Fast IE function
# Computes RR_all, RR_pos, RR_neg using M,
# and RR_pos.loss, RR_neg.loss using Mtilde.
#
# NEW: unknown ever-positivity
# Also computes:
#   RR_pos.unk, RR_neg.unk
#   RR_pos.unk.corr, RR_neg.unk.corr
##

IE_fast <- function(superpop, f1 = f_1, f2 = f_2) {
  
  ## Expect columns: G, D, M, Mtilde, R
  ## NEW: unknown ever-positivity
  ## Also expect column: U
  G      <- superpop[, "G"]
  D      <- superpop[, "D"]
  M      <- superpop[, "M"]
  Mtilde <- superpop[, "Mtilde"]
  R      <- superpop[, "R"]
  
  ## NEW: unknown ever-positivity
  if (!("U" %in% colnames(superpop))) {
    stop("IE_fast() expected column U in superpop, but U was not found.")
  }
  
  U <- superpop[, "U"]
  
  ## ---- small helper: safe mean in a subset ----
  mean_if <- function(x, idx) {
    if (!any(idx)) return(NA_real_)
    mean(x[idx], na.rm = TRUE)
  }
  
  ## ---- small helper: safe weighted risk ----
  weighted_risk <- function(D, idx, w) {
    if (!any(idx)) return(NA_real_)
    denom <- sum(w[idx], na.rm = TRUE)
    if (is.na(denom) || denom == 0) return(NA_real_)
    sum(w[idx] * D[idx], na.rm = TRUE) / denom
  }
  
  
  ## ==========================================================
  ## Gold-standard quantities using M
  ## ==========================================================
  
  idx_neg <- (M == 0)
  idx_pos <- (M == 1)
  
  risk_g1.neg <- mean_if(D == 1, (G == 1) & idx_neg)
  risk_g0.neg <- mean_if(D == 1, (G == 0) & idx_neg)
  
  risk_g1.pos <- mean_if(D == 1, (G == 1) & idx_pos)
  risk_g0.pos <- mean_if(D == 1, (G == 0) & idx_pos)
  
  risk_g1.all <- mean_if(D == 1, (G == 1))
  risk_g0.all <- mean_if(D == 1, (G == 0))
  
  RR_all_hat <- risk_g1.all / risk_g0.all
  RR_pos_hat <- risk_g1.pos / risk_g0.pos
  RR_neg_hat <- risk_g1.neg / risk_g0.neg
  
  
  ## ==========================================================
  ## Loss-of-signal naive quantities using Mtilde
  ## ==========================================================
  
  idx_neg.loss <- (Mtilde == 0)
  idx_pos.loss <- (Mtilde == 1)
  
  risk_g1.neg.loss <- mean_if(D == 1, (G == 1) & idx_neg.loss)
  risk_g0.neg.loss <- mean_if(D == 1, (G == 0) & idx_neg.loss)
  
  risk_g1.pos.loss <- mean_if(D == 1, (G == 1) & idx_pos.loss)
  risk_g0.pos.loss <- mean_if(D == 1, (G == 0) & idx_pos.loss)
  
  RR_pos_loss_hat <- risk_g1.pos.loss / risk_g0.pos.loss
  RR_neg_loss_hat <- risk_g1.neg.loss / risk_g0.neg.loss
  
  
  ## ==========================================================
  ## Bias-corrected RRpos under loss-of-signal
  ##
  ## Formula:
  ## P0(D+|M+) =
  ##   {1 + [P0(D-)/P0(D+)] *
  ##        [P0(Mtilde+|D-) / P1(Mtilde+|M+,D-)] /
  ##        [P0(Mtilde+|D+) / P1(Mtilde+|M+,D+)]
  ##   }^(-1)
  ##
  ## RR_pos.loss.corr = P1(D+|M+) / corrected P0(D+|M+)
  ##
  ## This uses the full-data version, not the sampling-weighted version.
  ## ==========================================================
  
  P0_Dpos <- mean_if(D == 1, (G == 0))
  P0_Dneg <- mean_if(D == 0, (G == 0))
  
  P0_Mt_pos_given_Dpos <- mean_if(Mtilde == 1, (G == 0) & (D == 1))
  P0_Mt_pos_given_Dneg <- mean_if(Mtilde == 1, (G == 0) & (D == 0))
  
  P1_Mt_pos_given_Mpos_Dpos <- mean_if(Mtilde == 1, (G == 1) & (M == 1) & (D == 1))
  P1_Mt_pos_given_Mpos_Dneg <- mean_if(Mtilde == 1, (G == 1) & (M == 1) & (D == 0))
  
  needed_corr <- c(
    P0_Dpos,
    P0_Dneg,
    P0_Mt_pos_given_Dpos,
    P0_Mt_pos_given_Dneg,
    P1_Mt_pos_given_Mpos_Dpos,
    P1_Mt_pos_given_Mpos_Dneg
  )
  
  if (any(is.na(needed_corr)) ||
      P0_Dpos == 0 ||
      P0_Mt_pos_given_Dpos == 0 ||
      P1_Mt_pos_given_Mpos_Dpos == 0 ||
      P1_Mt_pos_given_Mpos_Dneg == 0) {
    
    RR_pos_loss_corr_hat <- NA_real_
    
  } else {
    
    ratio_term <- (P0_Dneg / P0_Dpos) *
      (
        (P0_Mt_pos_given_Dneg / P1_Mt_pos_given_Mpos_Dneg) /
          (P0_Mt_pos_given_Dpos / P1_Mt_pos_given_Mpos_Dpos)
      )
    
    P0_Dpos_given_Mpos_corr <- 1 / (1 + ratio_term)
    
    P1_Dpos_given_Mpos <- risk_g1.pos
    
    RR_pos_loss_corr_hat <- P1_Dpos_given_Mpos / P0_Dpos_given_Mpos_corr
  }
  
  
  ## ==========================================================
  ## Stratified-sampling estimators in control arm
  ##
  ## In G=0, only use revealed observations R=1.
  ## Weight by inverse reveal fraction:
  ##   D=1 stratum: 1/f1
  ##   D=0 stratum: 1/f2
  ##
  ## RR_pos.samp:
  ##   numerator = P1(D+|M+) from fully observed G=1 screen arm
  ##   denominator = weighted P0(D+|M+) from sampled/revealed G=0 controls
  ##
  ## RR_neg.samp:
  ##   numerator = P1(D+|M-) from fully observed G=1 screen arm
  ##   denominator = weighted P0(D+|M-) from sampled/revealed G=0 controls
  ## ==========================================================
  
  if (f1 <= 0 || f2 <= 0) {
    
    RR_pos_samp_hat <- NA_real_
    RR_neg_samp_hat <- NA_real_
    
  } else {
    
    w0 <- rep(0, length(D))
    
    w0[G == 0 & R == 1 & D == 1] <- 1 / f1
    w0[G == 0 & R == 1 & D == 0] <- 1 / f2
    
    ## Weighted control-arm risks among true M+ and true M-
    risk_g0.pos.samp <- weighted_risk(
      D = D,
      idx = (G == 0) & (R == 1) & (M == 1),
      w = w0
    )
    
    risk_g0.neg.samp <- weighted_risk(
      D = D,
      idx = (G == 0) & (R == 1) & (M == 0),
      w = w0
    )
    
    ## Screen-arm numerator remains fully observed
    risk_g1.pos.samp <- risk_g1.pos
    risk_g1.neg.samp <- risk_g1.neg
    
    RR_pos_samp_hat <- risk_g1.pos.samp / risk_g0.pos.samp
    RR_neg_samp_hat <- risk_g1.neg.samp / risk_g0.neg.samp
  }
  
  
  ## ==========================================================
  ## NEW: unknown ever-positivity uncorrected quantities
  ##
  ## RR_pos.unk and RR_neg.unk are calculated among U = 0 only,
  ## i.e. among subjects whose M is observed.
  ##
  ## These can be biased for the gold-standard RR_pos and RR_neg
  ## if U depends on D and G.
  ## ==========================================================
  
  idx_known_M <- (U == 0)
  
  idx_neg.unk <- (M == 0) & idx_known_M
  idx_pos.unk <- (M == 1) & idx_known_M
  
  risk_g1.neg.unk <- mean_if(D == 1, (G == 1) & idx_neg.unk)
  risk_g0.neg.unk <- mean_if(D == 1, (G == 0) & idx_neg.unk)
  
  risk_g1.pos.unk <- mean_if(D == 1, (G == 1) & idx_pos.unk)
  risk_g0.pos.unk <- mean_if(D == 1, (G == 0) & idx_pos.unk)
  
  RR_pos_unk_hat <- risk_g1.pos.unk / risk_g0.pos.unk
  RR_neg_unk_hat <- risk_g1.neg.unk / risk_g0.neg.unk
  
  
  ## ==========================================================
  ## NEW: unknown ever-positivity corrected quantities
  ##
  ## Construct three D by G tables:
  ##
  ##   obs_never:
  ##     U = 0 and M = 0
  ##
  ##   obs_ever:
  ##     U = 0 and M = 1
  ##
  ##   unknown_pos:
  ##     U = 1
  ##
  ## The true M values among U = 1 subjects are not used in this
  ## corrected analysis.
  ## ==========================================================
  
  obs_never.unk <- make_D_by_G_table_for_unknown(
    superpop[(U == 0) & (M == 0), , drop = FALSE]
  )
  
  obs_ever.unk <- make_D_by_G_table_for_unknown(
    superpop[(U == 0) & (M == 1), , drop = FALSE]
  )
  
  unknown_pos.unk <- make_D_by_G_table_for_unknown(
    superpop[(U == 1), , drop = FALSE]
  )
  
  unk.corr <- correct_positive_tables_simple(
    obs_never   = obs_never.unk,
    obs_ever    = obs_ever.unk,
    unknown_pos = unknown_pos.unk,
    screen_col  = "screen",
    control_col = "control",
    round_output = FALSE
  )
  
  corrected_never.unk <- unk.corr$corrected_never_positive_table
  corrected_ever.unk  <- unk.corr$corrected_ever_positive_table
  
  RR_pos_unk_corr_hat <- RR_from_D_by_G_table_for_unknown(corrected_ever.unk)
  RR_neg_unk_corr_hat <- RR_from_D_by_G_table_for_unknown(corrected_never.unk)
  
  
  ## ==========================================================
  ## Return all quantities
  ## ==========================================================
  
  out <- c(
    RR_all = RR_all_hat,
    RR_pos = RR_pos_hat,
    RR_neg = RR_neg_hat,
    RR_pos.loss = RR_pos_loss_hat,
    RR_neg.loss = RR_neg_loss_hat,
    RR_pos.loss.corr = RR_pos_loss_corr_hat,
    RR_pos.samp = RR_pos_samp_hat,
    RR_neg.samp = RR_neg_samp_hat,
    
    ## NEW: unknown ever-positivity
    RR_pos.unk = RR_pos_unk_hat,
    RR_neg.unk = RR_neg_unk_hat,
    RR_pos.unk.corr = RR_pos_unk_corr_hat,
    RR_neg.unk.corr = RR_neg_unk_corr_hat
  )
  
  return(out)
}



##
# Set data-generating parameters
##

set.seed(1)

n <- 100e3              # total sample size across both arms/groups
p_Mplus <- 0.3          # P(M+)=P(M+|G), independent of arm/group
p_g0 <- 0.5             # P(G=0); code below uses exactly n/2 in each arm
pD_g0 <- 0.02           # P(D+|G=0)

RR <- 0.9               # P(D+|G=1) / P(D+|G=0)
RR_neg <- 1             # P(D+|G=1,M-) / P(D+|G=0,M-)
RR_pos <- 0.86667       # P(D+|G=1,M+) / P(D+|G=0,M+)


##
# Loss-of-signal parameters
#
# These are retention-of-signal probabilities among those with M = 1.
#
# If M = 1 and D = 1:
#   P(Mtilde = 1 | M = 1, D = 1, G = 0)
# = P(Mtilde = 1 | M = 1, D = 1, G = 1)
# = p_Mtilde_pos_given_Mpos_Dpos
#
# If M = 1 and D = 0:
#   P(Mtilde = 1 | M = 1, D = 0, G = 0)
# = P(Mtilde = 1 | M = 1, D = 0, G = 1)
# = p_Mtilde_pos_given_Mpos_Dneg
#
# If M = 0:
#   Mtilde = 0 with probability 1.
##

p_Mtilde_pos_given_Mpos_Dpos <- 0.90
p_Mtilde_pos_given_Mpos_Dneg <- 0.80

if (p_Mtilde_pos_given_Mpos_Dpos < 0 || p_Mtilde_pos_given_Mpos_Dpos > 1) {
  stop("p_Mtilde_pos_given_Mpos_Dpos must be between 0 and 1")
}

if (p_Mtilde_pos_given_Mpos_Dneg < 0 || p_Mtilde_pos_given_Mpos_Dneg > 1) {
  stop("p_Mtilde_pos_given_Mpos_Dneg must be between 0 and 1")
}


##
# NEW: unknown ever-positivity parameters
#
# U = 1 means unknown ever-positivity M.
# U = 0 means known/observed ever-positivity M.
#
# These are the 4 probabilities:
#
#   p_U_unknown_given_Dpos_G0 = P(U=1 | D=1, G=0)
#   p_U_unknown_given_Dpos_G1 = P(U=1 | D=1, G=1)
#   p_U_unknown_given_Dneg_G0 = P(U=1 | D=0, G=0)
#   p_U_unknown_given_Dneg_G1 = P(U=1 | D=0, G=1)
#
# Defaults correspond to more unknown ever-positivity in the control arm
# than in the screen arm, with the same unknown probability for D+ and D-
# within arm. Change these as desired.
##

p_U_unknown_given_Dpos_G0 <- 0.30
p_U_unknown_given_Dpos_G1 <- 0.20
p_U_unknown_given_Dneg_G0 <- 0.30
p_U_unknown_given_Dneg_G1 <- 0.20

p_U_check <- c(
  p_U_unknown_given_Dpos_G0 = p_U_unknown_given_Dpos_G0,
  p_U_unknown_given_Dpos_G1 = p_U_unknown_given_Dpos_G1,
  p_U_unknown_given_Dneg_G0 = p_U_unknown_given_Dneg_G0,
  p_U_unknown_given_Dneg_G1 = p_U_unknown_given_Dneg_G1
)

if (any(p_U_check < 0 | p_U_check > 1)) {
  print(p_U_check)
  stop("All unknown-ever-positivity probabilities must be between 0 and 1.")
}


##
# Stratified sampling of control-arm specimens to reveal M and Mtilde
# f_1 is for D=1 and f_2 is for D=0
##
f_1 <- 0.95
f_2 <- 0.50



##
# Simulation study settings
##

S <- 1000   # number of independently simulated datasets
B <- 100  # number of bootstrap replicates per dataset


## 
# No need to change anything below, just run it given the above
##


##
# Calculate the 4 probabilities P(D+|M,G)
##

pD_g1 <- RR * pD_g0
RD <- pD_g0 - pD_g1

## P(D+|M-,G=0)
pD_g0.neg <- ifelse(
  RR_pos - RR_neg == 0,
  pD_g0,
  (RR_pos * pD_g0 - pD_g1) / ((1 - p_Mplus) * (RR_pos - RR_neg))
)

## P(D+|M-,G=1)
pD_g1.neg <- RR_neg * pD_g0.neg

## P(D+|M+,G=0)
pD_g0.pos <- ifelse(
  RR_neg - RR_pos == 0,
  pD_g0,
  (RR_neg * pD_g0 - pD_g1) / (p_Mplus * (RR_neg - RR_pos))
)

## P(D+|M+,G=1)
pD_g1.pos <- RR_pos * pD_g0.pos


##
# Check probability validity
##

prob_check <- c(
  pD_g0.neg = pD_g0.neg,
  pD_g1.neg = pD_g1.neg,
  pD_g0.pos = pD_g0.pos,
  pD_g1.pos = pD_g1.pos
)

if (any(prob_check < 0 | prob_check > 1)) {
  print(prob_check)
  stop("At least one P(D+|M,G) probability is outside [0,1].")
}


##
# Create matrix of the 4 probabilities P(D+|M,G)
# Rows are G = 0, 1; columns are M = 0, 1
##

p_Dplus_given_GM <- matrix(
  c(
    pD_g0.neg, pD_g0.pos,   # G = 0, M = 0/1
    pD_g1.neg, pD_g1.pos    # G = 1, M = 0/1
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    G = c("0", "1"),
    M = c("0", "1")
  )
)


##
# Calculate true loss-of-signal parameters analytically
#
# These are the true values of the Mtilde-based estimands.
#
# For Mtilde = 1:
#
#   P(D=1 | G=g, Mtilde=1)
#   =
#   P(D=1, Mtilde=1 | G=g) / P(Mtilde=1 | G=g)
#
# Since Mtilde=1 only possible among M=1:
#
#   P(D=1, Mtilde=1 | G=g)
#   =
#   P(M=1) P(D=1 | G=g,M=1) P(Mtilde=1 | M=1,D=1)
#
# For Mtilde = 0, this includes all M=0 plus M=1 individuals
# whose signal is not retained.
##

true_loss_params <- function(
    p_Mplus,
    pD_g0.neg,
    pD_g1.neg,
    pD_g0.pos,
    pD_g1.pos,
    p_Mtilde_pos_given_Mpos_Dpos,
    p_Mtilde_pos_given_Mpos_Dneg
) {
  
  q_Dpos <- p_Mtilde_pos_given_Mpos_Dpos  # = P1(Mtilde+ | M+, D+)
  q_Dneg <- p_Mtilde_pos_given_Mpos_Dneg  # = P1(Mtilde+ | M+, D-)
  
  ## ---- Mtilde = 1 ("pos.loss") ----
  
  ## For G = 0, conditional on Mtilde=1
  numer_g0 <- p_Mplus * pD_g0.pos * q_Dpos
  denom_g0 <- p_Mplus * (pD_g0.pos * q_Dpos + (1 - pD_g0.pos) * q_Dneg)
  risk_g0_pos_loss <- numer_g0 / denom_g0
  
  ## For G = 1, conditional on Mtilde=1
  numer_g1 <- p_Mplus * pD_g1.pos * q_Dpos
  denom_g1 <- p_Mplus * (pD_g1.pos * q_Dpos + (1 - pD_g1.pos) * q_Dneg)
  risk_g1_pos_loss <- numer_g1 / denom_g1
  
  RR_pos.loss <- risk_g1_pos_loss / risk_g0_pos_loss
  
  
  ## ---- Mtilde = 0 ("neg.loss") ----
  
  ## For G = 0, conditional on Mtilde=0
  numer_g0 <- (1 - p_Mplus) * pD_g0.neg +
    p_Mplus * pD_g0.pos * (1 - q_Dpos)
  
  denom_g0 <- (1 - p_Mplus) +
    p_Mplus * (
      pD_g0.pos * (1 - q_Dpos) +
        (1 - pD_g0.pos) * (1 - q_Dneg)
    )
  
  risk_g0_neg_loss <- numer_g0 / denom_g0
  
  ## For G = 1, conditional on Mtilde=0
  numer_g1 <- (1 - p_Mplus) * pD_g1.neg +
    p_Mplus * pD_g1.pos * (1 - q_Dpos)
  
  denom_g1 <- (1 - p_Mplus) +
    p_Mplus * (
      pD_g1.pos * (1 - q_Dpos) +
        (1 - pD_g1.pos) * (1 - q_Dneg)
    )
  
  risk_g1_neg_loss <- numer_g1 / denom_g1
  
  RR_neg.loss <- risk_g1_neg_loss / risk_g0_neg_loss
  
  
  ## ---- TRUE RR_pos.loss.corr using your page-9 correction formula ----
  ## P0(D+|M+) = {1 + (P0(D-)/P0(D+)) *
  ##               ( P0(Mtilde+|D-) / P1(Mtilde+|M+,D-) ) /
  ##               ( P0(Mtilde+|D+) / P1(Mtilde+|M+,D+) )
  ##             }^(-1)
  ## RR_pos.loss.corr = P1(D+|M+) / P0(D+|M+)_(corr)
  
  ## Control-arm marginal outcome probs
  P0_Dpos <- p_Mplus * pD_g0.pos + (1 - p_Mplus) * pD_g0.neg
  P0_Dneg <- 1 - P0_Dpos
  
  ## Control-arm P(Mtilde+ | D+) and P(Mtilde+ | D-)
  ## Since Mtilde+ can only occur if M=1:
  ## P0(Mtilde+|D+) = P0(M=1|D+) * q_Dpos
  ## P0(Mtilde+|D-) = P0(M=1|D-) * q_Dneg
  
  P0_Mpos_given_Dpos <- (p_Mplus * pD_g0.pos) / P0_Dpos
  P0_Mpos_given_Dneg <- (p_Mplus * (1 - pD_g0.pos)) / P0_Dneg
  
  P0_Mt_pos_given_Dpos <- P0_Mpos_given_Dpos * q_Dpos
  P0_Mt_pos_given_Dneg <- P0_Mpos_given_Dneg * q_Dneg
  
  ## Screen-arm retest-positive fractions among true M+ by outcome
  ## Under your DGP these equal q_Dpos and q_Dneg:
  P1_Mt_pos_given_Mpos_Dpos <- q_Dpos
  P1_Mt_pos_given_Mpos_Dneg <- q_Dneg
  
  ratio_term <- (P0_Dneg / P0_Dpos) *
    ( (P0_Mt_pos_given_Dneg / P1_Mt_pos_given_Mpos_Dneg) /
        (P0_Mt_pos_given_Dpos / P1_Mt_pos_given_Mpos_Dpos) )
  
  P0_Dpos_given_Mpos_corr <- 1 / (1 + ratio_term)
  
  ## Numerator is the true screen-arm risk among true M+:
  P1_Dpos_given_Mpos <- pD_g1.pos
  
  RR_pos.loss.corr <- P1_Dpos_given_Mpos / P0_Dpos_given_Mpos_corr
  
  
  ## ---- Return as strictly named numeric vector ----
  out <- c(RR_pos.loss, RR_neg.loss, RR_pos.loss.corr)
  names(out) <- c("RR_pos.loss", "RR_neg.loss", "RR_pos.loss.corr")
  
  return(out)
}



true_loss <- true_loss_params(
  p_Mplus = p_Mplus,
  pD_g0.neg = pD_g0.neg,
  pD_g1.neg = pD_g1.neg,
  pD_g0.pos = pD_g0.pos,
  pD_g1.pos = pD_g1.pos,
  p_Mtilde_pos_given_Mpos_Dpos = p_Mtilde_pos_given_Mpos_Dpos,
  p_Mtilde_pos_given_Mpos_Dneg = p_Mtilde_pos_given_Mpos_Dneg
)


##
# True parameters corresponding exactly to columns returned by IE_fast()
##

true_par <- c(
  RR_all = RR,
  RR_pos = RR_pos,
  RR_neg = RR_neg,
  RR_pos.loss = unname(true_loss["RR_pos.loss"]),
  RR_neg.loss = unname(true_loss["RR_neg.loss"]),
  RR_pos.loss.corr = unname(true_loss["RR_pos.loss.corr"]),
  RR_pos.samp = RR_pos,
  RR_neg.samp = RR_neg,
  
  ## NEW: unknown ever-positivity
  ##
  ## These use the gold-standard targets in true_par.
  ## The actual coverage target is explicitly set below in true_for_coverage.
  RR_pos.unk = RR_pos,
  RR_neg.unk = RR_neg,
  RR_pos.unk.corr = RR_pos,
  RR_neg.unk.corr = RR_neg
)

print(true_par)

## 1) Create the truth vector USED FOR COVERAGE (bias-to-gold-standard check)
true_for_coverage <- true_par

## Force loss-based and sampling-based quantities to be judged vs gold-standard targets:
true_for_coverage["RR_pos.loss"]      <- true_par["RR_pos"]
true_for_coverage["RR_pos.loss.corr"] <- true_par["RR_pos"]
true_for_coverage["RR_pos.samp"]      <- true_par["RR_pos"]
true_for_coverage["RR_neg.loss"]      <- true_par["RR_neg"]
true_for_coverage["RR_neg.samp"]      <- true_par["RR_neg"]

## NEW: unknown ever-positivity
##
## For coverage:
##   RR_pos.unk and RR_pos.unk.corr use true RR_pos.
##   RR_neg.unk and RR_neg.unk.corr use true RR_neg.
true_for_coverage["RR_pos.unk"]       <- true_par["RR_pos"]
true_for_coverage["RR_pos.unk.corr"]  <- true_par["RR_pos"]
true_for_coverage["RR_neg.unk"]       <- true_par["RR_neg"]
true_for_coverage["RR_neg.unk.corr"]  <- true_par["RR_neg"]

stopifnot(all(names(true_for_coverage) == names(true_par)))



##
# Function to simulate one dataset superpop from the same distribution
# Now includes Mtilde.
#
# NEW: unknown ever-positivity
# Now also includes U.
##

simulate_superpop_fast <- function(
    n,
    p_Mplus,
    p_Dplus_given_GM,
    p_Mtilde_pos_given_Mpos_Dpos,
    p_Mtilde_pos_given_Mpos_Dneg,
    
    ## NEW: unknown ever-positivity
    p_U_unknown_given_Dpos_G0 = p_U_unknown_given_Dpos_G0,
    p_U_unknown_given_Dpos_G1 = p_U_unknown_given_Dpos_G1,
    p_U_unknown_given_Dneg_G0 = p_U_unknown_given_Dneg_G0,
    p_U_unknown_given_Dneg_G1 = p_U_unknown_given_Dneg_G1,
    
    f1 = f_1,
    f2 = f_2
) {
  
  if (n %% 2 != 0) {
    stop("n must be even because this function assigns exactly n/2 to each G arm.")
  }
  
  if (f1 < 0 || f1 > 1) {
    stop("f1 must be between 0 and 1.")
  }
  
  if (f2 < 0 || f2 > 1) {
    stop("f2 must be between 0 and 1.")
  }
  
  ## NEW: unknown ever-positivity
  p_U_check_local <- c(
    p_U_unknown_given_Dpos_G0 = p_U_unknown_given_Dpos_G0,
    p_U_unknown_given_Dpos_G1 = p_U_unknown_given_Dpos_G1,
    p_U_unknown_given_Dneg_G0 = p_U_unknown_given_Dneg_G0,
    p_U_unknown_given_Dneg_G1 = p_U_unknown_given_Dneg_G1
  )
  
  if (any(p_U_check_local < 0 | p_U_check_local > 1)) {
    print(p_U_check_local)
    stop("All unknown-ever-positivity probabilities must be between 0 and 1.")
  }
  
  ## Generate arm/group G with exactly n/2 in G=0 and n/2 in G=1
  G <- sample(c(rep(0, n / 2), rep(1, n / 2)), replace = FALSE)
  
  ## Generate gold-standard M based on P(M+)
  M <- rbinom(n, size = 1, prob = p_Mplus)
  
  ## Generate D based on P(D+|M,G)
  p_Dplus <- p_Dplus_given_GM[cbind(as.character(G), as.character(M))]
  D <- rbinom(n, size = 1, prob = p_Dplus)
  
  ## Generate Mtilde based on retention of signal among those with M = 1
  ##
  ## Everyone with M = 0 has Mtilde = 0.
  ## Among those with M = 1:
  ##   if D = 1, Mtilde ~ Bernoulli(p_Mtilde_pos_given_Mpos_Dpos)
  ##   if D = 0, Mtilde ~ Bernoulli(p_Mtilde_pos_given_Mpos_Dneg)
  
  p_Mtilde_pos <- rep(0, n)
  
  p_Mtilde_pos[M == 1 & D == 1] <- p_Mtilde_pos_given_Mpos_Dpos
  p_Mtilde_pos[M == 1 & D == 0] <- p_Mtilde_pos_given_Mpos_Dneg
  
  Mtilde <- rbinom(n, size = 1, prob = p_Mtilde_pos)
  
  ## Safety check: Mtilde+ can only occur among M+
  if (any(Mtilde == 1 & M == 0)) {
    stop("Error: Mtilde=1 occurred among observations with M=0.")
  }
  
  ## NEW: unknown ever-positivity
  ##
  ## Generate U using 4 probabilities depending on D and G.
  ##
  ## U = 1 means unknown ever-positivity M.
  ## U = 0 means known/observed ever-positivity M.
  
  p_U_unknown <- rep(NA_real_, n)
  
  p_U_unknown[D == 1 & G == 0] <- p_U_unknown_given_Dpos_G0
  p_U_unknown[D == 1 & G == 1] <- p_U_unknown_given_Dpos_G1
  p_U_unknown[D == 0 & G == 0] <- p_U_unknown_given_Dneg_G0
  p_U_unknown[D == 0 & G == 1] <- p_U_unknown_given_Dneg_G1
  
  if (any(is.na(p_U_unknown))) {
    stop("Error: at least one p_U_unknown value is NA.")
  }
  
  U <- rbinom(n, size = 1, prob = p_U_unknown)
  
  ## Generate reveal indicator R
  ##
  ## G=1 screen arm: M and Mtilde are observed/revealed for everyone.
  ## G=0 control arm:
  ##   D=1 stratum: reveal with probability f1
  ##   D=0 stratum: reveal with probability f2
  
  R <- integer(n)
  
  R[G == 1] <- 1
  
  idx_g0_d1 <- which(G == 0 & D == 1)
  idx_g0_d0 <- which(G == 0 & D == 0)
  
  R[idx_g0_d1] <- rbinom(length(idx_g0_d1), size = 1, prob = f1)
  R[idx_g0_d0] <- rbinom(length(idx_g0_d0), size = 1, prob = f2)
  
  ## Put dataset together
  superpop <- cbind(
    G = G,
    D = D,
    M = M,
    Mtilde = Mtilde,
    R = R,
    
    ## NEW: unknown ever-positivity
    U = U
  )
  
  return(superpop)
}


##
# Function to bootstrap one dataset and summarize IE output
##

bootstrap_IE_summary_fast <- function(superpop, B = 200, f1 = f_1, f2 = f_2) {
  
  ## Original-sample parameter estimates
  IE.orig <- IE_fast(superpop, f1 = f1, f2 = f2)
  
  par_names <- names(IE.orig)
  n_par <- length(IE.orig)
  
  ## Store bootstrap estimates in a numeric matrix
  IE.boot <- matrix(
    NA_real_,
    nrow = B,
    ncol = n_par,
    dimnames = list(NULL, par_names)
  )
  
  ## Precompute G-specific row indices once for this dataset
  idx_g0 <- which(superpop[, "G"] == 0)
  idx_g1 <- which(superpop[, "G"] == 1)
  
  n_g0 <- length(idx_g0)
  n_g1 <- length(idx_g1)
  
  for (b in seq_len(B)) {
    
    boot_idx_g0 <- sample(idx_g0, size = n_g0, replace = TRUE)
    boot_idx_g1 <- sample(idx_g1, size = n_g1, replace = TRUE)
    
    boot_idx <- c(boot_idx_g0, boot_idx_g1)
    
    superpop.boot <- superpop[boot_idx, , drop = FALSE]
    
    IE.boot[b, ] <- IE_fast(superpop.boot, f1 = f1, f2 = f2)
  }
  
  ## Bootstrap mean and SD for each parameter
  boot.mean <- colMeans(IE.boot, na.rm = TRUE)
  boot.sd <- apply(IE.boot, 2, sd, na.rm = TRUE)
  
  ## Normal-approximation 95% CI using bootstrap SD
  ci.lower <- IE.orig - 1.96 * boot.sd
  ci.upper <- IE.orig + 1.96 * boot.sd
  
  out <- data.frame(
    parameter = par_names,
    estimate = as.numeric(IE.orig),
    boot.mean = as.numeric(boot.mean),
    boot.sd = as.numeric(boot.sd),
    ci.lower = as.numeric(ci.lower),
    ci.upper = as.numeric(ci.upper),
    row.names = par_names
  )
  
  return(out)
}


# Names for all parameters
par_names <- names(true_par)
n_par <- length(true_par)


##
# Preallocate matrices for simulation results
##

estimate.mat <- matrix(
  NA_real_,
  nrow = S,
  ncol = n_par,
  dimnames = list(NULL, par_names)
)

boot.mean.mat <- matrix(
  NA_real_,
  nrow = S,
  ncol = n_par,
  dimnames = list(NULL, par_names)
)

boot.sd.mat <- matrix(
  NA_real_,
  nrow = S,
  ncol = n_par,
  dimnames = list(NULL, par_names)
)

ci.lower.mat <- matrix(
  NA_real_,
  nrow = S,
  ncol = n_par,
  dimnames = list(NULL, par_names)
)

ci.upper.mat <- matrix(
  NA_real_,
  nrow = S,
  ncol = n_par,
  dimnames = list(NULL, par_names)
)

covered.mat <- matrix(
  FALSE,
  nrow = S,
  ncol = n_par,
  dimnames = list(NULL, par_names)
)


##
# Run and time only the simulation loop
##

time.sim <- system.time({
  
  for (s in seq_len(S)) {
    
    superpop.s <- simulate_superpop_fast(
      n = n,
      p_Mplus = p_Mplus,
      p_Dplus_given_GM = p_Dplus_given_GM,
      p_Mtilde_pos_given_Mpos_Dpos = p_Mtilde_pos_given_Mpos_Dpos,
      p_Mtilde_pos_given_Mpos_Dneg = p_Mtilde_pos_given_Mpos_Dneg,
      
      ## NEW: unknown ever-positivity
      p_U_unknown_given_Dpos_G0 = p_U_unknown_given_Dpos_G0,
      p_U_unknown_given_Dpos_G1 = p_U_unknown_given_Dpos_G1,
      p_U_unknown_given_Dneg_G0 = p_U_unknown_given_Dneg_G0,
      p_U_unknown_given_Dneg_G1 = p_U_unknown_given_Dneg_G1,
      
      f1 = f_1,
      f2 = f_2
    )
    
    ## Optional sanity checks
    if (any(superpop.s[, "G"] == 1 & superpop.s[, "R"] != 1)) {
      stop("Unexpected result: Some G=1 observations have R != 1.")
    }
    
    if (mean(superpop.s[, "Mtilde"] == 1) > mean(superpop.s[, "M"] == 1)) {
      stop("Unexpected result: P(Mtilde+) > P(M+).")
    }
    
    ## NEW: unknown ever-positivity sanity check
    if (any(!(superpop.s[, "U"] %in% c(0, 1)))) {
      stop("Unexpected result: U contains values other than 0 and 1.")
    }
    
    boot.summary.s <- bootstrap_IE_summary_fast(
      superpop = superpop.s,
      B = B,
      f1 = f_1,
      f2 = f_2
    )
    
    ## Ensure returned parameter order matches preallocated matrices
    boot.summary.s <- boot.summary.s[par_names, , drop = FALSE]
    
    estimate.mat[s, ]  <- boot.summary.s$estimate
    boot.mean.mat[s, ] <- boot.summary.s$boot.mean
    boot.sd.mat[s, ]   <- boot.summary.s$boot.sd
    ci.lower.mat[s, ]  <- boot.summary.s$ci.lower
    ci.upper.mat[s, ]  <- boot.summary.s$ci.upper
    
    covered.mat[s, ] <- ci.lower.mat[s, ] <= true_for_coverage &
      true_for_coverage <= ci.upper.mat[s, ]
    
    if (s %% 10 == 0) {
      cat("Finished simulated dataset", s, "of", S, "\n")
    }
  }
})

system2("say", "Finished!", wait = FALSE)


IE.sim.summary <- data.frame(
  parameter = par_names,
  true = as.numeric(true_par),
  true.for.coverage = as.numeric(true_for_coverage),
  
  mean.estimate = colMeans(estimate.mat, na.rm = TRUE),
  sd.estimate = apply(estimate.mat, 2, sd, na.rm = TRUE),
  
  mean.boot.mean = colMeans(boot.mean.mat, na.rm = TRUE),
  mean.boot.sd = colMeans(boot.sd.mat, na.rm = TRUE),
  
  mean.ci.lower = colMeans(ci.lower.mat, na.rm = TRUE),
  mean.ci.upper = colMeans(ci.upper.mat, na.rm = TRUE),
  
  coverage = colMeans(covered.mat, na.rm = TRUE),
  n.sims = colSums(!is.na(covered.mat)),
  
  row.names = NULL
)

##
# Optional: create long-format simulation results
# This is not needed for IE.sim.summary, but can be useful for diagnostics.
##

sim.results <- do.call(
  rbind,
  lapply(seq_len(S), function(s) {
    data.frame(
      sim = s,
      parameter = par_names,
      true = as.numeric(true_par),
      estimate = estimate.mat[s, ],
      boot.mean = boot.mean.mat[s, ],
      boot.sd = boot.sd.mat[s, ],
      ci.lower = ci.lower.mat[s, ],
      ci.upper = ci.upper.mat[s, ],
      covered = covered.mat[s, ],
      row.names = NULL
    )
  })
)



##
# Print outputs
##

time.sim

true_par

true_for_coverage

IE.sim.summary

head(sim.results)












# #####
# ## OLD CODE BELOW all commented out
# #####
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##
# # Fast simulation + nonparametric bootstrap by G
# # Now includes loss of signal through Mtilde.
# #
# # The original/gold-standard quantities use M:
# #   RR_all, RR_pos, RR_neg
# #
# # The loss-of-signal quantities use Mtilde:
# #   RR_pos.loss, RR_neg.loss
# #
# # Mtilde is generated as follows:
# #   - If M = 0, then Mtilde = 0.
# #   - If M = 1 and D = 1, then Mtilde ~ Bernoulli(p_Mtilde_pos_given_Mpos_Dpos).
# #   - If M = 1 and D = 0, then Mtilde ~ Bernoulli(p_Mtilde_pos_given_Mpos_Dneg).
# #   - These retention probabilities do not depend on G.
# ##
# 
# ##
# # Bootstrap one dataset, preserving the original number of G=0 and G=1
# ##
# 
# bootstrap_superpop_by_G <- function(superpop) {
#   
#   idx_g0 <- which(superpop[, "G"] == 0)
#   idx_g1 <- which(superpop[, "G"] == 1)
#   
#   boot_idx_g0 <- sample(idx_g0, size = length(idx_g0), replace = TRUE)
#   boot_idx_g1 <- sample(idx_g1, size = length(idx_g1), replace = TRUE)
#   
#   boot_idx <- c(boot_idx_g0, boot_idx_g1)
#   
#   superpop.boot <- superpop[boot_idx, , drop = FALSE]
#   
#   return(superpop.boot)
# }
# 
# 
# ##
# # Fast IE function
# # Computes RR_all, RR_pos, RR_neg using M,
# # and RR_pos.loss, RR_neg.loss using Mtilde.
# ##
# 
# IE_fast <- function(superpop, f1 = f_1, f2 = f_2) {
#   
#   ## Expect columns: G, D, M, Mtilde, R
#   G      <- superpop[, "G"]
#   D      <- superpop[, "D"]
#   M      <- superpop[, "M"]
#   Mtilde <- superpop[, "Mtilde"]
#   R      <- superpop[, "R"]
#   
#   ## ---- small helper: safe mean in a subset ----
#   mean_if <- function(x, idx) {
#     if (!any(idx)) return(NA_real_)
#     mean(x[idx], na.rm = TRUE)
#   }
#   
#   ## ---- small helper: safe weighted risk ----
#   weighted_risk <- function(D, idx, w) {
#     if (!any(idx)) return(NA_real_)
#     denom <- sum(w[idx], na.rm = TRUE)
#     if (is.na(denom) || denom == 0) return(NA_real_)
#     sum(w[idx] * D[idx], na.rm = TRUE) / denom
#   }
#   
#   
#   ## ==========================================================
#   ## Gold-standard quantities using M
#   ## ==========================================================
#   
#   idx_neg <- (M == 0)
#   idx_pos <- (M == 1)
#   
#   risk_g1.neg <- mean_if(D == 1, (G == 1) & idx_neg)
#   risk_g0.neg <- mean_if(D == 1, (G == 0) & idx_neg)
#   
#   risk_g1.pos <- mean_if(D == 1, (G == 1) & idx_pos)
#   risk_g0.pos <- mean_if(D == 1, (G == 0) & idx_pos)
#   
#   risk_g1.all <- mean_if(D == 1, (G == 1))
#   risk_g0.all <- mean_if(D == 1, (G == 0))
#   
#   RR_all_hat <- risk_g1.all / risk_g0.all
#   RR_pos_hat <- risk_g1.pos / risk_g0.pos
#   RR_neg_hat <- risk_g1.neg / risk_g0.neg
#   
#   
#   ## ==========================================================
#   ## Loss-of-signal naive quantities using Mtilde
#   ## ==========================================================
#   
#   idx_neg.loss <- (Mtilde == 0)
#   idx_pos.loss <- (Mtilde == 1)
#   
#   risk_g1.neg.loss <- mean_if(D == 1, (G == 1) & idx_neg.loss)
#   risk_g0.neg.loss <- mean_if(D == 1, (G == 0) & idx_neg.loss)
#   
#   risk_g1.pos.loss <- mean_if(D == 1, (G == 1) & idx_pos.loss)
#   risk_g0.pos.loss <- mean_if(D == 1, (G == 0) & idx_pos.loss)
#   
#   RR_pos_loss_hat <- risk_g1.pos.loss / risk_g0.pos.loss
#   RR_neg_loss_hat <- risk_g1.neg.loss / risk_g0.neg.loss
#   
#   
#   ## ==========================================================
#   ## Bias-corrected RRpos under loss-of-signal
#   ##
#   ## Formula:
#   ## P0(D+|M+) =
#   ##   {1 + [P0(D-)/P0(D+)] *
#   ##        [P0(Mtilde+|D-) / P1(Mtilde+|M+,D-)] /
#   ##        [P0(Mtilde+|D+) / P1(Mtilde+|M+,D+)]
#   ##   }^(-1)
#   ##
#   ## RR_pos.loss.corr = P1(D+|M+) / corrected P0(D+|M+)
#   ##
#   ## This uses the full-data version, not the sampling-weighted version.
#   ## ==========================================================
#   
#   P0_Dpos <- mean_if(D == 1, (G == 0))
#   P0_Dneg <- mean_if(D == 0, (G == 0))
#   
#   P0_Mt_pos_given_Dpos <- mean_if(Mtilde == 1, (G == 0) & (D == 1))
#   P0_Mt_pos_given_Dneg <- mean_if(Mtilde == 1, (G == 0) & (D == 0))
#   
#   P1_Mt_pos_given_Mpos_Dpos <- mean_if(Mtilde == 1, (G == 1) & (M == 1) & (D == 1))
#   P1_Mt_pos_given_Mpos_Dneg <- mean_if(Mtilde == 1, (G == 1) & (M == 1) & (D == 0))
#   
#   needed_corr <- c(
#     P0_Dpos,
#     P0_Dneg,
#     P0_Mt_pos_given_Dpos,
#     P0_Mt_pos_given_Dneg,
#     P1_Mt_pos_given_Mpos_Dpos,
#     P1_Mt_pos_given_Mpos_Dneg
#   )
#   
#   if (any(is.na(needed_corr)) ||
#       P0_Dpos == 0 ||
#       P0_Mt_pos_given_Dpos == 0 ||
#       P1_Mt_pos_given_Mpos_Dpos == 0 ||
#       P1_Mt_pos_given_Mpos_Dneg == 0) {
#     
#     RR_pos_loss_corr_hat <- NA_real_
#     
#   } else {
#     
#     ratio_term <- (P0_Dneg / P0_Dpos) *
#       (
#         (P0_Mt_pos_given_Dneg / P1_Mt_pos_given_Mpos_Dneg) /
#           (P0_Mt_pos_given_Dpos / P1_Mt_pos_given_Mpos_Dpos)
#       )
#     
#     P0_Dpos_given_Mpos_corr <- 1 / (1 + ratio_term)
#     
#     P1_Dpos_given_Mpos <- risk_g1.pos
#     
#     RR_pos_loss_corr_hat <- P1_Dpos_given_Mpos / P0_Dpos_given_Mpos_corr
#   }
#   
#   
#   ## ==========================================================
#   ## Stratified-sampling estimators in control arm
#   ##
#   ## In G=0, only use revealed observations R=1.
#   ## Weight by inverse reveal fraction:
#   ##   D=1 stratum: 1/f1
#   ##   D=0 stratum: 1/f2
#   ##
#   ## RR_pos.samp:
#   ##   numerator = P1(D+|M+) from fully observed G=1 screen arm
#   ##   denominator = weighted P0(D+|M+) from sampled/revealed G=0 controls
#   ##
#   ## RR_neg.samp:
#   ##   numerator = P1(D+|M-) from fully observed G=1 screen arm
#   ##   denominator = weighted P0(D+|M-) from sampled/revealed G=0 controls
#   ## ==========================================================
#   
#   if (f1 <= 0 || f2 <= 0) {
#     
#     RR_pos_samp_hat <- NA_real_
#     RR_neg_samp_hat <- NA_real_
#     
#   } else {
#     
#     w0 <- rep(0, length(D))
#     
#     w0[G == 0 & R == 1 & D == 1] <- 1 / f1
#     w0[G == 0 & R == 1 & D == 0] <- 1 / f2
#     
#     ## Weighted control-arm risks among true M+ and true M-
#     risk_g0.pos.samp <- weighted_risk(
#       D = D,
#       idx = (G == 0) & (R == 1) & (M == 1),
#       w = w0
#     )
#     
#     risk_g0.neg.samp <- weighted_risk(
#       D = D,
#       idx = (G == 0) & (R == 1) & (M == 0),
#       w = w0
#     )
#     
#     ## Screen-arm numerator remains fully observed
#     risk_g1.pos.samp <- risk_g1.pos
#     risk_g1.neg.samp <- risk_g1.neg
#     
#     RR_pos_samp_hat <- risk_g1.pos.samp / risk_g0.pos.samp
#     RR_neg_samp_hat <- risk_g1.neg.samp / risk_g0.neg.samp
#   }
#   
#   
#   ## ==========================================================
#   ## Return all quantities
#   ## ==========================================================
#   
#   out <- c(
#     RR_all = RR_all_hat,
#     RR_pos = RR_pos_hat,
#     RR_neg = RR_neg_hat,
#     RR_pos.loss = RR_pos_loss_hat,
#     RR_neg.loss = RR_neg_loss_hat,
#     RR_pos.loss.corr = RR_pos_loss_corr_hat,
#     RR_pos.samp = RR_pos_samp_hat,
#     RR_neg.samp = RR_neg_samp_hat
#   )
#   
#   return(out)
# }
# 
# 
# 
# ##
# # Set data-generating parameters
# ##
# 
# set.seed(1)
# 
# n <- 100e3              # total sample size across both arms/groups
# p_Mplus <- 0.3          # P(M+)=P(M+|G), independent of arm/group
# p_g0 <- 0.5             # P(G=0); code below uses exactly n/2 in each arm
# pD_g0 <- 0.02           # P(D+|G=0)
# 
# RR <- 0.9               # P(D+|G=1) / P(D+|G=0)
# RR_neg <- 1             # P(D+|G=1,M-) / P(D+|G=0,M-)
# RR_pos <- 0.86667       # P(D+|G=1,M+) / P(D+|G=0,M+)
# 
# 
# ##
# # Loss-of-signal parameters
# #
# # These are retention-of-signal probabilities among those with M = 1.
# #
# # If M = 1 and D = 1:
# #   P(Mtilde = 1 | M = 1, D = 1, G = 0)
# # = P(Mtilde = 1 | M = 1, D = 1, G = 1)
# # = p_Mtilde_pos_given_Mpos_Dpos
# #
# # If M = 1 and D = 0:
# #   P(Mtilde = 1 | M = 1, D = 0, G = 0)
# # = P(Mtilde = 1 | M = 1, D = 0, G = 1)
# # = p_Mtilde_pos_given_Mpos_Dneg
# #
# # If M = 0:
# #   Mtilde = 0 with probability 1.
# ##
# 
# p_Mtilde_pos_given_Mpos_Dpos <- 1
# p_Mtilde_pos_given_Mpos_Dneg <- 0.50
# 
# if (p_Mtilde_pos_given_Mpos_Dpos < 0 || p_Mtilde_pos_given_Mpos_Dpos > 1) {
#   stop("p_Mtilde_pos_given_Mpos_Dpos must be between 0 and 1")
# }
# 
# if (p_Mtilde_pos_given_Mpos_Dneg < 0 || p_Mtilde_pos_given_Mpos_Dneg > 1) {
#   stop("p_Mtilde_pos_given_Mpos_Dneg must be between 0 and 1")
# }
# 
# ##
# # Stratified sampling of control-arm specimens to reveal M and Mtilde
# # f_1 is for D=1 and f_2 is for D=0
# ##
# f_1 <- 0.95
# f_2 <- 0.50
# 
# 
# ##
# # Calculate the 4 probabilities P(D+|M,G)
# ##
# 
# pD_g1 <- RR * pD_g0
# RD <- pD_g0 - pD_g1
# 
# ## P(D+|M-,G=0)
# pD_g0.neg <- ifelse(
#   RR_pos - RR_neg == 0,
#   pD_g0,
#   (RR_pos * pD_g0 - pD_g1) / ((1 - p_Mplus) * (RR_pos - RR_neg))
# )
# 
# ## P(D+|M-,G=1)
# pD_g1.neg <- RR_neg * pD_g0.neg
# 
# ## P(D+|M+,G=0)
# pD_g0.pos <- ifelse(
#   RR_neg - RR_pos == 0,
#   pD_g0,
#   (RR_neg * pD_g0 - pD_g1) / (p_Mplus * (RR_neg - RR_pos))
# )
# 
# ## P(D+|M+,G=1)
# pD_g1.pos <- RR_pos * pD_g0.pos
# 
# 
# ##
# # Check probability validity
# ##
# 
# prob_check <- c(
#   pD_g0.neg = pD_g0.neg,
#   pD_g1.neg = pD_g1.neg,
#   pD_g0.pos = pD_g0.pos,
#   pD_g1.pos = pD_g1.pos
# )
# 
# if (any(prob_check < 0 | prob_check > 1)) {
#   print(prob_check)
#   stop("At least one P(D+|M,G) probability is outside [0,1].")
# }
# 
# 
# ##
# # Create matrix of the 4 probabilities P(D+|M,G)
# # Rows are G = 0, 1; columns are M = 0, 1
# ##
# 
# p_Dplus_given_GM <- matrix(
#   c(
#     pD_g0.neg, pD_g0.pos,   # G = 0, M = 0/1
#     pD_g1.neg, pD_g1.pos    # G = 1, M = 0/1
#   ),
#   nrow = 2,
#   byrow = TRUE,
#   dimnames = list(
#     G = c("0", "1"),
#     M = c("0", "1")
#   )
# )
# 
# 
# ##
# # Calculate true loss-of-signal parameters analytically
# #
# # These are the true values of the Mtilde-based estimands.
# #
# # For Mtilde = 1:
# #
# #   P(D=1 | G=g, Mtilde=1)
# #   =
# #   P(D=1, Mtilde=1 | G=g) / P(Mtilde=1 | G=g)
# #
# # Since Mtilde=1 only possible among M=1:
# #
# #   P(D=1, Mtilde=1 | G=g)
# #   =
# #   P(M=1) P(D=1 | G=g,M=1) P(Mtilde=1 | M=1,D=1)
# #
# # For Mtilde = 0, this includes all M=0 plus M=1 individuals
# # whose signal is not retained.
# ##
# 
# true_loss_params <- function(
#     p_Mplus,
#     pD_g0.neg,
#     pD_g1.neg,
#     pD_g0.pos,
#     pD_g1.pos,
#     p_Mtilde_pos_given_Mpos_Dpos,
#     p_Mtilde_pos_given_Mpos_Dneg
# ) {
#   
#   q_Dpos <- p_Mtilde_pos_given_Mpos_Dpos  # = P1(Mtilde+ | M+, D+)
#   q_Dneg <- p_Mtilde_pos_given_Mpos_Dneg  # = P1(Mtilde+ | M+, D-)
#   
#   ## ---- Mtilde = 1 ("pos.loss") ----
#   
#   ## For G = 0, conditional on Mtilde=1
#   numer_g0 <- p_Mplus * pD_g0.pos * q_Dpos
#   denom_g0 <- p_Mplus * (pD_g0.pos * q_Dpos + (1 - pD_g0.pos) * q_Dneg)
#   risk_g0_pos_loss <- numer_g0 / denom_g0
#   
#   ## For G = 1, conditional on Mtilde=1
#   numer_g1 <- p_Mplus * pD_g1.pos * q_Dpos
#   denom_g1 <- p_Mplus * (pD_g1.pos * q_Dpos + (1 - pD_g1.pos) * q_Dneg)
#   risk_g1_pos_loss <- numer_g1 / denom_g1
#   
#   RR_pos.loss <- risk_g1_pos_loss / risk_g0_pos_loss
#   
#   
#   ## ---- Mtilde = 0 ("neg.loss") ----
#   
#   ## For G = 0, conditional on Mtilde=0
#   numer_g0 <- (1 - p_Mplus) * pD_g0.neg +
#     p_Mplus * pD_g0.pos * (1 - q_Dpos)
#   
#   denom_g0 <- (1 - p_Mplus) +
#     p_Mplus * (
#       pD_g0.pos * (1 - q_Dpos) +
#         (1 - pD_g0.pos) * (1 - q_Dneg)
#     )
#   
#   risk_g0_neg_loss <- numer_g0 / denom_g0
#   
#   ## For G = 1, conditional on Mtilde=0
#   numer_g1 <- (1 - p_Mplus) * pD_g1.neg +
#     p_Mplus * pD_g1.pos * (1 - q_Dpos)
#   
#   denom_g1 <- (1 - p_Mplus) +
#     p_Mplus * (
#       pD_g1.pos * (1 - q_Dpos) +
#         (1 - pD_g1.pos) * (1 - q_Dneg)
#     )
#   
#   risk_g1_neg_loss <- numer_g1 / denom_g1
#   
#   RR_neg.loss <- risk_g1_neg_loss / risk_g0_neg_loss
#   
#   
#   ## ---- TRUE RR_pos.loss.corr using your page-9 correction formula ----
#   ## P0(D+|M+) = {1 + (P0(D-)/P0(D+)) *
#   ##               ( P0(Mtilde+|D-) / P1(Mtilde+|M+,D-) ) /
#   ##               ( P0(Mtilde+|D+) / P1(Mtilde+|M+,D+) )
#   ##             }^(-1)
#   ## RR_pos.loss.corr = P1(D+|M+) / P0(D+|M+)_(corr)
#   
#   ## Control-arm marginal outcome probs
#   P0_Dpos <- p_Mplus * pD_g0.pos + (1 - p_Mplus) * pD_g0.neg
#   P0_Dneg <- 1 - P0_Dpos
#   
#   ## Control-arm P(Mtilde+ | D+) and P(Mtilde+ | D-)
#   ## Since Mtilde+ can only occur if M=1:
#   ## P0(Mtilde+|D+) = P0(M=1|D+) * q_Dpos
#   ## P0(Mtilde+|D-) = P0(M=1|D-) * q_Dneg
#   
#   P0_Mpos_given_Dpos <- (p_Mplus * pD_g0.pos) / P0_Dpos
#   P0_Mpos_given_Dneg <- (p_Mplus * (1 - pD_g0.pos)) / P0_Dneg
#   
#   P0_Mt_pos_given_Dpos <- P0_Mpos_given_Dpos * q_Dpos
#   P0_Mt_pos_given_Dneg <- P0_Mpos_given_Dneg * q_Dneg
#   
#   ## Screen-arm retest-positive fractions among true M+ by outcome
#   ## Under your DGP these equal q_Dpos and q_Dneg:
#   P1_Mt_pos_given_Mpos_Dpos <- q_Dpos
#   P1_Mt_pos_given_Mpos_Dneg <- q_Dneg
#   
#   ratio_term <- (P0_Dneg / P0_Dpos) *
#     ( (P0_Mt_pos_given_Dneg / P1_Mt_pos_given_Mpos_Dneg) /
#         (P0_Mt_pos_given_Dpos / P1_Mt_pos_given_Mpos_Dpos) )
#   
#   P0_Dpos_given_Mpos_corr <- 1 / (1 + ratio_term)
#   
#   ## Numerator is the true screen-arm risk among true M+:
#   P1_Dpos_given_Mpos <- pD_g1.pos
#   
#   RR_pos.loss.corr <- P1_Dpos_given_Mpos / P0_Dpos_given_Mpos_corr
#   
#   
#   ## ---- Return as strictly named numeric vector ----
#   out <- c(RR_pos.loss, RR_neg.loss, RR_pos.loss.corr)
#   names(out) <- c("RR_pos.loss", "RR_neg.loss", "RR_pos.loss.corr")
#   
#   return(out)
# }
# 
# 
# 
# true_loss <- true_loss_params(
#   p_Mplus = p_Mplus,
#   pD_g0.neg = pD_g0.neg,
#   pD_g1.neg = pD_g1.neg,
#   pD_g0.pos = pD_g0.pos,
#   pD_g1.pos = pD_g1.pos,
#   p_Mtilde_pos_given_Mpos_Dpos = p_Mtilde_pos_given_Mpos_Dpos,
#   p_Mtilde_pos_given_Mpos_Dneg = p_Mtilde_pos_given_Mpos_Dneg
# )
# 
# 
# ##
# # True parameters corresponding exactly to columns returned by IE_fast()
# ##
# 
# true_par <- c(
#   RR_all = RR,
#   RR_pos = RR_pos,
#   RR_neg = RR_neg,
#   RR_pos.loss = unname(true_loss["RR_pos.loss"]),
#   RR_neg.loss = unname(true_loss["RR_neg.loss"]),
#   RR_pos.loss.corr = unname(true_loss["RR_pos.loss.corr"]),
#   RR_pos.samp = RR_pos,
#   RR_neg.samp = RR_neg
# )
# 
# print(true_par)
# 
# ## 1) Create the truth vector USED FOR COVERAGE (bias-to-gold-standard check)
# true_for_coverage <- true_par
# 
# ## Force loss-based and sampling-based quantities to be judged vs gold-standard targets:
# true_for_coverage["RR_pos.loss"]      <- true_par["RR_pos"]
# true_for_coverage["RR_pos.loss.corr"] <- true_par["RR_pos"]
# true_for_coverage["RR_pos.samp"]      <- true_par["RR_pos"]
# true_for_coverage["RR_neg.loss"]      <- true_par["RR_neg"]
# true_for_coverage["RR_neg.samp"]      <- true_par["RR_neg"]
# 
# stopifnot(all(names(true_for_coverage) == names(true_par)))
# 
# 
# 
# ##
# # Function to simulate one dataset superpop from the same distribution
# # Now includes Mtilde.
# ##
# 
# simulate_superpop_fast <- function(
#     n,
#     p_Mplus,
#     p_Dplus_given_GM,
#     p_Mtilde_pos_given_Mpos_Dpos,
#     p_Mtilde_pos_given_Mpos_Dneg,
#     f1 = f_1,
#     f2 = f_2
# ) {
#   
#   if (n %% 2 != 0) {
#     stop("n must be even because this function assigns exactly n/2 to each G arm.")
#   }
#   
#   if (f1 < 0 || f1 > 1) {
#     stop("f1 must be between 0 and 1.")
#   }
#   
#   if (f2 < 0 || f2 > 1) {
#     stop("f2 must be between 0 and 1.")
#   }
#   
#   ## Generate arm/group G with exactly n/2 in G=0 and n/2 in G=1
#   G <- sample(c(rep(0, n / 2), rep(1, n / 2)), replace = FALSE)
#   
#   ## Generate gold-standard M based on P(M+)
#   M <- rbinom(n, size = 1, prob = p_Mplus)
#   
#   ## Generate D based on P(D+|M,G)
#   p_Dplus <- p_Dplus_given_GM[cbind(as.character(G), as.character(M))]
#   D <- rbinom(n, size = 1, prob = p_Dplus)
#   
#   ## Generate Mtilde based on retention of signal among those with M = 1
#   ##
#   ## Everyone with M = 0 has Mtilde = 0.
#   ## Among those with M = 1:
#   ##   if D = 1, Mtilde ~ Bernoulli(p_Mtilde_pos_given_Mpos_Dpos)
#   ##   if D = 0, Mtilde ~ Bernoulli(p_Mtilde_pos_given_Mpos_Dneg)
#   
#   p_Mtilde_pos <- rep(0, n)
#   
#   p_Mtilde_pos[M == 1 & D == 1] <- p_Mtilde_pos_given_Mpos_Dpos
#   p_Mtilde_pos[M == 1 & D == 0] <- p_Mtilde_pos_given_Mpos_Dneg
#   
#   Mtilde <- rbinom(n, size = 1, prob = p_Mtilde_pos)
#   
#   ## Safety check: Mtilde+ can only occur among M+
#   if (any(Mtilde == 1 & M == 0)) {
#     stop("Error: Mtilde=1 occurred among observations with M=0.")
#   }
#   
#   ## Generate reveal indicator R
#   ##
#   ## G=1 screen arm: M and Mtilde are observed/revealed for everyone.
#   ## G=0 control arm:
#   ##   D=1 stratum: reveal with probability f1
#   ##   D=0 stratum: reveal with probability f2
#   
#   R <- integer(n)
#   
#   R[G == 1] <- 1
#   
#   idx_g0_d1 <- which(G == 0 & D == 1)
#   idx_g0_d0 <- which(G == 0 & D == 0)
#   
#   R[idx_g0_d1] <- rbinom(length(idx_g0_d1), size = 1, prob = f1)
#   R[idx_g0_d0] <- rbinom(length(idx_g0_d0), size = 1, prob = f2)
#   
#   ## Put dataset together
#   superpop <- cbind(
#     G = G,
#     D = D,
#     M = M,
#     Mtilde = Mtilde,
#     R = R
#   )
#   
#   return(superpop)
# }
# 
# 
# ##
# # Function to bootstrap one dataset and summarize IE output
# ##
# 
# bootstrap_IE_summary_fast <- function(superpop, B = 200, f1 = f_1, f2 = f_2) {
#   
#   ## Original-sample parameter estimates
#   IE.orig <- IE_fast(superpop, f1 = f1, f2 = f2)
#   
#   par_names <- names(IE.orig)
#   n_par <- length(IE.orig)
#   
#   ## Store bootstrap estimates in a numeric matrix
#   IE.boot <- matrix(
#     NA_real_,
#     nrow = B,
#     ncol = n_par,
#     dimnames = list(NULL, par_names)
#   )
#   
#   ## Precompute G-specific row indices once for this dataset
#   idx_g0 <- which(superpop[, "G"] == 0)
#   idx_g1 <- which(superpop[, "G"] == 1)
#   
#   n_g0 <- length(idx_g0)
#   n_g1 <- length(idx_g1)
#   
#   for (b in seq_len(B)) {
#     
#     boot_idx_g0 <- sample(idx_g0, size = n_g0, replace = TRUE)
#     boot_idx_g1 <- sample(idx_g1, size = n_g1, replace = TRUE)
#     
#     boot_idx <- c(boot_idx_g0, boot_idx_g1)
#     
#     superpop.boot <- superpop[boot_idx, , drop = FALSE]
#     
#     IE.boot[b, ] <- IE_fast(superpop.boot, f1 = f1, f2 = f2)
#   }
#   
#   ## Bootstrap mean and SD for each parameter
#   boot.mean <- colMeans(IE.boot, na.rm = TRUE)
#   boot.sd <- apply(IE.boot, 2, sd, na.rm = TRUE)
#   
#   ## Normal-approximation 95% CI using bootstrap SD
#   ci.lower <- IE.orig - 1.96 * boot.sd
#   ci.upper <- IE.orig + 1.96 * boot.sd
#   
#   out <- data.frame(
#     parameter = par_names,
#     estimate = as.numeric(IE.orig),
#     boot.mean = as.numeric(boot.mean),
#     boot.sd = as.numeric(boot.sd),
#     ci.lower = as.numeric(ci.lower),
#     ci.upper = as.numeric(ci.upper),
#     row.names = par_names
#   )
#   
#   return(out)
# }
# 
# ##
# # Simulation study settings
# ##
# 
# S <- 10   # number of independently simulated datasets
# B <- 20  # number of bootstrap replicates per dataset
# 
# par_names <- names(true_par)
# n_par <- length(true_par)
# 
# 
# ##
# # Preallocate matrices for simulation results
# ##
# 
# estimate.mat <- matrix(
#   NA_real_,
#   nrow = S,
#   ncol = n_par,
#   dimnames = list(NULL, par_names)
# )
# 
# boot.mean.mat <- matrix(
#   NA_real_,
#   nrow = S,
#   ncol = n_par,
#   dimnames = list(NULL, par_names)
# )
# 
# boot.sd.mat <- matrix(
#   NA_real_,
#   nrow = S,
#   ncol = n_par,
#   dimnames = list(NULL, par_names)
# )
# 
# ci.lower.mat <- matrix(
#   NA_real_,
#   nrow = S,
#   ncol = n_par,
#   dimnames = list(NULL, par_names)
# )
# 
# ci.upper.mat <- matrix(
#   NA_real_,
#   nrow = S,
#   ncol = n_par,
#   dimnames = list(NULL, par_names)
# )
# 
# covered.mat <- matrix(
#   FALSE,
#   nrow = S,
#   ncol = n_par,
#   dimnames = list(NULL, par_names)
# )
# 
# 
# ##
# # Run and time only the simulation loop
# ##
# 
# time.sim <- system.time({
#   
#   for (s in seq_len(S)) {
#     
#     superpop.s <- simulate_superpop_fast(
#       n = n,
#       p_Mplus = p_Mplus,
#       p_Dplus_given_GM = p_Dplus_given_GM,
#       p_Mtilde_pos_given_Mpos_Dpos = p_Mtilde_pos_given_Mpos_Dpos,
#       p_Mtilde_pos_given_Mpos_Dneg = p_Mtilde_pos_given_Mpos_Dneg,
#       f1 = f_1,
#       f2 = f_2
#     )
#     
#     ## Optional sanity checks
#     if (any(superpop.s[, "G"] == 1 & superpop.s[, "R"] != 1)) {
#       stop("Unexpected result: Some G=1 observations have R != 1.")
#     }
#     
#     if (mean(superpop.s[, "Mtilde"] == 1) > mean(superpop.s[, "M"] == 1)) {
#       stop("Unexpected result: P(Mtilde+) > P(M+).")
#     }
#     
#     boot.summary.s <- bootstrap_IE_summary_fast(
#       superpop = superpop.s,
#       B = B,
#       f1 = f_1,
#       f2 = f_2
#     )
#     
#     ## Ensure returned parameter order matches preallocated matrices
#     boot.summary.s <- boot.summary.s[par_names, , drop = FALSE]
#     
#     estimate.mat[s, ]  <- boot.summary.s$estimate
#     boot.mean.mat[s, ] <- boot.summary.s$boot.mean
#     boot.sd.mat[s, ]   <- boot.summary.s$boot.sd
#     ci.lower.mat[s, ]  <- boot.summary.s$ci.lower
#     ci.upper.mat[s, ]  <- boot.summary.s$ci.upper
#     
#     covered.mat[s, ] <- ci.lower.mat[s, ] <= true_for_coverage &
#       true_for_coverage <= ci.upper.mat[s, ]
#     
#     if (s %% 10 == 0) {
#       cat("Finished simulated dataset", s, "of", S, "\n")
#     }
#   }
# })
# 
# system2("say", "Finished!", wait = FALSE)
# 
# 
# IE.sim.summary <- data.frame(
#   parameter = par_names,
#   true = as.numeric(true_par),
#   true.for.coverage = as.numeric(true_for_coverage),
#   
#   mean.estimate = colMeans(estimate.mat, na.rm = TRUE),
#   sd.estimate = apply(estimate.mat, 2, sd, na.rm = TRUE),
#   
#   mean.boot.mean = colMeans(boot.mean.mat, na.rm = TRUE),
#   mean.boot.sd = colMeans(boot.sd.mat, na.rm = TRUE),
#   
#   mean.ci.lower = colMeans(ci.lower.mat, na.rm = TRUE),
#   mean.ci.upper = colMeans(ci.upper.mat, na.rm = TRUE),
#   
#   coverage = colMeans(covered.mat, na.rm = TRUE),
#   n.sims = colSums(!is.na(covered.mat)),
#   
#   row.names = NULL
# )
# 
# ##
# # Optional: create long-format simulation results
# # This is not needed for IE.sim.summary, but can be useful for diagnostics.
# ##
# 
# sim.results <- do.call(
#   rbind,
#   lapply(seq_len(S), function(s) {
#     data.frame(
#       sim = s,
#       parameter = par_names,
#       true = as.numeric(true_par),
#       estimate = estimate.mat[s, ],
#       boot.mean = boot.mean.mat[s, ],
#       boot.sd = boot.sd.mat[s, ],
#       ci.lower = ci.lower.mat[s, ],
#       ci.upper = ci.upper.mat[s, ],
#       covered = covered.mat[s, ],
#       row.names = NULL
#     )
#   })
# )
# 
# 
# 
# ##
# # Print outputs
# ##
# 
# time.sim
# 
# true_par
# 
# IE.sim.summary
# 
# head(sim.results)
# 
# 
# 
# 
# 
# 
# 
# ## ------------------------------------------------------------
# ## Input tables from the PowerPoint slide
# ## ------------------------------------------------------------
# 
# observed_never_positive <- matrix(
#   c(
#     200,   175,
#     37800, 33075
#   ),
#   nrow = 2,
#   byrow = TRUE,
#   dimnames = list(
#     D = c("D+", "D-"),
#     arm = c("screen", "control")
#   )
# )
# 
# observed_ever_positive <- matrix(
#   c(
#     520,  525,
#     1480, 1225
#   ),
#   nrow = 2,
#   byrow = TRUE,
#   dimnames = list(
#     D = c("D+", "D-"),
#     arm = c("screen", "control")
#   )
# )
# 
# unknown_positive <- matrix(
#   c(
#     180,   300,
#     9820, 14700
#   ),
#   nrow = 2,
#   byrow = TRUE,
#   dimnames = list(
#     D = c("D+", "D-"),
#     arm = c("screen", "control")
#   )
# )
# 
# 
# ## ------------------------------------------------------------
# ## Corrected function
# ## ------------------------------------------------------------
# 
# correct_positive_tables_simple <- function(obs_never,
#                                            obs_ever,
#                                            unknown_pos,
#                                            screen_col = "screen",
#                                            control_col = "control",
#                                            round_output = TRUE) {
#   
#   obs_never   <- as.matrix(obs_never)
#   obs_ever    <- as.matrix(obs_ever)
#   unknown_pos <- as.matrix(unknown_pos)
#   
#   if (!identical(dim(obs_never), dim(obs_ever))) {
#     stop("obs_never and obs_ever must have the same dimensions.")
#   }
#   
#   if (!identical(dim(obs_never), dim(unknown_pos))) {
#     stop("obs_never and unknown_pos must have the same dimensions.")
#   }
#   
#   if (is.null(colnames(obs_never)) ||
#       is.null(colnames(obs_ever)) ||
#       is.null(colnames(unknown_pos))) {
#     stop("All three input tables must have column names.")
#   }
#   
#   required_cols <- c(screen_col, control_col)
#   
#   if (!all(required_cols %in% colnames(obs_never))) {
#     stop("obs_never must contain both screen_col and control_col.")
#   }
#   
#   if (!all(required_cols %in% colnames(obs_ever))) {
#     stop("obs_ever must contain both screen_col and control_col.")
#   }
#   
#   if (!all(required_cols %in% colnames(unknown_pos))) {
#     stop("unknown_pos must contain both screen_col and control_col.")
#   }
#   
#   ## Standard trial analysis table, reconstructed cellwise
#   ## from observed never-positive + observed ever-positive + unknown positive.
#   standard_table <- obs_never + obs_ever + unknown_pos
#   
#   ## Row-specific non-compliance fractions.
#   noncomp_screen <- unknown_pos[, screen_col] /
#     standard_table[, screen_col]
#   
#   noncomp_control <- unknown_pos[, control_col] /
#     standard_table[, control_col]
#   
#   ## Row-specific compliance fractions.
#   compliance_screen <- 1 - noncomp_screen
#   compliance_control <- 1 - noncomp_control
#   
#   ## Row-specific correction factor for the control column.
#   correction_factor <- compliance_screen / compliance_control
#   
#   ## Initialize corrected tables as copies of the observed tables.
#   corrected_never <- obs_never
#   corrected_ever  <- obs_ever
#   
#   ## Screen columns are unchanged.
#   corrected_never[, screen_col] <- obs_never[, screen_col]
#   corrected_ever[,  screen_col] <- obs_ever[,  screen_col]
#   
#   ## Control columns are corrected using the row-specific factor.
#   corrected_never[, control_col] <-
#     obs_never[, control_col] * correction_factor
#   
#   corrected_ever[, control_col] <-
#     obs_ever[, control_col] * correction_factor
#   
#   if (round_output) {
#     corrected_never <- round(corrected_never)
#     corrected_ever  <- round(corrected_ever)
#   }
#   
#   list(
#     standard_trial_analysis_table = standard_table,
#     
#     noncomp_screen = noncomp_screen,
#     noncomp_control = noncomp_control,
#     
#     compliance_screen = compliance_screen,
#     compliance_control = compliance_control,
#     
#     correction_factor = correction_factor,
#     
#     corrected_never_positive_table = corrected_never,
#     corrected_ever_positive_table = corrected_ever
#   )
# }
# 
# 
# ## ------------------------------------------------------------
# ## Run the function
# ## ------------------------------------------------------------
# 
#   res <- correct_positive_tables_simple(
#   obs_never   = observed_never_positive,
#   obs_ever    = observed_ever_positive,
#   unknown_pos = unknown_positive
# )
# 
# res$standard_trial_analysis_table
# res$correction_factor
# res$corrected_never_positive_table
# res$corrected_ever_positive_table
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 

# ###
# # OLD STUFF BELOW
# ###
# 
# ##
# # Function to conduct the simulation of the bootstrap variances for the IE analysis
# #
# # Sampling of which blood samples in control arm to retest:
# # frac.g0.plus: fraction of control arm (g0) D+ whose blood samples are tested
# # frac.g0.minus:fraction of control arm (g0) D- whose blood samples are tested
# #
# # Fractions of stored control-arm samples with retention of signal (complement of loss of signal)
# # frac.g0.plus.pos: fraction of control-arm (g0) stored samples that retain signal among D+
# # frac.g0.minus.pos:fraction of control-arm (g0) stored samples that retain signal among D-
# ##
# conceal.reveal <- function(nsim=1000, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02,
#                            RR=0.9, RR_neg=1, RR_pos=0.86667,
#                            frac.g0.plus=1, frac.g0.minus=1,
#                            frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,
#                            censor.screen.Dplus = 0, censor.control.Dplus = 0,
#                            censor.screen.Dminus = 0, censor.control.Dminus = 0,
#                            move.screen.Mplus.Dplus = 0, move.screen.Mminus.Dplus = 0,
#                            bootsims=1, sim=FALSE)
# {
#   # RR_neg = P(D+|M-,G=1) / P(D+|M-,G=0)
#   if (RR>RR_neg) stop("unrealistic for RR>RR_neg")
#   # RR_pos = P(D+|M+,G=1) / P(D+|M+,G=0)
#   if (RR<RR_pos) stop("unrealistic for RR<RR_pos")
#   
#   # Fix P(D+|G=0) as population disease mortality and expected mortality RR=P(D+|G=1)/P(D+|G=0)
#   # Calculate P(D+|G=1) = RR*P(D+|G=0), and true RD=P(D+|G=0)-P(D+|G=1)
#   pD_g1 <- RR*pD_g0
#   RD <- pD_g0 - pD_g1
#   
#   # P(D+|M-,G=0) = (RR_pos*P(D+|G=0) - P(D+|G=1)) / [P(M-)*(RR_pos-RR_neg)]
#   pD_g0.neg <- ifelse(RR_pos-RR_neg==0,
#                       pD_g0, # under the null, P(D+|M-,G=0)=P(D+|G=0)
#                       (RR_pos*pD_g0 - pD_g1) / ( (1-p_Mplus)*(RR_pos-RR_neg) ) )
#   
#   # P(D+|M-,G=1) = RR_neg * P(D+|M-,G=0)
#   pD_g1.neg <- RR_neg * pD_g0.neg
#   
#   # P(D+|M+,G=0) = (RR_neg*P(D+|G=0) - P(D+|G=1)) / [P(M+)(RR_neg-RR_pos)]
#   pD_g0.pos <- ifelse(RR_neg-RR_pos==0,
#                       pD_g0, # under the null, P(D+|M-,G=0)=P(D+|G=0)
#                       (RR_neg*pD_g0 - pD_g1) / ( p_Mplus*(RR_neg-RR_pos) ) )
#   
#   if (pD_g0.pos>1) stop(paste("pD_g0.pos=",pD_g0.pos,"is >1 because (RR_neg*pD_g0 - pD_g1) > p_Mplus*(RR_neg-RR_pos)"))
#   
#   # P(D+|M+,G=1) = RR_pos * P(D+|M+,G=0)
#   pD_g1.pos <- RR_pos * pD_g0.pos
#   
#   # P(D+|M-,G=1) = RR_neg * P(D+|M-,G=0) 
#   pD_g1.neg <- RR_neg * pD_g0.neg
#   
#   # RD_neg = P(D+|M-,G=0) - P(D+|M-,G=1)
#   RD_neg <- pD_g0.neg - pD_g1.neg
#   
#   # Given RD and RD_neg=P(D+|M-,G=0)-P(D+|M-,G=1), solve for RD_pos=P(D+|M+,G=0)-P(D+|M+,G=1)
#   # RD_pos <- (RD-RD_neg*(1-p_Mplus)) / p_Mplus
#   
#   # RD_pos = P(D+|M+,G=0) - P(D+|M+,G=1)
#   RD_pos <- pD_g0.pos - pD_g1.pos
#   
#   
#   # Deal with constraints
#   # stop if P(D+|M-,G) < 0, that is, we hit the boundary of no outcomes among never-positives
#   if (pD_g0.neg < 0 | pD_g1.neg < 0) stop("P(D+|M-,G) <0, meaning negative outcomes in never-positives. Check RR and RRpos.")
#   # stop if P(D+|M+,G) < 0, that is, we hit the boundary of no outcomes among ever-positives
#   if (pD_g0.pos < 0 | pD_g1.pos < 0) stop("P(D+|M+,G) <0, meaning negative outcomes in ever-positives. Check RR and RRpos.")
#   
#   # Specify PPV=P(D+|M+,G=0)
#   # Check for compatibilty with marginals by checking if sens=P(M+|D+,G=0)=PPV*P(M+)/P(D+|G=0) > 1
#   # I use 1.001 because sometimes there are numerical issues and sens is slightly > 1
#   # pD_g0.pos <- 0.3
#   if (pD_g0.pos*p_Mplus/pD_g0>1.001) stop("incompatible PPV conditional because sens>1")
#   
#   
#   # if (pD_g0.pos*p_Mplus/pD_g0<RD_pos) stop("PPV too small given large RD_pos")
#   
#   # # Now calculate the rest of the conditionals
#   # # P(D+|M+,G=1) = P(D+|M+,G=0) - RD_pos
#   # pD_g1.pos <- pD_g0.pos - RD_pos
#   # paste("conceal-reveal RR:", RR_pos <- pD_g1.pos/pD_g0.pos)
#   # # back out P(D+|M-,G=0) from P(D+|G=0) = P(D+|M-,G=0)P(M-) + P(D+|M+,G=0)P(M+)
#   # pD_g0.neg <- (pD_g0 - pD_g0.pos*p_Mplus) / (1-p_Mplus)
#   # # P(D+|M-,G=1) = P(D+|M-,G=0) - RD_pos
#   # pD_g1.neg <- pD_g0.neg - RD_neg
#   # RR_neg <- pD_g1.neg/pD_g0.neg # other conceal-reveal RR
#   
#   # Theoretical ratio of Z-statistics: RD_pos/RD * P(M+) * sqrt{P(M+)/(P(M+|D+)P(M+|D-)}, where:
#   # P(M+|D+) = P(D+|M+)P(M+)/P(D+)
#   # P(M+|D-) = P(D-|M+)(M+)/P(D-)
#   # P(D+|M+) = P(D+|M+,G=0)P(G=0) + P(D+|M+,G=1)P(G=1)
#   pD_pos <- pD_g0.pos*p_g0 + pD_g1.pos*(1-p_g0)
#   
#   # Make sure this expression agrees with the above, valid only for P(G=1)=0.5: equal sample size both arms
#   # P(D+|M+) = {0.5*(1+RR_pos)*P(D+|G=0)/P(M+)} * (RR_neg-RR)/(RR_neg-RR_pos)
#   pD_pos_check <- (0.5*(1+RR_pos)*pD_g0/p_Mplus) * (RR_neg-RR)/(RR_neg-RR_pos)
#   
#   # P(D+|M-) = P(D+|M-,G=0)P(G=0) + P(D+|M-,G=1)P(G=1)
#   pD_neg <- pD_g0.neg*p_g0 + pD_g1.neg*(1-p_g0)
#   # P(D+) = P(D+|M+)P(M+) + P(D+|M-)P(M-)
#   pD <- pD_pos*p_Mplus + pD_neg*(1-p_Mplus)
#   # P(M+|D+) = P(D+|M+)P(M+)/P(D+)
#   pMplus_Dplus <- pD_pos*p_Mplus/pD
#   # P(M+|D-) = P(D-|M+)P(M+)/P(D-)
#   pMplus_Dminus <- (1-pD_pos)*p_Mplus/(1-pD)
#   
#   # theoretical Z-ratio
#   # Z_ratio <- RD_pos/RD * p_Mplus * sqrt(p_Mplus/(pMplus_Dminus*pMplus_Dplus))
#   Z_ratio <- (1 - (RD_neg/RD)*(1-p_Mplus)) * sqrt(p_Mplus/(pMplus_Dminus*pMplus_Dplus))
#   
#   
#   ###
#   # Print out expected value of 2x2 tables
#   ###
#   # M+ table
#   # P(D+,G=1|M+) * #(M+)
#   a1 <- (pD_g1.pos*(1-p_g0) * n*p_Mplus) * (1-censor.screen.Dplus) 
#   move.a1.to.b1 <- (pD_g1.pos*(1-p_g0) * n*p_Mplus) * (censor.screen.Dplus)
#   # P(D+,G=0|M+) * #(M+)
#   a0 <- (pD_g0.pos*p_g0 * n*p_Mplus) * (1-censor.control.Dplus)
#   move.a0.to.b0 <- (pD_g0.pos*p_g0 * n*p_Mplus) * (censor.control.Dplus)
#   # P(D-,G=1|M+) * #(M+)
#   c1 <- ((1-pD_g1.pos)*(1-p_g0) * n*p_Mplus) * (1-censor.screen.Dminus) 
#   move.c1.to.d1 <- ((1-pD_g1.pos)*(1-p_g0) * n*p_Mplus) * (censor.screen.Dminus)
#   # P(D-,G=0|M+) * #(M+)
#   c0 <- ((1-pD_g0.pos)*p_g0 * n*p_Mplus) * (1-censor.control.Dminus)
#   move.c0.to.d0 <- ((1-pD_g0.pos)*p_g0 * n*p_Mplus) * (censor.control.Dminus)
#   
#   # M- table
#   # P(D+,G=1,M-) * #(M-)
#   b1 <- pD_g1.neg*(1-p_g0) * n*(1-p_Mplus) + move.a1.to.b1
#   # P(D+,G=0,M-) * #(M-)
#   b0 <- pD_g0.neg*p_g0 * n*(1-p_Mplus) + move.a0.to.b0
#   # P(D-,G=1,M-) * #(M-)
#   d1 <- (1-pD_g1.neg)*(1-p_g0) * n*(1-p_Mplus) + move.c1.to.d1
#   # P(D-,G=0,M-) * #(M-)
#   d0 <- (1-pD_g0.neg)*p_g0 * n*(1-p_Mplus) + move.c0.to.d0
#   
#   # Now account for those flipping from D+ to D- in screen arm due to non-adherence
#   # Do among ever-positives
#   move.c1.to.a1 <- c1 * move.screen.Mplus.Dplus
#   a1 <- a1 + move.c1.to.a1
#   c1 <- c1 - move.c1.to.a1
#   # Do among never-positives
#   move.d1.to.b1 <- d1 * move.screen.Mminus.Dplus
#   b1 <- b1 + move.d1.to.b1
#   d1 <- d1 - move.d1.to.b1
#   
#   # compute margins and final tables
#   n1pos <- a1 + c1 
#   n0pos <- a0 + c0 
#   Mpos <- matrix( c(a1,a0,c1,c0),nrow=2,byrow=TRUE,
#                   dimnames=list(c("D+","D-"),c("screen","control")))
#   n1neg <- b1 + d1 #n*(1-p_g0)*(1-p_Mplus)
#   n0neg <- b0 + d0 #n*p_g0*(1-p_Mplus)
#   # d1 <- n1neg - b1
#   # d0 <- n0neg - b0
#   Mneg <- matrix( c(b1,b0,d1,d0),nrow=2,byrow=TRUE,
#                   dimnames=list(c("D+","D-"),c("screen","control")))
#   
#   # Now include censoring of ever-test-positivity in the test-positive table
#   # These people will be misclassified into the test-negative table
#   print("Ever-positive table")
#   print(Mpos) ; Mpos.test <- prop.test(t(Mpos), correct=F)
#   print(paste("RR_pos:", Mpos.test$estimate[1]/Mpos.test$estimate[2],"RD_pos:", diff(Mpos.test$estimate), 
#               "p-value:",Mpos.test$p.val,
#               "per-arm N for 90% power:",round((n/2)*((qnorm(1-0.9)-1.96)/qnorm(Mpos.test$p.val/2))^2)))
#   
#   print("Never-positive table")
#   print(Mneg); Mneg.test <- prop.test(t(Mneg), correct=F)
#   print(paste("RR_neg:", Mneg.test$estimate[1]/Mneg.test$estimate[2],"RD_neg:", diff(Mneg.test$estimate), 
#               "p-value:",Mneg.test$p.val))
#   
#   print("Total table")
#   print(all<-Mneg+Mpos); all.test <- prop.test(t(all), correct=F)
#   print(paste("RR_all:", all.test$estimate[1]/all.test$estimate[2],"RD_all:", diff(all.test$estimate), 
#               "p-value:",all.test$p.val,
#               "per-arm N for 90% power:",round((n/2)*((qnorm(1-0.9)-1.96)/qnorm(all.test$p.val/2))^2)))
#   
#   # Do Peter Sasieni's "targeted IE"
#   targeted <- matrix( c(a1,a0,d1+c1+b1,d0+c0+b0),nrow=2,byrow=TRUE,
#                       dimnames=list(c("D+ and M+","otherwise"),c("screen","control")))
#   print(targeted) ; targeted.test <- prop.test(t(targeted), correct=F)
#   print(paste("Targeted RR_pos:", targeted.test$estimate[1]/targeted.test$estimate[2],
#               "RD_pos:", diff(targeted.test$estimate), 
#               "p-value:",targeted.test$p.val,
#               "per-arm N for 90% power:",round((n/2)*((qnorm(1-0.9)-1.96)/qnorm(targeted.test$p.val/2))^2)))
#   
#   # # Check if RRpos = RR * P(M+|D+,G=1)/P(M+|D+,G=0)
#   # pMplus_Dplus_g1 <- a1/(a1+b1)
#   # pMplus_Dplus_g0 <- a0/(a0+b0)
#   # RR_pos.check <- RR * pMplus_Dplus_g1/pMplus_Dplus_g0
#   # print(paste("RR_pos.check=",RR_pos.check,"RR_pos=",RR_pos))
#   # 
#   # # Check if RRneg = RR * P(M-|D+,G=1)/P(M-|D+,G=0)
#   # pMminus_Dplus_g1 <- 1-pMplus_Dplus_g1
#   # pMminus_Dplus_g0 <- 1-pMplus_Dplus_g0
#   # RR_neg.check <- RR * pMminus_Dplus_g1/pMminus_Dplus_g0
#   # print(paste("RR_neg.check=",RR_neg.check,"RR_neg=",RR_neg))
#   # # check if RRneg = (P(D+│G=0)⋅P(M-│D+,G=1)⋅RR)/(P(M-│G=1)-P(M-│D-,G=0)⋅P(D-|G=0))
#   # pMminus_Dminus_g0 <- d0/(d0+c0)
#   # print(paste("specificity in control arm=",pMminus_Dminus_g0))
#   # RR_neg.check1 <- (pD_g0*pMminus_Dplus_g1*RR) / ((1-p_Mplus)-pMminus_Dminus_g0*(1-pD_g0))
#   # print(paste("RR_neg.check1=",RR_neg.check1,"RR_neg=",RR_neg))
#   # pMminus_Dminus_g1 <- d1/(d1+c1)
#   # print(paste("specificity in screen arm=",pMminus_Dminus_g1))
#   # RR_neg.spec <- (pD_g0*pMminus_Dplus_g1*RR) / ((1-p_Mplus)-pMminus_Dminus_g1*(1-pD_g0))
#   # print(paste("RR_neg.spec=",RR_neg.spec,"RR_neg=",RR_neg))
#   # 
#   # # check if RRpos = (P(D+│G=0)⋅P(M+│D+,G=1)⋅RR)/(P(M+│G=1)-P(M+│D-,G=0)⋅P(D-|G=0))
#   # RR_pos.check1 <- (pD_g0*(1-pMminus_Dplus_g1)*RR) / (p_Mplus-(1-pMminus_Dminus_g0)*(1-pD_g0))
#   # print(paste("RR_pos.check1=",RR_pos.check1,"RR_pos=",RR_pos))
#   # # Inputting 1-spec from the screen arm (1-pMminus_Dminus_g1), is the approximation good?
#   # RR_pos.spec <- (pD_g0*(1-pMminus_Dplus_g1)*RR) / (p_Mplus-(1-pMminus_Dminus_g1)*(1-pD_g0))
#   # print(paste("RR_pos.spec=",RR_pos.spec,"RR_pos=",RR_pos))
#   
#   ##
#   # Output power for RDpos given fixed power for RD
#   # This is pnorm(Zratio * mean_Z - 1.96)
#   # mean Z is that from Goodman 1992 for replication probability
#   # i.e. power=0.5 is mean Z stat = 1.96, power=0.8 is mean Z stat = 2.802, power=0.9 is mean Z stat = 3.24
#   ##
#   if (sim==FALSE) { 
#     # RD_power <- 0.9
#     # return(pnorm(Z_ratio*(qnorm(1-0.05/2)+qnorm(RD_power))-1.96)) 
#     # return(Z_ratio)
#     tables <- NULL
#     tables$all <- all ; tables$Mpos <- Mpos ; tables$Mneg <- Mneg 
#     return(tables)
#   }
#   # return(data.frame(
#   #   RR=RR, RR_pos=RR_pos, RR_neg=RR_neg,
#   #   pD_g1.pos=pD_g1.pos,
#   #   pD_g0.pos=pD_g0.pos,
#   #   pD_g1.neg=pD_g1.neg,
#   #   pD_g0.neg=pD_g0.neg,
#   #   pMplus_Dminus=pMplus_Dminus,
#   #   pMplus_Dplus=pMplus_Dplus,
#   #   pD_pos=pD_pos,
#   #   pD_neg=pD_neg,
#   #   pD=pD,
#   #   RD=RD, RD_neg=RD_neg, RD_pos=RD_pos, Z_ratio=Z_ratio, 
#   #   power.pos=pnorm(Z_ratio*2.802-1.96)
#   # ))
#   
#   ##
#   # initialize simulation
#   ##
#   
#   # Simulate 2x2 tables
#   set.seed(406553735)
#   for(i in 1:nsim) {
#     
#     ##
#     # Order of data simulating:
#     # 1.  1st simulate arms assuming all assumptions hold, create tables for ever-positives and all-negatives
#     # 2.  2nd simulate non-compliance to blood collections, create table for the unknowns
#     # 3.  3rd simulate loss-of-signal in the control-arm ever-positives, reapportion them to negatives
#     # 4.  4th simulate subsampling, based on outcome, for control-arm blood testing
#     
#     
#     ##
#     # 1st: Generate M+ and M- sample sizes, fixing total sample sizes n*p_g0 and n*(1-p_g0) in each group
#     ##
#     
#     # #(M+,G=1) = n_1+
#     g1.pos <- rbinom(1, n*(1-p_g0), p_Mplus) 
#     
#     # #(M+,G=0) = n_0+
#     g0.pos <- rbinom(1, n*p_g0, p_Mplus) 
#     
#     # #(M-,G=1) = n_1-
#     g1.neg <- n*(1-p_g0)-g1.pos
#     
#     # #(M-,G=0) = n_0-
#     g0.neg <- n*p_g0-g0.pos
#     
#     ##
#     # generate the M- and M+ tables
#     ##
#     
#     # #(D+,G=1,M-) = b_1
#     Dplus_g1.neg <- rbinom(1, g1.neg, pD_g1.neg)
#     # #(D+,G=0,M-) = b_0
#     Dplus_g0.neg <- rbinom(1, g0.neg, pD_g0.neg)
#     # #(D+,G=1,M+) = a_1
#     Dplus_g1.pos <- rbinom(1, g1.pos, pD_g1.pos)
#     # #(D+,G=0,M+) = a_0
#     Dplus_g0.pos <- rbinom(1, g0.pos, pD_g0.pos)
#     
#     # Calculate the D- cells of the tables
#     Dminus_g1.pos <- g1.pos - Dplus_g1.pos # need to calculate true cell c_1
#     Dminus_g0.pos <- g0.pos - Dplus_g0.pos # need to calculate true cell c_0
#     Dminus_g1.neg <- g1.neg - Dplus_g1.neg # need to calculate true cell d_1
#     Dminus_g0.neg <- g0.neg - Dplus_g0.neg # need to calculate true cell d_0
#     
#     
#     
#     ##
#     # 2nd: Simulate non-adherence to blood draws in both group: dropout
#     # Now simulate dropout among the ever-positives in both screen and control arms
#     # They get redistributed to the never-positives in both screen and control arms
#     ##
#     
#     # First non-adherence among ever-positives in screen (G=1) and control (G=0) arms
#     # Calculate the non-censoring fractions, separately for D+ and D-, in each arm (although these are unobservable)
#     # The "move" variables are those who are unknown positivity but truly ever-positive - unobservable
#     # #(D+,G=1,M+) = a_1
#     move.a1.to.b1 <- rbinom(1, Dplus_g1.pos, censor.screen.Dplus)
#     Dplus_g1.pos <-  Dplus_g1.pos - move.a1.to.b1
#     # #(D+,G=0,M+) = a_0
#     move.a0.to.b0 <- rbinom(1, Dplus_g0.pos, censor.control.Dplus)
#     Dplus_g0.pos <- Dplus_g0.pos - move.a0.to.b0
#     # #(D-,G=1,M+) = c_1
#     move.c1.to.d1 <- rbinom(1, Dminus_g1.pos, censor.screen.Dminus)
#     Dminus_g1.pos <- Dminus_g1.pos - move.c1.to.d1
#     # #(D-,G=0,M+) = c_0
#     move.c0.to.d0 <- rbinom(1, Dminus_g0.pos, censor.control.Dminus)
#     Dminus_g0.pos <- Dminus_g0.pos - move.c0.to.d0
#     
#     # Second: non-adherence among never-positives in screen (G=1) and control (G=0) arms
#     # rename the "move" variables from above as ".unknown" to clarify ("move" is only if unknown is treated as negative)
#     # ".neg.unknown" is truly never-positive but observed as unknown (these are unobservable)
#     # #(D+,G=1,M-) = b_1
#     Dplus_g1.neg.unknown <- rbinom(1, Dplus_g1.neg, censor.screen.Dplus) 
#     Dplus_g1.neg <-  Dplus_g1.neg - Dplus_g1.neg.unknown
#     # #(D+,G=0,M-) = b_0
#     Dplus_g0.neg.unknown <- rbinom(1, Dplus_g0.neg, censor.control.Dplus)
#     Dplus_g0.neg <- Dplus_g0.neg - Dplus_g0.neg.unknown
#     # #(D-,G=1,M-) = d_1
#     Dminus_g1.neg.unknown <- rbinom(1, Dminus_g1.neg, censor.screen.Dminus)
#     Dminus_g1.neg <- Dminus_g1.neg - Dminus_g1.neg.unknown
#     # #(D-,G=0,M-) = d_0
#     Dminus_g0.neg.unknown <- rbinom(1, Dminus_g0.neg, censor.control.Dminus)
#     Dminus_g0.neg <- Dminus_g0.neg - Dminus_g0.neg.unknown
#     
#     # Third: Calculate the 2x2 table for group with unknown positivity ".unknown"
#     # Key point: each term in the sums below are unobservable, but the sum is observable
#     # First term is the 
#     Dplus_g1.unknown <- move.a1.to.b1 + Dplus_g1.neg.unknown 
#     Dplus_g0.unknown <- move.a0.to.b0 + Dplus_g0.neg.unknown
#     Dminus_g1.unknown <- move.c1.to.d1 + Dminus_g1.neg.unknown
#     Dminus_g0.unknown <- move.c0.to.d0 + Dminus_g0.neg.unknown
#     
#     # recalculate margins, where the unknown are separated from the never-positive
#     g1.pos <- Dplus_g1.pos + Dminus_g1.pos 
#     g0.pos <- Dplus_g0.pos + Dminus_g0.pos 
#     g1.neg <- Dplus_g1.neg + Dminus_g1.neg 
#     g0.neg <- Dplus_g0.neg + Dminus_g0.neg 
#     g1.unknown <- Dplus_g1.unknown + Dminus_g1.unknown
#     g0.unknown <- Dplus_g0.unknown + Dminus_g0.unknown
#     
#     # Define the 4 cells for the all table and its margins
#     Dplus_g1.all <- Dplus_g1.neg+Dplus_g1.pos+Dplus_g1.unknown
#     Dplus_g0.all <- Dplus_g0.neg+Dplus_g0.pos+Dplus_g0.unknown
#     Dminus_g1.all <- Dminus_g1.neg+Dminus_g1.pos+Dminus_g1.unknown
#     Dminus_g0.all <- Dminus_g0.neg+Dminus_g0.pos+Dminus_g0.unknown
#     g1.all <- g1.neg+g1.pos+g1.unknown
#     g0.all <- g0.neg+g0.pos+g0.unknown
#     
#     
#     ##
#     # 3rd: simulate loss-of-signal, based on outcome, in control-arm positives
#     # reapportion them into the negative table
#     # leave the unknown table unchanged
#     # ALSO: simulate the validation study in screen-arm to calculate retest-positive fractions
#     ##
#     
#     # Subsample among a_0=#(D+,M+,G=0): the fraction of these that test true positive M+
#     Dplus_g0.pos.loss <- rbinom(1, Dplus_g0.pos, frac.g0.plus.pos)
#     
#     # The a_0=#(D+,M+,G=0) that test false-negative, put them into b_0=#(D+,M-,G=0)
#     Dplus_g0.neg.loss <- Dplus_g0.neg + (Dplus_g0.pos - Dplus_g0.pos.loss)
#     
#     # Subsample among c_0=#(D-,M+,G=0): the fraction of these that test true positive M+
#     Dminus_g0.pos.loss <- rbinom(1, Dminus_g0.pos ,frac.g0.minus.pos)
#     
#     # The c_0=#(D-,M+,G=0) that test false-negative, put them into d_0=#(D-,M-,G=0)
#     Dminus_g0.neg.loss <- Dminus_g0.neg + (Dminus_g0.pos - Dminus_g0.pos.loss)
#     
#     # Simulate validation study in screened arm (G=1) that retests all M+ 
#     # Retest in G=1 all M+, and loss of signal depends on D+ vs. D- 
#     Dplus_g1.pos.retest <- rbinom(1, Dplus_g1.pos, frac.g0.plus.pos)
#     Dminus_g1.pos.retest <- rbinom(1, g1.pos - Dplus_g1.pos, frac.g0.minus.pos)
#     
#     ##
#     # 4th: simulate subsampling of controls for MCED testing, NOT incorporating loss-of-signal
#     #
#     # Now we only test a subsamples of G=0 for the marker M
#     # frac.g0.plus:  The fraction of #G=0,D+ that we will observe
#     # frac.g0.minus: The fraction of #G=0,D- that we will observe
#     # We will estimate these fractions and weight by them
#     # Note: we cannot estimate 4 sampling fractions (the above 2 by M+ and M-)
#     #       because we cannot know the total M+ or M- in the G=0 arm
#     ##
#     
#     # subsample among a_0=#(D+,M+,G=0)
#     # Dplus_g0.pos.subsamp <- rbinom(1,Dplus_g0.pos.loss,frac.g0.plus)
#     Dplus_g0.pos.subsamp <- rbinom(1,Dplus_g0.pos, frac.g0.plus)
#     
#     # subsample among c_0=#(D-,M+,G=0)
#     # Dminus_g0.pos.subsamp <- rbinom(1, g0.pos - Dplus_g0.pos.loss,frac.g0.minus)
#     Dminus_g0.pos.subsamp <- rbinom(1, g0.pos - Dplus_g0.pos, frac.g0.minus)
#     
#     # subsample among b_0=#(D+,M-,G=0)
#     # Dplus_g0.neg.subsamp <- rbinom(1,Dplus_g0.neg.loss, frac.g0.plus)
#     Dplus_g0.neg.subsamp <- rbinom(1, Dplus_g0.neg, frac.g0.plus)
#     
#     # subsample among d_0=#(D-,M-,G=0)
#     # Dminus_g0.neg.subsamp <- rbinom(1, g0.neg - Dplus_g0.neg.loss, frac.g0.minus)
#     Dminus_g0.neg.subsamp <- rbinom(1, g0.neg - Dplus_g0.neg, frac.g0.minus)
#     
#     # save parameter estimate for the simulation
#     sim_params <- ie_analysis(Dplus_g1.neg, Dplus_g0.neg, g1.neg, g0.neg, Dminus_g1.neg, Dminus_g0.neg,
#                               Dplus_g1.pos, Dplus_g0.pos, g1.pos, g0.pos, Dminus_g1.pos, Dminus_g0.pos,
#                               Dplus_g1.unknown, Dplus_g0.unknown, g1.unknown, g0.unknown, Dminus_g1.unknown, Dminus_g0.unknown,
#                               Dplus_g1.all, Dplus_g0.all, g1.all, g0.all, Dminus_g1.all, Dminus_g0.all,
#                               Dplus_g0.pos.loss, Dminus_g0.pos.loss, Dplus_g0.neg.loss, Dminus_g0.neg.loss,
#                               Dplus_g0.pos.subsamp, Dminus_g0.pos.subsamp, Dplus_g0.neg.subsamp, Dminus_g0.neg.subsamp,
#                               Dplus_g1.pos.retest, Dminus_g1.pos.retest,
#                               RR, RR_pos, RR_neg, RD_neg, RD_pos, RD, Z_ratio,
#                               frac.g0.plus, frac.g0.minus, frac.g0.plus.pos, frac.g0.minus.pos,
#                               pD_g1.pos, pD_g0.pos, pD_g1.neg, pD_g0.neg, pMplus_Dminus, pMplus_Dplus, pD_pos, pD_pos_check, pD_neg, pD)
#     
#     # Now create bootstrap samples and then analyze them
#     for (j in 1:bootsims) {
#       
#       # Bootstrap interior of tables, keep control (g0) and screen (g1) sample sizes the same
#       Dplus_g1.neg_sim <- rbinom(1, g1.neg, Dplus_g1.neg/g1.neg)
#       Dplus_g0.neg_sim <- rbinom(1, g0.neg, Dplus_g0.neg/g0.neg)
#       Dplus_g1.pos_sim <- rbinom(1, g1.pos, Dplus_g1.pos/g1.pos)
#       Dplus_g0.pos_sim <- rbinom(1, g0.pos, Dplus_g0.pos/g0.pos)
#       if (g1.unknown!=0) {
#         Dplus_g1.unknown_sim <- rbinom(1, g1.unknown, Dplus_g1.unknown/g1.unknown)
#         Dplus_g0.unknown_sim <- rbinom(1, g0.unknown, Dplus_g0.unknown/g0.unknown)
#       } else {
#         Dplus_g1.unknown_sim <- Dplus_g0.unknown_sim <- 0
#       }
#       
#       Dminus_g1.neg_sim <- g1.neg - Dplus_g1.neg_sim
#       Dminus_g0.neg_sim <- g0.neg - Dplus_g0.neg_sim
#       Dminus_g1.pos_sim <- g1.pos - Dplus_g1.pos_sim
#       Dminus_g0.pos_sim <- g0.pos - Dplus_g0.pos_sim
#       Dminus_g1.unknown_sim <- g1.unknown - Dplus_g1.unknown_sim
#       Dminus_g0.unknown_sim <- g0.unknown - Dplus_g0.unknown_sim
#       
#       Dplus_g1.all_sim <- Dplus_g1.neg_sim + Dplus_g1.pos_sim + Dplus_g1.unknown_sim
#       Dplus_g0.all_sim <- Dplus_g0.neg_sim + Dplus_g0.pos_sim + Dplus_g0.unknown_sim
#       Dminus_g1.all_sim <- g1.all - Dplus_g1.all_sim
#       Dminus_g0.all_sim <- g0.all - Dplus_g0.all_sim
#       
#       # Bootstrap the control-arm subsample (which incorporates loss-of-signal)
#       g0.pos.subsamp <- Dplus_g0.pos.subsamp + Dminus_g0.pos.subsamp
#       Dplus_g0.pos.subsamp_sim <- rbinom(1, g0.pos.subsamp, Dplus_g0.pos.subsamp/g0.pos.subsamp)
#       Dminus_g0.pos.subsamp_sim <- g0.pos.subsamp - Dplus_g0.pos.subsamp_sim
#       
#       g0.neg.subsamp <- Dplus_g0.neg.subsamp + Dminus_g0.neg.subsamp
#       Dplus_g0.neg.subsamp_sim <- rbinom(1, g0.neg.subsamp, Dplus_g0.neg.subsamp/g0.neg.subsamp)
#       Dminus_g0.neg.subsamp_sim <- g0.neg.subsamp - Dplus_g0.neg.subsamp_sim
#       
#       # Bootstrap control-arm subsample with loss-of-signal 
#       g0.pos.loss <- Dplus_g0.pos.loss + Dminus_g0.pos.loss
#       Dplus_g0.pos.loss_sim <- rbinom(1, g0.pos.loss, Dplus_g0.pos.loss/g0.pos.loss)
#       Dminus_g0.pos.loss_sim <- g0.pos.loss - Dplus_g0.pos.loss_sim
#       
#       g0.neg.loss <- Dplus_g0.neg.loss + Dminus_g0.neg.loss
#       Dplus_g0.neg.loss_sim <- rbinom(1, g0.neg.loss, Dplus_g0.neg.loss/g0.neg.loss)
#       Dminus_g0.neg.loss_sim <- g0.neg.loss - Dplus_g0.neg.loss_sim
#       
#       
#       # analyze bootstrap sim
#       sim_params[j+1,] <- ie_analysis(
#         Dplus_g1.neg_sim, Dplus_g0.neg_sim, g1.neg, g0.neg, Dminus_g1.neg_sim, Dminus_g0.neg_sim,
#         Dplus_g1.pos_sim, Dplus_g0.pos_sim, g1.pos, g0.pos, Dminus_g1.pos_sim, Dminus_g0.pos_sim,
#         Dplus_g1.unknown_sim, Dplus_g0.unknown_sim, g1.unknown, g0.unknown, Dminus_g1.unknown_sim, Dminus_g0.unknown_sim,
#         Dplus_g1.all_sim, Dplus_g0.all_sim, g1.all, g0.all, Dminus_g1.all_sim, Dminus_g0.all_sim,
#         Dplus_g0.pos.loss_sim, Dminus_g0.pos.loss_sim, Dplus_g0.neg.loss_sim, Dminus_g0.neg.loss_sim,
#         Dplus_g0.pos.subsamp_sim, Dminus_g0.pos.subsamp_sim, Dplus_g0.neg.subsamp_sim, Dminus_g0.neg.subsamp_sim,
#         Dplus_g1.pos.retest, Dminus_g1.pos.retest,
#         RR, RR_pos, RR_neg, RD_neg, RD_pos, RD, Z_ratio,
#         frac.g0.plus, frac.g0.minus, frac.g0.plus.pos, frac.g0.minus.pos,
#         pD_g1.pos, pD_g0.pos, pD_g1.neg, pD_g0.neg, pMplus_Dminus, pMplus_Dplus, pD_pos, pD_pos_check, pD_neg, pD)
#       
#     } # End bootstrap sims
#     
#     if (i == 1) { # initialize array if first superpopulation
#       # Define dimensions
#       rows <- rownames(sim_params)
#       cols <- colnames(sim_params)
#       # Initialize with dimnames
#       sim_outs <- array(NA, 
#                         dim = c(nrow(sim_params), ncol(sim_params), nsim),
#                         dimnames = list(rows, cols, NULL)) # Rows, Cols, Layers
#     }
#     
#     # Save bootstrap outputs for superpopulation i
#     sim_outs[,,i] <- as.matrix(sim_params)
#     
#   } # End simulating superpopulation
#   
#   return(sim_outs)
#   
# } # end function conceal.reveal()
# 
# 
# 
# 
# 
# ##
# # Function to do all IE analysis
# ##
# ie_analysis <- function(Dplus_g1.neg, Dplus_g0.neg, g1.neg, g0.neg, Dminus_g1.neg, Dminus_g0.neg,
#                         Dplus_g1.pos, Dplus_g0.pos, g1.pos, g0.pos, Dminus_g1.pos, Dminus_g0.pos,
#                         Dplus_g1.unknown, Dplus_g0.unknown, g1.unknown, g0.unknown, Dminus_g1.unknown, Dminus_g0.unknown,
#                         Dplus_g1.all, Dplus_g0.all, g1.all, g0.all, Dminus_g1.all, Dminus_g0.all,
#                         Dplus_g0.pos.loss, Dminus_g0.pos.loss, Dplus_g0.neg.loss, Dminus_g0.neg.loss,
#                         Dplus_g0.pos.subsamp, Dminus_g0.pos.subsamp, Dplus_g0.neg.subsamp, Dminus_g0.neg.subsamp,
#                         Dplus_g1.pos.retest, Dminus_g1.pos.retest,
#                         RR, RR_pos, RR_neg, RD_neg, RD_pos, RD, Z_ratio,
#                         frac.g0.plus, frac.g0.minus, frac.g0.plus.pos, frac.g0.minus.pos,
#                         pD_g1.pos, pD_g0.pos, pD_g1.neg, pD_g0.neg, pMplus_Dminus, pMplus_Dplus, pD_pos, pD_pos_check, pD_neg, pD)
# {
#   
#   ##
#   # initialize simulation
#   ##
#   Z_neg <- Z_pos <- Z_all <- p.val_neg <- p.val_pos <- p.val_all <- RD_neg.sim <- RD_pos.sim <- RD_all.sim <- RR_neg.sim <- RR_pos.sim <- RR_all.sim <- NA
#   Z_tar <- p.val_tar  <- RD_tar.sim  <- RR_tar.sim  <- RR_pos.noncomp.sim <- RR_tar.noncomp.sim <- RR_neg.noncomp.sim <- NA
#   p.val_neg.subsamp <- p.val_pos.subsamp <- RD_neg.subsamp.sim <- RD_pos.subsamp.sim <- p.val_pos.noncomp <- p.val_tar.noncomp <- p.val_neg.noncomp <- NA
#   RD_pos.loss.ss.sim <- RR_pos.loss.ss.sim <- RD_neg.loss.sim <- RD_pos.loss.sim <- RR_neg.loss.sim <- RR_pos.loss.sim <- p_Mplus.g0.hat <- p_Mplus.g1.hat <- check.RD_all.hat <- NA
#   RD_pos.retest <- RD_pos.BF <- RD_neg.BF <- RD_neg.BF.retest <- RD_pos.BF.retest <- RR_neg.BF.retest <- RR_pos.BF.retest <- NA
#   RD_pos_IE.sim <- RD_IE.sim <- RR_pos_IE.sim <- p.val.RD_pos_IE.sim <- NA
#   
#   ##
#   # Analysis 
#   ##
#   
#   ##### 
#   # Calculate naive statistics assuming all assumptions hold
#   # calculate statistics for M-, M+, and marginal tables
#   neg.out <- prop.test(c(Dplus_g1.neg, Dplus_g0.neg), c(g1.neg, g0.neg), correct=F)
#   pos.out <- prop.test(c(Dplus_g1.pos, Dplus_g0.pos), c(g1.pos, g0.pos), correct=F)
#   all.out <- prop.test(c(Dplus_g1.all, Dplus_g0.all), c(g1.all, g0.all), correct=F)
#   tar.out <- prop.test(c(Dplus_g1.pos, Dplus_g0.pos), c(g1.all, g0.all), correct=F)
#   
#   #####
#   # Correct for non-compliance to control-arm blood collections
#   # Correction for ever-positives (RRpos and RRtar) using observable data from the unknown-positive table
#   # The correction factor is the ratio of the noncensored fractions of cases in the screen vs control arms
#   Dplus_g0.pos.corrected <- Dplus_g0.pos * (1-Dplus_g1.unknown/Dplus_g1.all) / 
#     (1-Dplus_g0.unknown/Dplus_g0.all)
#   Dminus_g0.pos.corrected <- Dminus_g0.pos * (1-Dminus_g1.unknown/Dminus_g1.all) / 
#     (1-Dminus_g0.unknown/Dminus_g0.all)
#   g0.pos.corrected <- Dplus_g0.pos.corrected + Dminus_g0.pos.corrected
#   
#   # Correction for never-positives (RRneg) using observable data from the unknown-positive table
#   # The correction factor is the ratio of the noncensored fractions of cases in the screen vs control arms
#   Dplus_g0.neg.corrected <- Dplus_g0.neg * (1-Dplus_g1.unknown/Dplus_g1.all) / 
#     (1-Dplus_g0.unknown/Dplus_g0.all)
#   Dminus_g0.neg.corrected <- Dminus_g0.neg * (1-Dminus_g1.unknown/Dminus_g1.all) / 
#     (1-Dminus_g0.unknown/Dminus_g0.all)
#   g0.neg.corrected <- Dplus_g0.neg.corrected + Dminus_g0.neg.corrected
#   
#   # Calculate RRpos and RRneg with control-arm non-compliance corrected to equal that of the screen-arm
#   pos.noncomp.out <- prop.test(c(Dplus_g1.pos, Dplus_g0.pos.corrected), c(g1.pos, g0.pos.corrected), correct=F)
#   neg.noncomp.out <- prop.test(c(Dplus_g1.neg, Dplus_g0.neg.corrected), c(g1.neg, g0.neg.corrected), correct=F)
#   # Calculate RRtar with control-arm non-compliance corrected to equal that of the screen-arm
#   tar.noncomp.out <- prop.test(c(Dplus_g1.pos, Dplus_g0.pos.corrected), c(g1.all, g0.all), correct=F)
#   
#   
#   
#   # RDs and RRs
#   RD_neg.sim <- diff(neg.out$estimate)
#   RD_pos.sim <- diff(pos.out$estimate)
#   RD_all.sim <- diff(all.out$estimate)
#   RD_tar.sim <- diff(tar.out$estimate)
#   RR_neg.sim <- neg.out$estimate[1]/neg.out$estimate[2]
#   RR_pos.sim <- pos.out$estimate[1]/pos.out$estimate[2]
#   RR_all.sim <- all.out$estimate[1]/all.out$estimate[2]
#   RR_tar.sim <- tar.out$estimate[1]/tar.out$estimate[2]
#   RR_pos.noncomp.sim <- pos.noncomp.out$estimate[1]/pos.noncomp.out$estimate[2]
#   RR_neg.noncomp.sim <- neg.noncomp.out$estimate[1]/neg.noncomp.out$estimate[2]
#   RR_tar.noncomp.sim <- tar.noncomp.out$estimate[1]/tar.noncomp.out$estimate[2]
#   
#   # p-values for the RDs
#   p.val_neg <- neg.out$p.val
#   p.val_pos <- pos.out$p.val
#   p.val_all <- all.out$p.val
#   p.val_tar <- tar.out$p.val
#   p.val_pos.noncomp <- pos.noncomp.out$p.val
#   p.val_neg.noncomp <- neg.noncomp.out$p.val
#   p.val_tar.noncomp <- tar.noncomp.out$p.val
#   
#   # Z-stats: need to get the sign right because power<50% implies that mean Z-stat is negative
#   Z_neg  <- sign(RD_neg.sim)*sqrt(neg.out$stat)
#   Z_pos  <- sign(RD_pos.sim)*sqrt(pos.out$stat)
#   Z_all  <- sign(RD_all.sim)*sqrt(all.out$stat)
#   Z_tar  <- sign(RD_tar.sim)*sqrt(tar.out$stat)
#   
#   
#   
#   
#   #####
#   # Correct for control-arm subsampling of bloods to be tested
#   #
#   # Estimate sampling fractions for D+ and for D- within arm G=0
#   # These are P(subsampled|D+,G=0) and P(subsampled|D-,G=0)
#   # Estimating sampling fractions is important for D+ because it is the rare event
#   # Not important for D- because there are so many of them so true fraction is close 
#   # to the estimated fraction
#   ##
#   frac.g0.plus.hat <- (Dplus_g0.pos.subsamp + Dplus_g0.neg.subsamp) /
#     (Dplus_g0.pos + Dplus_g0.neg)
#   frac.g0.minus.hat<- (Dminus_g0.pos.subsamp + Dminus_g0.neg.subsamp) / 
#     ((g0.pos - Dplus_g0.pos) + (g0.neg - Dplus_g0.neg))
#   
#   ##
#   # calculate stats for M-,M+ with a subsampled table in control arm, 
#   # weight up each total
#   ##
#   
#   # Do for M- table: RD_neg = p_control - p_screen
#   # weight #(D+,M-,G=0,subsampled) * 1/P(subsampled|D+,G=0)
#   Dplus_g0.neg.subsamp.weighted <- Dplus_g0.neg.subsamp*1/frac.g0.plus.hat
#   # also weight #(D-,M-,G=0,subsampled) * 1/P(subsampled|D-,G=0)
#   g0.neg.weighted <- Dplus_g0.neg.subsamp.weighted + 
#     Dminus_g0.neg.subsamp*1/frac.g0.minus.hat
#   p_screen.neg  <- Dplus_g1.neg / g1.neg  # P(D+|M-,G=1)
#   p_control.neg <- Dplus_g0.neg.subsamp.weighted / g0.neg.weighted  # P(D+|M-,G=0)
#   RD_neg.subsamp.sim <- p_control.neg - p_screen.neg
#   RR_neg.subsamp.sim <- p_screen.neg/p_control.neg
#   
#   # Do for M+ table: RD_pos = p_control - p_screen
#   Dplus_g0.pos.subsamp.weighted <- Dplus_g0.pos.subsamp*1/frac.g0.plus.hat
#   g0.pos.weighted <- Dplus_g0.pos.subsamp.weighted + 
#     Dminus_g0.pos.subsamp*1/frac.g0.minus.hat
#   p_screen.pos  <- Dplus_g1.pos / g1.pos  # P(D+|M+,G=1)
#   p_control.pos <- Dplus_g0.pos.subsamp.weighted / g0.pos.weighted  # P(D+|M+,G=0)
#   RD_pos.subsamp.sim <- p_control.pos - p_screen.pos
#   RR_pos.subsamp.sim <- p_screen.pos/p_control.pos
#   
#   
#   
#   #####
#   # Account for loss-of-signal in control-arm bloods
#   ##
#   
#   # Estimate sampling fractions for D+ and for D- within arm G=0
#   frac.g0.plus.hat <- (Dplus_g0.pos.loss + Dplus_g0.neg.loss) /
#     (Dplus_g0.pos + Dplus_g0.neg)
#   frac.g0.minus.hat<- (Dminus_g0.pos.loss + Dminus_g0.neg.loss) /
#     ((g0.pos - Dplus_g0.pos) + (g0.neg - Dplus_g0.neg))
#   
#   # Do for M- table: RD_neg = p_control - p_screen
#   Dplus_g0.neg.loss.weighted <- Dplus_g0.neg.loss*1/frac.g0.plus.hat
#   g0.neg.weighted <- Dplus_g0.neg.loss.weighted + 
#     Dminus_g0.neg.loss*1/frac.g0.minus.hat
#   p_screen.neg  <- Dplus_g1.neg / g1.neg  # P(D+|M-,G=1)
#   p_control.neg <- Dplus_g0.neg.loss.weighted / g0.neg.weighted  # P(D+|M-,G=0)
#   RD_neg.loss.sim <- p_control.neg - p_screen.neg
#   RR_neg.loss.sim <- p_screen.neg / p_control.neg
#   
#   # Do for M+ table: RD_pos = p_control - p_screen
#   Dplus_g0.pos.loss.weighted <- Dplus_g0.pos.loss*1/frac.g0.plus.hat
#   g0.pos.weighted <- Dplus_g0.pos.loss.weighted + 
#     Dminus_g0.pos.loss*1/frac.g0.minus.hat
#   p_screen.pos  <- Dplus_g1.pos / g1.pos  # P(D+|M+,G=1)
#   p_control.pos <- Dplus_g0.pos.loss.weighted / g0.pos.weighted  # P(D+|M+,G=0)
#   RD_pos.loss.sim <- p_control.pos - p_screen.pos
#   RR_pos.loss.sim <- p_screen.pos / p_control.pos
#   
#   
#   # Check if estimated P(M+|G=0) is close to the true p_Mplus we input
#   # This will not be true if there is any dropout
#   p_Mplus.g0.hat <- g0.pos.weighted/(g0.pos+g0.neg)
#   p_Mplus.g1.hat <- g1.pos/(g1.pos+g1.neg)
#   
#   # Check if RD_pos and RD_neg are consistent with RD_all in the simulation
#   # when using p_Mplus.g0.hat as the (possibly biased) estimate for P(M+)
#   check.RD_all.hat <- RD_pos.loss.sim*p_Mplus.g0.hat + 
#     RD_neg.loss.sim*(1-p_Mplus.g0.hat)
#   
#   ##
#   # Calculate corrected RD_neg using prior odds and assuming randomization worked
#   # so can substitute P(M+|G=1) for P(M+|G=0)
#   # But this requires no loss of signal, so it cannot generalize
#   ##
#   
#   # prior odds = P(D+|G=0)/P(D-|G=0) = #(D+,G=0)/#(D-,G=0) 
#   # = {#(D+,M-,G=0)+#(D+,M+,G=0)} / {#(M-,G=0)+#(M+,G=0) - [#(D+,M-,G=0)+#(D+,M+,G=0)]}
#   Dplus_g0 <- Dplus_g0.neg + Dplus_g0.pos
#   Dminus_g0 <- g0.neg + g0.pos - Dplus_g0
#   odds.prior <- Dplus_g0 / Dminus_g0
#   
#   # Calculate P(M+|D+,G=0) = #(D+,M+,G=0)/#(D+,G=0), need to weight up numerator
#   p_Mplus_Dplus.g0.weighted <- Dplus_g0.pos.loss.weighted / Dplus_g0
#   
#   # Calculate P(M+|D-,G=0) = #(D-,M+,G=0)/#(D-,G=0)
#   p_Mplus_Dminus.g0.weighted <- (Dminus_g0.pos.loss*1/frac.g0.minus.hat) / Dminus_g0
#   
#   
#   #### From the validation study  
#   # Estimate retention of signal (i.e. 1-loss): the 2 retest-positive fractions
#   # P(~M+|M+,D+,G=1)
#   frac.g1.plus.pos.hat <- Dplus_g1.pos.retest / Dplus_g1.pos
#   # P(~M+|M+,D-,G=1)
#   frac.g1.minus.pos.hat <- Dminus_g1.pos.retest / (g1.pos - Dplus_g1.pos)
#   
#   # Now recalculate BF_neg, BF_pos, P(D+|M+/-,G=0), and RD_pos and RD_neg
#   BF_neg <- (1 - p_Mplus_Dplus.g0.weighted*1/frac.g1.plus.pos.hat) / 
#     (1 - p_Mplus_Dminus.g0.weighted*1/frac.g1.minus.pos.hat)
#   BF_pos <- (p_Mplus_Dplus.g0.weighted*1/frac.g1.plus.pos.hat) / 
#     (p_Mplus_Dminus.g0.weighted*1/frac.g1.minus.pos.hat)
#   pDplus_g0.neg <- 1/(1 + (odds.prior * BF_neg)^-1)
#   pDplus_g0.pos <- 1/(1 + (odds.prior * BF_pos)^-1)
#   RD_neg.BF.retest <- pDplus_g0.neg - Dplus_g1.neg / g1.neg
#   RD_pos.BF.retest <- pDplus_g0.pos - Dplus_g1.pos / g1.pos
#   RR_neg.BF.retest <- (Dplus_g1.neg / g1.neg) / pDplus_g0.neg
#   RR_pos.BF.retest <- (Dplus_g1.pos / g1.pos) / pDplus_g0.pos
#   
#   # Calculate estimated P(M+|G)
#   p_Mplus.g1.hat <- g1.pos/(g1.pos+g1.neg) # G=1 is easy because no loss of signal
#   # G=0 has to account for loss of signal:
#   # P(M+|G=0) = P(M+|D+,G=0)P(D+|G=0) + P(M+|D-,G=0)P(D-|G=0)
#   p_Mplus.g0.hat <- ( (p_Mplus_Dplus.g0.weighted*1/frac.g1.plus.pos.hat)*Dplus_g0 + 
#                         (p_Mplus_Dminus.g0.weighted*1/frac.g1.minus.pos.hat)*Dminus_g0 ) /
#     (Dplus_g0+Dminus_g0)
#   
#   # Use retest-positive fractions directly to correct P(D+|~M+,G=0) = p_control.pos
#   # so it estimates P(D+|M+,G=0). See 03_Intended-Effect-Stat.docx
#   # First calculate P(~M+|M+,G=1) directly
#   frac.g1.pos.hat <- (Dplus_g1.pos.retest + Dminus_g1.pos.retest) / g1.pos
#   # Now calculate P(D+|M+,G=0) = P(D+|~M+,G=0) * {P(~M+|M+,G=1)/P(~M+|M+,D+,G=1)}
#   # and plug into RDpos
#   RD_pos.retest <- p_control.pos*(frac.g1.pos.hat/frac.g1.plus.pos.hat) - Dplus_g1.pos / g1.pos
#   
#   
#   ##
#   # Pack up output for return
#   ##
#   
#   nsim <- 1 # hack since we are now doing 1 sim at a time but reusing old code
#   out <- data.frame(RR=RR, RR_pos=RR_pos, RR_neg=RR_neg, RD_neg=RD_neg, RD_pos=RD_pos, RD=RD,
#                     RD_neg.sim=mean(RD_neg.sim), RD_pos.sim=mean(RD_pos.sim), RD.sim=mean(RD_all.sim),
#                     Z_ratio=Z_ratio,
#                     Z_ratio.sim = mean(Z_pos)/mean(Z_all),
#                     # power.neg=mean(p.val_neg <= .05), # power for Fisher exact test not Wald test
#                     # power.pos=mean(p.val_pos <= .05),
#                     # power.all=mean(p.val_all <= .05),
#                     power.RDneg=mean(2*pnorm(-abs(RD_neg.sim)/sd(RD_neg.sim)) <= 0.05),
#                     power.RDpos=mean(2*pnorm(-abs(RD_pos.sim)/sd(RD_pos.sim)) <= 0.05),
#                     power.RDall=mean(2*pnorm(-abs(RD_all.sim)/sd(RD_all.sim)) <= 0.05),
#                     power.RRneg=mean(2*pnorm(-abs(log(RR_neg.sim))/sd(log(RR_neg.sim))) <= 0.05),
#                     power.RRpos=mean(2*pnorm(-abs(log(RR_pos.sim))/sd(log(RR_pos.sim))) <= 0.05),
#                     power.RRall=mean(2*pnorm(-abs(log(RR_all.sim))/sd(log(RR_all.sim))) <= 0.05),
#                     # p.neg=median(p.val_neg),
#                     # p.pos= median(p.val_pos),
#                     # p.all= median(p.val_all),
#                     # p.ratio=median(p.val_all)/median(p.val_pos),
#                     p.pos.lt.p.all=sum(p.val_pos<p.val_all)/nsim,
#                     RD_neg.subsamp=mean(RD_neg.subsamp.sim), 
#                     RD_pos.subsamp=mean(RD_pos.subsamp.sim), 
#                     # power.neg.subsamp=mean(p.val_neg.subsamp <= .05),
#                     # power.pos.subsamp=mean(p.val_pos.subsamp <= .05),
#                     # p.neg.subsamp=median(p.val_neg.subsamp),
#                     # p.pos.subsamp=median(p.val_pos.subsamp),
#                     # p.ratio.subsamp=median(p.val_all)/median(p.val_pos.subsamp),
#                     # p.pos.gt.p.all.subsamp=sum(p.val_pos.subsamp>p.val_all)/nsim,
#                     sd.RD_neg = sd(RD_neg.sim),
#                     sd.RD_pos = sd(RD_pos.sim),
#                     sd.RD_pos.subsamp = sd(RD_pos.subsamp.sim),
#                     sd.ratio = sd(RD_pos.subsamp.sim)/sd(RD_pos.sim),
#                     sd.RD_neg.subsamp = sd(RD_neg.subsamp.sim),
#                     # p.sim=median(2*pnorm(-abs(RD_all.sim)/sd(RD_all.sim))), # RD_all p-val agrees with simulation
#                     # p.pos.sim=median(2*pnorm(-abs(RD_pos.sim)/sd(RD_pos.sim))),# RD_pos p-val agrees with simulation
#                     # p.pos.subsamp.sim = median(2*pnorm(-abs(RD_pos.subsamp.sim)/sd(RD_pos.subsamp.sim))),
#                     power.pos.subsamp = mean(2*pnorm(-abs(RD_pos.subsamp.sim)/sd(RD_pos.subsamp.sim)) <= 0.05),
#                     # p.neg.subsamp.sim = median(2*pnorm(-abs(RD_neg.subsamp.sim)/sd(RD_neg.subsamp.sim))),
#                     power.neg.subsamp = mean(2*pnorm(-abs(RD_neg.subsamp.sim)/sd(RD_neg.subsamp.sim)) <= 0.05),
#                     RD_pos.loss = mean(RD_pos.loss.sim),
#                     RD_neg.loss = mean(RD_neg.loss.sim),
#                     sd.RD_pos.loss = sd(RD_pos.loss.sim),
#                     power.RD_pos.loss = mean(2*pnorm(-abs(RD_pos.loss.sim)/sd(RD_pos.loss.sim)) <= 0.05),
#                     power.RD_neg.loss = mean(2*pnorm(-abs(RD_neg.loss.sim)/sd(RD_neg.loss.sim)) <= 0.05),
#                     RR_pos.loss = mean(RR_pos.loss.sim),
#                     RR_neg.loss = mean(RR_neg.loss.sim),
#                     sd.RR_pos.loss = sd(RR_pos.loss.sim),
#                     power.RR_pos.loss = mean(2*pnorm(-abs(log(RR_pos.loss.sim))/sd(log(RR_pos.loss.sim))) <= 0.05),
#                     power.RR_neg.loss = mean(2*pnorm(-abs(log(RR_neg.loss.sim))/sd(log(RR_neg.loss.sim))) <= 0.05),
#                     p_Mplus.g0.hat = mean(p_Mplus.g0.hat),
#                     p_Mplus.g1.hat = mean(p_Mplus.g1.hat),
#                     check.RD_all.hat = mean(check.RD_all.hat),
#                     # RD_neg.BF = mean(RD_neg.BF),         
#                     # sd.RD_neg.BF = sd(RD_neg.BF),
#                     # power.neg.BF = mean(2*pnorm(-abs(RD_neg.BF)/sd(RD_neg.BF)) <= 0.05),
#                     # RD_pos.BF = mean(RD_pos.BF),         
#                     # sd.RD_pos.BF = sd(RD_pos.BF),
#                     # power.pos.BF = mean(2*pnorm(-abs(RD_pos.BF)/sd(RD_pos.BF)) <= 0.05),
#                     RD_neg.BF.retest=mean(RD_neg.BF.retest),
#                     RD_pos.BF.retest=mean(RD_pos.BF.retest),
#                     sd.RD_neg.BF.retest = sd(RD_neg.BF.retest),
#                     sd.RD_pos.BF.retest = sd(RD_pos.BF.retest),
#                     power.RD_neg.BF.retest=mean(2*pnorm(-abs(RD_neg.BF.retest)/sd(RD_neg.BF.retest)) <= 0.05),
#                     power.RD_pos.BF.retest=mean(2*pnorm(-abs(RD_pos.BF.retest)/sd(RD_pos.BF.retest)) <= 0.05),
#                     RR_neg.BF.retest=mean(RR_neg.BF.retest),
#                     RR_pos.BF.retest=mean(RR_pos.BF.retest),
#                     sd.RR_neg.BF.retest = sd(RR_neg.BF.retest),
#                     sd.RR_pos.BF.retest = sd(RR_pos.BF.retest),
#                     power.RR_neg.BF.retest=mean(2*pnorm(-abs(log(RR_neg.BF.retest))/sd(log(RR_neg.BF.retest))) <= 0.05),
#                     power.RR_pos.BF.retest=mean(2*pnorm(-abs(log(RR_pos.BF.retest))/sd(log(RR_pos.BF.retest))) <= 0.05),
#                     RD_pos_IE = mean(RD_pos_IE.sim),
#                     sd.RD_pos_IE = sd(RD_pos_IE.sim),
#                     power.RD_pos_IE = mean(2*pnorm(-abs(RD_pos_IE.sim)/sd(RD_pos_IE.sim)) <= 0.05),
#                     RD_IE = mean(RD_IE.sim),
#                     sd.RD_IE = sd(RD_IE.sim),
#                     sd.RD_all = sd(RD_all.sim),
#                     # RR_pos_IE = mean(RR_pos_IE.sim),
#                     # sd.RR_pos_IE = sd(RR_pos_IE.sim),
#                     # power.RR_pos_IE = mean(2*pnorm(-abs(log(RR_pos_IE.sim))/sd(log(RR_pos_IE.sim))) <= 0.05),
#                     medianp.RD_pos = median(2*pnorm(-abs(RD_pos.sim)/sd(RD_pos.sim))),
#                     medianp.RD_pos_IE = median(p.val.RD_pos_IE.sim),
#                     medianp.RD_all = median(2*pnorm(-abs(RD_all.sim)/sd(RD_all.sim))),
#                     p.RD_pos_IE_lt_p.RD_all = sum(p.val.RD_pos_IE.sim < p.val_all)/nsim,
#                     p.RD_pos_IE_lt_p.RD_pos = sum(p.val.RD_pos_IE.sim < p.val_pos)/nsim,
#                     RR_all = mean(RR_all.sim),
#                     RD_pos.retest = mean(RD_pos.retest),
#                     power.RD_pos.retest = mean(2*pnorm(-abs(RD_pos.retest)/sd(RD_pos.retest)) <= 0.05),
#                     RR_pos.sim = mean(RR_pos.sim),
#                     RD_pos.loss.ss = mean(RD_pos.loss.ss.sim),
#                     # RD_neg.loss = mean(RD_neg.loss.sim),
#                     # sd.RD_pos.loss = sd(RD_pos.loss.sim),
#                     power.RD_pos.loss.ss = mean(2*pnorm(-abs(RD_pos.loss.ss.sim)/sd(RD_pos.loss.ss.sim)) <= 0.05),
#                     # power.RD_neg.loss = mean(2*pnorm(-abs(RD_neg.loss.sim)/sd(RD_neg.loss.sim)) <= 0.05),
#                     RR_pos.loss.ss = mean(RR_pos.loss.ss.sim),
#                     # RR_neg.loss = mean(RR_neg.loss.sim),
#                     # sd.RR_pos.loss = sd(RR_pos.loss.sim),
#                     power.RR_pos.loss.ss = mean(2*pnorm(-abs(log(RR_pos.loss.ss.sim))/sd(log(RR_pos.loss.ss.sim))) <= 0.05),
#                     # power.RR_neg.loss = mean(2*pnorm(-abs(log(RR_neg.loss.sim))/sd(log(RR_neg.loss.sim))) <= 0.05),
#                     RR.tar=mean(RR_tar.sim),
#                     power.RRtar=mean(2*pnorm(-abs(log(RR_tar.sim))/sd(log(RR_tar.sim))) <= 0.05),
#                     RR_pos.noncomp = mean(RR_pos.noncomp.sim),
#                     RR_tar.noncomp = mean(RR_tar.noncomp.sim),
#                     power.RRpos.noncomp = mean(2*pnorm(-abs(log(RR_pos.noncomp.sim))/sd(log(RR_pos.noncomp.sim))) <= 0.05),
#                     power.RRtar.noncomp = mean(2*pnorm(-abs(log(RR_tar.noncomp.sim))/sd(log(RR_tar.noncomp.sim))) <= 0.05),
#                     RR_neg.noncomp = mean(RR_neg.noncomp.sim),
#                     power.RRneg.noncomp = mean(2*pnorm(-abs(log(RR_neg.noncomp.sim))/sd(log(RR_neg.noncomp.sim))) <= 0.05),
#                     RR_neg.sim = mean(RR_neg.sim),
#                     RR_pos.subsamp = mean(RR_pos.subsamp.sim),
#                     RR_neg.subsamp = mean(RR_neg.subsamp.sim),
#                     pD_g1.pos=pD_g1.pos,
#                     pD_g0.pos=pD_g0.pos,
#                     pD_g1.neg=pD_g1.neg,
#                     pD_g0.neg=pD_g0.neg,
#                     pMplus_Dminus=pMplus_Dminus,
#                     pMplus_Dplus=pMplus_Dplus,
#                     pD_pos=pD_pos,
#                     pD_pos_check,
#                     pD_neg=pD_neg,
#                     pD=pD
#   )
#   
#   # print(cbind(RD_pos_IE.sim,RD_pos.sim, p.val_pos, p.val.RD_pos_IE.sim, p.val_all))
#   
#   return(out)
#   
#   # return(list(RD_neg.sim=RD_neg.sim,RD_pos.sim=RD_pos.sim,RD_all.sim=RD_all.sim,
#   #              p.val_neg=p.val_neg,p.val_pos=p.val_pos,p.val_all=p.val_all,
#   #              Z_neg=Z_neg,Z_pos=Z_pos,Z_all=Z_all)
#   #        )
# }
# 
# 
# sd_summary <- function(x, ...) {
#   output <- c(summary(x, ...), "sd" = sd(x, na.rm = TRUE))
#   return(output)
# }
# 
# 
# ##
# # 1. Bootstrap for all assumptions hold: base case
# ##
# system.time({
#   basecase <- conceal.reveal(nsim=200, bootsims=100, sim=TRUE
#   )
#   res <- Sys.sleep(1) 
#   system('say "Finished!"') # MacOS only
# })
# mean.RRpos <- apply(basecase[-1,"RR_pos.sim",], 2, mean)
# se.RRpos <- apply(basecase[-1,"RR_pos.sim",], 2, sd)
# true.RRpos <- 0.867
# coverage <- mean(abs(mean.RRpos - true.RRpos) <= qnorm(0.975) * se.RRpos,  na.rm = TRUE)
# 
# # # Ignore RRneg
# # mean(basecase[-1,"RR_neg.sim",])
# # mean(apply(basecase[-1,"RR_neg.sim",], 2, sd))
# 
# 
# ##
# # 2. Bootstrap variance estimators that correct for loss of signal
# ##
# system.time({
#   basecase <- conceal.reveal(nsim=1, bootsims=1, sim=TRUE,
#                              frac.g0.plus.pos = 1, frac.g0.minus.pos = 1
#   )
#   res <- Sys.sleep(1) 
#   system('say "Finished!"') # MacOS only
# })
# # Τhis is the bias and lack of coverage when ignoring loss-of-signal
# est.RRpos <- basecase[1,"RR_pos.loss",]
# sd_summary(est.RRpos)
# se.RRpos <- apply(basecase[-1,"RR_pos.loss",], 2, sd)
# sd_summary(se.RRpos)
# true.RRpos <- 0.867
# mean(abs(est.RRpos - true.RRpos) <= qnorm(0.975) * se.RRpos,  na.rm = TRUE)
# # Τhis is the bias and lack of coverage when accounting loss-of-signal
# est.RRpos <- basecase[1,"RR_pos.BF.retest",]
# sd_summary(est.RRpos)
# se.RRpos <- apply(basecase[-1,"RR_pos.BF.retest",], 2, sd)
# sd_summary(se.RRpos)
# true.RRpos <- 0.867
# mean(abs(est.RRpos - true.RRpos) <= qnorm(0.975) * se.RRpos,  na.rm = TRUE)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## quantities for: base case all assumps hold
# # Observed means and stdevs across superpopulations
# mean(basecase[1,"RR_pos.sim",])
# sd(basecase[1,"RR_pos.sim",])
# mean(basecase[1,"RR_neg.sim",])
# sd(basecase[1,"RR_neg.sim",])
# # Boostrap means and mean of the bootstrap stdevs
# mean(basecase[-1,"RR_pos.sim",])
# mean(apply(basecase[-1,"RR_pos.sim",], 2, sd))
# mean(basecase[-1,"RR_neg.sim",])
# mean(apply(basecase[-1,"RR_neg.sim",], 2, sd))
# 
# ## quantities for: corrected with sampling weights
# # Observed means and stdevs across superpopulations
# mean(basecase[1,"RR_pos.subsamp",])
# sd(basecase[1,"RR_pos.subsamp",])
# mean(basecase[1,"RR_neg.subsamp",])
# sd(basecase[1,"RR_neg.subsamp",])
# # Boostrap means and mean of the bootstrap stdevs
# mean(basecase[-1,"RR_pos.subsamp",])
# mean(apply(basecase[-1,"RR_pos.subsamp",], 2, sd))
# mean(basecase[-1,"RR_neg.subsamp",])
# mean(apply(basecase[-1,"RR_neg.subsamp",], 2, sd))
# 
# ## quantities for: loss of signal corrected by validation study in screened arm 
# ## (G=1) that retests all M+ for loss of signal
# # Observed means and stdevs across superpopulations
# mean(basecase[1,"RR_pos.BF.retest",])
# sd(basecase[1,"RR_pos.BF.retest",])
# mean(basecase[1,"RR_neg.BF.retest",])
# sd(basecase[1,"RR_neg.BF.retest",])
# # Boostrap means and mean of the bootstrap stdevs
# mean(basecase[-1,"RR_pos.BF.retest",])
# mean(apply(basecase[-1,"RR_pos.BF.retest",], 2, sd))
# mean(basecase[-1,"RR_neg.BF.retest",])
# mean(apply(basecase[-1,"RR_neg.BF.retest",], 2, sd))
# 
# ## quantities for: corrected for non-compliance with control-arm blood collections
# # Observed means and stdevs across superpopulations
# mean(basecase[1,"RR_pos.noncomp",])
# sd(basecase[1,"RR_pos.noncomp",])
# mean(basecase[1,"RR_neg.noncomp",])
# sd(basecase[1,"RR_neg.noncomp",])
# # Boostrap means and mean of the bootstrap stdevs
# mean(basecase[-1,"RR_pos.noncomp",])
# mean(apply(basecase[-1,"RR_pos.noncomp",], 2, sd))
# mean(basecase[-1,"RR_neg.noncomp",])
# mean(apply(basecase[-1,"RR_neg.noncomp",], 2, sd))
# 
# 
# 
# 
# 
# 
# # the options to turn debugging on/off
# # options(error = recover) # debugging on
# # options(error = NULL) # debugging off
# 
# 
# 
# 
# 
# 
# 
# 
# # Below is old stuff from previous file
# #####
# ##### OLD STUFF FROM PREVIOUS FILE
# 
# 
# 
# ##
# # Table 3: effect of non-adherence to blood draws (dropout)
# # Compare IE to Targeted
# ##
# # Simulate and save the data
# basecase <- conceal.reveal(nsim=1, sim=TRUE)
# # n <- c(25e3,50e3,75e3,100e3,150e3,200e3) #seq(120e3,150e3,by=10e3)
# # RR_poss <- c(0.86667,0.8)
# # RR_negs <- c(1,1.05,0.95)
# # RR <- 0.9
# n <- 100e3 #c(10e3, 15e3, 20e3, 25e3,30e3,40e3,50e3) #seq(120e3,150e3,by=10e3)
# RR_poss <-0.86667 #c(0.8,0.86667)
# RR_negs <- 1 #c(1,1.05,0.95)
# RR <- 0.9
# censor <- matrix(c(0.00,	0.00,	0.00,	0.00, 0.00, 0.00,
#                    0.30,	0.30,	0.30,	0.30, 0.00, 0.00,
#                    0.20,	0.20,	0.30,	0.30, 0.00, 0.00,
#                    0.30,	0.10,	0.30,	0.10, 0.00, 0.00,
#                    0.30,	0.15,	0.15,	0.30, 0.00, 0.00,
#                    0.15,	0.30,	0.30,	0.15, 0.00, 0.00),ncol=6,byrow=TRUE)
# colnames(censor)<- c("screen D+","screen D-","control D+","control D-", "screen M+ to D-", "screen M- to D-")
# out <- matrix(NA,nrow=1,ncol=ncol(basecase))
# for (l in 1:nrow(censor)) {
#   for (k in 1:length(RR_negs)) {
#     for (j in 1:length(RR_poss)) {
#       for (i in 1:length(n)) {
#         out <- rbind(out,
#                      as.matrix(conceal.reveal(nsim=2e4, n=n, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02,
#                                               RR=RR,RR_neg=RR_negs[k], RR_pos=RR_poss[j],
#                                               frac.g0.plus = 1, 
#                                               frac.g0.minus = 1,
#                                               frac.g0.plus.pos = 1, 
#                                               frac.g0.minus.pos = 1,
#                                               censor.screen.Dplus = censor[l,1],
#                                               censor.control.Dplus = censor[l,3],
#                                               censor.screen.Dminus = censor[l,2],
#                                               censor.control.Dminus = censor[l,4],
#                                               move.screen.Mplus.Dplus = censor[l,5],
#                                               move.screen.Mminus.Dplus = censor[l,6],
#                                               sim=TRUE)) )
#       }}}}
# # basecase <- as.data.frame(cbind(censor, out[-1,4+c(38,63,43,60,48,34,35)]))
# # names(basecase)<-c(colnames(censor),"RDpos","RRpos","RRneg","RR","power RRpos","P(M+|G=0)","P(M+|G=1)")
# basecase <- as.data.frame(cbind(censor, out[-1,4+c(63,43,60,48,68,69,70,71,74,72,73,75,34,35)]))
# names(basecase)<-c(colnames(censor),"RRpos","RRneg","RR","power RRpos","RRtar","power RRtar","RRposFix",
#                    "RRtarFix","RRnegFix","power RRposFix","power RRtarFix","power RRnegFix","P(M+|G=0)","P(M+|G=1)")
# print(t(basecase), digits=3, row.names=FALSE)
# 
# 
# 
# 
# 
# 
# 
# 
# ##
# # Compare IE to targeted IE
# ##
# 
# # JNCI example, 20% vs 25% differential non-compliance in the ppt file Tiger-Team-6Sep2022.ppt
# conceal.reveal(nsim=1, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos=0.86667,
#                frac.g0.plus = 1, frac.g0.minus = 1,
#                censor.screen.Dplus = 0.2,
#                censor.control.Dplus = 0.2,
#                censor.screen.Dminus = 0.2,
#                censor.control.Dminus = 0.2,
#                move.screen.Mplus.Dplus = 0,
#                move.screen.Mminus.Dplus = 0,
#                sim=FALSE)
# 
# # JNCI example, strongly differential non-compliance in the ppt file Tiger-Team-6Sep2022.ppt
# conceal.reveal(nsim=60e3, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos=0.86667,
#                frac.g0.plus = 1, frac.g0.minus = 1,
#                censor.screen.Dplus = 0.4,
#                censor.control.Dplus = 0.8,
#                censor.screen.Dminus = 0.8,
#                censor.control.Dminus = 0.4,
#                move.screen.Mplus.Dplus = 0,
#                move.screen.Mminus.Dplus = 0,
#                sim=TRUE)
# 
# # Redo JNCI example with IE vs Targeted
# conceal.reveal(nsim=10e3, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos=0.86667,
#                frac.g0.plus = 1,
#                frac.g0.minus = 1,
#                sim=TRUE)
# 
# 
# 
# # JNCI example, strongly differential non-compliance *but null both overall and in ever/never-pos*
# conceal.reveal(nsim=30e3, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=1,RR_neg=1, RR_pos=1,
#                frac.g0.plus = 1, frac.g0.minus = 1,
#                censor.screen.Dplus = 0.4,
#                censor.control.Dplus = 0.8,
#                censor.screen.Dminus = 0.8,
#                censor.control.Dminus = 0.4,
#                move.screen.Mplus.Dplus = 0,
#                move.screen.Mminus.Dplus = 0,
#                sim=TRUE)
# 
# # Increase P(M+) to 50%, have slightly differential non-compliance
# conceal.reveal(nsim=10e3, n=100e3, p_Mplus=0.5, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos=0.86667,
#                frac.g0.plus = 1, frac.g0.minus = 1,
#                censor.screen.Dplus = 0.1,
#                censor.control.Dplus = 0.13,
#                censor.screen.Dminus = 0.12,
#                censor.control.Dminus = 0.15,
#                move.screen.Mplus.Dplus = 0,
#                move.screen.Mminus.Dplus = 0,
#                sim=TRUE)
# 
# # % reduction in test costs for IE vs targeted, 3 screens, $3000/person and $300/test
# 1-(53e3*(2*2000+ 2*3*200))/(73e3*(2*2000 + 1*3*200))
# # % reduction in test costs for IE vs standard, 3 screens, $3000/person and $300/test
# 1-(53e3*(2*2000+ 2*3*200))/(98e3*(2*2000 + 1*3*200))
# # Find the ratio of person/test cost so that IE and targeted are equal cost
# # Finds that the ratio is 742.5/300=2.475
# uniroot(function(x) 1-(53e3*(2*x+ 2*3*200))/(73e3*(2*x + 1*3*200)), c(0,3000))
# # Above but say we only need to test half the control arm
# uniroot(function(x) 1-(53e3*(2*x+ 1.5*3*200))/(73e3*(2*x + 1*3*200)), c(0,3000))
# 
# # Make sure I can reproduce Paul's figure 2 of IE sample sizes for 90% power in the JNCI paper
# # Run this line, take the p-value for the RRpos, in this case p=1.55680431351154e-05
# # Plug it into the line below to get the sample size for 90% power for the IE
# # Ex: for P(M+)=5%, RRpos=1-(1-RR)/(1-P_Evpos)=0.875, N=55165 per arm, which looks like it matches to the figure
# conceal.reveal(nsim=1, n=2*98e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos= 1-0.1/0.8, sim=FALSE)
# 
# conceal.reveal(nsim=1, n=2*50e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos= 0.86667, sim=FALSE)
# 
# 
# 
# ##
# # Loss of signal
# # Fix n=100k, but now vary differential loss of signal from 0.5 to 1 and 0.4 to 0.9
# ##
# basecase<-conceal.reveal(nsim=1e4,sim=TRUE)
# # fracs.plus <-c(seq(0.5,1,by=0.1),1)
# # fracs.minus <- c(seq(0.5,1,by=0.1),1.1) - 0.1 # allow for a base-case of 1.0 fraction
# fracs.plus <-c(1,0.45,0.9) #c(seq(0.5,1,by=0.1),1)
# fracs.minus <- c(1,0.45,0.45) # allow for a base-case of 1.0 fraction
# out <- matrix(NA,nrow=1,ncol=ncol(basecase))
# for (i in 1:length(fracs.plus)) {
#   # out <- rbind(out,
#   #              as.matrix(conceal.reveal(
#   #                nsim=10000, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
#   #                RR=0.9,RR_neg=1, RR_pos= 0.86667,
#   #                frac.g0.plus = 0.95, frac.g0.minus = 0.5, 
#   #                frac.g0.plus.pos = fracs.plus, frac.g0.minus.pos = fracs.minus, 
#   #                sim=TRUE))
#   # )
#   out <- rbind(out,
#                as.matrix(conceal.reveal(
#                  nsim=10000, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
#                  RR=0.9,RR_neg=1, RR_pos= 0.86667,
#                  frac.g0.plus = 1, frac.g0.minus = 1, 
#                  frac.g0.plus.pos = fracs.plus, frac.g0.minus.pos = fracs.minus, 
#                  sim=TRUE))
#   )
# }
# print(out[-1,])
# print(basecase <- as.data.frame(cbind(fracs.plus,fracs.minus,
#                                       out[-1,4+c(24,25,27,28,29,30,32,33)])), digits=2, row.names=FALSE)
# print(basecase <- as.data.frame(cbind(fracs.plus,fracs.minus,
#                                       out[-1,4+c(36,37,38,40,41,42,43,46,47,65-4,66-4)])), digits=2, row.names=FALSE)
# 
# 
# 
# 
# 
# 
# 
# # censor <- matrix(c(0.00,	0.00,	0.00,	0.00, 0.00, 0.00,
# #                    0.30,	0.30,	0.30,	0.30, 0.00, 0.00,
# #                    0.30,	0.30,	0.30,	0.30, 1e-3, 1e-3,
# #                    0.20,	0.20,	0.30,	0.30, 0.00, 0.00,
# #                    0.20,	0.20,	0.30,	0.30, 1e-3, 1e-3,
# #                    0.30,	0.10,	0.30,	0.10, 0.00, 0.00,
# #                    0.30,	0.10,	0.30,	0.10, 1e-3, 1e-3,
# #                    0.10,	0.05,	0.05,	0.10, 0.00, 0.00,
# #                    0.10,	0.05,	0.05,	0.10, 1e-3, 1e-3,
# #                    0.05,	0.10,	0.10,	0.05, 0.00, 0.00,
# #                    0.05,	0.10,	0.10,	0.05, 1e-3, 1e-3),ncol=6,byrow=TRUE)
# # censor <- matrix(c(0.30,	0.30,	0.30,	0.30, 0.00, 0.00,
# #                    0.30,	0.10,	0.30,	0.10, 0.00, 0.00),ncol=6,byrow=TRUE)
# 
# censor <- matrix(c(0.00,	0.00,	0.00,	0.00, 0.00, 0.00,
#                    0.10,	0.10,	0.10,	0.10, 0.00, 0.00,
#                    0.10,	0.10,	0.12,	0.12, 0.00, 0.00,
#                    0.10,	0.10,	0.13,	0.13, 0.00, 0.00,
#                    0.20,	0.20,	0.20,	0.20, 0.00, 0.00,
#                    0.20,	0.20,	0.21,	0.21, 0.00, 0.00,
#                    0.20,	0.20,	0.22,	0.22, 0.00, 0.00,
#                    0.20,	0.20,	0.23,	0.23, 0.00, 0.00),ncol=6,byrow=TRUE)
# 
# 
# ##
# # Table 2: effect of loss-of-signal in stored bloods in control arm
# ##
# # Simulate and save the data
# basecase <- conceal.reveal(nsim=1,sim=TRUE)
# # n <- c(25e3,50e3,75e3,100e3,150e3,200e3) #seq(120e3,150e3,by=10e3)
# # RR_poss <- c(0.86667,0.8)
# # RR_negs <- c(1,1.05,0.95)
# # RR <- 0.9
# n <- 100e3 #c(10e3, 15e3, 20e3, 25e3,30e3,40e3,50e3) #seq(120e3,150e3,by=10e3)
# RR_poss <-0.86667 #c(0.8,0.86667)
# RR_negs <- 1 #c(1,1.05,0.95)
# RR <- 0.9
# censor <- matrix(c(0.50,	0.20,	0.00,	0.00,
#                    0.50,  0.20, 0.50, 0.20),ncol=4,byrow=TRUE)
# colnames(censor)<- c("screen D+","screen D-","control D+","control D-")
# out <- matrix(NA,nrow=1,ncol=ncol(basecase))
# for (l in 1:nrow(censor)) {
#   for (k in 1:length(RR_negs)) {
#     for (j in 1:length(RR_poss)) {
#       for (i in 1:length(n)) {
#         out <- rbind(out,
#                      as.matrix(conceal.reveal(nsim=1e4, n=n, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02,
#                                               RR=RR,RR_neg=RR_negs[k], RR_pos=RR_poss[j],
#                                               # frac.g0.plus = 0.95, 
#                                               # frac.g0.minus = 0.5,
#                                               # frac.g0.plus.pos = 0.9, 
#                                               # frac.g0.minus.pos = 0.8,
#                                               frac.g0.plus = 1, 
#                                               frac.g0.minus = 1,
#                                               frac.g0.plus.pos = 0.9, 
#                                               frac.g0.minus.pos = 0.8,
#                                               censor.screen.Dplus = 0,
#                                               censor.control.Dplus = 0,
#                                               censor.screen.Dminus = 0,
#                                               censor.control.Dminus = 0,
#                                               move.screen.Mplus.Dplus = 0,
#                                               move.screen.Mminus.Dplus = 0,
#                                               sim=TRUE)) )
#       }}}}
# basecase <- as.data.frame(cbind(censor, out[-1,4+c(38,44,43,48,34,35)]))
# names(basecase)<-c(colnames(censor),"RDpos","RRpos","RRneg","power RRpos","P(M+|G=0)","P(M+|G=1)")
# print(basecase, digits=3, row.names=FALSE)
# 
# 
# 
# # Print 1 simulation of 2x2 tables with dropout
# conceal.reveal(nsim=1, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02,
#                RR=0.9,RR_neg=1, RR_pos=0.86667,
#                frac.g0.plus = 0.95,
#                frac.g0.minus = 0.5,
#                frac.g0.plus.pos = 0.9,
#                frac.g0.minus.pos = 0.8,
#                censor.screen.plus = 0.5,
#                censor.control.plus = 0,
#                censor.screen.minus = 0.2,
#                censor.control.minus = 0,
#                sim=TRUE)
# 
# 
# 
# 
# 
# ##
# # Figure 2, Stat Paper, power vs RRpos:
# #
# # Calculate Z_ratios (ratio increase in the noncentrality param) varying RRpos and P(M+)
# # fix RRneg=1, RR, P(D+|G=0), and P(G=0)=0.5
# basecase <- conceal.reveal(nsim=1,sim=TRUE)
# # Note that P(M+)=0.025, you get a non-monotonic Z-ratio at RRpos=RR
# # RR <- 0.8
# # n <- 50e3 #c(10e3, 15e3, 20e3, 25e3,30e3,40e3,50e3)
# # RR_poss <- c(0.8,0.7,0.6,0.3,0)
# # p_Mpluss <- c(0.025,0.05,0.10, 0.25, 0.50, 0.75, 0.9)
# RR <- 0.9
# n <- 50e3 #c(10e3, 15e3, 20e3, 25e3,30e3,40e3,50e3)
# RR_poss <- seq(0.9,0.55,by=-0.05) #c(0.9,0.8,0.7,0.5)
# p_Mpluss <- c(0.025,0.05,0.50,0.90)
# RR_negs <- c(1,1.05,0.95)
# out <- matrix(NA,nrow=1,ncol=ncol(basecase))
# for (l in 1:length(RR_negs)) {
#   for (k in 1:length(p_Mpluss)) {
#     for (j in 1:length(RR_poss)) {
#       for (i in 1:length(n)) {
#         out <- rbind(out,
#                      as.matrix(conceal.reveal(nsim=10000, n=n, p_Mplus=p_Mpluss[k], p_g0=0.5, pD_g0=0.02, 
#                                               RR=RR,RR_neg=RR_negs[l], RR_pos=RR_poss[j],
#                                               sim=TRUE))
#         )
#       }
#     }
#   }
# }
# out.fig1 <- as.data.frame(cbind(n,out[-1,c(2,10:11,13:14,73)],
#                                 rep(p_Mpluss,each=length(n)*length(RR_poss)),
#                                 rep(RR_negs,each=length(n)*length(RR_poss)*length(p_Mpluss))
# ))
# names(out.fig1)<-c("n","RR_pos","Zratio","Zratio_sim","IEpower","Stdpower","Targetedpower","everpos","RR_neg")
# # Replace everpos=0.9 with Targeted since they are essentially identical
# out.fig1[out.fig1$everpos==0.9,"IEpower"] <- out.fig1[out.fig1$everpos==0.9,"Targetedpower"]
# out.fig1[out.fig1$everpos==0.9,"everpos"] <- "Targeted"
# print(out.fig1, row.names=FALSE)
# 
# # Helper function to make plots
# plot.power.fig <- function(out.fig1) {
#   
#   # Don't think I need this
#   # plot_data <- reshape2::melt(out.fig1,id.var=c("everpos","RR_pos","RR_neg"))
#   # colnames(plot_data)[4:5] <- c("analysis","power")
#   # plot_data$power[plot_data$power==1] <- 0.9999 # cannot plot a power of 1 on probit scale
#   # print(plot_data)
#   
#   plot_data <- out.fig1
#   
#   # The color-blind palette with black as #3, move yellow to end:
#   cbbPalette <- c("#0072B2","#D55E00","#000000","#009E73","#56B4E9","#E69F00","#CC79A7", 
#                   "#F0E442")
#   
#   plotpowers <- 
#     ggplot(plot_data, aes(x=RR_pos,y=IEpower,color=as.factor(everpos))) +
#     scale_color_manual(values=cbbPalette) +
#     #scale_colour_manual(values=c("red","black","blue")) +
#     geom_point(size=3)+
#     geom_line()+#aes(lty=as.factor(everpos))) +
#     # scale_linetype_manual(values=c("dotted", "solid"))+ # rev so solid is the usual analysis
#     # scale_color_discrete(na.translate = F) +
#     scale_x_continuous(breaks=unique(out.fig1$RR_pos)) +
#     # scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1)) +
#     scale_y_continuous(trans = "probit", limits = c(0.1,0.9999), 
#                        breaks=c(seq(0.2,0.9,0.1),0.95,0.98,0.99,0.999)) +
#     geom_hline(yintercept = c(0.8,0.9),lty=3) +
#     geom_hline(yintercept = c(0.36),lty=1) + # std analysis has power 36%+
#     labs(x="RR_pos",y="power",color="ever positivity") +
#     #theme with white background
#     theme_bw() +
#     #eliminates background, gridlines, #and chart border
#     theme(
#       plot.background = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank()
#       # panel.border = element_blank()
#     ) +
#     NULL
# }
# 
# # plot the data
# plotpowers.rrneg.1 <- plot.power.fig(out.fig1[out.fig1$RR_neg==1,])
# plotpowers.rrneg.2 <- plot.power.fig(out.fig1[out.fig1$RR_neg==1.05,])
# plotpowers.rrneg.3 <- plot.power.fig(out.fig1[out.fig1$RR_neg==0.95,])
# 
# # Extract legend from one of the plots - they're all the same
# legend_b <- cowplot::get_legend(
#   plotpowers.rrneg.1 + 
#     guides(color = guide_legend(nrow = 1, byrow=TRUE)) +
#     theme(legend.position = "bottom")
# )
# 
# # Plot single analysis and subgroup analysis
# prow <- cowplot::plot_grid(plotpowers.rrneg.1 + theme(legend.position="none"), 
#                                  plotpowers.rrneg.2 + theme(legend.position="none"),
#                                  plotpowers.rrneg.3 + theme(legend.position="none"),
#                                  nrow=1, ncol=3, rel_heights = c(1, 1,1),
#                                  labels = c("RR_neg=1 (RR=0.9)","RR_neg=1.05 (RR=0.9)","RR_neg=0.95 (RR=0.9)"))
# # Include legen
# final_plot <- cowplot::plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))
# 
# # show the plot 
# final_plot
# 
# # Save the plot
# ggsave(filename="Fig-2-stat-paper.pdf", plot=final_plot, device="pdf",
#        units="in", width=15, height=8,  dpi=500)
# 
# 
# 
# 
# # Try examples, calculate Zratio
# conceal.reveal(nsim=1, n=50e3, p_Mplus=0.025, p_g0=0.5, pD_g0=0.02,
#                RR=0.9,RR_neg=0.95, RR_pos= 0.85,
#                frac.g0.plus = 1, frac.g0.minus = 1,
#                frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,sim=TRUE)
# 
# 
# 
# 
# 
# # Try examples where control arm gets MCD, so that RR and RRpos are higher
# conceal.reveal(nsim=1, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, RR=0.9,RR_neg=1, RR_pos=0.86667,
#                frac.g0.plus = 1,frac.g0.minus = 1,sim=FALSE) # std JNCI example
# 
# conceal.reveal(nsim=1, n=100e3, p_Mplus=0.03, p_g0=0.5, pD_g0=0.02, RR=0.95,RR_neg=1, RR_pos=0.5,
#                frac.g0.plus = 1,frac.g0.minus = 1,sim=FALSE) #
# 
# 
# 
# ##
# # Calculate 2x2 tables for tiger team presentation
# ##
# temp <- conceal.reveal(nsim=1e5, n=100e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
#                RR=0.9,RR_neg=1, RR_pos= 0.86667,
#                frac.g0.plus = 1, frac.g0.minus = 1, 
#                frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,sim=TRUE)
# print(temp)
# 2*(1-pnorm(1*(qnorm(1-0.05/2)+qnorm(temp$power.RDall))))
# 2*(1-pnorm(1*(qnorm(1-0.05/2)+qnorm(temp$power.RDpos))))
# # binomial SD for RD_all
# sd_RD_all <- 0.002/(qnorm(1-0.05/2)+qnorm(temp$power.RDall))
# # binomial SD for RD_IE
# sd_RD_all*(1/temp$Z_ratio)
# 
# # For RR_pos=0.86667 and 100k per arm
# Mpos <- matrix( c(1300,1500,3700,3500),nrow=2,byrow=TRUE,
#                 dimnames=list(c("D+","D-"),c("screen","control")))
# Mneg <- matrix( c(500,500,94500,94500),nrow=2,byrow=TRUE,
#                 dimnames=list(c("D+","D-"),c("screen","control")))
# 
# 
# conceal.reveal(nsim=1, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02, 
#                RR=0.9,RR_neg=1, RR_pos= 0.867,
#                frac.g0.plus = 1, frac.g0.minus = 1, 
#                frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,sim=TRUE)
# # # For RR_pos=0.8
# # Mpos <- matrix( c(800,1000,4200,4000),nrow=2,byrow=TRUE,
# #                 dimnames=list(c("D+","D-"),c("screen","control")))
# # Mneg <- matrix( c(1000,1000,94000,94000),nrow=2,byrow=TRUE,
# #                 dimnames=list(c("D+","D-"),c("screen","control")))
# 
# # Example with false-reassurance RR_neg=1.21
# conceal.reveal(nsim=1, n=200e3, p_Mplus=0.05, p_g0=0.5, pD_g0=0.02,
#                RR=0.9,RR_neg=1.21, RR_pos= 0.866667,
#                frac.g0.plus = 1, frac.g0.minus = 1,
#                frac.g0.plus.pos = 1, frac.g0.minus.pos = 1,sim=TRUE)
# Mpos <- matrix( c(1565,1806,3435,3194),nrow=2,byrow=TRUE,
#                 dimnames=list(c("D+","D-"),c("screen","control")))
# Mneg <- matrix( c(235,194,94765,94806),nrow=2,byrow=TRUE,
#                 dimnames=list(c("D+","D-"),c("screen","control")))
# 
# # Do 1/4 the sample size for the example
# # Mpos <- Mpos/4
# # Mneg <- Mneg/4
# print(Mpos)
# print(Mneg)
# print(all<-Mneg+Mpos)
# print(prop.test(t(all),correct=F))
# print(fisher.test(t(all)))
# prop.test(t(Mpos),correct=F)
# fisher.test(t(Mpos))
# prop.test(t(Mneg),correct=F)
# fisher.test(t(Mneg))
# 
# 
# 
