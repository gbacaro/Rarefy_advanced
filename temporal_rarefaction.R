# ==============================================================================
# temporal_rarefaction.R
# ==============================================================================
# R function to calculate temporal rarefaction profiles.
#
# Purpose
# -------
# temporal_rarefaction() calculates sample-based rarefaction curves independently
# for each time step, then measures temporal variability among those curves at
# each number of accumulated sampling units, m. The function is designed for
# replicated ecological surveys, for example permanent plots revisited through
# time.
#
# The output allows users to ask whether temporal variability is expressed mainly
# at the alpha-diversity scale (m = 1), at the gamma-diversity scale (m = N), or
# along the alpha-gamma accumulation profile between them.
#
# Input data organisation
# -----------------------
# Community data can be provided in either wide or long format.
#
# 1. Wide format, format = "wide"
#    Each row is one sampling unit observed at one time step. Species are columns.
#    Time, unit and optional group variables can either be columns in comm or
#    external vectors supplied through time, unit_id and group.
#
#    Example:
#      plot_id year treatment sp1 sp2 sp3 ... spK
#      P01     2001 Grazed     4   0   1      0
#      P02     2001 Grazed     0   3   2      1
#
# 2. Long format, format = "long"
#    Each row is one species observation in one sampling unit at one time step.
#    Required columns must be identified through time_col, unit_col, species_col
#    and abundance_col. group_col is optional.
#
#    Example:
#      plot_id year treatment species abundance
#      P01     2001 Grazed    sp1     4
#      P01     2001 Grazed    sp3     1
#
# Optional trait data
# -------------------
# If traits are supplied, rows must correspond to species and columns to traits.
# Species names can be row names, or can be stored in a column identified by
# trait_species_col. Traits should be numeric. By default they are standardised
# to mean 0 and standard deviation 1 before distances, PCA scores and CWM-based
# descriptors are calculated.
#
# Metrics implemented
# -------------------
# Taxonomic metrics:
#   richness  = number of species with positive abundance
#   shannon   = Shannon entropy, -sum(p log p)
#   simpson   = Gini-Simpson diversity, 1 - sum(p^2)
#
# Functional metrics, when traits are provided:
#   RaoQ      = Rao's quadratic entropy using Euclidean trait distances
#   FDis      = abundance-weighted distance from the community centroid
#   FRic2D    = convex-hull area in the first two axes of the trait ordination
#   CWMdist   = Euclidean distance between the rarefied CWM vector and a baseline
#
# Users can also pass custom metrics through custom_metrics. Each custom metric
# must be a named function with arguments:
#   function(abund, traits = NULL, trait_dist = NULL, trait_scores = NULL)
# where abund is a named numeric vector of pooled species abundances.
#
# Main outputs
# ------------
# A list with:
#   rarefaction_profiles   time-specific rarefaction curves
#   temporal_variability   temporal variance, SD and CV at each m
#   auc_summary            whole-profile summaries
#   alpha_gamma_summary    values at m = 1 and m = N for each time step
#   group_contrast         optional signed differences between two groups
#   settings               main settings used in the run
#
# Notes
# -----
# - The function evaluates all combinations of m sampling units when feasible.
#   Otherwise it uses random resampling without replacement.
# - For species richness, the expected value can optionally be replaced by the
#   exact sample-based richness expectation following the Shinozaki/Colwell
#   incidence formulation. Intervals still come from combinations/resampling.
# - This is a clear stand-alone implementation intended for GitHub distribution
#   and later integration into the Rarefy framework.
# ============================================================================== 

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

.tr_stop <- function(...) stop(paste0(...), call. = FALSE)

.tr_check_numeric_matrix <- function(x, object_name = "matrix") {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  x
}

.tr_n_choose_k <- function(n, k) {
  if (k < 0 || k > n) return(0)
  out <- suppressWarnings(choose(n, k))
  if (!is.finite(out)) Inf else out
}

.tr_make_subsets <- function(n, m, n_perm = 4999, use_exact = TRUE,
                             exact_threshold = 5000) {
  if (n < 1) .tr_stop("n must be >= 1.")
  if (m < 1 || m > n) .tr_stop("m must be between 1 and n.")

  n_comb <- .tr_n_choose_k(n, m)
  if (isTRUE(use_exact) && is.finite(n_comb) && n_comb <= exact_threshold) {
    out <- utils::combn(n, m, simplify = FALSE)
    attr(out, "subset_mode") <- "exact"
    return(out)
  }

  n_draws <- max(1L, as.integer(n_perm))
  out <- replicate(n_draws, sort(sample.int(n, m, replace = FALSE)), simplify = FALSE)
  attr(out, "subset_mode") <- "random"
  out
}

.tr_expected_richness_shinozaki <- function(mat, m) {
  # Exact expected sample-based richness for m sampling units.
  # A species contributes 1 minus the probability that all m selected units are
  # drawn from the N - incidence_i units where the species is absent.
  mat <- as.matrix(mat)
  N <- nrow(mat)
  incidence <- colSums(mat > 0, na.rm = TRUE)
  if (m < 1 || m > N) return(NA_real_)

  log_denom <- lchoose(N, m)
  prob_absent <- numeric(length(incidence))
  can_be_absent <- (N - incidence) >= m
  prob_absent[can_be_absent] <- exp(lchoose(N - incidence[can_be_absent], m) - log_denom)
  prob_absent[!can_be_absent] <- 0
  sum(1 - prob_absent)
}

.tr_taxonomic_metric <- function(abund, metric = "richness") {
  abund <- abund[is.finite(abund) & abund > 0]
  if (length(abund) == 0) return(0)
  p <- abund / sum(abund)

  switch(
    metric,
    richness = length(abund),
    shannon  = -sum(p * log(p)),
    simpson  = 1 - sum(p^2),
    .tr_stop("Unknown taxonomic metric: ", metric)
  )
}

.tr_weighted_cwm <- function(abund, trait_scaled) {
  if (is.null(names(abund))) .tr_stop("abund must be a named vector.")
  sp <- intersect(names(abund)[is.finite(abund) & abund > 0], rownames(trait_scaled))
  if (length(sp) == 0) return(rep(NA_real_, ncol(trait_scaled)))
  p <- abund[sp] / sum(abund[sp])
  out <- as.numeric(colSums(sweep(trait_scaled[sp, , drop = FALSE], 1, p, `*`)))
  names(out) <- colnames(trait_scaled)
  out
}

.tr_poly_area_2d <- function(xy) {
  xy <- as.matrix(xy)
  if (nrow(xy) < 3 || ncol(xy) < 2) return(0)
  xy <- xy[stats::complete.cases(xy), 1:2, drop = FALSE]
  if (nrow(xy) < 3) return(0)
  if (qr(xy)$rank < 2) return(0)
  hull <- grDevices::chull(xy[, 1], xy[, 2])
  pts <- xy[c(hull, hull[1]), , drop = FALSE]
  abs(sum(pts[-1, 1] * pts[-nrow(pts), 2] - pts[-nrow(pts), 1] * pts[-1, 2]) / 2)
}

.tr_prepare_traits <- function(traits, trait_cols = NULL, trait_species_col = NULL,
                               species_names, standardize_traits = TRUE) {
  if (is.null(traits)) return(NULL)

  traits <- as.data.frame(traits, stringsAsFactors = FALSE)
  if (!is.null(trait_species_col)) {
    if (!trait_species_col %in% names(traits)) {
      .tr_stop("trait_species_col not found in traits: ", trait_species_col)
    }
    rownames(traits) <- as.character(traits[[trait_species_col]])
    traits[[trait_species_col]] <- NULL
  }
  if (is.null(rownames(traits)) || any(rownames(traits) == "")) {
    .tr_stop("traits must have species row names or a valid trait_species_col.")
  }

  if (is.null(trait_cols)) {
    trait_cols <- names(traits)[vapply(traits, is.numeric, logical(1))]
  }
  if (length(trait_cols) == 0) .tr_stop("No numeric trait columns were found.")
  missing_trait_cols <- setdiff(trait_cols, names(traits))
  if (length(missing_trait_cols) > 0) {
    .tr_stop("These trait_cols are not present in traits: ", paste(missing_trait_cols, collapse = ", "))
  }

  X <- traits[, trait_cols, drop = FALSE]
  for (j in seq_along(X)) X[[j]] <- suppressWarnings(as.numeric(X[[j]]))
  X <- as.matrix(X)
  rownames(X) <- rownames(traits)
  storage.mode(X) <- "numeric"

  keep <- rownames(X) %in% species_names & stats::complete.cases(X)
  X <- X[keep, , drop = FALSE]
  if (nrow(X) < 2) {
    .tr_stop("Fewer than two species have complete trait data matching the community matrix.")
  }

  if (isTRUE(standardize_traits)) {
    mu <- colMeans(X, na.rm = TRUE)
    sig <- apply(X, 2, stats::sd, na.rm = TRUE)
    sig[!is.finite(sig) | sig == 0] <- 1
    X <- sweep(sweep(X, 2, mu, "-"), 2, sig, "/")
  }

  D <- as.matrix(stats::dist(X, method = "euclidean"))
  rownames(D) <- colnames(D) <- rownames(X)

  pca <- stats::prcomp(X, center = FALSE, scale. = FALSE)
  scores_all <- stats::predict(pca, newdata = X)
  if (is.null(dim(scores_all))) scores_all <- matrix(scores_all, ncol = 1)
  if (ncol(scores_all) < 2) {
    scores <- cbind(scores_all[, 1], 0)
  } else {
    scores <- scores_all[, 1:2, drop = FALSE]
  }
  colnames(scores) <- c("PC1", "PC2")
  rownames(scores) <- rownames(X)

  list(
    trait_scaled = X,
    trait_dist = D,
    trait_scores = scores,
    trait_cols = trait_cols,
    standardize_traits = standardize_traits
  )
}

.tr_make_community_from_long <- function(comm, time_col, unit_col, species_col,
                                         abundance_col, group_col = NULL) {
  needed <- c(time_col, unit_col, species_col, abundance_col)
  if (!is.null(group_col)) needed <- c(needed, group_col)
  missing <- setdiff(needed, names(comm))
  if (length(missing) > 0) .tr_stop("Missing columns in long comm: ", paste(missing, collapse = ", "))

  comm[[abundance_col]] <- suppressWarnings(as.numeric(comm[[abundance_col]]))
  comm[[abundance_col]][!is.finite(comm[[abundance_col]])] <- 0

  meta_cols <- c(unit_col, time_col, group_col)
  meta_cols <- meta_cols[!is.null(meta_cols)]
  meta <- unique(comm[, meta_cols, drop = FALSE])
  names(meta)[names(meta) == unit_col] <- ".unit"
  names(meta)[names(meta) == time_col] <- ".time"
  if (!is.null(group_col)) names(meta)[names(meta) == group_col] <- ".group"
  if (is.null(group_col)) meta$.group <- "All"

  # Aggregate duplicated records before reshaping.
  agg_cols <- c(unit_col, time_col, group_col, species_col)
  agg_cols <- agg_cols[!is.null(agg_cols)]
  agg <- stats::aggregate(
    comm[[abundance_col]],
    by = comm[, agg_cols, drop = FALSE],
    FUN = sum,
    na.rm = TRUE
  )
  names(agg)[ncol(agg)] <- ".abund"

  key <- paste(agg[[unit_col]], agg[[time_col]], if (!is.null(group_col)) agg[[group_col]] else "All", sep = "___")
  species <- as.character(agg[[species_col]])
  mat <- xtabs(.abund ~ key + species, data = data.frame(key = key, species = species, .abund = agg$.abund))
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  meta$.key <- paste(meta$.unit, meta$.time, meta$.group, sep = "___")
  meta <- meta[match(rownames(mat), meta$.key), c(".unit", ".time", ".group"), drop = FALSE]
  rownames(meta) <- rownames(mat)

  list(mat = mat, meta = meta)
}

.tr_make_community_from_wide <- function(comm, species_cols = NULL,
                                         time = NULL, unit_id = NULL, group = NULL,
                                         time_col = NULL, unit_col = NULL, group_col = NULL) {
  comm <- as.data.frame(comm, stringsAsFactors = FALSE)

  if (!is.null(time_col)) {
    if (!time_col %in% names(comm)) .tr_stop("time_col not found in comm: ", time_col)
    time <- comm[[time_col]]
  }
  if (!is.null(unit_col)) {
    if (!unit_col %in% names(comm)) .tr_stop("unit_col not found in comm: ", unit_col)
    unit_id <- comm[[unit_col]]
  }
  if (!is.null(group_col)) {
    if (!group_col %in% names(comm)) .tr_stop("group_col not found in comm: ", group_col)
    group <- comm[[group_col]]
  }

  n <- nrow(comm)
  if (is.null(time)) .tr_stop("time or time_col must be supplied.")
  if (length(time) != n) .tr_stop("time must have one value per row of comm.")
  if (is.null(unit_id)) unit_id <- seq_len(n)
  if (length(unit_id) != n) .tr_stop("unit_id must have one value per row of comm.")
  if (is.null(group)) group <- rep("All", n)
  if (length(group) != n) .tr_stop("group must have one value per row of comm.")

  if (is.null(species_cols)) {
    meta_names <- c(time_col, unit_col, group_col)
    meta_names <- meta_names[!is.null(meta_names)]
    species_cols <- setdiff(names(comm), meta_names)
    numeric_cols <- vapply(comm[, species_cols, drop = FALSE], function(z) {
      is.numeric(z) || all(is.na(suppressWarnings(as.numeric(as.character(z)))) == is.na(z))
    }, logical(1))
    species_cols <- species_cols[numeric_cols]
  }
  if (length(species_cols) == 0) .tr_stop("No species columns were identified.")
  missing_species_cols <- setdiff(species_cols, names(comm))
  if (length(missing_species_cols) > 0) {
    .tr_stop("These species_cols are not present in comm: ", paste(missing_species_cols, collapse = ", "))
  }

  mat <- comm[, species_cols, drop = FALSE]
  for (j in seq_along(mat)) mat[[j]] <- suppressWarnings(as.numeric(as.character(mat[[j]])))
  mat <- .tr_check_numeric_matrix(mat, "comm")
  rownames(mat) <- make.unique(as.character(unit_id))

  meta <- data.frame(
    .unit = as.character(unit_id),
    .time = time,
    .group = as.character(group),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- rownames(mat)
  list(mat = mat, meta = meta)
}

.tr_baseline_cwm <- function(rows, mat, trait_scaled) {
  abund <- colSums(mat[rows, , drop = FALSE], na.rm = TRUE)
  names(abund) <- colnames(mat)
  .tr_weighted_cwm(abund, trait_scaled)
}

.tr_functional_metrics <- function(abund, trait_object, min_trait_coverage = 0.50,
                                   baseline_cwm = NULL) {
  total_abund <- sum(abund[is.finite(abund) & abund > 0], na.rm = TRUE)
  if (total_abund <= 0) {
    return(c(RaoQ = NA_real_, FDis = NA_real_, FRic2D = NA_real_, CWMdist = NA_real_, trait_coverage = 0))
  }

  trait_scaled <- trait_object$trait_scaled
  trait_dist <- trait_object$trait_dist
  trait_scores <- trait_object$trait_scores

  present <- names(abund)[is.finite(abund) & abund > 0]
  matched <- Reduce(intersect, list(present, rownames(trait_scaled), rownames(trait_dist), rownames(trait_scores)))

  covered_abund <- sum(abund[matched], na.rm = TRUE)
  coverage <- covered_abund / total_abund

  if (length(matched) == 0 || !is.finite(coverage) || coverage < min_trait_coverage) {
    return(c(RaoQ = NA_real_, FDis = NA_real_, FRic2D = NA_real_, CWMdist = NA_real_, trait_coverage = coverage))
  }

  p <- abund[matched] / sum(abund[matched])

  D <- trait_dist[matched, matched, drop = FALSE]
  RaoQ <- as.numeric(t(p) %*% D %*% p)

  scores <- trait_scores[matched, , drop = FALSE]
  centroid <- colSums(sweep(scores, 1, p, `*`))
  dist_centroid <- sqrt(rowSums(sweep(scores, 2, centroid, "-")^2))
  FDis <- sum(p * dist_centroid)

  FRic2D <- .tr_poly_area_2d(scores)

  cwm <- .tr_weighted_cwm(abund, trait_scaled)
  if (is.null(baseline_cwm) || all(is.na(cwm)) || all(is.na(baseline_cwm))) {
    CWMdist <- NA_real_
  } else {
    CWMdist <- sqrt(sum((cwm - baseline_cwm)^2, na.rm = TRUE))
  }

  c(RaoQ = RaoQ, FDis = FDis, FRic2D = FRic2D, CWMdist = CWMdist, trait_coverage = coverage)
}

.tr_metric_values_for_subset <- function(abund, taxonomic_metrics,
                                         functional_metrics, trait_object,
                                         baseline_cwm, min_trait_coverage,
                                         custom_metrics) {
  out <- numeric(0)

  if (length(taxonomic_metrics) > 0) {
    tax <- vapply(taxonomic_metrics, function(metric) .tr_taxonomic_metric(abund, metric), numeric(1))
    out <- c(out, tax)
  }

  if (!is.null(trait_object) && length(functional_metrics) > 0) {
    fun_all <- .tr_functional_metrics(
      abund = abund,
      trait_object = trait_object,
      min_trait_coverage = min_trait_coverage,
      baseline_cwm = baseline_cwm
    )
    keep <- intersect(functional_metrics, names(fun_all))
    out <- c(out, fun_all[keep], trait_coverage = fun_all["trait_coverage"])
  }

  if (!is.null(custom_metrics) && length(custom_metrics) > 0) {
    for (nm in names(custom_metrics)) {
      fn <- custom_metrics[[nm]]
      val <- fn(
        abund = abund,
        traits = if (!is.null(trait_object)) trait_object$trait_scaled else NULL,
        trait_dist = if (!is.null(trait_object)) trait_object$trait_dist else NULL,
        trait_scores = if (!is.null(trait_object)) trait_object$trait_scores else NULL
      )
      out[nm] <- as.numeric(val)[1]
    }
  }

  out
}

.tr_auc <- function(y, x) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) == 0) return(NA_real_)
  if (length(x) == 1) return(y[1])
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  # Trapezoidal area. If x is 1:N, this is scale-aware. For exact compatibility
  # with simple sums used in some scripts, set auc_method = "sum" in the main function.
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}

# ------------------------------------------------------------------------------
# Main function
# ------------------------------------------------------------------------------

temporal_rarefaction <- function(
    comm,
    format = c("wide", "long"),
    time = NULL,
    unit_id = NULL,
    group = NULL,
    species_cols = NULL,
    time_col = NULL,
    unit_col = NULL,
    group_col = NULL,
    species_col = NULL,
    abundance_col = NULL,
    traits = NULL,
    trait_cols = NULL,
    trait_species_col = NULL,
    taxonomic_metrics = c("richness", "shannon"),
    functional_metrics = NULL,
    custom_metrics = NULL,
    min_trait_coverage = 0.50,
    standardize_traits = TRUE,
    cwm_baseline = c("first_time_by_group", "first_time_overall", "none"),
    n_perm = 4999,
    use_exact = TRUE,
    exact_threshold = 5000,
    richness_exact_shinozaki = TRUE,
    interval_probs = c(0.025, 0.975),
    auc_method = c("sum", "trapezoid"),
    contrast_pair = NULL,
    seed = NULL,
    verbose = TRUE
) {
  format <- match.arg(format)
  cwm_baseline <- match.arg(cwm_baseline)
  auc_method <- match.arg(auc_method)

  if (!is.null(seed)) set.seed(seed)

  if (format == "long") {
    if (is.null(time_col) || is.null(unit_col) || is.null(species_col) || is.null(abundance_col)) {
      .tr_stop("For format = 'long', time_col, unit_col, species_col and abundance_col must be supplied.")
    }
    prepared <- .tr_make_community_from_long(
      comm = comm,
      time_col = time_col,
      unit_col = unit_col,
      species_col = species_col,
      abundance_col = abundance_col,
      group_col = group_col
    )
  } else {
    prepared <- .tr_make_community_from_wide(
      comm = comm,
      species_cols = species_cols,
      time = time,
      unit_id = unit_id,
      group = group,
      time_col = time_col,
      unit_col = unit_col,
      group_col = group_col
    )
  }

  mat <- prepared$mat
  meta <- prepared$meta

  if (nrow(mat) != nrow(meta)) .tr_stop("Internal error: mat and meta have different row numbers.")
  if (anyDuplicated(paste(meta$.group, meta$.time, meta$.unit, sep = "___")) > 0) {
    warning("Some group-time-unit combinations are duplicated. They will be treated as separate rows unless aggregated before calling the function.")
  }

  taxonomic_metrics <- match.arg(taxonomic_metrics, choices = c("richness", "shannon", "simpson"), several.ok = TRUE)

  trait_object <- .tr_prepare_traits(
    traits = traits,
    trait_cols = trait_cols,
    trait_species_col = trait_species_col,
    species_names = colnames(mat),
    standardize_traits = standardize_traits
  )

  if (is.null(functional_metrics)) {
    functional_metrics <- character(0)
  } else {
    functional_metrics <- match.arg(functional_metrics,
                                    choices = c("RaoQ", "FDis", "FRic2D", "CWMdist"),
                                    several.ok = TRUE)
    if (is.null(trait_object)) .tr_stop("functional_metrics require a valid traits table.")
  }

  if (!is.null(custom_metrics)) {
    if (is.null(names(custom_metrics)) || any(names(custom_metrics) == "")) {
      .tr_stop("custom_metrics must be a named list of functions.")
    }
    if (!all(vapply(custom_metrics, is.function, logical(1)))) {
      .tr_stop("All elements of custom_metrics must be functions.")
    }
  }

  # Baseline CWM vectors for CWMdist.
  baseline_by_group <- list()
  if (!is.null(trait_object) && "CWMdist" %in% functional_metrics && cwm_baseline != "none") {
    if (cwm_baseline == "first_time_overall") {
      first_time <- sort(unique(meta$.time))[1]
      rows <- which(meta$.time == first_time)
      baseline_by_group[[".overall"]] <- .tr_baseline_cwm(rows, mat, trait_object$trait_scaled)
    } else {
      for (g in sort(unique(meta$.group))) {
        rows_g <- which(meta$.group == g)
        first_time <- sort(unique(meta$.time[rows_g]))[1]
        rows <- rows_g[meta$.time[rows_g] == first_time]
        baseline_by_group[[g]] <- .tr_baseline_cwm(rows, mat, trait_object$trait_scaled)
      }
    }
  }

  # Run rarefaction for each group x time.
  gt <- unique(data.frame(.group = meta$.group, .time = meta$.time, stringsAsFactors = FALSE))
  gt <- gt[order(gt$.group, gt$.time), , drop = FALSE]

  rare_list <- vector("list", nrow(gt))

  for (i in seq_len(nrow(gt))) {
    g <- gt$.group[i]
    tt <- gt$.time[i]
    rows <- which(meta$.group == g & meta$.time == tt)
    N <- length(rows)
    if (N < 1) next

    if (isTRUE(verbose)) message("Rarefying group = ", g, ", time = ", tt, ", N = ", N)

    baseline_cwm <- NULL
    if (!is.null(trait_object) && "CWMdist" %in% functional_metrics && cwm_baseline != "none") {
      baseline_cwm <- if (cwm_baseline == "first_time_overall") baseline_by_group[[".overall"]] else baseline_by_group[[g]]
    }

    rows_out <- list()
    out_counter <- 1L

    for (m in seq_len(N)) {
      subsets <- .tr_make_subsets(N, m, n_perm = n_perm, use_exact = use_exact, exact_threshold = exact_threshold)
      subset_mode <- attr(subsets, "subset_mode")

      metric_names <- c(taxonomic_metrics, functional_metrics)
      if (length(functional_metrics) > 0) metric_names <- c(metric_names, "trait_coverage")
      if (!is.null(custom_metrics)) metric_names <- c(metric_names, names(custom_metrics))
      metric_names <- unique(metric_names)

      vals <- matrix(NA_real_, nrow = length(subsets), ncol = length(metric_names))
      colnames(vals) <- metric_names

      for (j in seq_along(subsets)) {
        idx_local <- subsets[[j]]
        idx_global <- rows[idx_local]
        abund <- colSums(mat[idx_global, , drop = FALSE], na.rm = TRUE)
        names(abund) <- colnames(mat)

        tmp <- .tr_metric_values_for_subset(
          abund = abund,
          taxonomic_metrics = taxonomic_metrics,
          functional_metrics = functional_metrics,
          trait_object = trait_object,
          baseline_cwm = baseline_cwm,
          min_trait_coverage = min_trait_coverage,
          custom_metrics = custom_metrics
        )
        keep <- intersect(names(tmp), colnames(vals))
        vals[j, keep] <- tmp[keep]
      }

      for (metric in colnames(vals)) {
        value <- vals[, metric]
        mean_value <- mean(value, na.rm = TRUE)

        if (isTRUE(richness_exact_shinozaki) && metric == "richness") {
          mean_value <- .tr_expected_richness_shinozaki(mat[rows, , drop = FALSE], m)
        }

        rows_out[[out_counter]] <- data.frame(
          group = g,
          time = tt,
          N = N,
          m = m,
          metric = metric,
          n_subsets = length(subsets),
          subset_mode = subset_mode,
          n_valid = sum(is.finite(value)),
          mean = mean_value,
          sd = if (sum(is.finite(value)) > 1) stats::sd(value, na.rm = TRUE) else NA_real_,
          lwr = suppressWarnings(stats::quantile(value, probs = interval_probs[1], na.rm = TRUE, names = FALSE)),
          upr = suppressWarnings(stats::quantile(value, probs = interval_probs[2], na.rm = TRUE, names = FALSE)),
          stringsAsFactors = FALSE
        )
        out_counter <- out_counter + 1L
      }
    }

    rare_list[[i]] <- do.call(rbind, rows_out)
  }

  rarefaction_profiles <- do.call(rbind, rare_list)
  rownames(rarefaction_profiles) <- NULL
  rarefaction_profiles$time <- as.character(rarefaction_profiles$time)

  # Temporal variability profiles.
  rp_no_cov <- rarefaction_profiles[rarefaction_profiles$metric != "trait_coverage", , drop = FALSE]
  split_tv <- split(rp_no_cov, paste(rp_no_cov$group, rp_no_cov$metric, rp_no_cov$m, sep = "___"))
  temporal_variability <- do.call(rbind, lapply(split_tv, function(df) {
    y <- df$mean
    n_time <- sum(is.finite(y))
    mean_over_time <- mean(y, na.rm = TRUE)
    temporal_sd <- if (n_time > 1) stats::sd(y, na.rm = TRUE) else NA_real_
    temporal_variance <- if (n_time > 1) stats::var(y, na.rm = TRUE) else NA_real_
    temporal_cv <- if (is.finite(mean_over_time) && mean_over_time != 0) temporal_sd / abs(mean_over_time) else NA_real_
    data.frame(
      group = df$group[1],
      metric = df$metric[1],
      m = df$m[1],
      n_time = n_time,
      mean_over_time = mean_over_time,
      temporal_variance = temporal_variance,
      temporal_sd = temporal_sd,
      temporal_cv = temporal_cv,
      temporal_lwr = suppressWarnings(stats::quantile(y, probs = interval_probs[1], na.rm = TRUE, names = FALSE)),
      temporal_upr = suppressWarnings(stats::quantile(y, probs = interval_probs[2], na.rm = TRUE, names = FALSE)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(temporal_variability) <- NULL
  temporal_variability <- temporal_variability[order(temporal_variability$group, temporal_variability$metric, temporal_variability$m), ]

  # AUC and alpha-gamma summaries.
  split_auc <- split(temporal_variability, paste(temporal_variability$group, temporal_variability$metric, sep = "___"))
  auc_summary <- do.call(rbind, lapply(split_auc, function(df) {
    df <- df[order(df$m), ]
    max_m <- max(df$m, na.rm = TRUE)
    if (auc_method == "sum") {
      auc_var <- sum(df$temporal_variance, na.rm = TRUE)
      auc_sd  <- sum(df$temporal_sd, na.rm = TRUE)
      auc_cv  <- sum(df$temporal_cv, na.rm = TRUE)
    } else {
      auc_var <- .tr_auc(df$temporal_variance, df$m)
      auc_sd  <- .tr_auc(df$temporal_sd, df$m)
      auc_cv  <- .tr_auc(df$temporal_cv, df$m)
    }
    data.frame(
      group = df$group[1],
      metric = df$metric[1],
      max_m = max_m,
      auc_temporal_variance = auc_var,
      auc_temporal_sd = auc_sd,
      auc_temporal_cv = auc_cv,
      alpha_scale_cv = df$temporal_cv[df$m == min(df$m, na.rm = TRUE)][1],
      gamma_scale_cv = df$temporal_cv[df$m == max_m][1],
      alpha_scale_sd = df$temporal_sd[df$m == min(df$m, na.rm = TRUE)][1],
      gamma_scale_sd = df$temporal_sd[df$m == max_m][1],
      stringsAsFactors = FALSE
    )
  }))
  rownames(auc_summary) <- NULL

  split_ag <- split(rp_no_cov, paste(rp_no_cov$group, rp_no_cov$time, rp_no_cov$metric, sep = "___"))
  alpha_gamma_summary <- do.call(rbind, lapply(split_ag, function(df) {
    df <- df[order(df$m), ]
    min_m <- min(df$m, na.rm = TRUE)
    max_m <- max(df$m, na.rm = TRUE)
    a <- df[df$m == min_m, ][1, ]
    gg <- df[df$m == max_m, ][1, ]
    data.frame(
      group = df$group[1],
      time = df$time[1],
      metric = df$metric[1],
      N = max_m,
      alpha_scale_value = a$mean,
      gamma_scale_value = gg$mean,
      alpha_scale_lwr = a$lwr,
      alpha_scale_upr = a$upr,
      gamma_scale_lwr = gg$lwr,
      gamma_scale_upr = gg$upr,
      gamma_minus_alpha = gg$mean - a$mean,
      stringsAsFactors = FALSE
    )
  }))
  rownames(alpha_gamma_summary) <- NULL
  alpha_gamma_summary <- alpha_gamma_summary[order(alpha_gamma_summary$group, alpha_gamma_summary$time, alpha_gamma_summary$metric), ]

  # Optional group contrast on temporal CV.
  group_contrast <- NULL
  groups <- sort(unique(temporal_variability$group))
  if (is.null(contrast_pair) && length(groups) == 2) contrast_pair <- groups
  if (!is.null(contrast_pair)) {
    if (length(contrast_pair) != 2) .tr_stop("contrast_pair must have two group names.")
    if (!all(contrast_pair %in% groups)) {
      .tr_stop("contrast_pair groups not found. Available groups: ", paste(groups, collapse = ", "))
    }
    a <- temporal_variability[temporal_variability$group == contrast_pair[1], c("metric", "m", "temporal_cv")]
    b <- temporal_variability[temporal_variability$group == contrast_pair[2], c("metric", "m", "temporal_cv")]
    names(a)[3] <- paste0("temporal_cv_", contrast_pair[1])
    names(b)[3] <- paste0("temporal_cv_", contrast_pair[2])
    group_contrast <- merge(a, b, by = c("metric", "m"), all = FALSE)
    group_contrast$delta_cv <- group_contrast[[paste0("temporal_cv_", contrast_pair[2])]] - group_contrast[[paste0("temporal_cv_", contrast_pair[1])]]
    group_contrast$direction <- ifelse(
      is.na(group_contrast$delta_cv), NA_character_,
      ifelse(group_contrast$delta_cv > 0,
             paste(contrast_pair[2], ">", contrast_pair[1]),
             ifelse(group_contrast$delta_cv < 0,
                    paste(contrast_pair[1], ">", contrast_pair[2]),
                    "No difference"))
    )
    group_contrast <- group_contrast[order(group_contrast$metric, group_contrast$m), ]
  }

  settings <- list(
    format = format,
    taxonomic_metrics = taxonomic_metrics,
    functional_metrics = functional_metrics,
    custom_metrics = if (!is.null(custom_metrics)) names(custom_metrics) else character(0),
    min_trait_coverage = min_trait_coverage,
    standardize_traits = standardize_traits,
    cwm_baseline = cwm_baseline,
    n_perm = n_perm,
    use_exact = use_exact,
    exact_threshold = exact_threshold,
    richness_exact_shinozaki = richness_exact_shinozaki,
    interval_probs = interval_probs,
    auc_method = auc_method,
    contrast_pair = contrast_pair
  )

  out <- list(
    rarefaction_profiles = rarefaction_profiles,
    temporal_variability = temporal_variability,
    auc_summary = auc_summary,
    alpha_gamma_summary = alpha_gamma_summary,
    group_contrast = group_contrast,
    settings = settings
  )
  class(out) <- c("temporal_rarefaction", class(out))
  out
}

# Alias
temporal_rarefy <- temporal_rarefaction

# ------------------------------------------------------------------------------
# Convenience methods
# ------------------------------------------------------------------------------

print.temporal_rarefaction <- function(x, ...) {
  cat("Temporal rarefaction result\n")
  cat("  Rarefaction rows:       ", nrow(x$rarefaction_profiles), "\n", sep = "")
  cat("  Temporal profile rows:  ", nrow(x$temporal_variability), "\n", sep = "")
  cat("  Groups:                 ", paste(unique(x$rarefaction_profiles$group), collapse = ", "), "\n", sep = "")
  cat("  Metrics:                ", paste(unique(x$rarefaction_profiles$metric), collapse = ", "), "\n", sep = "")
  invisible(x)
}

# ------------------------------------------------------------------------------
# Minimal examples
# ------------------------------------------------------------------------------
# Example 1: taxonomic temporal rarefaction from wide data
#
# set.seed(1)
# comm_wide <- data.frame(
#   plot = rep(paste0("P", 1:6), times = 3),
#   year = rep(2001:2003, each = 6),
#   group = rep(c("A", "A", "A", "B", "B", "B"), times = 3),
#   sp1 = rpois(18, 3),
#   sp2 = rpois(18, 2),
#   sp3 = rpois(18, 1),
#   sp4 = rpois(18, 1)
# )
# ans <- temporal_rarefaction(
#   comm = comm_wide,
#   format = "wide",
#   time_col = "year",
#   unit_col = "plot",
#   group_col = "group",
#   taxonomic_metrics = c("richness", "shannon"),
#   n_perm = 999,
#   use_exact = TRUE,
#   seed = 1
# )
# ans$temporal_variability
# ans$auc_summary
#
# Example 2: add functional metrics
#
# traits <- data.frame(
#   species = c("sp1", "sp2", "sp3", "sp4"),
#   trait1 = c(1.1, 2.3, 0.7, 3.2),
#   trait2 = c(4.0, 2.1, 3.5, 1.2)
# )
# ans_fun <- temporal_rarefaction(
#   comm = comm_wide,
#   format = "wide",
#   time_col = "year",
#   unit_col = "plot",
#   group_col = "group",
#   traits = traits,
#   trait_species_col = "species",
#   taxonomic_metrics = c("richness", "shannon"),
#   functional_metrics = c("RaoQ", "FDis", "FRic2D", "CWMdist"),
#   n_perm = 999,
#   use_exact = TRUE,
#   seed = 1
# )
# ans_fun$group_contrast
#
# Example 3: custom metric
#
# evenness_fun <- function(abund, traits = NULL, trait_dist = NULL, trait_scores = NULL) {
#   abund <- abund[abund > 0]
#   if (length(abund) <= 1) return(NA_real_)
#   p <- abund / sum(abund)
#   H <- -sum(p * log(p))
#   H / log(length(p))
# }
# ans_custom <- temporal_rarefaction(
#   comm = comm_wide,
#   format = "wide",
#   time_col = "year",
#   unit_col = "plot",
#   group_col = "group",
#   custom_metrics = list(Pielou_evenness = evenness_fun),
#   n_perm = 999,
#   seed = 1
# )
