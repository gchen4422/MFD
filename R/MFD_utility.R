#' MFD: Decision-guided multi-ancestry fine-mapping
#'
#' Utilities for decision-guided multi-ancestry fine-mapping using
#' SuSiE post-hoc and MESuSiE.
#'
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"

.datatable.aware <- TRUE

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.check_required_cols <- function(x, cols, arg_name) {
  missing_cols <- setdiff(cols, colnames(x))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "%s is missing required columns: %s",
        arg_name,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

.check_pop_names <- function(pop_names) {
  if (!is.character(pop_names) || length(pop_names) != 2) {
    stop("pop_names must be a character vector of length 2.", call. = FALSE)
  }
}

.standardize_sumstats <- function(x, arg_name) {
  .check_required_cols(x, c("SNP", "PVAL"), arg_name)
  x <- data.table::copy(data.table::as.data.table(x))
  x[, SNP := as.character(SNP)]
  x
}

.standardize_gwas <- function(x, arg_name) {
  .check_required_cols(x, c("SNP", "CHR", "POS", "Z", "Beta", "Se", "N"), arg_name)
  x <- data.table::copy(data.table::as.data.table(x))
  x[, SNP := as.character(SNP)]
  x
}

.check_and_convert_cor <- function(mat) {
  if (!is.matrix(mat)) {
    stop("LD input must be a matrix.", call. = FALSE)
  }
  if (nrow(mat) != ncol(mat)) {
    stop("LD matrix must be square.", call. = FALSE)
  }
  if (any(abs(diag(mat) - 1) > 0.1, na.rm = TRUE)) {
    stats::cov2cor(mat)
  } else {
    mat
  }
}

.set_ld_dimnames <- function(ld, snps) {
  if (!is.matrix(ld)) {
    stop("LD input must be a matrix.", call. = FALSE)
  }
  if (nrow(ld) != length(snps) || ncol(ld) != length(snps)) {
    stop(
      "The LD matrix dimensions do not match the number of SNPs in the GWAS table.",
      call. = FALSE
    )
  }
  rownames(ld) <- snps
  colnames(ld) <- snps
  ld
}

.merge_snp_coordinates <- function(gwas_1, gwas_2) {
  out <- dplyr::full_join(
    dplyr::select(as.data.frame(gwas_1), dplyr::all_of(c("SNP", "CHR", "POS"))),
    dplyr::select(as.data.frame(gwas_2), dplyr::all_of(c("SNP", "CHR", "POS"))),
    by = "SNP"
  )
  
  out <- dplyr::mutate(
    out,
    CHR = dplyr::coalesce(.data$CHR.x, .data$CHR.y),
    POS = dplyr::coalesce(.data$POS.x, .data$POS.y)
  )
  
  dplyr::select(out, dplyr::all_of(c("SNP", "CHR", "POS")))
}

.replace_missing_except_snp <- function(df) {
  cols_to_fill <- setdiff(colnames(df), "SNP")
  dplyr::mutate(
    df,
    dplyr::across(dplyr::all_of(cols_to_fill), ~ tidyr::replace_na(., 0))
  )
}

.finalize_mf_result <- function(mf_result) {
  colnames(mf_result) <- c(
    "SNP", "CHR", "POS", "PIP_Either", "PIP_Shared",
    "PIP_Ancestry_1", "PIP_Ancestry_2", "CS",
    "CS_Ancestry_1", "CS_Ancestry_2"
  )
  dplyr::arrange(mf_result, .data$CHR, .data$POS)
}

.get_pop_cs_vec <- function(mesusie_obj, pop_name) {
  vec <- rep(0L, length(mesusie_obj$pip))
  pop_sets <- mesusie_obj$cs$cs[mesusie_obj$cs$cs_category == pop_name]
  
  if (length(pop_sets) > 0) {
    for (j in seq_along(pop_sets)) {
      vec[pop_sets[[j]]] <- j
    }
  }
  
  vec
}

# local replacement for susieR:::in_CS_x so package code does not rely on :::
.in_cs_x <- function(x, coverage = 0.95) {
  if (!is.numeric(x)) {
    stop("Input to .in_cs_x must be numeric.", call. = FALSE)
  }
  if (length(x) == 0) {
    return(integer(0))
  }
  if (all(is.na(x))) {
    return(rep(0L, length(x)))
  }
  
  x[is.na(x)] <- 0
  o <- order(x, decreasing = TRUE)
  cx <- cumsum(x[o])
  
  k <- which(cx >= coverage)[1]
  if (is.na(k)) {
    k <- length(x)
  }
  
  out <- rep(0L, length(x))
  out[o[seq_len(k)]] <- 1L
  out
}

# -------------------------------------------------------------------------
# Exported functions
# -------------------------------------------------------------------------

#' Decide between MESuSiE and SuSiE post-hoc
#'
#' Uses ancestry-specific significant variants and LD relationships to choose
#' between joint and post-hoc fine-mapping strategies.
#'
#' @param sum_eur Summary statistics for ancestry 1. Must contain `SNP` and `PVAL`.
#' @param sum_afr Summary statistics for ancestry 2. Must contain `SNP` and `PVAL`.
#' @param ld_eur LD matrix for ancestry 1.
#' @param ld_afr LD matrix for ancestry 2.
#' @param p_thresh P-value threshold used to define significant variants.
#' @param r2_thresh Squared-correlation threshold used to define high LD.
#'
#' @return A list with elements `method`, `reason_code`, and `reason`.
#' @export
decide_finemapping_method <- function(
    sum_eur,
    sum_afr,
    ld_eur,
    ld_afr,
    p_thresh = 5e-8,
    r2_thresh = 0.6
) {
  sum_eur <- .standardize_sumstats(sum_eur, "sum_eur")
  sum_afr <- .standardize_sumstats(sum_afr, "sum_afr")
  
  R_eur <- .check_and_convert_cor(ld_eur)
  R_afr <- .check_and_convert_cor(ld_afr)
  
  shared_snps <- intersect(sum_eur$SNP, sum_afr$SNP)
  asv_eur <- setdiff(sum_eur$SNP, sum_afr$SNP)
  asv_afr <- setdiff(sum_afr$SNP, sum_eur$SNP)
  
  sig_asv_eur <- sum_eur[SNP %in% asv_eur & PVAL < p_thresh, SNP]
  sig_asv_afr <- sum_afr[SNP %in% asv_afr & PVAL < p_thresh, SNP]
  has_sig_asv <- length(sig_asv_eur) > 0 || length(sig_asv_afr) > 0
  
  if (!has_sig_asv) {
    return(list(
      method = "MESuSiE",
      reason_code = 1L,
      reason = "AS-V not significant"
    ))
  }
  
  sig_shared_eur <- sum_eur[SNP %in% shared_snps & PVAL < p_thresh, SNP]
  sig_shared_afr <- sum_afr[SNP %in% shared_snps & PVAL < p_thresh, SNP]
  all_sig_shared <- unique(c(sig_shared_eur, sig_shared_afr))
  
  if (length(all_sig_shared) == 0) {
    return(list(
      method = "SuSiE post-hoc",
      reason_code = 2L,
      reason = "AS-V significant, but not in high LD with shared signals"
    ))
  }
  
  high_ld_found <- FALSE
  high_ld_partners <- character()
  
  if (length(sig_asv_eur) > 0) {
    valid_asv <- intersect(sig_asv_eur, rownames(R_eur))
    valid_shared <- intersect(all_sig_shared, colnames(R_eur))
    
    if (length(valid_asv) > 0 && length(valid_shared) > 0) {
      ld_sub <- R_eur[valid_asv, valid_shared, drop = FALSE]
      if (max(ld_sub^2, na.rm = TRUE) > r2_thresh) {
        high_ld_found <- TRUE
        high_ld_partners <- c(
          high_ld_partners,
          colnames(ld_sub)[apply(ld_sub^2, 2, max, na.rm = TRUE) > r2_thresh]
        )
      }
    }
  }
  
  if (length(sig_asv_afr) > 0) {
    valid_asv <- intersect(sig_asv_afr, rownames(R_afr))
    valid_shared <- intersect(all_sig_shared, colnames(R_afr))
    
    if (length(valid_asv) > 0 && length(valid_shared) > 0) {
      ld_sub <- R_afr[valid_asv, valid_shared, drop = FALSE]
      if (max(ld_sub^2, na.rm = TRUE) > r2_thresh) {
        high_ld_found <- TRUE
        high_ld_partners <- c(
          high_ld_partners,
          colnames(ld_sub)[apply(ld_sub^2, 2, max, na.rm = TRUE) > r2_thresh]
        )
      }
    }
  }
  
  high_ld_partners <- unique(high_ld_partners)
  
  if (!high_ld_found) {
    return(list(
      method = "SuSiE post-hoc",
      reason_code = 2L,
      reason = "AS-V significant, but not in high LD with shared signals"
    ))
  }
  
  shared_check <- dplyr::left_join(
    data.frame(SNP = high_ld_partners, stringsAsFactors = FALSE),
    data.frame(SNP = sum_eur$SNP, PVAL_eur = sum_eur$PVAL, stringsAsFactors = FALSE),
    by = "SNP"
  )
  shared_check <- dplyr::left_join(
    shared_check,
    data.frame(SNP = sum_afr$SNP, PVAL_afr = sum_afr$PVAL, stringsAsFactors = FALSE),
    by = "SNP"
  )
  
  is_sig_both <- shared_check$PVAL_eur < p_thresh & shared_check$PVAL_afr < p_thresh
  is_sig_both[is.na(is_sig_both)] <- FALSE
  
  if (any(is_sig_both)) {
    return(list(
      method = "MESuSiE",
      reason_code = 3L,
      reason = "AS-V significant and in high LD with a shared SNP significant in all ancestries"
    ))
  }
  
  list(
    method = "SuSiE post-hoc",
    reason_code = 4L,
    reason = "AS-V significant and in high LD with a shared SNP significant in only specific ancestry"
  )
}

#' Compute credible set membership
#'
#' Applies a credible set membership rule row-wise to a matrix of posterior
#' probabilities.
#'
#' @param res Numeric matrix of posterior probabilities.
#' @param coverage Credible set coverage threshold.
#'
#' @return A matrix of 0/1 indicators.
#' @export
in_CS <- function(res, coverage = 0.95) {
  if (!is.matrix(res)) {
    stop("res must be a matrix.", call. = FALSE)
  }
  t(apply(res, 1, function(x) .in_cs_x(x, coverage = coverage)))
}

#' Compute purity statistics for a credible set
#'
#' @param pos Integer vector of SNP indices.
#' @param Xcorr List of correlation matrices.
#'
#' @return A numeric vector containing minimum, mean, and median absolute
#'   correlation across the supplied matrices.
#' @export
get_purity <- function(pos, Xcorr) {
  if (!is.list(Xcorr) || length(Xcorr) == 0) {
    stop("Xcorr must be a non-empty list of correlation matrices.", call. = FALSE)
  }
  
  if (length(pos) == 1) {
    return(c(1, 1, 1))
  }
  
  value_list <- lapply(Xcorr, function(x) c(abs(x[pos, pos])))
  value_matrix <- Reduce(cbind, value_list)
  value_max <- do.call(pmax, data.frame(value_matrix))
  
  c(
    min(value_max, na.rm = TRUE),
    mean(value_max, na.rm = TRUE),
    median(value_max, na.rm = TRUE)
  )
}

#' Extract ancestry-specific credible sets from a MESuSiE result
#'
#' @param res A MESuSiE result object.
#' @param Xcorr List of correlation matrices.
#' @param target_idx Target configuration index.
#' @param coverage Credible set coverage threshold.
#' @param prior_tol Threshold for retaining effects.
#' @param cor_method Purity metric to threshold. One of
#'   `"min.abs.corr"`, `"mean.abs.corr"`, or `"median.abs.corr"`.
#' @param cor_threshold Minimum purity threshold.
#'
#' @return A list with elements `cs` and `cs_index`.
#' @export
meSuSie_get_cs_specific <- function(
    res,
    Xcorr,
    target_idx,
    coverage = 0.95,
    prior_tol = 1e-9,
    cor_method = "min.abs.corr",
    cor_threshold = 0.5
) {
  include_idx <- unlist(lapply(res$V, function(x) max(diag(x)) > prior_tol))
  
  alpha_specific <- t(Reduce(cbind, lapply(res$alpha, function(x) {
    target_vec <- x[, target_idx]
    if (sum(target_vec) < 1e-10) {
      return(rep(0, length(target_vec)))
    }
    target_vec / sum(target_vec)
  })))
  
  status <- in_CS(alpha_specific, coverage = coverage)
  cs <- lapply(seq_len(nrow(status)), function(i) which(status[i, ] != 0))
  
  include_idx <- include_idx * (lengths(cs) > 0)
  include_idx <- include_idx * (!duplicated(cs))
  include_idx <- as.logical(include_idx)
  
  if (sum(include_idx) == 0) {
    return(list(cs = NULL, cs_index = NULL))
  }
  
  cs <- cs[include_idx]
  
  purity <- data.frame(do.call(rbind, lapply(seq_along(cs), function(i) {
    get_purity(cs[[i]], Xcorr)
  })))
  colnames(purity) <- c("min.abs.corr", "mean.abs.corr", "median.abs.corr")
  
  if (!cor_method %in% colnames(purity)) {
    stop("Invalid cor_method.", call. = FALSE)
  }
  
  is_pure <- which(purity[, cor_method] >= cor_threshold)
  
  if (length(is_pure) == 0) {
    return(list(cs = NULL, cs_index = NULL))
  }
  
  cs <- cs[is_pure]
  cs_index <- which(include_idx)[is_pure]
  names(cs) <- paste0("L", cs_index)
  
  list(cs = cs, cs_index = cs_index)
}

#' Map credible set labels to a full SNP vector
#'
#' @param cs_res Result object containing `cs` and optionally `cs_index`.
#' @param n_snps Total number of SNPs.
#' @param renumber If `TRUE`, renumber sets sequentially.
#'
#' @return Integer vector of length `n_snps`.
#' @export
get_cs_index_vector <- function(cs_res, n_snps, renumber = TRUE) {
  cs_vec <- rep(0L, n_snps)
  
  if (!is.null(cs_res) && !is.null(cs_res$cs)) {
    for (j in seq_along(cs_res$cs)) {
      assign_id <- if (renumber) j else cs_res$cs_index[j]
      snp_indices <- cs_res$cs[[j]]
      cs_vec[snp_indices] <- assign_id
    }
  }
  
  cs_vec
}

#' Run decision-guided multi-ancestry fine-mapping
#'
#' Chooses between SuSiE post-hoc and MESuSiE and returns harmonized results.
#'
#' @param gwas_1 GWAS summary statistics for ancestry 1. Must contain
#'   `SNP`, `CHR`, `POS`, `Beta`, `Se`, `Z`, and `N`.
#' @param gwas_2 GWAS summary statistics for ancestry 2. Must contain
#'   `SNP`, `CHR`, `POS`, `Beta`, `Se`, `Z`, and `N`.
#' @param ld_1 LD matrix matched to `gwas_1`.
#' @param ld_2 LD matrix matched to `gwas_2`.
#' @param pop_names Character vector of length 2 giving ancestry labels.
#'
#' @return A data frame with harmonized PIPs and credible sets.
#' @export
run_mf_decision_fm <- function(
    gwas_1,
    gwas_2,
    ld_1,
    ld_2,
    pop_names = c("Pop1", "Pop2")
) {
  .check_pop_names(pop_names)
  
  gwas_1 <- .standardize_gwas(gwas_1, "gwas_1")
  gwas_2 <- .standardize_gwas(gwas_2, "gwas_2")
  
  ld_1 <- .set_ld_dimnames(ld_1, gwas_1$SNP)
  ld_2 <- .set_ld_dimnames(ld_2, gwas_2$SNP)
  
  decision <- decide_finemapping_method(
    sum_eur = gwas_1,
    sum_afr = gwas_2,
    ld_eur = ld_1,
    ld_afr = ld_2
  )
  
  if (decision$method == "SuSiE post-hoc") {
    susie_1 <- susieR::susie_rss(gwas_1$Z, ld_1, check_prior = FALSE)
    susie_2 <- susieR::susie_rss(gwas_2$Z, ld_2, check_prior = FALSE)
    
    cs_1 <- get_cs_index_vector(susie_1$sets, nrow(gwas_1), renumber = TRUE)
    cs_2 <- get_cs_index_vector(susie_2$sets, nrow(gwas_2), renumber = TRUE)
    
    res1_df <- data.frame(SNP = gwas_1$SNP, PIP1 = susie_1$pip, CS1 = cs_1)
    res2_df <- data.frame(SNP = gwas_2$SNP, PIP2 = susie_2$pip, CS2 = cs_2)
    
    mf_result <- .merge_snp_coordinates(gwas_1, gwas_2)
    mf_result <- dplyr::left_join(mf_result, res1_df, by = "SNP")
    mf_result <- dplyr::left_join(mf_result, res2_df, by = "SNP")
    mf_result <- .replace_missing_except_snp(mf_result)
    mf_result <- dplyr::mutate(
      mf_result,
      PIP_Either = pmax(.data$PIP1, .data$PIP2),
      PIP_Shared = pmin(.data$PIP1, .data$PIP2),
      CS = ifelse(.data$CS1 + .data$CS2 == 0, 0, 1)
    )
    mf_result <- dplyr::select(
      mf_result,
      dplyr::all_of(c(
        "SNP", "CHR", "POS", "PIP_Either", "PIP_Shared",
        "PIP1", "PIP2", "CS", "CS1", "CS2"
      ))
    )
    
    return(.finalize_mf_result(mf_result))
  }
  
  common_snps <- intersect(gwas_1$SNP, gwas_2$SNP)
  if (length(common_snps) == 0) {
    stop("No shared SNPs found across the two GWAS inputs.", call. = FALSE)
  }
  
  g1_sub <- gwas_1[SNP %in% common_snps][order(match(SNP, common_snps))]
  g2_sub <- gwas_2[SNP %in% common_snps][order(match(SNP, common_snps))]
  
  ld1_sub <- ld_1[common_snps, common_snps, drop = FALSE]
  ld2_sub <- ld_2[common_snps, common_snps, drop = FALSE]
  
  summary_stat_list <- list(as.data.frame(g1_sub), as.data.frame(g2_sub))
  names(summary_stat_list) <- pop_names
  
  R_mat_list <- list(ld1_sub, ld2_sub)
  names(R_mat_list) <- pop_names
  
  mesusie_res <- MESuSiE::meSuSie_core(
    R_mat_list = R_mat_list,
    summary_stat_list = summary_stat_list,
    L = 10
  )
  
  cs_either <- get_cs_index_vector(mesusie_res$cs, length(mesusie_res$pip), renumber = TRUE)
  cs_p1 <- .get_pop_cs_vec(mesusie_res, pop_names[1])
  cs_p2 <- .get_pop_cs_vec(mesusie_res, pop_names[2])
  
  mesusie_df <- data.frame(
    SNP = common_snps,
    PIP_Either = mesusie_res$pip,
    PIP1 = mesusie_res$pip_config[, 1],
    PIP2 = mesusie_res$pip_config[, 2],
    PIP_Shared = mesusie_res$pip_config[, 3],
    CS = cs_either,
    CS1 = cs_p1,
    CS2 = cs_p2
  )
  
  mf_result <- .merge_snp_coordinates(gwas_1, gwas_2)
  mf_result <- dplyr::left_join(mf_result, mesusie_df, by = "SNP")
  mf_result <- .replace_missing_except_snp(mf_result)
  mf_result <- dplyr::select(
    mf_result,
    dplyr::all_of(c(
      "SNP", "CHR", "POS", "PIP_Either", "PIP_Shared",
      "PIP1", "PIP2", "CS", "CS1", "CS2"
    ))
  )
  
  .finalize_mf_result(mf_result)
}

#' Run decision-guided multi-ancestry fine-mapping with raw outputs
#'
#' This version returns the model choice, harmonized results, and raw fitted
#' objects from either SuSiE or MESuSiE.
#'
#' @param gwas_1 GWAS summary statistics for ancestry 1. Must contain
#'   `SNP`, `CHR`, `POS`, `Beta`, `Se`, `Z`, and `N`.
#' @param gwas_2 GWAS summary statistics for ancestry 2. Must contain
#'   `SNP`, `CHR`, `POS`, `Beta`, `Se`, `Z`, and `N`.
#' @param ld_1 LD matrix matched to `gwas_1`.
#' @param ld_2 LD matrix matched to `gwas_2`.
#' @param pop_names Character vector of length 2 giving ancestry labels.
#' @param L Number of single-effect components.
#' @param prior_weights Optional prior weights passed to SuSiE or MESuSiE.
#' @param ancestry_weight Optional ancestry weights passed to MESuSiE.
#' @param p_thresh P-value threshold used in the decision rule.
#' @param r2_thresh LD threshold used in the decision rule.
#'
#' @return A list with elements `decision`, `results`, and `raw_objects`.
#' @export
run_mf_decision <- function(
    gwas_1,
    gwas_2,
    ld_1,
    ld_2,
    pop_names = c("Pop1", "Pop2"),
    L = 10,
    prior_weights = NULL,
    ancestry_weight = NULL,
    p_thresh = 5e-8,
    r2_thresh = 0.6
) {
  .check_pop_names(pop_names)
  
  gwas_1 <- .standardize_gwas(gwas_1, "gwas_1")
  gwas_2 <- .standardize_gwas(gwas_2, "gwas_2")
  
  ld_1 <- .set_ld_dimnames(ld_1, gwas_1$SNP)
  ld_2 <- .set_ld_dimnames(ld_2, gwas_2$SNP)
  
  decision <- decide_finemapping_method(
    sum_eur = gwas_1,
    sum_afr = gwas_2,
    ld_eur = ld_1,
    ld_afr = ld_2,
    p_thresh = p_thresh,
    r2_thresh = r2_thresh
  )
  
  raw_objects <- list()
  
  if (decision$method == "SuSiE post-hoc") {
    susie_1 <- susieR::susie_rss(
      z = gwas_1$Z,
      R = ld_1,
      L = L,
      prior_weights = prior_weights,
      check_prior = FALSE
    )
    susie_2 <- susieR::susie_rss(
      z = gwas_2$Z,
      R = ld_2,
      L = L,
      prior_weights = prior_weights,
      check_prior = FALSE
    )
    
    raw_objects$susie_1 <- susie_1
    raw_objects$susie_2 <- susie_2
    
    cs_1 <- get_cs_index_vector(susie_1$sets, nrow(gwas_1))
    cs_2 <- get_cs_index_vector(susie_2$sets, nrow(gwas_2))
    
    res1_df <- data.frame(
      SNP = gwas_1$SNP,
      PIP1 = susie_1$pip,
      CS1 = cs_1,
      stringsAsFactors = FALSE
    )
    res2_df <- data.frame(
      SNP = gwas_2$SNP,
      PIP2 = susie_2$pip,
      CS2 = cs_2,
      stringsAsFactors = FALSE
    )
    
    temp_res <- .merge_snp_coordinates(gwas_1, gwas_2)
    temp_res <- dplyr::left_join(temp_res, res1_df, by = "SNP")
    temp_res <- dplyr::left_join(temp_res, res2_df, by = "SNP")
    temp_res <- dplyr::mutate(
      temp_res,
      dplyr::across(
        dplyr::all_of(c("PIP1", "PIP2", "CS1", "CS2")),
        ~ tidyr::replace_na(., 0)
      )
    )
    temp_res <- dplyr::mutate(
      temp_res,
      PIP_Either = pmax(.data$PIP1, .data$PIP2)
    )
    
    cs_snps_idx <- which(temp_res$CS1 > 0 | temp_res$CS2 > 0)
    
    if (length(cs_snps_idx) > 0) {
      active_snps <- temp_res$SNP[cs_snps_idx]
      
      idx1 <- match(active_snps, colnames(ld_1))
      idx2 <- match(active_snps, colnames(ld_2))
      
      sub_ld1 <- ld_1[idx1, idx1, drop = FALSE]
      sub_ld2 <- ld_2[idx2, idx2, drop = FALSE]
      sub_ld1[is.na(sub_ld1)] <- 0
      sub_ld2[is.na(sub_ld2)] <- 0
      
      ld_either <- pmax(abs(sub_ld1), abs(sub_ld2))
      
      if (length(active_snps) > 1) {
        h_dynamic <- 1 - sqrt(r2_thresh)
        dist_mat <- stats::as.dist(1 - ld_either)
        hc <- stats::hclust(dist_mat, method = "complete")
        clusters <- stats::cutree(hc, h = h_dynamic)
      } else {
        clusters <- 1L
      }
      
      temp_res$CS <- 0L
      temp_res$CS[cs_snps_idx] <- clusters
    } else {
      temp_res$CS <- 0L
    }
    
    mf_result <- dplyr::mutate(
      temp_res,
      PIP_Shared = pmin(.data$PIP1, .data$PIP2)
    )
    mf_result <- dplyr::select(
      mf_result,
      dplyr::all_of(c(
        "SNP", "CHR", "POS", "PIP_Either", "PIP_Shared",
        "PIP1", "PIP2", "CS", "CS1", "CS2"
      ))
    )
    
    return(list(
      decision = decision,
      results = .finalize_mf_result(mf_result),
      raw_objects = raw_objects
    ))
  }
  
  common_snps <- intersect(gwas_1$SNP, gwas_2$SNP)
  if (length(common_snps) == 0) {
    stop("No shared SNPs found across the two GWAS inputs.", call. = FALSE)
  }
  
  g1_sub <- gwas_1[SNP %in% common_snps][order(match(SNP, common_snps))]
  g2_sub <- gwas_2[SNP %in% common_snps][order(match(SNP, common_snps))]
  
  ld1_sub <- ld_1[common_snps, common_snps, drop = FALSE]
  ld2_sub <- ld_2[common_snps, common_snps, drop = FALSE]
  
  summary_stat_list <- list(as.data.frame(g1_sub), as.data.frame(g2_sub))
  names(summary_stat_list) <- pop_names
  
  R_mat_list <- list(ld1_sub, ld2_sub)
  names(R_mat_list) <- pop_names
  
  mesusie_res <- MESuSiE::meSuSie_core(
    R_mat_list = R_mat_list,
    summary_stat_list = summary_stat_list,
    L = L,
    prior_weights = prior_weights,
    ancestry_weight = ancestry_weight
  )
  
  raw_objects$mesusie_res <- mesusie_res
  
  cs_either <- get_cs_index_vector(mesusie_res$cs, length(mesusie_res$pip), renumber = TRUE)
  cs_p1 <- .get_pop_cs_vec(mesusie_res, pop_names[1])
  cs_p2 <- .get_pop_cs_vec(mesusie_res, pop_names[2])
  
  mesusie_df <- data.frame(
    SNP = common_snps,
    PIP_Either = mesusie_res$pip,
    PIP1 = mesusie_res$pip_config[, 1],
    PIP2 = mesusie_res$pip_config[, 2],
    PIP_Shared = mesusie_res$pip_config[, 3],
    CS = cs_either,
    CS1 = cs_p1,
    CS2 = cs_p2
  )
  
  mf_result <- .merge_snp_coordinates(gwas_1, gwas_2)
  mf_result <- dplyr::left_join(mf_result, mesusie_df, by = "SNP")
  mf_result <- .replace_missing_except_snp(mf_result)
  mf_result <- dplyr::select(
    mf_result,
    dplyr::all_of(c(
      "SNP", "CHR", "POS", "PIP_Either", "PIP_Shared",
      "PIP1", "PIP2", "CS", "CS1", "CS2"
    ))
  )
  
  list(
    decision = decision,
    results = .finalize_mf_result(mf_result),
    raw_objects = raw_objects
  )
}