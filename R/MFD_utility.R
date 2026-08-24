#' MFD: Decision-guided multi-ancestry fine-mapping for K ancestries
#'
#' Utilities for decision-guided multi-ancestry fine-mapping using
#' SuSiE post-hoc and MESuSiE, generalized from 2 to K (>= 2) ancestry groups.
#'
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"

.datatable.aware <- TRUE
utils::globalVariables(c(":=", "ALT", "Beta", "MAF", "PVAL", "REF", "SNP", "Z"))

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
  if (!is.character(pop_names) || length(pop_names) < 2) {
    stop("pop_names must be a character vector of length >= 2.", call. = FALSE)
  }
}

.standardize_sumstats <- function(x, arg_name) {
  .check_required_cols(x, c("SNP", "PVAL"), arg_name)
  x <- data.table::copy(data.table::as.data.table(x))
  x[, SNP := as.character(SNP)]
  x
}

.standardize_gwas <- function(x, arg_name) {
  .check_required_cols(x, c("SNP", "CHR", "POS", "Z", "Beta", "Se", "N", "PVAL"), arg_name)
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

.harmonize_alleles <- function(g_subs, ld_subs, common_snps) {
  ref <- g_subs[[1]]
  K <- length(g_subs)

  for (i in 2:K) {
    tgt <- g_subs[[i]]
    same_al <- (tgt$ALT == ref$ALT) & (tgt$REF == ref$REF)
    flip_al <- (tgt$ALT == ref$REF) & (tgt$REF == ref$ALT)
    drop_al <- !same_al & !flip_al

    if (any(drop_al)) {
      keep <- !drop_al
      common_snps <- common_snps[keep]
      g_subs <- lapply(g_subs, function(g) g[keep])
      ld_subs <- lapply(ld_subs, function(ld) ld[keep, keep, drop = FALSE])
      flip_al <- flip_al[keep]
      ref <- g_subs[[1]]
      tgt <- g_subs[[i]]
    }

    flip_idx <- which(flip_al)
    if (length(flip_idx) > 0) {
      tgt[flip_idx, Z    := -Z]
      tgt[flip_idx, Beta := -Beta]
      tgt[flip_idx, MAF  := 1 - MAF]
      old_alt <- tgt$ALT[flip_idx]
      tgt[flip_idx, ALT := REF]
      tgt[flip_idx, REF := old_alt]
      g_subs[[i]] <- tgt

      ld_i <- ld_subs[[i]]
      ld_i[flip_idx, ] <- -ld_i[flip_idx, ]
      ld_i[, flip_idx] <- -ld_i[, flip_idx]
      ld_subs[[i]] <- ld_i
    }
  }

  list(g_subs = g_subs, ld_subs = ld_subs, common_snps = common_snps)
}

.merge_snp_coordinates <- function(gwas_list) {
  coord_list <- lapply(gwas_list, function(g) {
    dplyr::select(as.data.frame(g), dplyr::all_of(c("SNP", "CHR", "POS")))
  })
  out <- coord_list[[1]]
  if (length(coord_list) > 1) {
    for (i in 2:length(coord_list)) {
      out <- dplyr::full_join(out, coord_list[[i]], by = "SNP",
                              suffix = c("", paste0(".", i)))
      chr_cols <- grep("^CHR", colnames(out), value = TRUE)
      pos_cols <- grep("^POS", colnames(out), value = TRUE)
      out$CHR <- do.call(dplyr::coalesce, out[chr_cols])
      out$POS <- do.call(dplyr::coalesce, out[pos_cols])
      out <- dplyr::select(out, dplyr::all_of(c("SNP", "CHR", "POS")))
    }
  }
  out
}

.replace_missing_except_snp <- function(df) {
  cols_to_fill <- setdiff(colnames(df), "SNP")
  dplyr::mutate(
    df,
    dplyr::across(dplyr::all_of(cols_to_fill), ~ tidyr::replace_na(., 0))
  )
}

.finalize_mf_result <- function(mf_result, pop_names) {
  K <- length(pop_names)
  old_pip <- paste0("PIP_", pop_names)
  old_cs  <- paste0("CS_", pop_names)
  new_pip <- paste0("PIP_Ancestry_", seq_len(K))
  new_cs  <- paste0("CS_Ancestry_", seq_len(K))
  for (i in seq_len(K)) {
    colnames(mf_result)[colnames(mf_result) == old_pip[i]] <- new_pip[i]
    colnames(mf_result)[colnames(mf_result) == old_cs[i]]  <- new_cs[i]
  }
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

#' Decide between MESuSiE and SuSiE post-hoc for K ancestries
#'
#' Uses ancestry-specific significant variants and LD relationships to choose
#' between joint and post-hoc fine-mapping strategies.  Each significant AS-SV
#' is evaluated jointly across all its signal-bearing ancestries: a single
#' K-way shared proxy must be in high LD in **every** ancestry where the AS-SV
#' is significant.
#'
#' @param gwas_list Named list of K GWAS summary statistics data frames.
#'   Each must contain `SNP` and `PVAL`.
#' @param ld_list Named list of K LD matrices corresponding to `gwas_list`.
#' @param p_thresh P-value threshold used to define significant variants.
#' @param r2_thresh Squared-correlation threshold used to define high LD.
#'
#' @return A list with elements:
#' \describe{
#'   \item{method}{`"MESuSiE"` or `"SuSiE post-hoc"`.}
#'   \item{reason_code}{Integer code (1--4).}
#'   \item{reason}{Human-readable reason.}
#'   \item{diagnostics}{A data frame (`NULL` when no significant AS-SVs exist)
#'     with one row per (AS-SV, qualifying proxy) pair.  Columns:
#'     `AS_SV`, `signal_ancestries`, `proxy_SNP`, `ld_values`,
#'     `proxy_pvals`, `proxy_globally_sig`, `classification`.
#'     AS-SVs with no qualifying proxy have a single row with `proxy_SNP = NA`
#'     and `classification = "untagged"`.}
#' }
#' @export
decide_finemapping_method <- function(
    gwas_list,
    ld_list,
    p_thresh = 5e-8,
    r2_thresh = 0.6
) {
  K <- length(gwas_list)
  pop_names <- names(gwas_list)
  if (is.null(pop_names)) pop_names <- paste0("Ancestry_", seq_len(K))

  gwas_list <- lapply(seq_len(K), function(i) {
    .standardize_sumstats(gwas_list[[i]], paste0("gwas_", i))
  })
  ld_list <- lapply(ld_list, .check_and_convert_cor)

  all_snp_sets <- lapply(gwas_list, function(g) g$SNP)
  shared_snps  <- Reduce(intersect, all_snp_sets)

  # AS-SV candidates per ancestry: SNPs present in that ancestry but absent from the K-way shared set
  asv_per_pop <- lapply(all_snp_sets, function(s) setdiff(s, shared_snps))

  # Significant AS-SV per ancestry
  sig_asv_per_pop <- mapply(function(gwas, asv) {
    gwas[SNP %in% asv & PVAL < p_thresh, SNP]
  }, gwas_list, asv_per_pop, SIMPLIFY = FALSE)

  has_sig_asv <- any(vapply(sig_asv_per_pop, length, integer(1)) > 0)
  if (!has_sig_asv) {
    return(list(
      method = "MESuSiE",
      reason_code = 1L,
      reason = "AS-SV not significant",
      diagnostics = NULL
    ))
  }

  all_sig_asvs <- unique(unlist(sig_asv_per_pop))

  # Shared SNPs significant in at least one ancestry
  sig_shared_per_pop <- lapply(gwas_list, function(g) {
    g[SNP %in% shared_snps & PVAL < p_thresh, SNP]
  })
  all_sig_shared <- unique(unlist(sig_shared_per_pop))

  # Pre-compute p-values for shared significant SNPs across all K ancestries
  if (length(all_sig_shared) > 0) {
    shared_pvals <- do.call(cbind, lapply(gwas_list, function(g) {
      m <- match(all_sig_shared, g$SNP)
      ifelse(is.na(m), NA_real_, g$PVAL[m])
    }))
    if (!is.matrix(shared_pvals)) {
      shared_pvals <- matrix(shared_pvals, nrow = 1)
    }
    rownames(shared_pvals) <- all_sig_shared
    colnames(shared_pvals) <- pop_names
    is_global <- apply(shared_pvals < p_thresh, 1, all, na.rm = FALSE)
    is_global[is.na(is_global)] <- FALSE
    globally_sig_shared <- all_sig_shared[is_global]
  } else {
    shared_pvals <- matrix(nrow = 0, ncol = K,
                           dimnames = list(NULL, pop_names))
    globally_sig_shared <- character(0)
  }

  # Evaluate each unique significant AS-SV jointly across ALL its
  # signal-bearing ancestries: a single proxy must be in high LD in every
  # ancestry where the AS-SV is significant.
  any_high_ld <- FALSE
  all_explained <- TRUE
  diag_list <- list()

  for (snp in all_sig_asvs) {
    bearing_idx   <- which(vapply(sig_asv_per_pop,
                                  function(s) snp %in% s, logical(1)))
    bearing_names <- pop_names[bearing_idx]

    # Progressively intersect high-LD proxies across signal-bearing ancestries
    candidate_proxies <- all_sig_shared

    for (ai in bearing_idx) {
      if (length(candidate_proxies) == 0) break
      R_i <- ld_list[[ai]]
      if (!(snp %in% rownames(R_i))) {
        candidate_proxies <- character(0)
        break
      }
      valid <- intersect(candidate_proxies, colnames(R_i))
      if (length(valid) == 0) {
        candidate_proxies <- character(0)
        break
      }
      r2_vals <- R_i[snp, valid]^2
      candidate_proxies <- valid[r2_vals >= r2_thresh]
    }

    # No qualifying proxy — AS-SV is untagged
    if (length(candidate_proxies) == 0) {
      diag_list[[length(diag_list) + 1]] <- data.frame(
        AS_SV              = snp,
        signal_ancestries  = paste(bearing_names, collapse = ","),
        proxy_SNP          = NA_character_,
        ld_values          = NA_character_,
        proxy_pvals        = NA_character_,
        proxy_globally_sig = NA,
        classification     = "untagged",
        stringsAsFactors   = FALSE
      )
      all_explained <- FALSE
      next
    }

    any_high_ld <- TRUE
    has_global <- any(candidate_proxies %in% globally_sig_shared)
    cls <- if (has_global) "explained" else "partial"
    if (!has_global) all_explained <- FALSE

    # One diagnostic row per qualifying proxy
    for (proxy in candidate_proxies) {
      ld_per_anc <- vapply(bearing_idx, function(ai) {
        ld_list[[ai]][snp, proxy]^2
      }, numeric(1))
      ld_str <- paste(sprintf("%s=%.4f", bearing_names, ld_per_anc),
                      collapse = ";")
      pval_str <- paste(sprintf("%s=%.2e", pop_names,
                                shared_pvals[proxy, ]),
                        collapse = ";")

      diag_list[[length(diag_list) + 1]] <- data.frame(
        AS_SV              = snp,
        signal_ancestries  = paste(bearing_names, collapse = ","),
        proxy_SNP          = proxy,
        ld_values          = ld_str,
        proxy_pvals        = pval_str,
        proxy_globally_sig = proxy %in% globally_sig_shared,
        classification     = cls,
        stringsAsFactors   = FALSE
      )
    }
  }

  diagnostics <- if (length(diag_list) > 0) do.call(rbind, diag_list) else NULL

  if (!any_high_ld) {
    return(list(
      method      = "SuSiE post-hoc",
      reason_code = 2L,
      reason      = "AS-SV significant, but no high-LD K-way shared proxy found",
      diagnostics = diagnostics
    ))
  }

  if (all_explained) {
    return(list(
      method      = "MESuSiE",
      reason_code = 3L,
      reason      = "All significant AS-SVs explained by shared SNPs significant in all ancestries",
      diagnostics = diagnostics
    ))
  }

  list(
    method      = "SuSiE post-hoc",
    reason_code = 4L,
    reason      = "Some significant AS-SVs lack a globally significant shared proxy",
    diagnostics = diagnostics
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
    stats::median(value_max, na.rm = TRUE)
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

#' Run decision-guided multi-ancestry fine-mapping (simplified, K ancestries)
#'
#' Chooses between SuSiE post-hoc and MESuSiE and returns harmonized results.
#' This is the simplified version without raw object returns.
#'
#' @param gwas_list Named list of K GWAS summary statistics data frames.
#'   Each must contain `SNP`, `CHR`, `POS`, `Beta`, `Se`, `Z`, `N`, and `PVAL`.
#' @param ld_list Named list of K LD matrices corresponding to `gwas_list`.
#' @param pop_names Character vector of length K giving ancestry labels.
#'   Defaults to `names(gwas_list)`.
#' @param L Number of single-effect components.
#'
#' @return A data frame with harmonized PIPs and credible sets.
#' @export
run_mf_decision_fm <- function(
    gwas_list,
    ld_list,
    pop_names = names(gwas_list),
    L = 10
) {
  K <- length(gwas_list)
  stopifnot(is.list(gwas_list), is.list(ld_list), length(ld_list) == K, K >= 2)
  if (is.null(pop_names)) pop_names <- paste0("Pop", seq_len(K))
  .check_pop_names(pop_names)

  gwas_list <- stats::setNames(
    lapply(seq_len(K), function(i) .standardize_gwas(gwas_list[[i]], paste0("gwas_", i))),
    pop_names
  )
  ld_list <- stats::setNames(
    lapply(seq_len(K), function(i) .set_ld_dimnames(ld_list[[i]], gwas_list[[i]]$SNP)),
    pop_names
  )

  decision <- decide_finemapping_method(
    gwas_list = gwas_list,
    ld_list   = ld_list
  )

  pip_cols <- paste0("PIP_", pop_names)
  cs_cols  <- paste0("CS_", pop_names)

  if (decision$method == "SuSiE post-hoc") {
    susie_list <- lapply(seq_len(K), function(i) {
      susieR::susie_rss(gwas_list[[i]]$Z, ld_list[[i]], L = L,
                        n = stats::median(gwas_list[[i]]$N), check_prior = FALSE)
    })
    names(susie_list) <- pop_names

    cs_list <- lapply(seq_len(K), function(i) {
      get_cs_index_vector(susie_list[[i]]$sets, nrow(gwas_list[[i]]), renumber = TRUE)
    })

    res_dfs <- lapply(seq_len(K), function(i) {
      df <- data.frame(
        SNP = gwas_list[[i]]$SNP,
        PIP = susie_list[[i]]$pip,
        CS  = cs_list[[i]],
        stringsAsFactors = FALSE
      )
      colnames(df) <- c("SNP", pip_cols[i], cs_cols[i])
      df
    })

    mf_result <- .merge_snp_coordinates(gwas_list)
    for (i in seq_len(K)) {
      mf_result <- dplyr::left_join(mf_result, res_dfs[[i]], by = "SNP")
    }
    mf_result <- .replace_missing_except_snp(mf_result)

    mf_result$PIP_Either <- do.call(pmax, mf_result[pip_cols])
    mf_result$PIP_Shared <- do.call(pmin, mf_result[pip_cols])
    mf_result$CS <- ifelse(rowSums(mf_result[cs_cols] > 0) == 0, 0L, 1L)

    select_cols <- c("SNP", "CHR", "POS", "PIP_Either", "PIP_Shared",
                     pip_cols, "CS", cs_cols)
    mf_result <- dplyr::select(mf_result, dplyr::all_of(select_cols))

    return(.finalize_mf_result(mf_result, pop_names))
  }

  # MESuSiE path
  common_snps <- Reduce(intersect, lapply(gwas_list, function(g) g$SNP))
  if (length(common_snps) == 0) {
    stop("No shared SNPs found across the GWAS inputs.", call. = FALSE)
  }

  g_subs <- lapply(gwas_list, function(g) {
    g[SNP %in% common_snps][order(match(SNP, common_snps))]
  })
  ld_subs <- lapply(ld_list, function(ld) {
    ld[common_snps, common_snps, drop = FALSE]
  })

  # Harmonize alleles to first ancestry, drop inconsistent SNPs
  harm <- .harmonize_alleles(g_subs, ld_subs, common_snps)
  g_subs      <- harm$g_subs
  ld_subs     <- harm$ld_subs
  common_snps <- harm$common_snps

  summary_stat_sd_list <- stats::setNames(lapply(g_subs, as.data.frame), pop_names)
  R_mat_list           <- stats::setNames(ld_subs, pop_names)

  mesusie_res <- MESuSiE::meSuSie_core(
    R_mat_list        = R_mat_list,
    summary_stat_list = summary_stat_sd_list,
    L = L
  )

  # Extract per-ancestry PIPs from pip_config directly
  # pip_config columns use combinatorial (combn) ordering:
  #   columns 1..K = single-ancestry; last column = all shared
  #   k=2: col 1=pop1, col 2=pop2, col 3=shared
  #   k=3: col 1=pop1, col 2=pop2, col 3=pop3, col 4=pop1_pop2, col 5=pop1_pop3, col 6=pop2_pop3, col 7=shared
  pip_per_pop <- sapply(seq_len(K), function(i) {
    mesusie_res$pip_config[, i]
  })
  colnames(pip_per_pop) <- pop_names

  pip_shared <- mesusie_res$pip_config[, ncol(mesusie_res$pip_config)]

  cs_either  <- get_cs_index_vector(mesusie_res$cs, length(mesusie_res$pip), renumber = TRUE)
  cs_per_pop <- sapply(pop_names, function(pn) .get_pop_cs_vec(mesusie_res, pn))

  mesusie_df <- data.frame(
    SNP        = common_snps,
    PIP_Either = mesusie_res$pip,
    PIP_Shared = pip_shared,
    pip_per_pop,
    CS         = cs_either,
    cs_per_pop,
    stringsAsFactors = FALSE
  )
  colnames(mesusie_df) <- c("SNP", "PIP_Either", "PIP_Shared",
                            pip_cols, "CS", cs_cols)

  mf_result <- .merge_snp_coordinates(gwas_list)
  mf_result <- dplyr::left_join(mf_result, mesusie_df, by = "SNP")
  mf_result <- .replace_missing_except_snp(mf_result)

  select_cols <- c("SNP", "CHR", "POS", "PIP_Either", "PIP_Shared",
                   pip_cols, "CS", cs_cols)
  mf_result <- dplyr::select(mf_result, dplyr::all_of(select_cols))

  .finalize_mf_result(mf_result, pop_names)
}

#' Run decision-guided multi-ancestry fine-mapping with raw outputs (K ancestries)
#'
#' This version returns the model choice, harmonized results, and raw fitted
#' objects from either SuSiE or MESuSiE.
#'
#' @param gwas_list Named list of K GWAS summary statistics data frames.
#'   Each must contain `SNP`, `CHR`, `POS`, `Beta`, `Se`, `Z`, `N`, and `PVAL`.
#' @param ld_list Named list of K LD matrices corresponding to `gwas_list`.
#' @param pop_names Character vector of length K giving ancestry labels.
#'   Defaults to `names(gwas_list)`.
#' @param L Number of single-effect components.
#' @param prior_weights Optional prior weights passed to SuSiE or MESuSiE.
#' @param ancestry_weight Optional ancestry weights passed to MESuSiE.
#' @param p_thresh P-value threshold used in the decision rule.
#' @param r2_thresh LD threshold used in the decision rule.
#'
#' @return A list with elements `decision`, `results`, and `raw_objects`.
#' @export
run_mf_decision <- function(
    gwas_list,
    ld_list,
    pop_names = names(gwas_list),
    L = 10,
    prior_weights = NULL,
    ancestry_weight = NULL,
    p_thresh = 5e-8,
    r2_thresh = 0.6
) {
  K <- length(gwas_list)
  stopifnot(is.list(gwas_list), is.list(ld_list), length(ld_list) == K, K >= 2)
  if (is.null(pop_names)) pop_names <- paste0("Pop", seq_len(K))
  .check_pop_names(pop_names)

  gwas_list <- stats::setNames(
    lapply(seq_len(K), function(i) .standardize_gwas(gwas_list[[i]], paste0("gwas_", i))),
    pop_names
  )
  ld_list <- stats::setNames(
    lapply(seq_len(K), function(i) .set_ld_dimnames(ld_list[[i]], gwas_list[[i]]$SNP)),
    pop_names
  )

  decision <- decide_finemapping_method(
    gwas_list = gwas_list,
    ld_list   = ld_list,
    p_thresh  = p_thresh,
    r2_thresh = r2_thresh
  )

  raw_objects <- list()

  pip_cols <- paste0("PIP_", pop_names)
  cs_cols  <- paste0("CS_", pop_names)

  # ----- SuSiE post-hoc branch -----
  if (decision$method == "SuSiE post-hoc") {
    susie_list <- lapply(seq_len(K), function(i) {
      susieR::susie_rss(
        z = gwas_list[[i]]$Z,
        R = ld_list[[i]],
        L = L,
        n = stats::median(gwas_list[[i]]$N),
        prior_weights = prior_weights,
        check_prior = FALSE
      )
    })
    names(susie_list) <- pop_names
    raw_objects$susie_list <- susie_list

    # Per-ancestry PIP and CS
    cs_list <- lapply(seq_len(K), function(i) {
      get_cs_index_vector(susie_list[[i]]$sets, nrow(gwas_list[[i]]))
    })

    res_dfs <- lapply(seq_len(K), function(i) {
      df <- data.frame(
        SNP = gwas_list[[i]]$SNP,
        PIP = susie_list[[i]]$pip,
        CS  = cs_list[[i]],
        stringsAsFactors = FALSE
      )
      colnames(df) <- c("SNP", pip_cols[i], cs_cols[i])
      df
    })

    # Merge SNP coordinates across all k GWAS
    temp_res <- .merge_snp_coordinates(gwas_list)
    for (i in seq_len(K)) {
      temp_res <- dplyr::left_join(temp_res, res_dfs[[i]], by = "SNP")
    }

    all_pip_cs <- c(pip_cols, cs_cols)
    temp_res <- dplyr::mutate(
      temp_res,
      dplyr::across(dplyr::all_of(all_pip_cs), ~ tidyr::replace_na(., 0))
    )

    # PIP_Either = max across ancestries
    temp_res$PIP_Either <- do.call(pmax, temp_res[pip_cols])

    # Cluster CS SNPs using max |LD| across all ancestries
    cs_snps_idx <- which(rowSums(temp_res[cs_cols] > 0) > 0)

    if (length(cs_snps_idx) > 0) {
      active_snps <- temp_res$SNP[cs_snps_idx]

      # Element-wise max of |LD| across k ancestries
      sub_lds <- lapply(seq_len(K), function(i) {
        idx <- match(active_snps, colnames(ld_list[[i]]))
        m <- ld_list[[i]][idx, idx, drop = FALSE]
        m[is.na(m)] <- 0
        abs(m)
      })
      ld_either <- Reduce(pmax, sub_lds)

      if (length(active_snps) > 1) {
        h_dynamic <- 1 - sqrt(r2_thresh)
        dist_mat  <- stats::as.dist(1 - ld_either)
        hc        <- stats::hclust(dist_mat, method = "complete")
        clusters  <- stats::cutree(hc, h = h_dynamic)
      } else {
        clusters <- 1L
      }

      temp_res$CS <- 0L
      temp_res$CS[cs_snps_idx] <- clusters
    } else {
      temp_res$CS <- 0L
    }

    mf_result <- temp_res
    mf_result$PIP_Shared <- do.call(pmin, mf_result[pip_cols])

    select_cols <- c("SNP", "CHR", "POS", "PIP_Either", "PIP_Shared",
                     pip_cols, "CS", cs_cols)
    mf_result <- dplyr::select(mf_result, dplyr::all_of(select_cols))

    return(list(
      decision    = decision,
      results     = .finalize_mf_result(mf_result, pop_names),
      raw_objects = raw_objects
    ))
  }

  # ----- MESuSiE branch -----
  common_snps <- Reduce(intersect, lapply(gwas_list, function(g) g$SNP))
  if (length(common_snps) == 0) {
    stop("No shared SNPs found across the GWAS inputs.", call. = FALSE)
  }

  g_subs <- lapply(gwas_list, function(g) {
    g[SNP %in% common_snps][order(match(SNP, common_snps))]
  })
  ld_subs <- lapply(ld_list, function(ld) {
    ld[common_snps, common_snps, drop = FALSE]
  })

  # Harmonize alleles to first ancestry, drop inconsistent SNPs
  harm <- .harmonize_alleles(g_subs, ld_subs, common_snps)
  g_subs      <- harm$g_subs
  ld_subs     <- harm$ld_subs
  common_snps <- harm$common_snps

  summary_stat_sd_list <- stats::setNames(lapply(g_subs, as.data.frame), pop_names)
  R_mat_list           <- stats::setNames(ld_subs, pop_names)

  mesusie_res <- MESuSiE::meSuSie_core(
    R_mat_list        = R_mat_list,
    summary_stat_list = summary_stat_sd_list,
    L = L,
    prior_weights     = prior_weights,
    ancestry_weight   = ancestry_weight
  )
  raw_objects$mesusie_res <- mesusie_res

  # Extract per-ancestry PIPs from pip_config directly
  # pip_config columns use combinatorial (combn) ordering:
  #   columns 1..K = single-ancestry; last column = all shared
  pip_per_pop <- sapply(seq_len(K), function(i) {
    mesusie_res$pip_config[, i]
  })
  colnames(pip_per_pop) <- pop_names

  pip_shared <- mesusie_res$pip_config[, ncol(mesusie_res$pip_config)]

  # CS extraction
  cs_either  <- get_cs_index_vector(mesusie_res$cs, length(mesusie_res$pip), renumber = TRUE)
  cs_per_pop <- sapply(pop_names, function(pn) .get_pop_cs_vec(mesusie_res, pn))

  mesusie_df <- data.frame(
    SNP        = common_snps,
    PIP_Either = mesusie_res$pip,
    PIP_Shared = pip_shared,
    pip_per_pop,
    CS         = cs_either,
    cs_per_pop,
    stringsAsFactors = FALSE
  )
  colnames(mesusie_df) <- c("SNP", "PIP_Either", "PIP_Shared",
                            pip_cols, "CS", cs_cols)

  mf_result <- .merge_snp_coordinates(gwas_list)
  mf_result <- dplyr::left_join(mf_result, mesusie_df, by = "SNP")
  mf_result <- .replace_missing_except_snp(mf_result)

  select_cols <- c("SNP", "CHR", "POS", "PIP_Either", "PIP_Shared",
                   pip_cols, "CS", cs_cols)
  mf_result <- dplyr::select(mf_result, dplyr::all_of(select_cols))

  list(
    decision    = decision,
    results     = .finalize_mf_result(mf_result, pop_names),
    raw_objects = raw_objects
  )
}

# -------------------------------------------------------------------------
# Backward-compatible wrappers for the old 2-ancestry signature
# -------------------------------------------------------------------------

#' Run decision-guided fine-mapping (2-ancestry wrapper)
#'
#' Thin wrapper around [run_mf_decision()] that accepts the original
#' positional `gwas_1, gwas_2, ld_1, ld_2` signature.
#'
#' @param gwas_1 GWAS summary statistics for ancestry 1.
#' @param gwas_2 GWAS summary statistics for ancestry 2.
#' @param ld_1 LD matrix matched to `gwas_1`.
#' @param ld_2 LD matrix matched to `gwas_2`.
#' @param pop_names Character vector of length 2 giving ancestry labels.
#' @param ... Additional arguments passed to [run_mf_decision()].
#'
#' @return See [run_mf_decision()].
#' @export
run_mf_decision_2pop <- function(
    gwas_1, gwas_2, ld_1, ld_2,
    pop_names = c("Pop1", "Pop2"),
    ...
) {
  run_mf_decision(
    gwas_list = stats::setNames(list(gwas_1, gwas_2), pop_names),
    ld_list   = stats::setNames(list(ld_1, ld_2), pop_names),
    pop_names = pop_names,
    ...
  )
}

#' Run decision-guided fine-mapping simplified (2-ancestry wrapper)
#'
#' Thin wrapper around [run_mf_decision_fm()] that accepts the original
#' positional `gwas_1, gwas_2, ld_1, ld_2` signature.
#'
#' @param gwas_1 GWAS summary statistics for ancestry 1.
#' @param gwas_2 GWAS summary statistics for ancestry 2.
#' @param ld_1 LD matrix matched to `gwas_1`.
#' @param ld_2 LD matrix matched to `gwas_2`.
#' @param pop_names Character vector of length 2 giving ancestry labels.
#'
#' @return See [run_mf_decision_fm()].
#' @export
run_mf_decision_fm_2pop <- function(
    gwas_1, gwas_2, ld_1, ld_2,
    pop_names = c("Pop1", "Pop2")
) {
  run_mf_decision_fm(
    gwas_list = stats::setNames(list(gwas_1, gwas_2), pop_names),
    ld_list   = stats::setNames(list(ld_1, ld_2), pop_names),
    pop_names = pop_names
  )
}
