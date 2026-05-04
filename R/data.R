#' GWAS summary statistics for ancestry 1 (European) — example data
#'
#' Simulated GWAS summary statistics for a single genomic locus (chromosome 6)
#' representing the European ancestry cohort. Used to demonstrate
#' \code{\link{run_mf_decision}}.
#'
#' @format A data.table with 1,933 rows and 9 columns:
#' \describe{
#'   \item{SNP}{Variant ID in \code{chr:pos} format}
#'   \item{CHR}{Chromosome number}
#'   \item{POS}{Base-pair position (GRCh38)}
#'   \item{Signal}{Simulated causal signal indicator (0 = non-causal)}
#'   \item{Beta}{Effect size estimate}
#'   \item{Se}{Standard error of the effect size}
#'   \item{Z}{Z-score (Beta / Se)}
#'   \item{N}{Sample size}
#'   \item{PVAL}{Two-sided p-value}
#' }
#' @seealso \code{\link{summary_stat_2}}, \code{\link{susie_EU_cov}}, \code{\link{susie_BB_cov}}
#' @examples
#' data(summary_stat_1)
#' head(summary_stat_1)
"summary_stat_1"

#' GWAS summary statistics for ancestry 2 (African) — example data
#'
#' Simulated GWAS summary statistics for a single genomic locus (chromosome 6)
#' representing the African ancestry cohort. Used to demonstrate
#' \code{\link{run_mf_decision}}.
#'
#' @format A data.table with 2,762 rows and 9 columns:
#' \describe{
#'   \item{SNP}{Variant ID in \code{chr:pos} format}
#'   \item{CHR}{Chromosome number}
#'   \item{POS}{Base-pair position (GRCh38)}
#'   \item{Signal}{Simulated causal signal indicator (0 = non-causal)}
#'   \item{Beta}{Effect size estimate}
#'   \item{Se}{Standard error of the effect size}
#'   \item{Z}{Z-score (Beta / Se)}
#'   \item{N}{Sample size}
#'   \item{PVAL}{Two-sided p-value}
#' }
#' @seealso \code{\link{summary_stat_1}}, \code{\link{susie_EU_cov}}, \code{\link{susie_BB_cov}}
#' @examples
#' data(summary_stat_2)
#' head(summary_stat_2)
"summary_stat_2"

#' LD correlation matrix for ancestry 1 (European) — example data
#'
#' A 1,933 x 1,933 LD correlation matrix for the European ancestry cohort,
#' corresponding to the variants in \code{\link{summary_stat_1}}.
#' Row and column names are SNP IDs matching the \code{SNP} column.
#'
#' @format A numeric matrix with 1,933 rows and 1,933 columns.
#' @seealso \code{\link{summary_stat_1}}, \code{\link{susie_BB_cov}}
#' @examples
#' data(susie_EU_cov)
#' dim(susie_EU_cov)
"susie_EU_cov"

#' LD correlation matrix for ancestry 2 (African) — example data
#'
#' A 2,762 x 2,762 LD correlation matrix for the African ancestry cohort,
#' corresponding to the variants in \code{\link{summary_stat_2}}.
#' Row and column names are SNP IDs matching the \code{SNP} column.
#'
#' @format A numeric matrix with 2,762 rows and 2,762 columns.
#' @seealso \code{\link{summary_stat_2}}, \code{\link{susie_EU_cov}}
#' @examples
#' data(susie_BB_cov)
#' dim(susie_BB_cov)
"susie_BB_cov"
