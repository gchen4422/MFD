# MFD

**MFD** (Multi-ancestry Fine-mapping Decision) is an R package that adaptively selects between joint and post-hoc multi-ancestry fine-mapping strategies based on locus-level genetic architecture.

> Analysis code for the benchmarking study is available at [gchen4422/MFD_analysis](https://github.com/gchen4422/MFD_analysis/tree/main).

> **Detailed tutorials** (installation, data preparation, worked examples, and visualization) will be available at the MFD website (coming soon).

## How It Works

MFD chooses between two fine-mapping paths using a three-step decision rule:

1. **No significant ancestry-specific variants (AS-Vs)?** → Use **MESuSiE** (joint modeling leverages shared LD structure for higher resolution)
2. **Significant AS-Vs with low LD to shared signals (r² < 0.6)?** → Use **SuSiE post-hoc** (runs SuSiE-RSS independently per ancestry, then merges with a consensus LD matrix)
3. **Significant AS-Vs in high LD with shared signals?** → Use **MESuSiE** if shared signal is genome-wide significant in all ancestries; otherwise **SuSiE post-hoc**

![MFD flowchart](MFD_flowchart.png)

## Installation

```r
# Install dependencies first
devtools::install_github("borangao/MESuSiE")

# Install MFD
devtools::install_github("gchen4422/MFD")
```

## Quick Start

```r
library(MFD)

# gwas_1, gwas_2: data frames with columns SNP, CHR, POS, Z, Beta, Se, N
# ld_1, ld_2:     LD correlation matrices with row/colnames matching SNP IDs

result <- run_mf_decision(
  gwas_1, gwas_2,
  ld_1,   ld_2,
  pop_names = c("EUR", "AFR")
)

result$decision   # which method was used and why
result$results    # per-SNP PIPs and credible sets
result$raw_objects  # raw SuSiE or MESuSiE model objects
```

### Input format (`gwas_1` / `gwas_2`)

| Column | Required | Description |
|--------|----------|-------------|
| `SNP`  | Yes | Variant ID (rsID or chr:pos); must match LD matrix names |
| `CHR`  | Yes | Chromosome |
| `POS`  | Yes | Base-pair position |
| `Z`    | Yes | Z-score (Beta / Se) |
| `Beta` | Optional | Effect size |
| `Se`   | Optional | Standard error |
| `N`    | Optional | Sample size |

### Output (`result$results`)

| Column | Description |
|--------|-------------|
| `PIP_Either` | Probability of being causal in **at least one** ancestry (primary discovery metric) |
| `PIP_Shared` | Probability of being causal in **both** ancestries |
| `PIP_Ancestry_1/2` | Ancestry-specific causal probability |
| `CS` | "Either" 95% credible set membership (0 = not in any CS) |
| `CS_Ancestry_1/2` | Ancestry-specific credible set membership |

A common threshold to declare a fine-mapped signal is `PIP_Either > 0.5`.

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `L` | `10` | Maximum number of causal effects per region |
| `p_thresh` | `5e-8` | P-value threshold to define significant AS-Vs and shared variants |
| `r2_thresh` | `0.6` | LD threshold for the decision rule and post-hoc CS merging |
| `prior_weights` | `NULL` | Per-SNP prior probability (e.g., from functional annotations) |
| `ancestry_weight` | `NULL` | Ancestry weights passed to MESuSiE |
| `pop_names` | `c("Pop1","Pop2")` | Labels for the two ancestries |

## Dependencies

`MESuSiE`, `susieR`, `data.table`, `dplyr`, `tidyr`
