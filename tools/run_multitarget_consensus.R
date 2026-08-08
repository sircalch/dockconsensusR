#!/usr/bin/env Rscript

# Usage: Rscript tools/run_multitarget_consensus.R <derived-data-dir> <output-dir>
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: Rscript tools/run_multitarget_consensus.R <derived-data-dir> <output-dir>", call. = FALSE)
}

script_argument <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
if (length(script_argument) != 1L) stop("Could not determine the path of this script.", call. = FALSE)
project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_argument)), ".."), mustWork = TRUE)
source(file.path(project_root, "R", "consensus.R"), local = TRUE)

derived_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
output_dir <- arguments[[2L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

targets <- c("ada", "ampc", "comt")
agreement_rows <- list()
consensus_rows <- list()

for (target in targets) {
  input_path <- file.path(derived_dir, paste0("dude-", target), "pilot_runs_normalized.csv")
  runs <- utils::read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
  input <- runs[, c("ligand_id", "scenario_id", "score"), drop = FALSE]
  names(input)[names(input) == "scenario_id"] <- "method_id"
  validate_consensus_inputs(input)
  if (length(unique(input$method_id)) != 3L || any(table(input$method_id) != 24L)) {
    stop("Expected three complete 24-ligand seed rankings for ", target, ".", call. = FALSE)
  }

  agreement <- summarise_rank_agreement(input, score_direction = "lower")
  agreement$target <- target
  agreement_rows[[length(agreement_rows) + 1L]] <- agreement[, c(
    "target", "method_first", "method_second", "ligands_shared", "spearman_rho"
  )]

  consensus <- summarise_consensus(input, score_direction = "lower")
  consensus$target <- target
  consensus_rows[[length(consensus_rows) + 1L]] <- consensus[, c(
    "target", "ligand_id", "methods_observed", "methods_total", "method_coverage",
    "mean_rank", "mean_normalized_rank", "borda_score", "interpretation"
  )]
}

agreements <- do.call(rbind, agreement_rows)
consensus <- do.call(rbind, consensus_rows)
if (!all(consensus$method_coverage == 1)) stop("An input ligand lacked full seed coverage.", call. = FALSE)
utils::write.csv(agreements, file.path(output_dir, "multitarget_rank_agreement.csv"), row.names = FALSE)
utils::write.csv(consensus, file.path(output_dir, "multitarget_consensus.csv"), row.names = FALSE)
message("Wrote consensus summaries to: ", normalizePath(output_dir))
