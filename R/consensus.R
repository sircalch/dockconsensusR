#' Validate precomputed rankings for consensus analysis
#'
#' @param data A data frame with `ligand_id`, `method_id`, and numeric `score`.
#'
#' @return Invisibly returns `data` when it satisfies the input contract.
#' @export
validate_consensus_inputs <- function(data) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  required <- c("ligand_id", "method_id", "score")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("`data` is missing required column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(data) == 0L) stop("`data` must contain at least one result.", call. = FALSE)
  if (anyNA(data$ligand_id) || anyNA(data$method_id) ||
      any(!nzchar(as.character(data$ligand_id))) || any(!nzchar(as.character(data$method_id)))) {
    stop("`ligand_id` and `method_id` must be complete non-empty values.", call. = FALSE)
  }
  if (!is.numeric(data$score) || any(!is.finite(data$score))) {
    stop("`score` must be a finite numeric column.", call. = FALSE)
  }
  key <- paste(data$ligand_id, data$method_id, sep = "\r")
  if (anyDuplicated(key)) {
    stop("Each `ligand_id` and `method_id` pair must occur once.", call. = FALSE)
  }
  invisible(data)
}

#' Rank inputs within each method
#'
#' @param data A validated consensus-input table.
#' @param score_direction Whether lower or higher scores are preferred.
#'
#' @return A copy of `data` with `rank_within_method` and `method_size`.
#' @export
rank_consensus_inputs <- function(data, score_direction = c("lower", "higher")) {
  validate_consensus_inputs(data)
  score_direction <- match.arg(score_direction)
  ranking_score <- if (identical(score_direction, "lower")) data$score else -data$score
  method <- as.character(data$method_id)
  data$rank_within_method <- stats::ave(ranking_score, method, FUN = function(x) rank(x, ties.method = "min"))
  data$method_size <- as.integer(stats::ave(ranking_score, method, FUN = length))
  data
}

#' Summarise a rank consensus across methods
#'
#' Scores are never pooled. Each method contributes a within-method rank, then
#' the returned mean rank and Borda score are calculated only over the methods
#' that actually evaluated that ligand.
#'
#' @param data A validated consensus-input table.
#' @param score_direction Whether lower or higher scores are preferred.
#'
#' @return One row per ligand with method coverage, mean rank, normalized mean
#'   rank, and Borda score. These are consensus summaries only.
#' @export
summarise_consensus <- function(data, score_direction = c("lower", "higher")) {
  ranked <- rank_consensus_inputs(data, match.arg(score_direction))
  all_methods <- length(unique(ranked$method_id))
  pieces <- split(ranked, as.character(ranked$ligand_id), drop = TRUE)
  output <- lapply(pieces, function(x) {
    normalized_rank <- x$rank_within_method / x$method_size
    data.frame(
      ligand_id = as.character(x$ligand_id[[1L]]),
      methods_observed = nrow(x),
      methods_total = all_methods,
      method_coverage = nrow(x) / all_methods,
      mean_rank = mean(x$rank_within_method),
      mean_normalized_rank = mean(normalized_rank),
      borda_score = sum((x$method_size + 1) - x$rank_within_method),
      interpretation = "Rank-consensus summary only; not binding or activity evidence.",
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, output)
  result[order(result$mean_normalized_rank, -result$method_coverage, -result$borda_score, result$ligand_id), , drop = FALSE]
}

#' Summarise pairwise rank agreement between methods
#'
#' @param data A validated consensus-input table.
#' @param score_direction Whether lower or higher scores are preferred.
#'
#' @return One row for each method pair, including shared ligand count and
#'   Spearman rank correlation. Missing ligands are not treated as poor ranks.
#' @export
summarise_rank_agreement <- function(data, score_direction = c("lower", "higher")) {
  ranked <- rank_consensus_inputs(data, match.arg(score_direction))
  methods <- sort(unique(as.character(ranked$method_id)))
  if (length(methods) < 2L) stop("At least two methods are required.", call. = FALSE)
  pairs <- utils::combn(methods, 2L, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    first <- ranked[ranked$method_id == pair[[1L]], c("ligand_id", "rank_within_method"), drop = FALSE]
    second <- ranked[ranked$method_id == pair[[2L]], c("ligand_id", "rank_within_method"), drop = FALSE]
    shared <- merge(first, second, by = "ligand_id", suffixes = c("_first", "_second"))
    rho <- if (nrow(shared) < 2L) NA_real_ else stats::cor(
      shared$rank_within_method_first, shared$rank_within_method_second, method = "spearman"
    )
    data.frame(
      method_first = pair[[1L]], method_second = pair[[2L]],
      ligands_shared = nrow(shared), spearman_rho = rho,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
