rankings <- data.frame(
  ligand_id = rep(c("A", "B", "C"), 2L),
  method_id = rep(c("first", "second"), each = 3L),
  score = c(-9, -8, -7, -7, -9, -8)
)

test_that("method-local ranks underpin consensus", {
  ranked <- rank_consensus_inputs(rankings)
  expect_equal(ranked$rank_within_method, c(1, 2, 3, 3, 1, 2))
  result <- summarise_consensus(rankings)
  expect_identical(result$ligand_id, c("B", "A", "C"))
  expect_equal(result$methods_observed, rep(2L, 3L))
  expect_match(result$interpretation[1], "consensus")
})

test_that("consensus preserves missing-method coverage", {
  incomplete <- rankings[-6, ]
  result <- summarise_consensus(incomplete)
  c_row <- result[result$ligand_id == "C", ]
  expect_equal(c_row$methods_observed, 1L)
  expect_equal(c_row$methods_total, 2L)
  expect_equal(c_row$method_coverage, 0.5)
})

test_that("pairwise agreement uses shared ligands only", {
  agreement <- summarise_rank_agreement(rankings)
  expect_equal(agreement$ligands_shared, 3L)
  expect_equal(agreement$spearman_rho, -0.5)
  expect_error(summarise_rank_agreement(rankings[rankings$method_id == "first", ]), "two methods")
})

test_that("invalid inputs fail clearly", {
  duplicate <- rbind(rankings, rankings[1, ])
  expect_error(validate_consensus_inputs(duplicate), "must occur once")
  invalid <- rankings
  invalid$score[1] <- NA_real_
  expect_error(validate_consensus_inputs(invalid), "finite")
})
