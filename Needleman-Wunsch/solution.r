needleman_wunsch <- function(seq1, seq2, match = 1, mismatch = -1, gap = -1) {
  len1 <- nchar(seq1)
  len2 <- nchar(seq2)

  seq1_chars <- c("-", unlist(strsplit(seq1, split = "")))
  seq2_chars <- c("-", unlist(strsplit(seq2, split = "")))

  scoring_matrix <- matrix(0, nrow = len1 + 1, ncol = len2 + 1)
  rownames(scoring_matrix) <- seq1_chars
  colnames(scoring_matrix) <- seq2_chars

  for (i in 1:(len1 + 1)) {
    scoring_matrix[i, 1] <- (i - 1) * gap
  }
  for (j in 1:(len2 + 1)) {
    scoring_matrix[1, j] <- (j - 1) * gap
  }

  for (i in 2:(len1 + 1)) {
    for (j in 2:(len2 + 1)) {
      char1 <- seq1_chars[i]
      char2 <- seq2_chars[j]
      match_score <- ifelse(char1 == char2, match, mismatch)
      score_diag <- scoring_matrix[i - 1, j - 1] + match_score
      score_up   <- scoring_matrix[i - 1, j] + gap
      score_left <- scoring_matrix[i, j - 1] + gap

      scoring_matrix[i, j] <- max(score_diag, score_up, score_left)
    }
  }

  return(scoring_matrix)
}

print(needleman_wunsch("GATTACA", "GCATGCU"))