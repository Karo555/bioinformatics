verify_sequence <- function(seq) {
  if (nchar(seq) == 0) {
    stop("Input sequence is empty.")
  }
  seq <- toupper(seq)
  if (!grepl("^[ACGTU]*$", seq)) {
    stop("Invalid sequence: only A, C, G, T, U characters are allowed.")
  }
  return(seq)
}

initialize_scoring_matrix <- function(seq1, seq2, gap = -1) {
  verify_sequence(seq1)
  verify_sequence(seq2)
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
  list(
    scoring_matrix = scoring_matrix,
    seq1_chars = seq1_chars,
    seq2_chars = seq2_chars
  )
}

compute_scoring_matrix <- function(
  scoring_matrix, seq1_chars, seq2_chars, match, mismatch, gap
) {
  len1 <- nrow(scoring_matrix) - 1
  len2 <- ncol(scoring_matrix) - 1

  for (i in 2:(len1 + 1)) {
    for (j in 2:(len2 + 1)) {
      char1 <- seq1_chars[i]
      char2 <- seq2_chars[j]
      score_diag <- scoring_matrix[i - 1, j - 1] + 
        ifelse(char1 == char2, match, mismatch)
      score_up <- scoring_matrix[i - 1, j] + gap
      score_left <- scoring_matrix[i, j - 1] + gap
      scoring_matrix[i, j] <- max(score_diag, score_up, score_left)
    }
  }
  scoring_matrix
}

traceback_alignment <- function(scoring_matrix, seq1_chars, seq2_chars, match, mismatch, gap) {
  i <- nrow(scoring_matrix)
  j <- ncol(scoring_matrix)
  align1 <- ""
  align2 <- ""

  while (i > 1 || j > 1) {
    current <- scoring_matrix[i, j]
    diag <- if (i > 1 && j > 1) scoring_matrix[i - 1, j - 1] else -Inf
    up <- if (i > 1) scoring_matrix[i - 1, j] else -Inf
    left <- if (j > 1) scoring_matrix[i, j - 1] else -Inf

    char1 <- seq1_chars[i]
    char2 <- seq2_chars[j]

    if (i > 1 && j > 1 &&
          current == diag + ifelse(char1 == char2, match, mismatch)) {
      align1 <- paste0(char1, align1)
      align2 <- paste0(char2, align2)
      i <- i - 1
      j <- j - 1
    } else if (i > 1 && current == up + gap) {
      align1 <- paste0(seq1_chars[i], align1)
      align2 <- paste0("-", align2)
      i <- i - 1
    } else {
      align1 <- paste0("-", align1)
      align2 <- paste0(seq2_chars[j], align2)
      j <- j - 1
    }
  }
  list(align1 = align1, align2 = align2)
}

compute_alignment_stats <- function(align1, align2) {
  alignment_length <- nchar(align1)
  identity <- 0
  gaps <- 0

  for (k in 1:alignment_length) {
    a1 <- substr(align1, k, k)
    a2 <- substr(align2, k, k)
    if (a1 == "-" || a2 == "-") {
      gaps <- gaps + 1
    } else if (a1 == a2) {
      identity <- identity + 1
    }
  }

  identity_pct <- round((identity / alignment_length) * 100, 2)
  gaps_pct <- round((gaps / alignment_length) * 100, 2)
  list(identity_pct = identity_pct, gaps_pct = gaps_pct)
}

needleman_wunsch <- function(seq1, seq2, match = 1, mismatch = -1, gap = -1) {
  init <- initialize_scoring_matrix(seq1, seq2, gap)
  scoring_matrix <- compute_scoring_matrix(
    init$scoring_matrix, init$seq1_chars, init$seq2_chars, match, mismatch, gap
  )
  traceback_result <- traceback_alignment(
    scoring_matrix, init$seq1_chars, init$seq2_chars, match, mismatch, gap
  )
  stats <- compute_alignment_stats(
    traceback_result$align1, traceback_result$align2
  )

  cat("\nScoring matrix:\n")
  print(scoring_matrix)
  cat("\nAlignment:\n")
  cat(traceback_result$align1, "\n")
  cat(traceback_result$align2, "\n")
  cat("\nAlignment stats:\n")
  cat("Identity:", stats$identity_pct, "%\n")
  cat("Gaps:", stats$gaps_pct, "%\n")
  alignment_score <- scoring_matrix[nrow(scoring_matrix), ncol(scoring_matrix)]
  cat("Alignment score:", alignment_score, "\n")
  }