#!/usr/bin/env Rscript
# TEST VERSION - v3.0 Distance Analysis
# Tests on ONE event pair and ONE structure

library(tidyverse)

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   TEST: v3.0 DISTANCE ANALYSIS                                ║\n")
cat("║   Testing on A5SS_gain + A3SS_gain, 5-exon structures         ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

output_dir <- "/Users/petecastaldi/claude_projects/nmd/results/isoform_transitions/v3.0_reference_based"

# ============================================================================
# Load Data
# ============================================================================

cat("Loading v3.0 filtered event data...\n")
transitions <- readRDS(file.path(output_dir, "reference_event_vectors_v3.0_filtered.rds"))
cat("Loaded:", nrow(transitions), "transitions\n\n")

# ============================================================================
# Helper Functions (same as main script)
# ============================================================================

classify_topology_a5ss_a3ss <- function(a5ss_exon, a3ss_exon) {
  distance <- abs(a5ss_exon - a3ss_exon)
  if (distance == 0) {
    return("same_exon")
  } else if (a3ss_exon < a5ss_exon) {
    return(paste0("dist_", distance, "_f2f"))
  } else {
    return(paste0("dist_", distance, "_b2b"))
  }
}

classify_topology_with_se <- function(exon_a, exon_b) {
  distance <- abs(exon_a - exon_b)
  if (distance == 0) return("same_exon")
  else return(paste0("dist_", distance))
}

extract_event_positions <- function(event_vector, event_type) {
  if (is.null(event_vector) || length(event_vector) == 0) {
    return(integer(0))
  }
  events <- event_vector %>% filter(event_type == !!event_type)
  if (nrow(events) == 0) return(integer(0))
  return(events$exon_number)
}

get_valid_positions <- function(event_type, n_exons) {
  base_type <- str_remove(event_type, "_(gain|loss)$")
  if (base_type == "A5SS") {
    # A5SS at downstream boundaries: all exons EXCEPT last (E1, E2, ..., En-1)
    return(1:(n_exons - 1))
  } else if (base_type == "A3SS") {
    # A3SS at upstream boundaries: all exons EXCEPT first (E2, E3, ..., En)
    return(2:n_exons)
  } else if (base_type == "SE") {
    # SE at internal exons only (E2, ..., En-1)
    return(2:(n_exons - 1))
  } else {
    return(integer(0))
  }
}

# ============================================================================
# Test Parameters
# ============================================================================

TEST_EVENT_A <- "A5SS_gain"
TEST_EVENT_B <- "A3SS_gain"
TEST_N_EXONS <- 5
N_PERMUTATIONS <- 100  # Reduced for testing

cat("Test parameters:\n")
cat("  Event pair:", TEST_EVENT_A, "+", TEST_EVENT_B, "\n")
cat("  Structure:", TEST_N_EXONS, "exons\n")
cat("  Permutations:", N_PERMUTATIONS, "\n\n")

# ============================================================================
# Filter to Test Case
# ============================================================================

cat("Filtering data...\n")

# Filter to transitions with both events
has_A <- transitions[[paste0("has_", TEST_EVENT_A)]]
has_B <- transitions[[paste0("has_", TEST_EVENT_B)]]

trans_subset <- transitions %>%
  filter(has_A & has_B)

cat("  Transitions with both events:", nrow(trans_subset), "\n")

# Add n_exons (CORRECTED: exon_number means different things for different events)
trans_subset <- trans_subset %>%
  mutate(n_exons = map_int(event_vector, function(ev) {
    if (is.null(ev) || nrow(ev) == 0) return(NA_integer_)

    # If Alt_TES exists, use its exon_number (actual last exon)
    tes_events <- ev %>% filter(event_type == "Alt_TES")
    if (nrow(tes_events) > 0) {
      return(as.integer(tes_events$exon_number[1]))
    }

    # Otherwise, use max junction index + 1
    junction_events <- ev %>% filter(event_type %in% c("A5SS", "A3SS", "CONST", "SE"))
    if (nrow(junction_events) > 0) {
      max_junction <- max(junction_events$exon_number, na.rm = TRUE)
      if (is.infinite(max_junction)) return(NA_integer_)
      return(as.integer(max_junction + 1))
    }

    return(1L)
  })) %>%
  filter(!is.na(n_exons))

cat("  After extracting n_exons:", nrow(trans_subset), "\n")

# Filter to test structure
obs_struct <- trans_subset %>%
  filter(n_exons == TEST_N_EXONS)

cat("  With", TEST_N_EXONS, "exons:", nrow(obs_struct), "\n\n")

if (nrow(obs_struct) == 0) {
  stop("No observations found for test case!")
}

# ============================================================================
# Extract Observed Topologies
# ============================================================================

cat("Extracting observed topologies...\n")

event_A_base <- str_remove(TEST_EVENT_A, "_(gain|loss)$")
event_B_base <- str_remove(TEST_EVENT_B, "_(gain|loss)$")

cat("  Event A base:", event_A_base, "\n")
cat("  Event B base:", event_B_base, "\n\n")

# Check first few event vectors
cat("Sample event vectors (first 3 transitions):\n")
for (i in 1:min(3, nrow(obs_struct))) {
  ev <- obs_struct$event_vector[[i]]
  cat("  Transition", i, ":\n")
  if (!is.null(ev) && nrow(ev) > 0) {
    cat("    Columns:", paste(colnames(ev), collapse = ", "), "\n")
    cat("    Rows:", nrow(ev), "\n")
    # Show event types present
    event_types <- unique(ev$event_type)
    cat("    Event types:", paste(event_types, collapse = ", "), "\n")
    # Show exon numbers
    if ("exon_number" %in% colnames(ev)) {
      cat("    Exon numbers:", paste(unique(ev$exon_number), collapse = ", "), "\n")
    }
  } else {
    cat("    Empty or NULL\n")
  }
}
cat("\n")

# Extract all observed topologies
observed_topologies <- map(obs_struct$event_vector, function(ev) {
  pos_A <- extract_event_positions(ev, event_A_base)
  pos_B <- extract_event_positions(ev, event_B_base)

  if (length(pos_A) == 0 || length(pos_B) == 0) {
    return(character(0))
  }

  n_combinations <- length(pos_A) * length(pos_B)
  topologies <- character(n_combinations)
  idx <- 1

  for (i in seq_along(pos_A)) {
    for (j in seq_along(pos_B)) {
      if (event_A_base == "A5SS" && event_B_base == "A3SS") {
        topo <- classify_topology_a5ss_a3ss(pos_A[i], pos_B[j])
      } else if (event_A_base == "A3SS" && event_B_base == "A5SS") {
        topo <- classify_topology_a5ss_a3ss(pos_B[j], pos_A[i])
      } else {
        topo <- classify_topology_with_se(pos_A[i], pos_B[j])
      }
      topologies[idx] <- topo
      idx <- idx + 1
    }
  }

  return(topologies)
})

observed_topologies <- unlist(observed_topologies)

if (length(observed_topologies) == 0) {
  stop("No topologies extracted! Check event_vector structure.")
}

obs_counts <- table(observed_topologies)
n_obs_actual <- length(observed_topologies)

cat("Observed topologies:\n")
print(obs_counts)
cat("\nTotal observations:", n_obs_actual, "from", nrow(obs_struct), "transitions\n\n")

# ============================================================================
# Run Permutation Test
# ============================================================================

cat("Running permutation test (", N_PERMUTATIONS, "iterations)...\n")

valid_A <- get_valid_positions(TEST_EVENT_A, TEST_N_EXONS)
valid_B <- get_valid_positions(TEST_EVENT_B, TEST_N_EXONS)

cat("  Valid positions for", event_A_base, ":", paste(valid_A, collapse = ", "), "\n")
cat("  Valid positions for", event_B_base, ":", paste(valid_B, collapse = ", "), "\n")

# Determine which events are terminal
is_A_terminal <- length(valid_A) == 0
is_B_terminal <- length(valid_B) == 0

if (is_A_terminal && is_B_terminal) {
  stop("Cannot permute - both events are terminal")
}

cat("  Terminal events: A =", is_A_terminal, ", B =", is_B_terminal, "\n\n")

null_distribution <- replicate(N_PERMUTATIONS, {
  perm_topologies <- map(obs_struct$event_vector, function(ev) {
    pos_A <- extract_event_positions(ev, event_A_base)
    pos_B <- extract_event_positions(ev, event_B_base)

    if (length(pos_A) == 0 || length(pos_B) == 0) {
      return(character(0))
    }

    # For terminal events, use observed positions (fixed)
    # For internal events, sample random positions (permuted)
    if (is_A_terminal) {
      perm_pos_A <- pos_A  # Keep terminal event fixed
    } else {
      perm_pos_A <- sample(valid_A, length(pos_A), replace = TRUE)
    }

    if (is_B_terminal) {
      perm_pos_B <- pos_B  # Keep terminal event fixed
    } else {
      perm_pos_B <- sample(valid_B, length(pos_B), replace = TRUE)
    }

    topologies <- character(length(perm_pos_A) * length(perm_pos_B))
    idx <- 1
    for (i in seq_along(perm_pos_A)) {
      for (j in seq_along(perm_pos_B)) {
        if (event_A_base == "A5SS" && event_B_base == "A3SS") {
          topo <- classify_topology_a5ss_a3ss(perm_pos_A[i], perm_pos_B[j])
        } else if (event_A_base == "A3SS" && event_B_base == "A5SS") {
          topo <- classify_topology_a5ss_a3ss(perm_pos_B[j], perm_pos_A[i])
        } else {
          topo <- classify_topology_with_se(perm_pos_A[i], perm_pos_B[j])
        }
        topologies[idx] <- topo
        idx <- idx + 1
      }
    }
    return(topologies)
  })

  table(unlist(perm_topologies))
}, simplify = FALSE)

cat("Permutation test complete.\n\n")

# ============================================================================
# Calculate Statistics
# ============================================================================

cat("Calculating statistics...\n\n")

all_topology_names <- unique(c(names(obs_counts),
                               unlist(map(null_distribution, names))))

expected_counts <- map_dbl(all_topology_names, function(topo) {
  counts <- map_dbl(null_distribution, function(perm) {
    ifelse(topo %in% names(perm), perm[topo], 0)
  })
  mean(counts)
})
names(expected_counts) <- all_topology_names

# Calculate OR and p-values
results <- map_dfr(all_topology_names, function(topo) {
  obs_count <- ifelse(topo %in% names(obs_counts), obs_counts[topo], 0)
  exp_count <- expected_counts[topo]

  obs_other <- sum(obs_counts) - obs_count
  exp_other <- sum(expected_counts) - exp_count

  if (obs_other == 0 || exp_other == 0 || obs_count == 0 || exp_count == 0) {
    return(tibble(
      topology = topo,
      obs_count = obs_count,
      obs_pct = 100 * obs_count / sum(obs_counts),
      exp_count = exp_count,
      exp_pct = 100 * exp_count / sum(expected_counts),
      odds_ratio = NA_real_,
      p_value = NA_real_
    ))
  }

  contingency <- matrix(c(obs_count, obs_other, exp_count, exp_other), nrow = 2)

  fisher_result <- tryCatch({
    fisher.test(contingency)
  }, error = function(e) {
    list(estimate = NA_real_, p.value = NA_real_)
  })

  tibble(
    topology = topo,
    obs_count = obs_count,
    obs_pct = 100 * obs_count / sum(obs_counts),
    exp_count = exp_count,
    exp_pct = 100 * exp_count / sum(expected_counts),
    odds_ratio = as.numeric(fisher_result$estimate),
    p_value = fisher_result$p.value
  )
}) %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(desc(obs_count))

# ============================================================================
# Display Results
# ============================================================================

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   TEST RESULTS                                                 ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("Event Pair:", TEST_EVENT_A, "+", TEST_EVENT_B, "\n")
cat("Structure:", TEST_N_EXONS, "exons\n")
cat("Transitions:", nrow(obs_struct), "\n")
cat("Observations:", n_obs_actual, "\n")
cat("Permutations:", N_PERMUTATIONS, "\n\n")

cat("RESULTS:\n")
print(results, n = Inf)

cat("\n\n")
cat("INTERPRETATION:\n")
significant_enriched <- results %>% filter(p_value < 0.05, odds_ratio > 1)
significant_depleted <- results %>% filter(p_value < 0.05, odds_ratio < 1)

if (nrow(significant_enriched) > 0) {
  cat("ENRICHED topologies (OR > 1, p < 0.05):\n")
  for (i in seq_len(nrow(significant_enriched))) {
    cat(sprintf("  - %s: OR=%.2f, p=%.4f %s\n",
                significant_enriched$topology[i],
                significant_enriched$odds_ratio[i],
                significant_enriched$p_value[i],
                significant_enriched$significance[i]))
  }
  cat("\n")
}

if (nrow(significant_depleted) > 0) {
  cat("DEPLETED topologies (OR < 1, p < 0.05):\n")
  for (i in seq_len(nrow(significant_depleted))) {
    cat(sprintf("  - %s: OR=%.2f, p=%.4f %s\n",
                significant_depleted$topology[i],
                significant_depleted$odds_ratio[i],
                significant_depleted$p_value[i],
                significant_depleted$significance[i]))
  }
  cat("\n")
}

if (nrow(significant_enriched) == 0 && nrow(significant_depleted) == 0) {
  cat("No significant enrichments or depletions found.\n\n")
}

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   TEST COMPLETE                                                ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
