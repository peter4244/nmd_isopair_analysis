# Adaptive Permutation Testing Function
# Progressively increases permutations for significant results to enable precise FDR calculation

adaptive_permutation_test <- function(obs_topologies,
                                      perm_function,  # Function to generate one permutation
                                      initial_n = 100,
                                      stages = list(
                                        list(threshold = 0.1, max_perms = 100),      # Stage 1: screening
                                        list(threshold = 0.01, max_perms = 1000),    # Stage 2: significant
                                        list(threshold = 0.001, max_perms = 10000),  # Stage 3: highly significant
                                        list(threshold = 0, max_perms = 100000)      # Stage 4: very highly sig (max)
                                      )) {

  obs_counts <- table(obs_topologies)
  all_topology_names <- names(obs_counts)

  # Storage for permutation results
  perm_counts_list <- list()
  n_perms_done <- 0

  # Progressive stages
  for (stage in stages) {
    n_to_run <- stage$max_perms - n_perms_done

    if (n_to_run > 0) {
      # Run additional permutations
      cat(sprintf("  Running %d permutations (total: %d)...\n", n_to_run, stage$max_perms))

      new_perms <- replicate(n_to_run, {
        perm_topologies <- perm_function()
        table(perm_topologies)
      }, simplify = FALSE)

      perm_counts_list <- c(perm_counts_list, new_perms)
      n_perms_done <- stage$max_perms

      # Calculate current p-values
      results <- calculate_pvalues(obs_counts, perm_counts_list, all_topology_names)

      # Check if we can stop early
      min_p <- min(results$p_value, na.rm = TRUE)

      cat(sprintf("    Min p-value: %.6f (after %d permutations)\n", min_p, n_perms_done))

      # If minimum p-value is above threshold, stop
      if (min_p > stage$threshold) {
        cat(sprintf("    All p-values > %.2f, stopping at %d permutations\n",
                   stage$threshold, n_perms_done))
        break
      }
    }
  }

  # Final results
  cat(sprintf("  Complete. %d permutations total.\n", n_perms_done))

  results$n_permutations <- n_perms_done
  return(results)
}

# Helper function to calculate p-values from permutation counts
calculate_pvalues <- function(obs_counts, perm_counts_list, all_topology_names) {

  # Calculate expected counts from permutations
  expected_counts <- map_dbl(all_topology_names, function(topo) {
    counts <- map_dbl(perm_counts_list, function(perm) {
      ifelse(topo %in% names(perm), perm[topo], 0)
    })
    mean(counts)
  })
  names(expected_counts) <- all_topology_names

  # Calculate empirical permutation p-values for each topology
  results <- map_dfr(all_topology_names, function(topo) {
    obs_count <- ifelse(topo %in% names(obs_counts), obs_counts[topo], 0)
    exp_count <- expected_counts[topo]

    # Extract counts for this topology from all permutations
    perm_counts_for_topo <- map_dbl(perm_counts_list, function(perm) {
      ifelse(topo %in% names(perm), perm[topo], 0)
    })

    n_perms <- length(perm_counts_for_topo)

    # Empirical permutation p-value
    # Enrichment: proportion of permutations with count >= observed (ties count as extreme)
    n_greater_or_equal <- sum(perm_counts_for_topo >= obs_count)
    p_enrichment <- n_greater_or_equal / n_perms

    # Depletion: proportion of permutations with count <= observed
    n_less_or_equal <- sum(perm_counts_for_topo <= obs_count)
    p_depletion <- n_less_or_equal / n_perms

    # Two-tailed p-value
    p_value <- 2 * min(p_enrichment, p_depletion)
    p_value <- min(p_value, 1.0)  # Cap at 1.0

    # Odds ratio: observed odds / expected odds
    obs_other <- sum(obs_counts) - obs_count
    exp_other <- sum(expected_counts) - exp_count

    if (obs_count == 0 || exp_count == 0 || obs_other == 0 || exp_other == 0) {
      odds_ratio <- NA_real_
    } else {
      odds_ratio <- (obs_count / obs_other) / (exp_count / exp_other)
    }

    tibble(
      topology = topo,
      obs_count = obs_count,
      obs_pct = 100 * obs_count / sum(obs_counts),
      exp_count = exp_count,
      exp_pct = 100 * exp_count / sum(expected_counts),
      odds_ratio = odds_ratio,
      p_value = p_value
    )
  })

  return(results)
}
