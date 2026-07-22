# Exact Topological Test using Mathematical Null Distribution
# Replaces permutation testing with exact probability calculations

exact_topology_test <- function(obs_topologies,
                                obs_struct,
                                event_A_base,
                                event_B_base,
                                valid_A,
                                valid_B,
                                extract_event_positions_func,
                                classify_topology_func) {

  # Calculate observed counts
  obs_counts <- table(obs_topologies)
  n_obs_total <- length(obs_topologies)

  # Calculate exact expected distribution mathematically
  # For each transition, calculate probability of each topology under uniform null

  is_A_terminal <- length(valid_A) == 0
  is_B_terminal <- length(valid_B) == 0

  # Build expected distribution from probability calculations
  expected_topology_probs <- list()

  for (i in seq_len(nrow(obs_struct))) {
    ev <- obs_struct$event_vector[[i]]

    pos_A <- extract_event_positions_func(ev, event_A_base)
    pos_B <- extract_event_positions_func(ev, event_B_base)

    if (length(pos_A) == 0 || length(pos_B) == 0) next

    # For terminal events, positions are fixed (no probability distribution)
    # For internal events, uniform distribution across valid positions

    # Calculate probability for each possible topology configuration
    # under uniform distribution

    if (is_A_terminal && is_B_terminal) {
      # Both fixed - only one configuration possible
      # This matches observed (no randomness)
      topology_probs <- tibble(topology = character(0), prob = numeric(0))
      for (ia in seq_along(pos_A)) {
        for (ib in seq_along(pos_B)) {
          topo <- classify_topology_func(pos_A[ia], pos_B[ib])
          topology_probs <- bind_rows(topology_probs,
                                     tibble(topology = topo, prob = 1.0))
        }
      }
    } else if (is_A_terminal) {
      # A is fixed, B is uniform over valid_B
      topology_probs <- tibble(topology = character(0), prob = numeric(0))
      for (ia in seq_along(pos_A)) {
        for (valid_b in valid_B) {
          topo <- classify_topology_func(pos_A[ia], valid_b)
          # Probability: (1 / length(valid_B)) for each B position
          prob <- 1.0 / length(valid_B)
          topology_probs <- bind_rows(topology_probs,
                                     tibble(topology = topo, prob = prob))
        }
      }
    } else if (is_B_terminal) {
      # B is fixed, A is uniform over valid_A
      topology_probs <- tibble(topology = character(0), prob = numeric(0))
      for (ib in seq_along(pos_B)) {
        for (valid_a in valid_A) {
          topo <- classify_topology_func(valid_a, pos_B[ib])
          prob <- 1.0 / length(valid_A)
          topology_probs <- bind_rows(topology_probs,
                                     tibble(topology = topo, prob = prob))
        }
      }
    } else {
      # Both are uniform distributions
      topology_probs <- tibble(topology = character(0), prob = numeric(0))
      for (valid_a in valid_A) {
        for (valid_b in valid_B) {
          topo <- classify_topology_func(valid_a, valid_b)
          # Probability: (1/m_A) × (1/m_B)
          prob <- (1.0 / length(valid_A)) * (1.0 / length(valid_B))
          topology_probs <- bind_rows(topology_probs,
                                     tibble(topology = topo, prob = prob))
        }
      }
    }

    # Aggregate probabilities for this transition
    topology_probs_agg <- topology_probs %>%
      group_by(topology) %>%
      summarise(prob = sum(prob), .groups = "drop")

    expected_topology_probs[[i]] <- topology_probs_agg
  }

  # Sum expected probabilities across all transitions to get expected counts
  all_expected <- bind_rows(expected_topology_probs)
  expected_counts_df <- all_expected %>%
    group_by(topology) %>%
    summarise(expected_count = sum(prob), .groups = "drop")

  expected_counts <- setNames(expected_counts_df$expected_count,
                              expected_counts_df$topology)

  # Get all topology names
  all_topology_names <- unique(c(names(obs_counts), names(expected_counts)))

  # Run binomial test for each topology
  results <- map_dfr(all_topology_names, function(topo) {
    obs_count <- ifelse(topo %in% names(obs_counts), obs_counts[topo], 0)
    exp_count <- ifelse(topo %in% names(expected_counts), expected_counts[topo], 0)

    # Expected proportion under null
    exp_prop <- exp_count / sum(expected_counts)

    # Binomial test
    if (obs_count == 0 && exp_count == 0) {
      p_value <- NA_real_
    } else if (exp_prop == 0 || exp_prop == 1) {
      p_value <- NA_real_
    } else {
      binom_result <- tryCatch({
        binom.test(obs_count, n_obs_total, exp_prop, alternative = "two.sided")
      }, error = function(e) {
        list(p.value = NA_real_)
      })
      p_value <- binom_result$p.value
    }

    # Odds ratio
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

  # Chi-square omnibus test for overall distribution
  # Build vectors of observed and expected counts (aligned by topology)
  obs_vec <- sapply(all_topology_names, function(t) {
    ifelse(t %in% names(obs_counts), obs_counts[t], 0)
  })
  exp_vec <- sapply(all_topology_names, function(t) {
    ifelse(t %in% names(expected_counts), expected_counts[t], 0)
  })

  # Chi-square test
  chisq_result <- tryCatch({
    chisq.test(x = obs_vec, p = exp_vec / sum(exp_vec), rescale.p = TRUE)
  }, error = function(e) {
    list(statistic = NA_real_, parameter = NA_real_, p.value = NA_real_,
         method = "Chi-square test (failed)")
  }, warning = function(w) {
    # Run test but capture result even with warnings
    suppressWarnings(
      chisq.test(x = obs_vec, p = exp_vec / sum(exp_vec), rescale.p = TRUE)
    )
  })

  chisq_summary <- tibble(
    chi_squared = as.numeric(chisq_result$statistic),
    df = as.numeric(chisq_result$parameter),
    p_value_chisq = chisq_result$p.value,
    n_topologies = length(all_topology_names),
    n_observations = n_obs_total
  )

  # Return both detailed results and chi-square summary
  return(list(
    topology_tests = results,
    chisq_test = chisq_summary
  ))
}
