**SF38 | Attention distribution across the model's five candidate open reading frames, NMD susceptible vs non-NMD.**

The deep-learning model evaluates each transcript at up to five candidate open reading frames and learns an attention weight for each. The candidate ORFs are ranked by a priority rule that prefers the annotated coding sequence first, then the standard caller's choice, then the open reading frame with the strongest Kozak start context. Most attention concentrates at the priority ORF in both classes.

**(A)** Attention weight per candidate-ORF rank, by class (NMD-susceptible in coral, non-NMD controls in blue). Most attention concentrates at the priority ORF (rank 0) in both classes (median ~0.74–0.76). The median rank-0 attention is slightly lower in NMD, with a corresponding shift of attention onto non-priority ORFs at ranks 1–4 and a long tail of NMD isoforms placing more than 25% of attention on a non-priority candidate.

**(B)** Shannon entropy of the per-isoform attention vector across the five candidate ORFs (higher entropy = more broadly distributed attention). The NMD density is shifted toward higher entropy than the Control density, quantifying the broader spread shown in Panel A.

**Abbreviations.** NMD, nonsense-mediated decay; ORF, open reading frame.
