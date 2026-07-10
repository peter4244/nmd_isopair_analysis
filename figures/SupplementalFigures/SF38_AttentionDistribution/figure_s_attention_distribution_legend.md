**SF38 | Attention Distribution Across the Model's Five Candidate Open Reading Frames, NMD Susceptible vs Non-NMD.** The deep-learning model evaluates each transcript at up to five candidate open reading frames and learns an attention weight for each. The candidates are ranked by a priority rule that prefers the annotated reference CDS first (projected from the gene's dominant non-NMD isoform when available), then the TransDecoder2-called CDS, then the open reading frame with the strongest Kozak start context. Most attention concentrates at the main (rank-0) open reading frame in both classes.

(A) Attention weight per candidate open reading frame rank, by class (NMD susceptible in coral, non-NMD in blue). Most attention concentrates at the rank-0 open reading frame in both classes (median ~0.74–0.76). The median rank-0 attention is slightly lower in NMD, with a corresponding shift of attention onto ranks 1–4 and a long tail of NMD isoforms placing more than 25% of attention on a non-rank-0 candidate.

(B) Shannon entropy of the per-isoform attention vector across the five candidate open reading frames (higher entropy, more broadly distributed attention). The NMD density is shifted toward higher entropy than the Control density, quantifying the broader spread shown in Panel A.

Abbreviations. NMD, nonsense-mediated decay.
