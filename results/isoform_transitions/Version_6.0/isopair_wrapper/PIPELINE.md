# Isopair Wrapper Pipeline — Execution Order

## Script Dependency DAG

```
01_prepare_data_mashr.R          ← isocall count matrix, mashr DE CSVs
  │
  ├→ data_mashr/expression_data.rds
  ├→ data_mashr/sample_metadata.rds
  ├→ data_mashr/gene_map.rds
  ├→ data_mashr/nmd_classification.rds
  ├→ data_mashr/dmso_samples.rds
  └→ data_mashr/smg1i_samples.rds

02_build_profiles_mashr.R        ← 01 outputs + Isopair infrastructure
  │
  ├→ data_mashr/structures.rds (copy)
  ├→ data_mashr/cds.rds (copy)
  ├→ data_mashr/ptc.rds (copy)
  ├→ data_mashr/profiles_c2_*.rds
  └→ data_mashr/profiles_c4_*.rds

04_run_analyses_mashr.R          ← 02 outputs (profiles)
  │
  └→ data_mashr/analysis_cache/
       ├→ ptc_c2_*.rds, ptc_c4_*.rds
       ├→ fw_c2_*.rds, fw_c4_*.rds
       ├→ div_c2_*.rds, div_c4_*.rds
       └→ fc_c2_*.rds, fc_c4_*.rds

05k_utr5_all_isoforms.R         ← structures, cds, nmd_classification, FASTA
  │
  └→ data_mashr/analysis_cache/utr5_features_all.rds

05l_unified_model.R              ← nmd_classification, cds, ptc, utr5_features_all
  │
  └→ data_mashr/analysis_cache/unified_model.rds

05r_ref_atg_analysis.R           ← cds, structures, profiles, ptc, FASTA
  │
  └→ data_mashr/analysis_cache/ref_atg_analysis.rds

05s_orfik_scan.R                 ← structures, cds, nmd_classification, FASTA
  │
  └→ data_mashr/analysis_cache/orfik_scan.rds

05t_ref_cds_features.R           ← expression_data, sample_metadata, nmd_classification,
  │                                 cds, structures, gene_map, FASTA
  └→ data_mashr/analysis_cache/ref_cds_features_all.rds

05u_paralog_annotation.R         ← structures, nmd_classification (+ biomaRt query)
  │
  └→ data_mashr/analysis_cache/paralog_genes.rds

05v_model_comparison.R           ← unified_model, ref_cds_features_all, paralog_genes
  │
  └→ data_mashr/analysis_cache/model_comparison.rds

05_final_report_mashr.Rmd        ← ALL of the above
  │
  └→ 05_final_report_mashr.html
```

## Execution Order

Scripts must be run in this order (dependencies are strict):

```
1.  01_prepare_data_mashr.R
2.  02_build_profiles_mashr.R
3.  04_run_analyses_mashr.R
4.  05k_utr5_all_isoforms.R     (independent of 04)
5.  05l_unified_model.R          (depends on 05k)
6.  05r_ref_atg_analysis.R       (independent of 05k/05l)
7.  05s_orfik_scan.R             (independent of 05k/05l/05r)
8.  05t_ref_cds_features.R       (independent of 05l/05r/05s)
9.  05u_paralog_annotation.R     (independent of 05l-05t)
10. 05v_model_comparison.R       (depends on 05l, 05t, 05u)
11. 05_final_report_mashr.Rmd    (depends on ALL above)
```

Steps 4-9 can run in parallel after step 3.
Step 10 requires steps 5, 8, and 9.
Step 11 requires all previous steps.

## External Dependencies

- **SQANTI corrected FASTA**: Used by 05k, 05r, 05s, 05t
- **mashr DE CSVs**: Used by 01
- **limma DE CSVs**: Used by 05_final_report (dose-response analysis)
- **biomaRt (Ensembl)**: Used by 05u (requires internet)
- **Isopair package**: Used by 02, 04, 05k, 05_final_report
