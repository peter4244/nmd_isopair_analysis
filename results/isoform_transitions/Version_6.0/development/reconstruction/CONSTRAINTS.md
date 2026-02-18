# Development Folder Constraints

## Primary Constraint

**ONLY USE SCRIPTS FROM THIS DIRECTORY (`development/reconstruction/`)**

- ❌ Do NOT use scripts from `../../scripts/`
- ❌ Do NOT use scripts from `../../testing/`
- ✅ ONLY use scripts listed below in "Available Resources"

## Available Resources

### Pipeline Scripts (in this directory)

1. **prepare_test_data.R**
   - Validates GTF files
   - Randomly selects dominant isoforms
   - Creates `dominant_isoforms.tsv`
   - Usage: `Rscript prepare_test_data.R <input.gtf> <output_dominant_mapping.tsv>`

2. **build_atomic_union_exons.R**
   - Creates atomic union exons from GTF
   - Creates tabix-indexed union exon file
   - Usage: `Rscript build_atomic_union_exons.R <input.gtf> <output_file>`
   - Outputs: `<output_file>.gz` and `<output_file>.gz.tbi` index

3. **generate_pairs_from_dominant.R**
   - Generates pairs.tsv from dominant mapping
   - Usage: `Rscript generate_pairs_from_dominant.R <gtf> <dominant_mapping.tsv> <output_pairs.tsv>`

4. **detect_and_save_events.R**
   - Detects splicing events between isoform pairs
   - Sources: `../../scripts/event_detection_functions.R` (ALLOWED - required dependency)
   - Usage: `Rscript detect_and_save_events.R <gtf> <pairs.tsv> <output_events.tsv>`
   - Must be run from reconstruction/ directory (relative path dependency)

5. **reconstruct_dominant_isoforms.R**
   - Reconstructs dominant isoforms from comparator + events
   - Sources: `reconstruction_functions.R` (local)
   - Usage: `Rscript reconstruct_dominant_isoforms.R <comparator.gtf> <events.tsv> <union_exons.tsv.gz> <output.gtf> <output_log.tsv>`

6. **verify_reconstruction.R**
   - Verifies reconstructed isoforms match originals
   - Usage: `Rscript verify_reconstruction.R <original.gtf> <reconstructed.gtf> <events.tsv> <output_verification.tsv>`

7. **reconstruction_functions.R**
   - Library of reconstruction functions
   - Sourced by reconstruct_dominant_isoforms.R
   - Not called directly

8. **visualize_verification.R**
   - Visualizes verification failures
   - Usage: `Rscript visualize_verification.R <verification.tsv> <original.gtf> <reconstructed.gtf> <events.tsv> <output_dir>`

9. **visualization_functions.R**
   - Library of visualization functions
   - Sourced by visualize_verification.R
   - Not called directly

10. **visualize_all_three.R**
   - Visualizes all three isoforms (dominant, reconstructed, comparator) for each gene
   - Usage: `Rscript visualize_all_three.R <original.gtf> <reconstructed.gtf> <comparator.gtf> <events.tsv> <output.pdf>`

## Complete Validation Workflow

Starting from a GTF file in a test data folder (e.g., `synthetic_data/`):

```bash
# Set working directory
cd development/reconstruction/

# Step 1: Validate GTF and select dominants
Rscript prepare_test_data.R ...

# Step 2: Build atomic union exons
Rscript build_atomic_union_exons.R ...

# Step 3: Generate pairs
Rscript generate_pairs_from_dominant.R ...

# Step 4: Detect events
Rscript detect_and_save_events.R ...

# Step 5: Extract comparator GTF (bash command)
tail -n +2 <data_dir>/events.tsv | cut -f3 | sort -u > /tmp/comparator_ids.txt
grep -Ff /tmp/comparator_ids.txt <data_dir>/base_events.gtf > <data_dir>/comparator.gtf

# Step 6: Reconstruct
Rscript reconstruct_dominant_isoforms.R <data_dir>/comparator.gtf <data_dir>/events.tsv <data_dir>/union_exons.tsv.gz <data_dir>/reconstructed.gtf <data_dir>/reconstruction_log.tsv

# Step 7: Verify
Rscript verify_reconstruction.R <data_dir>/base_events.gtf <data_dir>/reconstructed.gtf <data_dir>/events.tsv <data_dir>/verification.tsv

# Step 8 (Optional): Visualize all three isoforms
Rscript visualize_all_three.R <data_dir>/base_events.gtf <data_dir>/reconstructed.gtf <data_dir>/comparator.gtf <data_dir>/events.tsv <data_dir>/three_isoform_viz.pdf
```

## Working Directories

- **synthetic_data/** - Clean synthetic test data (44 test cases)
- **real_data/** - Real gene data testing

## Pre-Execution Checklist

Before running any script, verify:

- [ ] Script is in `development/reconstruction/` directory
- [ ] All input files exist
- [ ] Output paths are specified correctly
- [ ] Working directory is set appropriately (usually `reconstruction/` for scripts with relative imports)

## When in Doubt

If you're about to use a script:
1. Check: Is it in the list above?
2. If NO → STOP and find the correct script in this directory
3. If YES → Proceed

**Remember: The constraint exists to ensure reproducibility and test the complete self-contained workflow.**
