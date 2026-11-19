# ISO1 → BL54591: guide provenance and target-site reporting

This document describes the **new outputs**, not a change in the Bowtie or
Doench CFD scoring algorithm. The published `conserved_pams.csv`,
`failed_tolerance_pams.csv`, and `conservation_by_tolerance.tsv` continue to use
the earlier **one selected pair per protospacer** approach and should be
explicitly labeled *legacy guide-level summaries*. They are not the new
site-level conservation estimate.

## Guide discovery

`workflow/scripts/pam_discovery.py` now writes a third synchronized output:

- `<prefix>.fa` — unchanged original 23-mer candidates, named `PAM_1`, ...
- `<prefix>.bed` — unchanged six-column BED, local to the extracted euchromatin FASTA
- `<prefix>.tsv` — guide provenance manifest (one source occurrence per `PAM_n`)

Source manifest columns: `protospacerID`, `source_chr`,
`source_start_local`, `source_end_local`, `source_start`, `source_end`,
`source_strand`, `guide_23mer`, `guide_20mer`, `source_pam`.

**Coordinate warning:** `extract_and_fix_euchromatin` generates a new FASTA from
`config/euchromatin_boundaries.bed` and strips interval offsets from FASTA
headers. Consequently the existing reference CFD `ref_start`/`ref_end` and
PAM discovery BED refer to the **local, extracted FASTA**. In contrast,
`source_start`/`source_end` in the TSV are zero-based, half-open positions on
the original ISO1 chromosome, recovered by adding the start of the extraction
interval. Query CFD coordinates are on query scaffolds. Never compare raw
reference and query positions without synteny projection.

The offset conversion requires one euchromatin extraction interval per
chromosome, as in the current `config/euchromatin_boundaries.bed`; the discovery
script raises an error if the provided BED includes multiple intervals for the
same chromosome rather than silently producing incorrect positions.

The source PAM is derived from the 23-mer (`guide_23mer[21:23]`, usually `GG`).
The CFD CSV `pam` remains the **genomic PAM dinucleotide of each mapped hit**,
not the source PAM. The scoring code has not been changed.

## New synteny-mapper outputs

All are in `results/07_summary/ISO1_release6_BL54591/`.

| File | Unit of observation | Use |
| --- | --- | --- |
| `target_sites.tsv.gz` | One retained Bowtie/CFD alignment in one assembly | Canonical site/CFD dataset; includes source provenance and site class |
| `site_relationships.tsv.gz` | One reference × query alignment pairing sharing `protospacerID` | Evidence table, including repeated/repetitive pairs |
| `guide_targeting_summary.tsv.gz` | One guide detected in at least one raw CFD table | `n_ref`, `n_qry`, `delta_hits`, matched, unmatched, unresolved alignment counts, source CFD, strongest unmatched-query CFD |
| `relationship_multiplicity_summary.tsv` | One `(n_ref, n_qry)` bin among *synteny-assigned* hits in both genomes | Ref × query multiplicity matrix; pairing counts |

The legacy summary TSV and figures continue to use one pair per guide.
**Do not combine their denominators with the new site counts.**

### Relationship classes (one row per potential pair)

- `syntenic_within_tolerance`: same chromosome and syntenic block, and
  projected query midpoint within `--tol` bp (default 1,000) of the projected
  reference midpoint.
- `same_block_outside_tolerance`: same chromosome/block, projection distance
  greater than the threshold.
- `same_chr_different_block`: same chromosome, different blocks.
- `different_chr`: different chromosome arms/scaffolds.
- `unresolved`: no comparable coordinate. In this pipeline, alignments without
  unique syntenic blocks are excluded **before** the pairing step, so this
  relationship category normally has zero entries.

The original raw alignment tables can have synteny-unassigned hits; those are
included in `target_sites.tsv.gz` with the site-level status `unresolved`.
Projecting numeric coordinates between different chromosomes is meaningless;
`comparable_proj_dist` is null for those pairs. Same-chromosome/different-block
numeric distances are preserved as *diagnostics*, but not treated as evidence
of an orthologous locus without further validation.

### Site classes (one row per genomic alignment)

**Reference / ISO1**

- `design_source`: the exact local source chromosome, interval and oriented
  23-mer identified in the discovery manifest.
- `conserved_alternative`: a non-source reference hit with ≥1 qualifying
  syntenic query relationship.
- `reference_unmatched`: a non-source, synteny-assigned reference hit without
  a qualifying query relationship.
- `unresolved`: a non-source reference hit without a unique syntenic block.

**Query / BL54591**

- `corresponding_source`: query hit with ≥1 qualifying relationship to the
  *reference design-source* alignment.
- `conserved_alternative`: other query hit with ≥1 qualifying relationship to
  a reference alternative alignment.
- `query_unmatched`: synteny-assigned query hit without a qualifying reference
  partner, including within-guide `1×1` nonsyntenic cases.
- `unresolved`: query hit without a unique syntenic block.

A source reference hit lacking a qualifying query match stays labeled
`design_source`; use its `has_syntenic_partner` boolean or the per-guide counts
to distinguish matched from unmatched source sites.

A repeated hit may match multiple sites in the other genome. We preserve every
qualifying edge, but count each alignment **once** in site and guide summaries.
Thus 10×10 can produce 100 qualifying relationships and only 10 reference plus
10 query target sites.

These are **detectable alignments under Bowtie's `-v 3 -m 10` restrictions**,
not exhaustive guides discovered independently in BL54591. A query-unmatched
hit is a *candidate query-specific detectable target* under the matching
criteria, **not** demonstrated evolutionary gain or absolute absence from ISO1.
Particularly cautious interpretation is required when a guide has unresolved
reference hits.

`design_source_CFD_ref` and `corresponding_source_CFD_query` in the site table
are repeated per-guide context values. If multiple possible query source
correspondences exist, the latter uses the maximum CFD for a convenient summary;
inspect `site_relationships.tsv.gz` to retain all correspondence-level CFD scores.

### Examples

```bash
# Inspect new site-level classes without decompressing the whole file to disk.
python3 - <<'PY'
import pandas as pd
p = 'results/07_summary/ISO1_release6_BL54591/target_sites.tsv.gz'
print(pd.read_csv(p, sep='\t', nrows=20)[['protospacerID','name','assembly','site_class','CFD']])
PY

# Reconstruct the reference × query guide multiplicity matrix:
python3 - <<'PY'
import pandas as pd
p = 'results/07_summary/ISO1_release6_BL54591/relationship_multiplicity_summary.tsv'
x = pd.read_csv(p, sep='\t')
print(x.pivot(index='n_ref', columns='n_qry', values='protospacers').fillna(0).astype(int))
PY
```

### Testing and compute considerations

```bash
python3 -m unittest discover -s test -p 'test_provenance_relationships.py' -v
snakemake -s workflow/Snakefile -n --use-singularity
```

The tests exercise source orientation, offset arithmetic, 1×1 nonsyntenic,
unequal multiplicity, unresolved mappings, a 2×2 fully qualifying many-to-many
case, and same-block/outside-tolerance classification. SVMU2 and full-scale
Snakemake were **not available** in the editing environment; run the complete
pipeline's dry run and verify resource allocations on your system before
production. The mapper already holds a large Cartesian DataFrame; the
additional compressed tables incur substantial disk I/O and require memory
benchmarking. The current `mem_mb: 16000` and `runtime: 60` settings are
unchanged from the supplied archive and may not suffice for new reports.

**Reproducibility:** Back up existing `results/07_summary/...` outputs if you
want a direct comparison. The revised mapper writes files with the same legacy
names and will replace them when run.
