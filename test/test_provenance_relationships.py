"""Small regression checks for source provenance and site-level reporting.

Run from the repository root:
    python3 -m unittest discover -s test -p 'test_provenance_relationships.py'

These tests exercise reporting logic without requiring the large CFD data or
SVMU2; its interval-tree assignment is intentionally tested separately.
"""
import os
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from pam_discovery import discover_pams, source_offsets, write_outputs
from target_site_reporting import mark_design_source, write_analysis_tables


class ProvenanceTests(unittest.TestCase):
    def test_forward_reverse_and_offset(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            # Both orientations, with a spacer that does not itself contain an NGG.
            fwd = "A" * 20 + "AGG"
            rev = "CC" + "A" * 21
            fasta = tmp / "masked.fa"
            fasta.write_text(f">X\n{fwd}{'N'*6}{rev}\n")
            bed = tmp / "intervals.bed"
            bed.write_text("X\t277911\t278011\n")
            count = write_outputs(discover_pams(fasta), str(tmp / "out"), source_offsets(bed))
            self.assertGreaterEqual(count, 2)
            tsv = pd.read_csv(tmp / "out.tsv", sep="\t")
            self.assertEqual(tsv.loc[0, "source_start"], 277911)
            self.assertEqual(tsv.loc[0, "source_end_local"], 23)
            self.assertEqual(tsv.loc[0, "source_pam"], "GG")
            self.assertIn("-", set(tsv["source_strand"]))
            self.assertEqual((tmp / "out.bed").read_text().splitlines()[0].split("\t")[:4],
                             ["X", "0", "23", "PAM_1"])

    def test_ref_query_sites_and_repetitive_pairs(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            manifest = pd.DataFrame([
                ("PAM_1", "X", 10, 33, 277921, 277944, "+", "A"*21+"GG", "A"*20, "GG"),
                ("PAM_2", "X", 50, 73, 277961, 277984, "+", "C"*21+"GG", "C"*20, "GG"),
                ("PAM_3", "2L", 20, 43, 82475, 82498, "+", "G"*21+"GG", "G"*20, "GG"),
                ("PAM_4", "X", 300, 323, 278211, 278234, "+", "T"*21+"GG", "T"*20, "GG"),
                ("PAM_5", "X", 350, 373, 278261, 278284, "+", "T"*20+"AGG", "T"*20, "GG"),
                ("PAM_6", "X", 400, 423, 278311, 278334, "+", "A"*20+"TGG", "A"*20, "GG"),
                ("PAM_7", "X", 600, 623, 278511, 278534, "+", "A"*20+"CGG", "A"*20, "GG"),
                ("PAM_8", "X", 800, 823, 278711, 278734, "+", "C"*20+"AGG", "C"*20, "GG"),
            ], columns=["protospacerID", "source_chr", "source_start_local", "source_end_local",
                         "source_start", "source_end", "source_strand", "guide_23mer",
                         "guide_20mer", "source_pam"])

            def hit(guide, name, chrom, start, block, cfd, yhat=None, source=False):
                guide_seq = manifest.set_index('protospacerID').loc[guide, 'guide_23mer']
                return dict(protospacerID=guide, name=name, ref_chr=chrom,
                            ref_start=start, ref_end=start+23, midpoint=start+11.5,
                            syntenic_block=block, yhat=yhat, CFD=cfd, pam="GG",
                            genome_target=guide_seq if source else 'A'*20 + 'TGG',
                            wt=guide_seq[:20], sg=guide_seq[:20], mismatches=0,
                            protospacer=guide_seq)

            ref = pd.DataFrame([
                hit('PAM_1','PAM_1_1','X',10,0,1.0,110.5,True),
                hit('PAM_1','PAM_1_2','3R',200,1,.06,500.5),
                hit('PAM_2','PAM_2_1','X',50,0,1.0,150.5,True),
                hit('PAM_3','PAM_3_1','2L',20,5,1.0,220.5,True),
                hit('PAM_4','PAM_4_1','X',300,3,1.0,400.5,True),
                hit('PAM_6','PAM_6_1','X',400,float('nan'),1.0,float('nan'),True),
                hit('PAM_7','PAM_7_1','X',600,9,1.0,700.5,True),
                hit('PAM_7','PAM_7_2','X',650,9,.2,750.5),
                hit('PAM_8','PAM_8_1','X',800,11,1.0,900.5,True),
            ])
            qry = pd.DataFrame([
                hit('PAM_1','PAM_1_1','X',99,0,.9),
                hit('PAM_1','PAM_1_2','3R',489,1,.06),
                hit('PAM_1','PAM_1_3','2L',800,3,.75), # additional query off-target
                hit('PAM_2','PAM_2_1','3R',140,1,.45), # 1x1 but nonsyntenic
                hit('PAM_3','PAM_3_1','2L',210,6,.55), # same chromosome, different block
                hit('PAM_5','PAM_5_1','X',700,7,.7),  # no reference CFD hit
                hit('PAM_6','PAM_6_1','X',500,4,.7),  # source unresolved
                hit('PAM_7','PAM_7_1','X',690,9,.8),
                hit('PAM_7','PAM_7_2','X',740,9,.2), # all four pairs qualify
                hit('PAM_8','PAM_8_1','X',5000,11,.7), # >1 kb within block
            ]).drop(columns=['yhat'])
            source = mark_design_source(ref, manifest)
            qry["is_design_source"] = False
            self.assertEqual(ref["is_design_source"].sum(), 7)
            mapped_ref = ref.dropna(subset=['syntenic_block','yhat'])
            mapped_qry = qry.dropna(subset=['syntenic_block'])
            merged = mapped_ref.merge(mapped_qry, on='protospacerID',suffixes=('_ref','_qry'))
            merged['same_chr'] = merged['ref_chr_ref'].eq(merged['ref_chr_qry'])
            merged['same_block'] = (merged['same_chr']
                                    & merged['syntenic_block_ref'].eq(merged['syntenic_block_qry']))
            merged['comparable_proj_dist'] = (merged['midpoint_qry'] - merged['yhat']).abs().where(merged['same_chr'])
            import numpy as np
            tol = 1000
            merged['relationship'] = np.select(
                [merged['same_block'] & merged['comparable_proj_dist'].le(tol),
                 merged['same_block'] & merged['comparable_proj_dist'].gt(tol),
                 merged['same_chr'] & ~merged['same_block'], ~merged['same_chr']],
                ['syntenic_within_tolerance','same_block_outside_tolerance',
                 'same_chr_different_block','different_chr'], default='unresolved')
            write_analysis_tables(ref, qry, merged, source, out, 'ISO1_release6', 'BL54591')
            sites = pd.read_csv(out / 'target_sites.tsv.gz',sep='\t')
            ps = pd.read_csv(out / 'guide_targeting_summary.tsv.gz',sep='\t').set_index('protospacerID')
            edges = pd.read_csv(out / 'site_relationships.tsv.gz',sep='\t')
            self.assertEqual(len(sites), len(ref)+len(qry))
            self.assertEqual(len(edges), 13)  # 2x3 + 1x1 + 1x1 + 2x2 + 1x1
            self.assertEqual(ps.loc['PAM_1','n_qry_unmatched'],1)
            self.assertEqual(ps.loc['PAM_1','n_qry_matched'],2)
            self.assertEqual(ps.loc['PAM_1','n_ref_matched'],2)
            self.assertEqual(ps.loc['PAM_1','delta_hits'],1)
            self.assertAlmostEqual(ps.loc['PAM_1','max_unmatched_qry_CFD'],.75)
            self.assertEqual(ps.loc['PAM_2','n_qry_unmatched'],1)
            self.assertEqual(ps.loc['PAM_3','n_qry_unmatched'],1)
            self.assertEqual(ps.loc['PAM_4','n_qry'],0)
            self.assertEqual(ps.loc['PAM_5','n_ref'],0)
            self.assertEqual(ps.loc['PAM_6','n_ref_unresolved'],1)
            self.assertEqual(ps.loc['PAM_7','n_ref_matched'],2)
            self.assertEqual(ps.loc['PAM_7','n_qry_matched'],2)
            self.assertEqual((edges.query("protospacerID=='PAM_7'")['relationship']
                              == 'syntenic_within_tolerance').sum(), 4)
            self.assertEqual(edges.query("protospacerID=='PAM_8'")['relationship'].iloc[0],
                             'same_block_outside_tolerance')
            self.assertEqual(ps.loc['PAM_8','n_qry_unmatched'],1)
            self.assertEqual(sites.query("protospacerID=='PAM_1' and assembly=='ISO1_release6' and ref_chr=='X'")['site_class'].iloc[0], 'design_source')
            self.assertEqual(sites.query("protospacerID=='PAM_1' and assembly=='BL54591' and ref_chr=='X'")['site_class'].iloc[0], 'corresponding_source')
            self.assertEqual(sites.query("protospacerID=='PAM_1' and assembly=='BL54591' and ref_chr=='2L'")['site_class'].iloc[0], 'query_unmatched')
            self.assertIn('different_chr',set(edges.relationship))
            self.assertIn('same_chr_different_block',set(edges.relationship))
            self.assertAlmostEqual(ps.loc['PAM_1', 'corresponding_source_CFD_query'],.9)

            # Exercise the real mapper's aggregation and legacy output path as
            # well. Stub only unavailable SVMU2 imports (not the mapper code).
            import types
            from unittest.mock import patch
            fake_parse = types.ModuleType("svmu2.orchestration.parse")
            fake_parse.parse = lambda path: (None, {})
            fake_syn = types.ModuleType("svmu2.orchestration.synteny")
            fake_syn.resolve_synteny = lambda **kwargs: None
            fake_root = types.ModuleType("svmu2")
            fake_orch = types.ModuleType("svmu2.orchestration")
            stubs = {"svmu2": fake_root, "svmu2.orchestration": fake_orch,
                     "svmu2.orchestration.parse": fake_parse,
                     "svmu2.orchestration.synteny": fake_syn}
            with patch.dict(sys.modules, stubs):
                import importlib
                mapper = importlib.import_module("synteny_mapper")
            counts = mapper.candidate_counts(ref, qry)
            mapper.merge_and_filter(ref, qry, out, 1000, [1000], counts,
                                    source, 'ISO1_release6', 'BL54591')
            grouped = pd.read_csv(out / 'relationship_multiplicity_summary.tsv', sep='\t')
            self.assertEqual(grouped.total_pairs.sum(), 13)
            self.assertEqual(grouped.syntenic_pairs.sum(), 6)  # PAM1=2, PAM7=4
            self.assertEqual(grouped.same_block_outside_tolerance_pairs.sum(), 1)
            self.assertTrue((out / 'conservation_by_tolerance.tsv').exists())


if __name__ == '__main__':
    unittest.main()
