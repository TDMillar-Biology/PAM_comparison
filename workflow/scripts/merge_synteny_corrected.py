#!/usr/bin/env python3
import sys
import pandas as pd

# ------------------------------------------------------------
# Helper: collapse duplicate hits INSIDE a genome
# ------------------------------------------------------------
def collapse_by_best_cfd(df, label):
    print(f"[INFO] Collapsing {label}: starting with {len(df):,} rows")

    # Drop rows missing block_id (they cannot be merged but we keep them)
    no_block = df[df["block_id"].isna()].copy()
    with_block = df[df["block_id"].notna()].copy()

    # Group by protospacerID + block_id and pick highest CFD
    collapsed = (
        with_block.sort_values("CFD", ascending=False)
        .groupby(["protospacerID", "block_id"], as_index=False)
        .first()
    )

    print(f"[INFO] Collapsed {label}: {len(collapsed):,} rows remain with blocks")
    print(f"[INFO] {len(no_block):,} rows had no block_id (kept separately)")

    return collapsed, no_block


# ------------------------------------------------------------
# Main merge logic
# ------------------------------------------------------------
def merge_synteny(iso1_csv, bl_csv, out_csv):
    print("[INFO] Loading input CSVs...")
    iso1 = pd.read_csv(iso1_csv)
    bl   = pd.read_csv(bl_csv)

    # --- Collapse duplicates by best CFD ---
    iso1_collapsed, iso1_noblock = collapse_by_best_cfd(iso1, "ISO1")
    bl_collapsed,   bl_noblock   = collapse_by_best_cfd(bl, "BL54591")

    print("[INFO] Merging ISO1 & BL54591 by (protospacerID, block_id)...")

    merged = iso1_collapsed.merge(
        bl_collapsed,
        on=["protospacerID", "block_id"],
        how="outer",
        suffixes=("_iso1", "_bl")
    )

    print(f"[INFO] Merge produced {len(merged):,} rows")

    # ------------------------------------------------------------
    # Add lost flags
    # ------------------------------------------------------------
    merged["lost_iso1"] = merged["CFD_iso1"].isna()
    merged["lost_bl"]   = merged["CFD_bl"].isna()

    # ------------------------------------------------------------
    # Append rows that had no block_id in either genome
    # ------------------------------------------------------------
    final = pd.concat([merged, iso1_noblock, bl_noblock], ignore_index=True)
    print(f"[INFO] Final table after adding noblock rows: {len(final):,}")

    # ------------------------------------------------------------
    # Optional: sort by protospacerID for nicer readability
    # ------------------------------------------------------------
    final = final.sort_values(["protospacerID", "block_id"], na_position="last")

    # ------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------
    final.to_csv(out_csv, index=False)
    print(f"[✔] Wrote final merged table → {out_csv}")

# ------------------------------------------------------------
# CLI entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python merge_synteny_corrected.py iso1.csv bl.csv output.csv")
        sys.exit(1)

    merge_synteny(sys.argv[1], sys.argv[2], sys.argv[3])
