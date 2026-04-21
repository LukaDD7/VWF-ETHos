#!/usr/bin/env python3
"""
VWF Alpha-Matrix Data Integration
==================================
Phase 1: Merge AlphaGenome + AF3 + Type-2 Labels

Input:
- results/03_inference_results.csv (AlphaGenome 1198 variants)
- VWF_Type2_variants.csv (100 known Type-2 labels)
- Proteo-Structure-Pipeline/output/af3_batches_type2/AF3_Results/extracted/ (62 AF3 structures)

Output:
- VWF_Alpha_Matrix.parquet (unified feature matrix)

Key Design (per Alphagenome-DFR-Phase3.md):
- ag_rna_delta: 区分 Type 1 vs Type 2A (表达量断崖→Type1)
- ag_splice_delta: 强制介入 Type 2A/2M (剪接异常)
- D4 domain: 用ag_rna_delta解决Type1 vs 2A竞争

Author: Claude Code
Date: 2026-04-19
"""

import pandas as pd
import numpy as np
import re
import json
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import warnings

# ============================================================================
# CONFIG
# ============================================================================

BASE_DIR = Path("/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan")
AG_CSV = BASE_DIR / "results" / "03_inference_results.csv"
TYPE2_CSV = BASE_DIR / "VWF_Type2_variants.csv"
AF3_EXTRACTED = BASE_DIR / "Proteo-Structure-Pipeline" / "output" / "af3_batches_type2" / "AF3_Results" / "extracted"
OUTPUT_PARQUET = BASE_DIR / "VWF_Alpha_Matrix.parquet"

# VWF Domain Architecture (from literature)
VWF_DOMAINS = {
    'D1': (1, 233),
    'D2': (234, 480),
    'D3': (481, 728),
    "D'": (729, 763),
    'D3_extended': (764, 1233),
    'A1': (1271, 1492),
    'A2': (1493, 1684),
    'A3': (1685, 1873),
    'D4': (1874, 2255),
    'C1': (2256, 2299),
    'C2': (2300, 2363),
    'C3': (2364, 2411),
    'C4': (2412, 2450),
    'C5': (2451, 2459),
    'C6': (2460, 2527),
    'CK': (2528, 2813),
}

# Type-2 subdomain associations (from Blood 2024/2026 literature)
TYPE2_DOMAIN_RULES = {
    '2A': ['A2', 'D4', 'C1-C6', 'CK'],
    '2B': ['A1'],
    '2M': ['A1', 'A3'],
    '2N': ["D'", 'D3'],
}

# ============================================================================
# PARSE UTILITIES
# ============================================================================

def parse_cdna_to_protein_position(cdna: str) -> Optional[int]:
    """
    Parse cDNA notation like 'c.8324C>G' to protein position.
    Extracts the numeric position from cDNA change.
    Example: 'c.8324C>G' -> 2775 (Serine position from S2775C)
    """
    if pd.isna(cdna):
        return None

    # Pattern: c.(\d+)[A-Z]>.*
    match = re.search(r'c\.(\d+)', str(cdna))
    if match:
        dna_pos = int(match.group(1))
        # VWF has 8667 bp coding sequence, ~2813 aa
        # Approximate: protein_pos ≈ dna_pos / 3
        return (dna_pos - 1) // 3 + 1

    return None


def parse_aa_change(aa_change: str) -> Optional[Tuple[str, int, str]]:
    """
    Parse amino acid change like 'S2775C' -> (ref='S', pos=2775, alt='C')
    """
    if pd.isna(aa_change):
        return None

    match = re.match(r'^([A-Zaiz*])(\d+)([A-Zaiz*])$', str(aa_change))
    if match:
        return (match.group(1), int(match.group(2)), match.group(3))

    return None


def get_domain_from_position(pos: int) -> Optional[str]:
    """Get VWF domain from protein position."""
    for domain, (start, end) in VWF_DOMAINS.items():
        if start <= pos <= end:
            return domain
    return None


def extract_af3_features_from_folder(folder_path: Path) -> Dict:
    """
    Extract AF3 features from a fold_vwf_XXX folder.
    Reads full_data JSON for atom_plddts and pae data.
    """
    features = {
        'af3_has_structure': False,
        'af3_plddt_mean': np.nan,
        'af3_plddt_min': np.nan,
        'af3_local_rmsd': np.nan,
        'af3_pae_interface': np.nan,  # PAE at interface (for 2N/2B binding)
    }

    if not folder_path.exists():
        return features

    # Check for model files
    cif_files = list(folder_path.glob("fold_vwf_*_model_*.cif"))
    if not cif_files:
        return features

    features['af3_has_structure'] = True

    # Read pLDDT from full_data JSON (atom_plddts array)
    full_data_file = folder_path / f"{folder_path.name}_full_data_0.json"
    if full_data_file.exists():
        try:
            with open(full_data_file) as f:
                data = json.load(f)

            # atom_plddts is a list of pLDDT scores per atom
            if 'atom_plddts' in data and len(data['atom_plddts']) > 0:
                plddts = data['atom_plddts']
                features['af3_plddt_mean'] = np.mean(plddts)
                features['af3_plddt_min'] = np.min(plddts)

            # PAE matrix for interface analysis (first token pair)
            if 'pae' in data and len(data['pae']) > 0:
                # PAE is a 2D matrix, get the first off-diagonal element as interface indicator
                pae = data['pae']
                if isinstance(pae, list) and len(pae) > 0 and len(pae[0]) > 1:
                    features['af3_pae_interface'] = pae[0][1]  # Proxy for interface confidence

        except Exception as e:
            print(f"    Warning: Could not parse {full_data_file}: {e}")

    return features


def parse_af3_folder_name(folder_name: str) -> Optional[Tuple[str, int, str]]:
    """
    Parse AF3 folder name like 'fold_vwf_r816q' -> (ref='R', pos=816, alt='Q')
    Case-insensitive - always returns uppercase for amino acids.
    """
    match = re.match(r'fold_vwf_([a-zA-Z*])(\d+)([a-zA-Z*])', folder_name)
    if match:
        return (match.group(1).upper(), int(match.group(2)), match.group(3).upper())
    return None


# ============================================================================
# MAIN INTEGRATION
# ============================================================================

def merge_alpha_features():
    """
    Main integration function.
    Creates VWF_Alpha_Matrix by merging:
    1. AlphaGenome inference results (AG)
    2. Type-2 variant labels (ground truth)
    3. AF3 structural features
    """

    print("=" * 60)
    print("VWF Alpha-Matrix Data Integration (Phase 1)")
    print("=" * 60)

    # ========================================================================
    # STEP 1: Load AlphaGenome Results
    # ========================================================================
    print("\n[Step 1] Loading AlphaGenome results...")
    ag_df = pd.read_csv(AG_CSV)
    print(f"  - Loaded {len(ag_df)} variants from {AG_CSV.name}")
    print(f"  - Columns: {ag_df.columns.tolist()}")

    # Clean up AG data
    ag_df = ag_df.rename(columns={
        'rna_seq_delta_max': 'ag_rna_delta',
        'splice_site_delta_max': 'ag_splice_delta',
        'delta_score': 'ag_delta_score',
    })

    # ========================================================================
    # STEP 2: Load Type-2 Labels
    # ========================================================================
    print("\n[Step 2] Loading Type-2 variant labels...")
    # Type-2 CSV structure:
    # Row 0: Multi-line header labels (Transcript\nConsequence, etc.)
    # Row 1: First data row (c.8324C>G, etc.)
    # Data cols: 0=cDNA, 2=ProteinChange, 14=Type2Label, 19=AA_Change
    # Use header=None, then skip row 0 (header labels)
    type2_df_raw = pd.read_csv(TYPE2_CSV, header=None)
    # Skip row 0 (header labels), use row 1+ as data
    type2_data = type2_df_raw.iloc[1:].copy()
    type2_data = type2_data.reset_index(drop=True)

    type2_cols = type2_data.iloc[:, [0, 2, 14, 19]].copy()
    type2_cols.columns = ['cdna_change', 'protein_change', 'type2_label', 'aa_change']

    # Parse protein position from AA_Change
    type2_cols['protein_pos'] = type2_cols['aa_change'].apply(
        lambda x: parse_aa_change(x)[1] if parse_aa_change(x) else None
    )
    type2_cols['ref_aa'] = type2_cols['aa_change'].apply(
        lambda x: parse_aa_change(x)[0] if parse_aa_change(x) else None
    )
    type2_cols['alt_aa'] = type2_cols['aa_change'].apply(
        lambda x: parse_aa_change(x)[2] if parse_aa_change(x) else None
    )

    # Normalize type2_label (e.g., "type2A" -> "2A", "type2M" -> "2M")
    type2_cols['type2_subtype'] = type2_cols['type2_label'].apply(
        lambda x: str(x).replace('type', '').strip() if pd.notna(x) else None
    )

    print(f"  - Type-2 subtypes: {type2_cols['type2_subtype'].value_counts().to_dict()}")
    print(f"  - Protein positions: {type2_cols['protein_pos'].min()} - {type2_cols['protein_pos'].max()}")

    # ========================================================================
    # STEP 3: Extract AF3 Features
    # ========================================================================
    print("\n[Step 3] Extracting AF3 structural features...")

    af3_folders = [f for f in AF3_EXTRACTED.iterdir() if f.is_dir() and f.name.startswith('fold_vwf_')]
    print(f"  - Found {len(af3_folders)} AF3 structure folders")

    af3_records = []
    for folder in af3_folders:
        parsed = parse_af3_folder_name(folder.name)
        if parsed:
            ref_aa, pos, alt_aa = parsed
            features = extract_af3_features_from_folder(folder)
            af3_records.append({
                'protein_pos': pos,
                'ref_aa': ref_aa,
                'alt_aa': alt_aa,
                **features
            })

    af3_df = pd.DataFrame(af3_records)
    print(f"  - Extracted features from {len(af3_df)} AF3 structures")
    print(f"  - Structures with pLDDT: {af3_df['af3_has_structure'].sum()}")

    # ========================================================================
    # STEP 4: Merge AlphaGenome with Type-2 Labels (CORRECT METHOD)
    # ========================================================================
    print("\n[Step 4] Merging AlphaGenome with Type-2 labels...")
    print("  - Using cDNA position extraction from AG consequences (NOT naive offset)")

    # Build a lookup: cDNA position -> AG features
    # Extract cDNA position directly from AG consequence strings
    def extract_cdna_variants(ag_df):
        """Extract cDNA position and ref/alt from AG consequence strings."""
        ag_lookup = {}  # cdna_pos -> {ref, alt, ag_rna_delta, ag_splice_delta, ag_delta_score, genomic_pos}

        for idx, row in ag_df.iterrows():
            cons_str = str(row['consequences'])
            # Pattern: c.\d+[A-Z]>[A-Z] (e.g., c.100C>T)
            matches = re.findall(r'c\.(\d+)([A-Z])>([A-Z])', cons_str)
            for match in matches:
                cdna_pos = int(match[0])
                ref = match[1]
                alt = match[2]
                if cdna_pos not in ag_lookup:
                    ag_lookup[cdna_pos] = {
                        'ref': ref,
                        'alt': alt,
                        'ag_rna_delta': row['ag_rna_delta'],
                        'ag_splice_delta': row['ag_splice_delta'],
                        'ag_delta_score': row['ag_delta_score'],
                        'genomic_pos': row['position']
                    }
        return ag_lookup

    ag_lookup = extract_cdna_variants(ag_df)
    print(f"  - Built AG lookup with {len(ag_lookup)} cDNA positions")

    # Parse cDNA position from Type-2 cdna_change
    def parse_cdna_pos(cdna_str):
        match = re.search(r'c\.(\d+)', str(cdna_str))
        return int(match.group(1)) if match else None

    # Match Type-2 variants to AG via cDNA position
    merged = type2_cols.copy()
    merged['cdna_pos'] = merged['cdna_change'].apply(parse_cdna_pos)

    for col in ['ag_rna_delta', 'ag_splice_delta', 'ag_delta_score', 'ag_genomic_pos']:
        merged[col] = np.nan

    ag_matches = 0
    for idx, row in merged.iterrows():
        cdna_pos = row['cdna_pos']
        if cdna_pos in ag_lookup:
            ag_data = ag_lookup[cdna_pos]
            merged.at[idx, 'ag_rna_delta'] = ag_data['ag_rna_delta']
            merged.at[idx, 'ag_splice_delta'] = ag_data['ag_splice_delta']
            merged.at[idx, 'ag_delta_score'] = ag_data['ag_delta_score']
            merged.at[idx, 'ag_genomic_pos'] = ag_data['genomic_pos']
            ag_matches += 1

    print(f"  - Matched {ag_matches}/{len(merged)} Type-2 variants with AG features via cDNA")
    print(f"  - Match rate: {ag_matches/len(merged)*100:.1f}%")

    # Add domain information
    merged['domain'] = merged['protein_pos'].apply(get_domain_from_position)

    # ========================================================================
    # STEP 5: Merge AF3 Features
    # ========================================================================
    print("\n[Step 5] Merging AF3 structural features...")

    af3_merged = merged.merge(
        af3_df[['protein_pos', 'ref_aa', 'alt_aa', 'af3_has_structure', 'af3_plddt_mean', 'af3_plddt_min']],
        on=['protein_pos', 'ref_aa', 'alt_aa'],
        how='left'
    )

    af3_match_count = af3_merged['af3_has_structure'].sum()
    print(f"  - Matched {af3_match_count}/{len(af3_merged)} Type-2 variants with AF3 structures")

    # ========================================================================
    # STEP 6: Add Derived Features (per Phase3 architecture)
    # ========================================================================
    print("\n[Step 6] Adding derived features per Phase3 architecture...")

    # Type-1 signal: ag_rna_delta极低 (表达量断崖)
    af3_merged['is_rna_drop'] = af3_merged['ag_rna_delta'] < 0.5

    # Splice override: ag_splice_delta极高 → 强制2A/2M
    af3_merged['is_splice_override'] = af3_merged['ag_splice_delta'] > 0.5

    # D4 domain Type1/2A competition signal
    af3_merged['is_d4_domain'] = af3_merged['domain'] == 'D4'

    # Clean up and select final columns
    final_df = af3_merged[[
        'cdna_change', 'aa_change', 'protein_pos', 'ref_aa', 'alt_aa',
        'domain', 'type2_subtype', 'type2_label',
        'ag_rna_delta', 'ag_splice_delta', 'ag_delta_score',
        'is_rna_drop', 'is_splice_override', 'is_d4_domain',
        'af3_has_structure', 'af3_plddt_mean', 'af3_plddt_min'
    ]].copy()

    # ========================================================================
    # STEP 7: Summary Statistics
    # ========================================================================
    print("\n" + "=" * 60)
    print("INTEGRATION SUMMARY")
    print("=" * 60)
    print(f"Total variants: {len(final_df)}")
    print(f"  - With AG features: {final_df['ag_rna_delta'].notna().sum()}")
    print(f"  - With AF3 structures: {final_df['af3_has_structure'].sum()}")
    print(f"\nType-2 subtype distribution:")
    for subtype, count in final_df['type2_subtype'].value_counts().items():
        print(f"  - {subtype}: {count}")
    print(f"\nDomain distribution:")
    for domain, count in final_df['domain'].value_counts().items():
        print(f"  - {domain}: {count}")
    print(f"\nD4 domain variants (Type1 vs 2A competition): {final_df['is_d4_domain'].sum()}")

    # ========================================================================
    # STEP 8: Save Output
    # ========================================================================
    print("\n[Step 8] Saving VWF_Alpha_Matrix.parquet...")
    final_df.to_parquet(OUTPUT_PARQUET, index=False)
    print(f"  - Saved to: {OUTPUT_PARQUET}")
    print(f"  - Size: {OUTPUT_PARQUET.stat().st_size / 1024:.1f} KB")

    return final_df


if __name__ == "__main__":
    df = merge_alpha_features()
    print("\n✓ Phase 1 (Data Integration) complete!")
    print(f"  Output: {OUTPUT_PARQUET}")