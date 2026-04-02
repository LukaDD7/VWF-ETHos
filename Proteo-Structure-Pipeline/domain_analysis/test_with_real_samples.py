#!/usr/bin/env python3
"""
Test VWF Type-2 Pipeline with Real Samples
Uses actual AF3 structures from /output/af3_batches_type2
"""

import pandas as pd
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent / "domain_analysis"))

from vwf_type2_literature_based_classifier import VWFType2Classifier, VWFVariant

def load_type2_variants():
    """Load Type-2 variants from the reference table."""
    df = pd.read_csv("/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/VWF_Type2_AF3_Reference_Table.csv")
    # Filter out WT control
    df = df[df['Type2_Subtype'] != 'WT_Control']
    return df

def test_classifier():
    """Test the classifier with real Type-2 variants."""
    print("=" * 80)
    print("Testing VWF Type-2 Classifier with Real Samples")
    print("=" * 80)
    print()

    # Load variants
    df = load_type2_variants()
    print(f"Loaded {len(df)} Type-2 variants")
    print()

    # Initialize classifier
    classifier = VWFType2Classifier()

    # Test a subset of variants
    test_variants = [
        # Type 2A (A2 domain)
        ("G1609R", 1609, "G", "R", "2A"),
        ("D1614N", 1614, "D", "N", "2A"),
        ("R1597W", 1597, "R", "W", "2A"),
        ("I1628T", 1628, "I", "T", "2A"),

        # Type 2M (A3 domain - collagen binding)
        ("S1783A", 1783, "S", "A", "2M"),
        ("W1745C", 1745, "W", "C", "2M"),
        ("S1731L", 1731, "S", "L", "2M"),

        # Type 2A (other)
        ("S2775C", 2775, "S", "C", "2A"),
    ]

    print("Detailed Classification Tests:")
    print("-" * 80)

    correct = 0
    total = 0

    for var_id, pos, ref, alt, known_type in test_variants:
        # Get classification by position
        result = classifier.classify_by_position(pos)

        # Create variant and predict
        variant = VWFVariant(
            variant_id=var_id,
            protein_change=f"p.{ref}{pos}{alt}",
            position=pos,
            ref_aa=ref,
            alt_aa=alt,
            acmg_classification="Pathogenic",
            type2_subtype=known_type
        )

        predicted, confidence = classifier.predict_subtype(variant)

        # Check if prediction matches known type
        match = (predicted.upper() == known_type.upper())
        if match:
            correct += 1
        total += 1

        print(f"\nVariant: {var_id} (Position {pos})")
        print(f"  Known Type:     {known_type}")
        print(f"  Predicted Type: {predicted} (confidence: {confidence:.2f})")
        print(f"  Domain:         {result['domain']}")
        print(f"  Match:          {'✓ YES' if match else '✗ NO'}")

        if result['associations']:
            print(f"  Associations:")
            for assoc in result['associations']:
                print(f"    - Type {assoc['type']}: {assoc['mechanism'][:60]}...")

    print()
    print("-" * 80)
    print(f"Accuracy: {correct}/{total} = {correct/total*100:.1f}%")
    print()

    # Full analysis of all variants
    print("=" * 80)
    print("Full Dataset Analysis")
    print("=" * 80)
    print()

    # Count by domain
    print("Variants by VWF Domain:")
    domain_counts = df['VWF_Domain'].value_counts()
    for domain, count in domain_counts.items():
        print(f"  {domain:40s}: {count:3d} variants")
    print()

    # Count by Type-2 subtype
    print("Variants by Type-2 Subtype:")
    type_counts = df['Type2_Subtype'].value_counts()
    for subtype, count in type_counts.items():
        print(f"  Type {subtype:6s}: {count:3d} variants")
    print()

    # Count by ACMG classification
    print("Variants by ACMG Classification:")
    acmg_counts = df['ACMG_Classification'].value_counts()
    for acmg, count in acmg_counts.items():
        print(f"  {acmg:20s}: {count:3d} variants")
    print()

    # Check if structures exist
    print("=" * 80)
    print("AF3 Structure Availability Check")
    print("=" * 80)
    print()

    structures_dir = Path("/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/Proteo-Structure-Pipeline/output/af3_batches_type2/AF3_Results/extracted")

    found = 0
    not_found = []

    for _, row in df.head(20).iterrows():  # Check first 20
        aa_change = row['AA_Change']
        var_id = f"fold_vwf_{aa_change.lower()}"
        var_dir = structures_dir / var_id

        if var_dir.exists():
            found += 1
            # Count CIF files
            cif_files = list(var_dir.glob("*.cif"))
            print(f"  {aa_change:10s}: ✓ Found ({len(cif_files)} models)")
        else:
            not_found.append(aa_change)
            print(f"  {aa_change:10s}: ✗ Not found")

    print()
    print(f"Structures found: {found}/20")
    if not_found:
        print(f"Missing: {', '.join(not_found[:5])}...")
    print()

    print("=" * 80)
    print("Test Complete!")
    print("=" * 80)

if __name__ == "__main__":
    test_classifier()
