#!/usr/bin/env python3
"""
Full Pipeline Test with AF3 Structures
Tests feature extraction and classification with real CIF files
"""

import pandas as pd
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "domain_analysis"))

from vwf_type2_literature_based_classifier import VWFType2Classifier, VWFVariant
from vwf_structure_feature_extractor import VWFFeatureExtractor

def test_feature_extraction():
    """Test feature extraction with actual AF3 structures."""
    print("=" * 80)
    print("Testing Feature Extraction with AF3 Structures")
    print("=" * 80)
    print()

    # Path to structures
    structures_dir = Path("/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/Proteo-Structure-Pipeline/output/af3_batches_type2/AF3_Results/extracted")

    # Initialize extractor (without WT for now)
    extractor = VWFFeatureExtractor()

    # Test variants with different characteristics
    test_cases = [
        ("D1614N", 1614, "D", "N", "2A"),
        ("R1597W", 1597, "R", "W", "2A"),
        ("S1783A", 1783, "S", "A", "2M"),
        ("W1745C", 1745, "W", "C", "2M"),
    ]

    print("Extracting structural features from AF3 predictions:")
    print("-" * 80)

    for var_id, pos, ref, alt, known_type in test_cases:
        print(f"\n{var_id} (Position {pos}, Type {known_type}):")

        # Find structure directory
        var_dir = structures_dir / f"fold_vwf_{var_id.lower()}"

        if not var_dir.exists():
            print(f"  ✗ Structure directory not found")
            continue

        # Find CIF files
        cif_files = list(var_dir.glob("*_model_*.cif"))
        if not cif_files:
            print(f"  ✗ No CIF files found")
            continue

        # Use first model
        cif_file = cif_files[0]
        print(f"  Structure: {cif_file.name}")

        # Extract domain features
        try:
            domain_features = extractor.extract_domain_features(str(cif_file))

            # Find domain for this position
            for domain_name, features in domain_features.items():
                if features.start_residue <= pos <= features.end_residue:
                    print(f"  Domain: {domain_name}")
                    print(f"    pLDDT mean: {features.plddt_mean:.2f}")
                    print(f"    pLDDT range: [{features.plddt_min:.2f}, {features.plddt_max:.2f}]")
                    print(f"    Low confidence regions: {len(features.low_confidence_regions)}")
                    print(f"    Contains functional site: {features.contains_functional_site}")
                    if features.functional_site_names:
                        print(f"    Functional sites: {', '.join(features.functional_site_names)}")
                    break

            # Calculate mutation impact features
            mut_features = extractor.calculate_mutation_impact(
                variant_id=var_id,
                position=pos,
                ref_aa=ref,
                alt_aa=alt,
                mut_structure_path=str(cif_file)
            )

            print(f"  Mutation impact:")
            print(f"    Hydrophobicity change: {mut_features.hydrophobicity_change:+.2f}")
            print(f"    Charge change: {mut_features.charge_change:+d}")
            print(f"    Size change: {mut_features.size_change:+.1f} Å³")
            print(f"    Is functional site: {mut_features.is_functional_site}")
            print(f"    Is hotspot: {mut_features.is_hotspot}")

        except Exception as e:
            print(f"  ✗ Error extracting features: {e}")

    print()
    print("=" * 80)
    print("Feature Extraction Test Complete")
    print("=" * 80)


def test_full_pipeline():
    """Test the complete pipeline with a subset of variants."""
    print()
    print("=" * 80)
    print("Full Pipeline Integration Test")
    print("=" * 80)
    print()

    # Load variant data
    df = pd.read_csv("/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/VWF_Type2_AF3_Reference_Table.csv")
    df = df[df['Type2_Subtype'] != 'WT_Control']

    classifier = VWFType2Classifier()
    structures_dir = Path("/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/Proteo-Structure-Pipeline/output/af3_batches_type2/AF3_Results/extracted")

    # Test on first 10 variants
    test_df = df.head(10)

    results = []

    print("Classifying variants with literature-based rules:")
    print("-" * 80)

    for _, row in test_df.iterrows():
        aa_change = row['AA_Change']
        pos = int(row['Position'])  # Convert to int
        known_type = row['Type2_Subtype']
        domain = row['VWF_Domain']

        # Parse AA change
        ref_aa = aa_change[0]
        alt_aa = aa_change[-1]

        # Create variant
        variant = VWFVariant(
            variant_id=aa_change,
            protein_change=f"p.{ref_aa}{pos}{alt_aa}",
            position=pos,
            ref_aa=ref_aa,
            alt_aa=alt_aa,
            acmg_classification=row['ACMG_Classification'],
            type2_subtype=known_type
        )

        # Predict
        predicted, confidence = classifier.predict_subtype(variant)

        # Check structure availability
        var_dir = structures_dir / f"fold_vwf_{aa_change.lower()}"
        has_structure = var_dir.exists()

        # Normalize type formats for comparison
        # known_type: "type2A" -> "2A", predicted: "2M" -> "2M"
        known_normalized = known_type.upper().replace("TYPE", "")
        predicted_normalized = predicted.upper()
        match = (predicted_normalized == known_normalized)

        results.append({
            'variant': aa_change,
            'position': pos,
            'domain': domain,
            'known': known_type,
            'predicted': predicted,
            'confidence': confidence,
            'match': match,
            'has_structure': has_structure
        })

        status = "✓" if match else "✗"
        struct_status = "📦" if has_structure else "❌"
        print(f"  {status} {struct_status} {aa_change:10s} (Pos {pos:4d}): Known={known_type:6s}, Predicted={predicted:6s} ({confidence:.2f})")

    # Summary
    print()
    print("Summary:")
    print(f"  Total tested: {len(results)}")
    print(f"  Correct predictions: {sum(r['match'] for r in results)}")
    print(f"  Accuracy: {sum(r['match'] for r in results)/len(results)*100:.1f}%")
    print(f"  Structures available: {sum(r['has_structure'] for r in results)}/{len(results)}")

    print()
    print("=" * 80)


if __name__ == "__main__":
    test_feature_extraction()
    test_full_pipeline()
