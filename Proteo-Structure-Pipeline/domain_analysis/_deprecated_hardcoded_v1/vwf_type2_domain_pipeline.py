#!/usr/bin/env python3
"""
VWF Type-2 Domain-Specific Analysis Pipeline
Main analysis script integrating literature knowledge, AF3 structures, and classification

This pipeline performs:
1. Loads Type-2 VWF variants with known classifications
2. Extracts structural features from AlphaFold3 predictions
3. Classifies variants using literature-based rules
4. Generates comprehensive analysis report

Usage:
    python vwf_type2_domain_pipeline.py \
        --variants /path/to/type2_variants.csv \
        --structures /path/to/af3_structures/ \
        --output /path/to/output/

Author: Claude Code
Date: 2026-04-02
"""

import argparse
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from vwf_type2_literature_based_classifier import (
    VWFType2Classifier,
    VWFVariant,
    VWF_DOMAIN_ARCHITECTURE,
    FUNCTIONAL_REGIONS,
    KNOWN_TYPE2_HOTSPOTS
)

from vwf_structure_feature_extractor import (
    VWFFeatureExtractor,
    DomainFeatures,
    MutationImpactFeatures
)


class VWFType2AnalysisPipeline:
    """
    Main pipeline for VWF Type-2 variant analysis.
    Combines literature knowledge with structural data.
    """

    def __init__(self, wt_structure_path: Optional[str] = None):
        self.classifier = VWFType2Classifier()
        self.extractor = VWFFeatureExtractor(wt_structure_path)

        # Results storage
        self.variants: List[VWFVariant] = []
        self.domain_features: Dict[str, DomainFeatures] = {}
        self.analysis_results: List[Dict] = []

    def load_variants(self, variants_csv: str) -> pd.DataFrame:
        """Load Type-2 variants from CSV."""
        df = pd.read_csv(variants_csv)
        print(f"Loaded {len(df)} variants from {variants_csv}")
        return df

    def process_variant(
        self,
        row: pd.Series,
        structures_dir: Path
    ) -> Optional[Dict]:
        """Process a single variant through the pipeline."""

        # Create variant object
        variant = VWFVariant(
            variant_id=row.get("variant_id", row.get("Name", "")),
            protein_change=row.get("protein_change", row.get("Protein.change", "")),
            position=int(row.get("position", row.get("AA_Position", 0))),
            ref_aa=row.get("ref_aa", row.get("REF_AA", "")),
            alt_aa=row.get("alt_aa", row.get("ALT_AA", "")),
            acmg_classification=row.get("acmg_classification", row.get("INFO_acmg_classification_base", "VUS")),
            type2_subtype=row.get("type2_subtype", row.get("VWD_Type", None))
        )

        # Find structure file
        structure_name = f"fold_vwf_{variant.variant_id.lower()}"
        structure_file = structures_dir / structure_name / f"{variant.variant_id}.cif"

        if not structure_file.exists():
            # Try alternative naming
            structure_file = structures_dir / structure_name / "model.cif"

        if structure_file.exists():
            # Extract structural features
            mut_features = self.extractor.calculate_mutation_impact(
                variant_id=variant.variant_id,
                position=variant.position,
                ref_aa=variant.ref_aa,
                alt_aa=variant.alt_aa,
                mut_structure_path=str(structure_file)
            )

            # Update variant with structural data
            variant.plddt_wt = mut_features.pldtt_delta_site  # Need WT reference
            variant.plddt_mut = mut_features.pldtt_delta_site
            variant.plddt_delta = mut_features.pldtt_delta_site
            variant.local_rmsd = mut_features.local_rmsd
            variant.domain = mut_features.domain_name

        # Classify using literature-based rules
        predicted_subtype, confidence = self.classifier.predict_subtype(variant)
        variant.predicted_subtype = predicted_subtype
        variant.confidence = confidence

        self.variants.append(variant)

        # Compile result
        result = {
            "variant_id": variant.variant_id,
            "protein_change": variant.protein_change,
            "position": variant.position,
            "ref_aa": variant.ref_aa,
            "alt_aa": variant.alt_aa,
            "domain": variant.domain,
            "known_subtype": variant.type2_subtype,
            "predicted_subtype": variant.predicted_subtype,
            "prediction_confidence": variant.confidence,
            "acmg_classification": variant.acmg_classification,
            "is_functional_site": mut_features.is_functional_site if structure_file.exists() else None,
            "is_hotspot": mut_features.is_hotspot if structure_file.exists() else None,
            "plddt_delta": variant.plddt_delta,
            "local_rmsd": variant.local_rmsd,
            "hydrophobicity_change": mut_features.hydrophobicity_change if structure_file.exists() else None,
            "charge_change": mut_features.charge_change if structure_file.exists() else None,
            "size_change": mut_features.size_change if structure_file.exists() else None,
        }

        return result

    def run_analysis(
        self,
        variants_csv: str,
        structures_dir: str,
        output_dir: str
    ) -> pd.DataFrame:
        """Run the complete analysis pipeline."""

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        structures_dir = Path(structures_dir)

        # Load variants
        df = self.load_variants(variants_csv)

        print("\n" + "=" * 80)
        print("Processing Variants")
        print("=" * 80)

        # Process each variant
        results = []
        for idx, row in df.iterrows():
            if idx % 10 == 0:
                print(f"  Processing variant {idx+1}/{len(df)}...")

            result = self.process_variant(row, structures_dir)
            if result:
                results.append(result)

        # Create results DataFrame
        results_df = pd.DataFrame(results)

        # Generate analysis report
        self._generate_report(results_df, output_dir)

        # Save results
        results_df.to_csv(output_dir / "type2_classification_results.csv", index=False)

        return results_df

    def _generate_report(self, results_df: pd.DataFrame, output_dir: Path):
        """Generate comprehensive analysis report."""

        report = []

        report.append("=" * 80)
        report.append("VWF Type-2 Variant Analysis Report")
        report.append("Based on Literature Knowledge and AlphaFold3 Structures")
        report.append("=" * 80)
        report.append("")

        # Summary statistics
        report.append("SUMMARY STATISTICS")
        report.append("-" * 80)
        report.append(f"Total variants analyzed: {len(results_df)}")
        report.append("")

        # Domain distribution
        report.append("Domain Distribution:")
        domain_counts = results_df['domain'].value_counts()
        for domain, count in domain_counts.items():
            report.append(f"  {domain:20s}: {count:3d} variants")
        report.append("")

        # Known vs predicted subtypes
        if 'known_subtype' in results_df.columns and results_df['known_subtype'].notna().any():
            report.append("Known Type-2 Subtypes:")
            known_counts = results_df['known_subtype'].value_counts()
            for subtype, count in known_counts.items():
                report.append(f"  Type {subtype:3s}: {count:3d} variants")
            report.append("")

        report.append("Predicted Type-2 Subtypes:")
        pred_counts = results_df['predicted_subtype'].value_counts()
        for subtype, count in pred_counts.items():
            report.append(f"  Type {subtype:10s}: {count:3d} variants")
        report.append("")

        # High confidence predictions
        high_conf = results_df[results_df['prediction_confidence'] > 0.7]
        report.append(f"High confidence predictions (>70%): {len(high_conf)}")
        report.append("")

        # Domain-specific analysis
        report.append("=" * 80)
        report.append("DOMAIN-SPECIFIC ANALYSIS")
        report.append("=" * 80)
        report.append("")

        # A1 Domain (Type 2B/2M)
        a1_variants = results_df[results_df['domain'] == 'A1']
        if len(a1_variants) > 0:
            report.append("A1 Domain (GPIbα binding) - Associated with Type 2B/2M:")
            report.append(f"  Total variants: {len(a1_variants)}")
            report.append(f"  Functional site variants: {a1_variants['is_functional_site'].sum()}")
            report.append(f"  Hotspot variants: {a1_variants['is_hotspot'].sum()}")
            report.append("")

        # A2 Domain (Type 2A)
        a2_variants = results_df[results_df['domain'] == 'A2']
        if len(a2_variants) > 0:
            report.append("A2 Domain (ADAMTS13 cleavage) - Associated with Type 2A:")
            report.append(f"  Total variants: {len(a2_variants)}")
            report.append(f"  Functional site variants: {a2_variants['is_functional_site'].sum()}")
            report.append(f"  Hotspot variants: {a2_variants['is_hotspot'].sum()}")
            report.append("")

        # A3 Domain (Type 2M)
        a3_variants = results_df[results_df['domain'] == 'A3']
        if len(a3_variants) > 0:
            report.append("A3 Domain (Collagen binding) - Associated with Type 2M:")
            report.append(f"  Total variants: {len(a3_variants)}")
            report.append(f"  Functional site variants: {a3_variants['is_functional_site'].sum()}")
            report.append(f"  Hotspot variants: {a3_variants['is_hotspot'].sum()}")
            report.append("")

        # D'D3 Domain (Type 2N)
        d3_variants = results_df[results_df['domain'].isin(['D_prime', 'D3'])]
        if len(d3_variants) > 0:
            report.append("D'/D3 Domain (FVIII binding) - Associated with Type 2N:")
            report.append(f"  Total variants: {len(d3_variants)}")
            report.append(f"  Functional site variants: {d3_variants['is_functional_site'].sum()}")
            report.append(f"  Hotspot variants: {d3_variants['is_hotspot'].sum()}")
            report.append("")

        # Structural impact analysis
        report.append("=" * 80)
        report.append("STRUCTURAL IMPACT ANALYSIS")
        report.append("=" * 80)
        report.append("")

        if 'plddt_delta' in results_df.columns and results_df['plddt_delta'].notna().any():
            report.append("pLDDT Changes (structural confidence):")
            report.append(f"  Mean pLDDT delta: {results_df['plddt_delta'].mean():.2f}")
            report.append(f"  Variants with pLDDT < -10: {(results_df['plddt_delta'] < -10).sum()}")
            report.append("")

        if 'local_rmsd' in results_df.columns and results_df['local_rmsd'].notna().any():
            report.append("Local RMSD (structural deviation):")
            report.append(f"  Mean RMSD: {results_df['local_rmsd'].mean():.2f} Å")
            report.append(f"  High RMSD (>2Å): {(results_df['local_rmsd'] > 2).sum()}")
            report.append("")

        # Detailed variant list
        report.append("=" * 80)
        report.append("DETAILED CLASSIFICATION RESULTS")
        report.append("=" * 80)
        report.append("")

        for _, row in results_df.iterrows():
            report.append(f"Variant: {row['variant_id']} ({row['protein_change']})")
            report.append(f"  Position: {row['position']} | Domain: {row['domain']}")

            known = row.get('known_subtype', 'N/A')
            pred = row.get('predicted_subtype', 'N/A')
            conf = row.get('prediction_confidence', 0)

            report.append(f"  Known Type: {known} | Predicted: {pred} (confidence: {conf:.2f})")

            if row.get('is_functional_site'):
                report.append(f"  ⚠️  Located in functional site")
            if row.get('is_hotspot'):
                report.append(f"  ⚠️  Located in known Type-2 hotspot")

            report.append("")

        # Save report
        report_text = "\n".join(report)
        with open(output_dir / "analysis_report.txt", "w") as f:
            f.write(report_text)

        print(f"\nReport saved to: {output_dir / 'analysis_report.txt'}")


def main():
    parser = argparse.ArgumentParser(
        description="VWF Type-2 Domain-Specific Analysis Pipeline"
    )
    parser.add_argument(
        "--variants",
        required=True,
        help="CSV file with Type-2 variants"
    )
    parser.add_argument(
        "--structures",
        required=True,
        help="Directory containing AlphaFold3 structures"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        "--wt-structure",
        default=None,
        help="Path to WT VWF structure for comparison (optional)"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("VWF Type-2 Domain-Specific Analysis Pipeline")
    print("=" * 80)
    print(f"Variants: {args.variants}")
    print(f"Structures: {args.structures}")
    print(f"Output: {args.output}")
    print("")

    # Run pipeline
    pipeline = VWFType2AnalysisPipeline(args.wt_structure)
    results = pipeline.run_analysis(
        variants_csv=args.variants,
        structures_dir=args.structures,
        output_dir=args.output
    )

    print("\n" + "=" * 80)
    print("Analysis Complete!")
    print(f"Results saved to: {args.output}")
    print("=" * 80)


if __name__ == "__main__":
    main()
