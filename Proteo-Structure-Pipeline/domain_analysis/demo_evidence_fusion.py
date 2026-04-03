#!/usr/bin/env python3
"""
Quick Demo: Evidence-Based Classification for VWF Type-2 Variants
Shows how literature evidence enhances classification robustness and interpretability
"""

import sys
sys.path.insert(0, '/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/Proteo-Structure-Pipeline/domain_analysis')

from vwf_evidence_based_classifier import EvidenceBasedClassifier


def print_section(title):
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80 + "\n")


def demo_comparison():
    """Compare basic vs evidence-based classification."""

    print_section("Comparison: Basic vs Evidence-Based Classification")

    classifier = EvidenceBasedClassifier()

    # Test variants from different domains
    test_cases = [
        ("V89A", 89, "V", "A", "D1-D2 (Propeptide)"),
        ("D1614N", 1614, "D", "N", "A2 (ADAMTS13 cleavage)"),
        ("R1306W", 1306, "R", "W", "A1 (GPIbα binding)"),
        ("P1888L", 1888, "P", "L", "D4 (Multimerization)"),
        ("R816W", 816, "R", "W", "D3 (FVIII binding)"),
    ]

    print(f"{'Variant':<12} {'Domain':<25} {'Basic Conf':<12} {'Evidence':<10} {'Enhanced':<12} {'#Papers':<8}")
    print("-" * 80)

    for var_id, pos, ref, alt, domain_desc in test_cases:
        predicted, confidence, evidence = classifier.classify_with_evidence(
            var_id, pos, ref, alt
        )

        # Calculate evidence boost
        evidence_score = evidence.get_evidence_score()
        evidence_boost = min(evidence_score * 0.1, 0.2)
        enhanced = min(confidence + evidence_boost, 1.0)

        n_papers = len(evidence.supporting_papers)

        print(f"{var_id:<12} {domain_desc:<25} {confidence*100:>6.1f}%     {evidence_score:>5.2f}     {enhanced*100:>6.1f}%     {n_papers:>3}")


def demo_explainability():
    """Demonstrate explainable classification output."""

    print_section("Explainable Classification: V89A (Propeptide Variant)")

    classifier = EvidenceBasedClassifier()

    # Generate full report
    report = classifier.generate_explainable_report("V89A", 89, "V", "A")
    print(report)


def demo_evidence_quality():
    """Show evidence quality metrics."""

    print_section("Evidence Quality Analysis by Domain")

    classifier = EvidenceBasedClassifier()

    domains = [
        ("propeptide_D1", "D1 Propeptide"),
        ("D4", "D4 Domain"),
        ("A2", "A2 Domain"),
        ("A1", "A1 Domain"),
        ("CK", "Cystine Knot"),
    ]

    print(f"{'Domain':<20} {'#Papers':<10} {'Avg Year':<12} {'Avg Citations':<15} {'Evidence Score':<15}")
    print("-" * 80)

    for domain_key, domain_name in domains:
        if domain_key in classifier.evidence_db:
            papers = classifier.evidence_db[domain_key]
            n_papers = len(papers)
            avg_year = sum(p.year for p in papers) / n_papers
            avg_citations = sum(p.citation_count for p in papers) / n_papers
            avg_confidence = sum(p.confidence for p in papers) / n_papers

            print(f"{domain_name:<20} {n_papers:<10} {avg_year:>6.0f}      {avg_citations:>8.0f}        {avg_confidence*100:>6.1f}%")


def demo_clinical_recommendations():
    """Show domain-specific clinical recommendations."""

    print_section("Clinical Recommendations by Domain")

    classifier = EvidenceBasedClassifier()

    # Map domains to recommendations
    recommendations = {
        "propeptide_D1": [
            "Order VWF multimer analysis (expect: pronounced dimer band)",
            "Check VWFpp/Ag ratio",
            "Consider propeptide supplementation therapy (experimental)",
        ],
        "D4": [
            "Check VWFpp/Ag ratio (expect: 1-<2 for secretion defect)",
            "Multimer analysis (may be normal or smeary)",
            "Rule out Type 1 misdiagnosis",
        ],
        "A2": [
            "Multimer analysis (expect: loss of large multimers)",
            "ADAMTS13 activity and inhibitor screen",
        ],
        "A1": [
            "RIPA test (expect: increased aggregation at low ristocetin)",
            "Platelet count (may have thrombocytopenia)",
        ],
    }

    for domain, recs in recommendations.items():
        print(f"\n{domain.upper()} Domain:")
        for i, rec in enumerate(recs, 1):
            print(f"  {i}. {rec}")


def demo_uncertainty_quantification():
    """Show uncertainty quantification for different scenarios."""

    print_section("Uncertainty Quantification")

    classifier = EvidenceBasedClassifier()

    scenarios = [
        ("V89A", 89, "V", "A", "Strong evidence (5 papers)"),
        ("P1888L", 1888, "P", "L", "Moderate evidence (3 papers)"),
        ("R1986C", 1986, "R", "C", "D4 variant (same evidence)"),
    ]

    print(f"{'Variant':<12} {'Evidence Level':<30} {'Confidence':<12} {'Uncertainty':<12}")
    print("-" * 70)

    for var_id, pos, ref, alt, desc in scenarios:
        predicted, confidence, evidence = classifier.classify_with_evidence(
            var_id, pos, ref, alt
        )

        # Calculate uncertainty
        n_papers = len(evidence.supporting_papers)
        evidence_coverage = min(n_papers / 5, 1.0)  # Normalize to 5 papers
        uncertainty = 1 - (confidence * evidence_coverage)

        print(f"{var_id:<12} {desc:<30} {confidence*100:>6.1f}%      {uncertainty*100:>6.1f}%")


def main():
    """Run all demos."""

    print("\n" + "=" * 80)
    print("VWF Type-2 Evidence-Based Classification Pipeline - Demo")
    print("=" * 80)
    print("\nThis demo shows how literature evidence enhances:")
    print("  1. Classification robustness (multi-source evidence)")
    print("  2. Interpretability (mechanism explanations)")
    print("  3. Clinical utility (domain-specific recommendations)")
    print("  4. Uncertainty quantification (evidence coverage)")

    demo_comparison()
    demo_explainability()
    demo_evidence_quality()
    demo_clinical_recommendations()
    demo_uncertainty_quantification()

    print_section("Demo Complete")
    print("Files created:")
    print("  - vwf_evidence_based_classifier.py  (Evidence-based classifier)")
    print("  - INTEGRATION_GUIDE.md              (Integration documentation)")
    print("\nTo integrate into main pipeline, see INTEGRATION_GUIDE.md")


if __name__ == "__main__":
    main()
