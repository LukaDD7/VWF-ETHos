# Domain-Ligand Complex Generator for Boltz-2 / AlphaFold 3

This directory contains the standardized pipeline for converting VWF/VWD mutations (like Type 2B) into highly specific, domain-and-ligand-focused complex structures.

## Motivation
Previous workflows included the full 2813 amino acid sequence of VWF in structural predictions. This is computationally expensive, reduces structural accuracy for specific interaction interfaces, and might introduce noise into AlphaFold3 / Boltz-2 attention layers.
This pipeline resolves that by dynamically slicing targeted domains and pairing them with specific ligand sequences.

## 1. Standardized Input Specification

To generate structured JSON files, you **must** provide a standardized CSV file with the following exact columns:

| Column Name | Description | Example |
| :--- | :--- | :--- |
| `Variant_ID` | A unique identifier or name for the mutation | `VWF_Pro1266Leu` |
| `WT_AA` | Wild-type amino acid (1-letter) | `P` |
| `Position` | Full VWF sequence position (1-indexed) | `1266` |
| `Mut_AA` | Mutated amino acid (1-letter) | `L` |
| `Subtype` | VWF mutation subtype | `2B` |
| `Domain_Start`| The 1-indexed start position of the target domain | `1260` (A1 domain start) |
| `Domain_End` | The 1-indexed end position of the target domain | `1479` (A1 domain end) |
| `Ligand_Name` | Name of the interacting protein/ligand | `GPIb_alpha` |
| `Ligand_Seq` | The amino acid sequence of the ligand | `MPLLL...` |

*(Optional columns can be left blank, but the column headers must exist. If `Ligand_Seq` is blank, it generates a monomer.)*

## 2. Pipeline Scripts

### Phase A: Standardization (`standardize_inputs.py`)
Since the project contains messy exported tables (like `2B型突变.xlsx` with varying column names like `Amino acid\nchange`), this script parses those loose formats and creates the strictly defined standard CSV.
**Usage:**
```bash
python standardize_inputs.py --input ../../2B型突变.xlsx --format type2b_excel --output ../../output/standardized_2b.csv
```

### Phase B: Boltz-2 / AF3 JSON Generation (`generate_domain_complex.py`)
This script reads the standardized CSV, applies the mutation to the full WT sequence, slices out the specific domain, and pairs it with the ligand sequence to generate Boltz-2 / AF3 compatible batch files.
**Usage:**
```bash
python generate_domain_complex.py --csv ../../output/standardized_2b.csv --wt-fasta ../../structures/wt/VWF_P04275_WT.fasta --format boltz2
```

## Features
- **Precise Domain Slicing**: Extracts e.g. `[1260:1479]` of the mutated sequence.
- **Support for Boltz-2 Format**: Native support for the `{ "jobs": [ { "name": ..., "sequences": ... } ] }` spec containing `"protein"` objects.
- **Support for AF3 Web Format**: Optionally export the `[ { "name": ..., "sequences": [{"proteinChain": ...}] } ]` array format.
- **Residue Math Safety**: Validates that the mutation position falls within the requested `Domain_Start` and `Domain_End`, and ensures the wild-type residue at `Position` correctly matches the FASTA reference before mutating.
