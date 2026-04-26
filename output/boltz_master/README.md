# VWF Type-1 & Type-2 Master Structural Batch

This directory contains the aggregated, domain-sliced, and structurally normalized Boltz-2 input file for all discovered Type-1 and Type-2 mutations within the `VWF-ETHos` workspace. 

## 🏗 Architectural Paradigm (Crucial for Downstream Agents)

Unlike naive whole-protein analysis, this pipeline incorporates a **Physics-Aware Domain-Ligand Mapping Strategy**. VWF is a massive 2813-AA multimeric protein, and its mutations affect distinct physiological mechanisms. Thus, this pipeline parses each mutation via a **Dual-Key Route (Subtype × Local Domain)** to generate the correct interaction complex:

| Clinical Subtype | Molecular Consequence | Extracted Domain | Target Ligand Complex Assigned | Boltz-2 Config |
| :--- | :--- | :--- | :--- | :--- | 
| **Type 2B** | Spontaneous Gain-of-Function (High Affinity) | **A1** (1242-1479) | **GPIb_alpha** (Mature LRR, no signal peptide) | Complex (A+B) |
| **Type 2M** | Loss-of-Function (Decreased Affinity) | **A1** | **GPIb_alpha** | Complex (A+B) |
| **Type 2M** | Loss-of-Function (Decreased Collagen Bind) | **A3** (1673-1874) | **COL3A1_TH** (Triple-Helix Peptide) | **Homo-Trimer** (A+B+C+D) |
| **Type 2N** | Loss-of-Function (Deficient FVIII Bind) | **D'-D3** (764-1241) | **Factor VIII a3 acidic region** | Complex (A+B) * |
| **Type 2A** | Cleavage by ADAMTS13 / Unfolding / Shear | **A2** (1480-1672) | *Skipped by Boltz-2 (No Ligand Assigned)* | **FoldX** (Monomer ddG) |
| **Type 1 / 3** | General Secretion / Clearance defects | **Any** | *Skipped by Boltz-2 (No Ligand Assigned)* | **FoldX** (Monomer ddG) |

> ⚠️ *Note on Type 2N / FVIII*: Boltz-2 cannot process Sulfation post-translational modifications out of the box (e.g., FVIII Tyr1680). The absolute `Kd` affinity metric will be artificially low, therefore downstream agents must strictly use the **relative ΔΔG (Mutant vs WT)** to draw scientific conclusions.

## 📝 Execution Protocol

### 1. Engine Core & Generation
If further variations need to be included, they must be formatted into standard `data/raw_tables/` sheets, followed by running:
```bash
# 1. Update routing matrices & extract mutations
python3 scripts/aggregate_all_mutations.py

# 2. Build explicit Boltz-2 topology JSON
python3 Proteo-Structure-Pipeline/domain_ligand_generator/generate_domain_complex.py \
  --csv data/processed/master_type1_type2.csv \
  --wt-fasta structures/wt/VWF_P04275_WT.fasta \
  --format boltz2 --batch-size 9999 --output-dir output/boltz_master
```

### 2. Server Inference (Boltz-2)
Copy the latest batch JSON and execute inference directly:
```bash
scp output/boltz_master/boltz2_batch_01.json user@server:/path/to/workspace/
boltz predict boltz2_batch_01.json --out_dir ./vwf_boltz_results/ --use_msa_server
```

## 🧩 Structure of Generated JSON payload
Structural elements requiring homo-trimers (Collagen COL3A1) will dynamically inherit 3 additional replicated receptor chains, while the VWF segment acts as Chain `A`. An `affinity` property is uniformly injected instructing Boltz-2 models to score the interfacial binding confidence around `A`.

```json
{
  "name": "VWF_P1266L",
  "sequences": [
    { "protein": { "id": "A", "sequence": "..." } },
    { "protein": { "id": "B", "sequence": "..." } },
    ...
  ],
  "properties": [
    { "affinity": { "binder": "A" } }
  ]
}
```
