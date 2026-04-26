# VWF Type-1 & Type-2 Master Boltz-2 Batch

This directory contains the aggregated and normalized Boltz-2 input file for all discovered Type-1 and Type-2 mutations in the `VWF-ETHos` workspace.

## Details
- **Total Valid Mutations**: 88
- **Strategy**: Instead of utilizing the full 2813 amino acid sequence of VWF, this batch uses dynamically generated *Bounding Boxes* (Domain slices) based on the mutation's position. This ensures Boltz-2 focuses exclusively on the local structural domain (e.g. A1, A2, A3, D1-D4 domains) where the mutation occurs, vastly reducing inference time and potential noise.
- **Complex Formation**: For relevant subtypes like Type 2B (which heavily modifies interactions in the A1 domain), the `GPIb_alpha` ligand has been automatically paired as Chain B in the complex to simulate structural dynamics correctly.

## How to execute on the server

1. **Transfer the batch file**
Copy `boltz2_batch_01.json` to your inference server where Boltz-2 is installed.
```bash
scp output/boltz_master/boltz2_batch_01.json user@server:/path/to/boltz/workspace/
```

2. **Run Boltz-2**
Simply invoke the Boltz-2 prediction pipeline passing the generated standard JSON:
```bash
boltz predict boltz2_batch_01.json --out_dir ./vwf_boltz_results/ --use_msa_server
```

3. **Check Results**
Boltz-2 will create a directory for each mutation (e.g., `VWF_P1266L`) inside `--out_dir`, containing the predicted `.cif` (or `.pdb`) structures, PAE matrices, and confidence scores (pLDDT).

## Structure
All variants in the JSON are defined strictly as:
```json
{
  "name": "VWF_...",
  "sequences": [
    {
      "chain_name": "A",
      "protein": { "sequence": "...", ... }
    },
    // Optional ligand for complex tasks
    {
      "chain_name": "B",
      "protein": { "sequence": "...", ... }
    }
  ]
}
```

## Directory Organization Changes
The `VWF-ETHos` structure has also been refactored following standard SE rules to ensure a clean workspace for future operations:
- `data/raw_tables/`: Contains old exported mutation Excel sheets.
- `data/processed/`: Contains normalized CSVs like the `master_type1_type2.csv` used for this build.
- `scripts/`: Holds quick python operations.
- `archived/`: Deprecated or faulty historical JSONs/tables are securely shelved here.
