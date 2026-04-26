#!/usr/bin/env python3
import argparse
import json
import os
import pandas as pd

def generate_boltz2_batch(mutations_df, wt_seq, batch_size=20):
    """
    Generates batches in Boltz-2 / AF3 complex format.
    Format: {"name": "batch_name", "jobs": [ { "name": ..., "sequences": [...] } ] }
    """
    batch_data_list = []
    
    current_batch_index = 1
    current_jobs = []
    
    for idx, row in mutations_df.iterrows():
        variant_id = row['Variant_ID']
        wt_aa = row['WT_AA']
        pos = int(row['Position'])
        mut_aa = row['Mut_AA']
        domain_start = int(row['Domain_Start'])
        domain_end = int(row['Domain_End'])
        # Extract newly placed variables
        ligand_name = row.get('Ligand_Name', '')
        ligand_seq = row.get('Ligand_Seq', '')
        n_chains = row.get('Ligand_NChains', 0)
        
        # Validation bounds
        if pd.isna(n_chains): n_chains = 0
        n_chains = int(n_chains)
        
        # Bypass monomers if we are strictly producing Boltz-2 complex batches
        # Wait, if `n_chains` == 0, it is a FoldX/Stability type, we can still generate the monomer JSON for Boltz-2 if we want, or we can skip. The user said type 2A and type 1 should bypass.
        if n_chains == 0:
            print(f"Skipping {variant_id} for Boltz2 complex generation (marked for FoldX stability analysis).")
            continue

        # Validation
        if len(wt_seq) < pos:
            print(f"Warning: Pos {pos} out of bounds for {variant_id}. Skipping.")
            continue
            
        actual_wt = wt_seq[pos - 1]
        if actual_wt != wt_aa:
            print(f"Warning: WT mismatch for {variant_id}. Expected {wt_aa}, found {actual_wt}. Skipping.")
            continue
            
        if pos < domain_start or pos > domain_end:
            print(f"Warning: Mutation {variant_id} is outside the specified domain slice ({domain_start}-{domain_end}). Skipping.")
            continue
            
        # Apply mutation to full sequence
        mutated_full_seq = wt_seq[:pos-1] + mut_aa + wt_seq[pos:]
        
        # Slice domain
        # Domain bounds are 1-indexed, so start-1 to end
        sliced_domain = mutated_full_seq[domain_start-1:domain_end]
        
        # Construct Sequences array
        sequences = []
        
        # Add primary domain
        sequences.append({
            "protein": {
                "id": "A",
                "sequence": sliced_domain,
            }
        })
        
        # Add ligand chains if present (B, C, D...)
        current_chain_idx = ord('B')
        for _ in range(n_chains):
            chain_id = chr(current_chain_idx)
            sequences.append({
                "protein": {
                    "id": chain_id,
                    "sequence": str(ligand_seq).strip(),
                }
            })
            current_chain_idx += 1
            
        job = {
            "name": variant_id,
            "sequences": sequences,
            "properties": [{"affinity": {"binder": "A"}}]
        }
        
        current_jobs.append(job)
        
        # Split logic
        if len(current_jobs) >= batch_size:
            batch_data_list.append({
                "name": f"boltz2_batch_{current_batch_index:02d}",
                "jobs": current_jobs
            })
            current_batch_index += 1
            current_jobs = []
            
    if current_jobs:
        batch_data_list.append({
            "name": f"boltz2_batch_{current_batch_index:02d}",
            "jobs": current_jobs
        })
        
    return batch_data_list

def generate_af3_web_batch(mutations_df, wt_seq, batch_size=20):
    """
    Generates batches in AF3 Web Server Batch format.
    Format: [ { "name": ..., "sequences": [{"proteinChain": {"sequence": "...", "count": 1}}] } ]
    """
    batch_data_list = []
    current_batch_index = 1
    current_jobs = []
    
    for idx, row in mutations_df.iterrows():
        variant_id = row['Variant_ID']
        wt_aa = row['WT_AA']
        pos = int(row['Position'])
        mut_aa = row['Mut_AA']
        domain_start = int(row['Domain_Start'])
        domain_end = int(row['Domain_End'])
        ligand_seq = row.get('Ligand_Seq', '')
        
        # Apply bounds & mutations
        if len(wt_seq) < pos or wt_seq[pos - 1] != wt_aa or pos < domain_start or pos > domain_end:
            continue
            
        mutated_full_seq = wt_seq[:pos-1] + mut_aa + wt_seq[pos:]
        sliced_domain = mutated_full_seq[domain_start-1:domain_end]
        
        sequences = []
        sequences.append({
            "proteinChain": {
                "sequence": sliced_domain,
                "count": 1
            }
        })
        
        if pd.notna(ligand_seq) and str(ligand_seq).strip():
            sequences.append({
                "proteinChain": {
                    "sequence": str(ligand_seq).strip(),
                    "count": 1
                }
            })
            
        job = {
            "name": variant_id,
            "sequences": sequences,
            "properties": [
                {
                    "affinity": { "binder": "A" }
                }
            ]
        }
        
        current_jobs.append(job)
        
        if len(current_jobs) >= batch_size:
            batch_data_list.append(current_jobs)
            current_batch_index += 1
            current_jobs = []
            
    if current_jobs:
        batch_data_list.append(current_jobs)
        
    return batch_data_list

def read_fasta(file_path):
    seq = ""
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq += line.strip()
    return seq

def main():
    parser = argparse.ArgumentParser(description="Generate domain-and-ligand focused structure inputs.")
    parser.add_argument('--csv', required=True, help='Path to the standardized CSV file')
    parser.add_argument('--wt-fasta', required=True, help='Path to the Wild-Type FASTA record for VWF')
    parser.add_argument('--format', choices=['boltz2', 'af3_web'], default='boltz2', 
                        help='Target JSON format (boltz2 handles job arrays, af3_web handles simple lists)')
    parser.add_argument('--output-dir', default='./generated_batches', help='Directory to save JSON batches')
    parser.add_argument('--batch-size', type=int, default=30, help='Max number of complex models per batch file')
    
    args = parser.parse_args()
    
    print(f"Loading Wild Type sequence from {args.wt_fasta}")
    try:
        wt_sequence = read_fasta(args.wt_fasta)
    except Exception as e:
        print(f"Error loading FASTA: {e}")
        return
        
    print(f"Loading standardized input CSV from {args.csv}")
    df = pd.read_csv(args.csv)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.format == 'boltz2':
        batches = generate_boltz2_batch(df, wt_sequence, args.batch_size)
    else:
        batches = generate_af3_web_batch(df, wt_sequence, args.batch_size)
        
    total_jobs = 0
    for i, batch_data in enumerate(batches):
        out_path = os.path.join(args.output_dir, f"{args.format}_batch_{i+1:02d}.json")
        with open(out_path, 'w') as f:
            json.dump(batch_data, f, indent=2)
        
        num_jobs = len(batch_data['jobs']) if args.format == 'boltz2' else len(batch_data)
        total_jobs += num_jobs
        print(f"Saved {out_path} [{num_jobs} jobs]")
        
    print(f"\nSuccessfully generated {len(batches)} batches containing {total_jobs} total structural predictions.")
    print("Optimization: Included only specified structure domains & ligands to improve prediction quality.")

if __name__ == '__main__':
    main()
