#!/usr/bin/env python3
import argparse
import pandas as pd
import re
import os

# Dictionary to convert 3-letter amino acid codes to 1-letter codes
AA_MAP = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
}

def parse_mutation_string(mut_str):
    """
    Parses a string like 'Pro1266Leu' or 'P1266L' into (WT, Pos, Mut).
    """
    if pd.isna(mut_str):
        return None, None, None
    
    mut_str = str(mut_str).strip()
    
    # Try 3-letter code match: e.g., Pro1266Leu
    match_3 = re.match(r'^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', mut_str)
    if match_3:
        wt3, pos, mut3 = match_3.groups()
        if wt3 in AA_MAP and mut3 in AA_MAP:
            return AA_MAP[wt3], int(pos), AA_MAP[mut3]

    # Try 1-letter code match: e.g., P1266L
    match_1 = re.match(r'^([A-Z])(\d+)([A-Z])$', mut_str)
    if match_1:
        wt, pos, mut = match_1.groups()
        return wt, int(pos), mut
        
    return None, None, None

def main():
    parser = argparse.ArgumentParser(description="Standardize chaotic VWF mutation tables into pipeline format")
    parser.add_argument('--input', required=True, help='Path to the input unstructured Excel or CSV file')
    parser.add_argument('--output', required=True, help='Path to save the standardized CSV file')
    parser.add_argument('--format', choices=['type2b_excel', 'generic_csv'], default='generic_csv',
                        help='Pre-defined logic for parsing legacy tables')
    parser.add_argument('--subtype', default='Type2', help='Subtype label to attach to entries')
    parser.add_argument('--domain-start', type=int, default=1260, help='Default 1-indexed domain start')
    parser.add_argument('--domain-end', type=int, default=1479, help='Default 1-indexed domain end')
    parser.add_argument('--ligand-name', type=str, default='GPIb_alpha', help='Name of generic ligand')
    parser.add_argument('--ligand-seq', type=str, 
                        default='MPLLLLLLLLPSPLHPHPICEVSKVASHLEVNCDKRNLTALPPDLPKDTTILHLSENLLYTFSLATLMPYTRLTQLNLDRCELTKLQVDGTLPVLGTLDLSHNQLQSLPLLGQTLPALTVLDVSFNRLTSLPLGALRGLGELQELYLKGNELKTLPPGLLTPTPKLEKLSLANNNLTELPAGLLNGLENLDTLLLQENSLYTIPKGFFGSHLLPFAFLHGNPWLCNCEILYFRRWLQDNAENVYVWKQGVDVKAMTSNVASVQCDNSDKFPVYKYPGKGCPTLGDEGDTDLYD',
                        help='Sequence of the generic ligand (e.g. GPIb alpha for 2B)')
    
    args = parser.parse_args()
    
    print(f"Reading input file: {args.input}")
    
    try:
        if args.input.endswith('.xlsx'):
            df_in = pd.read_excel(args.input, header=1) # Many legacy tables have title in row 0
        else:
            df_in = pd.read_csv(args.input)
    except Exception as e:
        print(f"Failed to read file. Trying without header skip. Error: {e}")
        if args.input.endswith('.xlsx'):
            df_in = pd.read_excel(args.input)
        else:
            df_in = pd.read_csv(args.input)

    output_rows = []
    
    # Logic for parsing Type 2B specific legacy excel format
    if args.format == 'type2b_excel':
        # Find the column containing 'Amino acid change' or similar
        mut_col = None
        for col in df_in.columns:
            if 'amino acid' in str(col).lower() and 'change' in str(col).lower():
                mut_col = col
                break
        
        if not mut_col:
            # Fallback
            mut_col = df_in.columns[0]
            print(f"Could not explicitly identify mutation column, defaulting to first column: {mut_col}")
            
        for _, row in df_in.iterrows():
            mut_str = row[mut_col]
            wt, pos, mut = parse_mutation_string(mut_str)
            
            if pos is not None:
                output_rows.append({
                    'Variant_ID': f"VWF_{wt}{pos}{mut}",
                    'WT_AA': wt,
                    'Position': pos,
                    'Mut_AA': mut,
                    'Subtype': args.subtype,
                    'Domain_Start': args.domain_start,
                    'Domain_End': args.domain_end,
                    'Ligand_Name': args.ligand_name,
                    'Ligand_Seq': args.ligand_seq
                })
    
    # Generic CSV handling assumes there's an 'AA_change', 'Mutation' or 'Variant' column
    elif args.format == 'generic_csv':
        potential_cols = ['AA_change', 'Mutation', 'Variant', 'mut']
        mut_col = None
        for col in df_in.columns:
            if any(p.lower() in str(col).lower() for p in potential_cols):
                mut_col = col
                break
                
        if not mut_col:
            print(f"Available columns: {df_in.columns.tolist()}")
            raise ValueError("Could not auto-identify mutation column and --format generic_csv was specified.")
            
        for _, row in df_in.iterrows():
            mut_str = row[mut_col]
            wt, pos, mut = parse_mutation_string(mut_str)
            
            if pos is not None:
                output_rows.append({
                    'Variant_ID': f"VWF_{wt}{pos}{mut}",
                    'WT_AA': wt,
                    'Position': pos,
                    'Mut_AA': mut,
                    'Subtype': row.get('Subtype', args.subtype),
                    'Domain_Start': row.get('Domain_Start', args.domain_start),
                    'Domain_End': row.get('Domain_End', args.domain_end),
                    'Ligand_Name': row.get('Ligand_Name', args.ligand_name),
                    'Ligand_Seq': row.get('Ligand_Seq', args.ligand_seq)
                })

    df_out = pd.DataFrame(output_rows)
    # Deduplicate
    df_out.drop_duplicates(subset=['Position', 'Mut_AA'], inplace=True)
    
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    df_out.to_csv(args.output, index=False)
    
    print(f"Successfully processed {len(df_out)} unique mutations.")
    print(f"Standardized mapping saved to: {args.output}")

if __name__ == '__main__':
    main()
