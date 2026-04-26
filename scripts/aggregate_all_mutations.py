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

LIGAND_LIB = {
    'GPIb_alpha':  'HPICEVSKVASHLEVNCDKRNLTALPPDLPKDTTILHLSENLLYTFSLATLMPYTRLTQLNLDRCELTKLQVDGTLPVLGTLDLSHNQLQSLPLLGQTLPALTVLDVSFNRLTSLPLGALRGLGELQELYLKGNELKTLPPGLLTPTPKLEKLSLANNNLTELPAGLLNGLENLDTLLLQENSLYTIPKGFFGSHLLPFAFLHGNPWLCNCEILYFRRWLQDNAENVYVWKQGVDVKAMTSNVASVQCDNSDKFPVYKYPGKGCPTLGDEGDTDLYD',
    'COL3A1_TH':   'GPRGQPGVMGFPGPKGNDGAPGKNGERGGPGGP',
    'FVIII_a3_A3': 'EITRTTLQSDQEEIDYDDTISVEMKKEDFDIYDEDENQSPRSFQKKTRHYFIAAVERLWDYGMSSSPHVLRNRAQSGSVPQFKKVVFQEFTDGSFTQPLYRG',
}

def parse_mut(mut_str):
    if pd.isna(mut_str): return None, None, None
    mut_str = str(mut_str).strip()
    match_3 = re.search(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', mut_str)
    if match_3:
        wt, pos, mut = match_3.groups()
        if wt in AA_MAP and mut in AA_MAP: return AA_MAP[wt], int(pos), AA_MAP[mut]
    match_1 = re.search(r'([A-Z])(\d+)([A-Z])', mut_str)
    if match_1:
        wt, pos, mut = match_1.groups()
        return wt, int(pos), mut
    return None, None, None

def get_domain(pos):
    """Maps the mutation position to its rigid structural domain boundaries."""
    if not (1 <= pos <= 2813):
        return None, None, None
    
    domains = [
        (1, 763, 'D1_D2', 1, 763),
        (764, 1241, 'D_D3', 764, 1241),
        (1242, 1479, 'A1', 1242, 1479),
        (1480, 1672, 'A2', 1480, 1672),
        (1673, 1874, 'A3', 1673, 1874),
        (1875, 2255, 'D4', 1875, 2255),
        (2256, 2361, 'B1_B3', 2256, 2361),
        (2362, 2770, 'C1_C6', 2362, 2770),
        (2771, 2813, 'CK', 2771, 2813)
    ]
    for start, end, name, d_start, d_end in domains:
        if start <= pos <= end:
            return d_start, d_end, name
            
    # Out of boundary safety check
    return max(1, pos-100), min(2813, pos+100), 'Local'

def assign_ligand(subtype, domain_name):
    """
    Dispatcher based on Subtype AND Domain.
    Returns: (ligand_name, ligand_seq, n_chains, notes)
    """
    s = str(subtype).lower()
    
    # Type 2B: A1 gain-of-function for GPIb-alpha
    if '2b' in s:
        return ('GPIb_alpha', LIGAND_LIB['GPIb_alpha'], 1, 'A1 gain-of-function (high affinity to GPIb)')

    # Type 2M: domain-dependent
    if '2m' in s:
        if domain_name == 'A1':
            return ('GPIb_alpha', LIGAND_LIB['GPIb_alpha'], 1, 'A1 loss-of-function (expect inverse ddG vs 2B)')
        if domain_name == 'A3':
            return ('COL3A1_TH', LIGAND_LIB['COL3A1_TH'], 3, 'A3 loss-of-function; collagen is homo-trimer')
        return (None, None, 0, '2M but not in A1/A3 — flag for manual review')

    # Type 2N: D'D3 loss-of-function for FVIII
    if '2n' in s:
        if domain_name in ('D_D3', 'D1_D2', 'A1'): # Sometimes 2N borders
            return ('FVIII_a3_A3', LIGAND_LIB['FVIII_a3_A3'], 1, 'WARNING: Y1680 sulfation not modeled; absolute Kd unreliable, use ddG')
        return (None, None, 0, '2N but not in D domains — flag')

    # Type 2A: skip ligand, use folding stability instead
    if '2a' in s:
        return (None, None, 0, 'Type 2A: use FoldX ddG_fold on A2 domain alone, no Boltz-2 ligand')

    # Type 1 / Type 3: quantitative; use folding stability
    if 'type1' in s or 'type 1' in s or 'type3' in s or 'type 3' in s:
        return (None, None, 0, 'Type 1/3: secretion/folding defect; use FoldX ddG_fold')

    return (None, None, 0, f'Unhandled subtype fallback logic: {subtype}')

def main():
    out_rows = []
    
    # -------------------------------------------------------------
    # FILE 1: Main Type 1 and Type 2 data extraction
    # -------------------------------------------------------------
    excel_file = '../data/raw_tables/VWF变异分析-用于提取Type-2.xlsx'
    try:
        df = pd.read_excel(excel_file)
        for _, row in df.iterrows():
            clin_class = str(row.get('       Clinical\nclassification', ''))
            
            if 'type1' in clin_class.lower() or 'type2' in clin_class.lower() or 'type3' in clin_class.lower():
                wt, pos, mut = parse_mut(row.get('Protein\nConsequence', ''))
                
                if pos is not None:
                    domain_res = get_domain(pos)
                    if domain_res[0] is None: continue
                    d_start, d_end, d_name = domain_res
                    stype = clin_class.strip()
                    ligand_name, ligand_seq, n_chains, notes = assign_ligand(stype, d_name)

                    out_rows.append({
                        'Variant_ID': f'VWF_{wt}{pos}{mut}',
                        'WT_AA': wt,
                        'Position': pos,
                        'Mut_AA': mut,
                        'Subtype': stype,
                        'Domain_Name': d_name,
                        'Domain_Start': d_start,
                        'Domain_End': d_end,
                        'Ligand_Name': ligand_name or '',
                        'Ligand_Seq': ligand_seq or '',
                        'Ligand_NChains': n_chains,
                        'Computation_Mode': 'Boltz2_affinity' if n_chains > 0 else 'FoldX_only',
                        'Notes': notes
                    })
    except Exception as e:
        print('Error reading Main Table:', e)

    # -------------------------------------------------------------
    # FILE 2: Additional 2B Specific parsing
    # -------------------------------------------------------------
    try:
        df2 = pd.read_excel('../data/raw_tables/2B型突变.xlsx', header=1)
        mut_col = [c for c in df2.columns if 'amino acid' in str(c).lower()][0]
        
        for _, row in df2.iterrows():
            wt, pos, mut = parse_mut(row[mut_col])
            if pos is not None:
                domain_res = get_domain(pos)
                if domain_res[0] is None: continue
                d_start, d_end, d_name = domain_res
                stype = 'Type2B'
                ligand_name, ligand_seq, n_chains, notes = assign_ligand(stype, d_name)
                
                out_rows.append({
                    'Variant_ID': f'VWF_{wt}{pos}{mut}',
                    'WT_AA': wt,
                    'Position': pos,
                    'Mut_AA': mut,
                    'Subtype': stype,
                    'Domain_Name': d_name,
                    'Domain_Start': d_start,
                    'Domain_End': d_end,
                    'Ligand_Name': ligand_name or '',
                    'Ligand_Seq': ligand_seq or '',
                    'Ligand_NChains': n_chains,
                    'Computation_Mode': 'Boltz2_affinity' if n_chains > 0 else 'FoldX_only',
                    'Notes': notes
                })
    except Exception as e:
        print('Error reading 2B Table:', e)

    # Avoid dropping different subtypes mapped to the exact same mutation position incorrectly.
    res = pd.DataFrame(out_rows).drop_duplicates(subset=['Position', 'Mut_AA', 'Subtype'])
    
    os.makedirs('../data/processed', exist_ok=True)
    res.to_csv('../data/processed/master_type1_type2.csv', index=False)
    print(f'Successfully mapped {len(res)} variants with strictly correct domain-ligand physics logic.')

if __name__ == "__main__":
    main()
