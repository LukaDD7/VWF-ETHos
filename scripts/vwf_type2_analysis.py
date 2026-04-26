#!/usr/bin/env python3
"""
VWF Type-2 可解释分诊工具

目标不是直接给出不可解释的终诊，而是输出：
DNA变异/位置 → 机制分数 → Type-2亚型倾向 → 推荐验证实验

示例：
  python vwf_type2_analysis.py --position 1437 --aa-change p.Arg1437Trp
  python vwf_type2_analysis.py --batch-predict variants.csv --output predictions.csv
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from typing import Any

import pandas as pd


DOMAINS = {
    'D1-D2': (1, 272),
    "D'D3": (273, 761),
    'A1': (762, 1035),
    'A2': (1036, 1227),
    'A3': (1228, 1451),
    'D4': (1452, 1670),
    'C1-C2': (1671, 2051),
    'CT': (2052, 2813),
}

SUBTYPES = ['type2A', 'type2B', 'type2M', 'type2N', 'uncertain']

# 机制 -> 亚型的映射强度
MECHANISM_TO_SUBTYPE = {
    'multimer_assembly_defect': {'type2A': 1.0},
    'adamts13_susceptibility': {'type2A': 0.9},
    'gp1b_gain_of_function': {'type2B': 1.0},
    'gp1b_loss_of_function': {'type2M': 0.8},
    'collagen_binding_loss': {'type2M': 0.9},
    'fviii_binding_loss': {'type2N': 1.0},
    'secretion_defect': {'type2A': 0.7, 'uncertain': 0.3},
    'splice_regulatory_risk': {'uncertain': 0.8, 'type2A': 0.2},
}

DOMAIN_PRIORS = {
    "D'D3": {
        'fviii_binding_loss': 0.80,
        'multimer_assembly_defect': 0.35,
        'secretion_defect': 0.20,
    },
    'A1': {
        'gp1b_gain_of_function': 0.62,
        'gp1b_loss_of_function': 0.48,
        'fviii_binding_loss': 0.18,
    },
    'A2': {
        'adamts13_susceptibility': 0.72,
        'multimer_assembly_defect': 0.35,
    },
    'A3': {
        'collagen_binding_loss': 0.72,
        'multimer_assembly_defect': 0.18,
    },
    'D4': {
        'multimer_assembly_defect': 0.62,
        'secretion_defect': 0.30,
    },
    'D1-D2': {
        'secretion_defect': 0.60,
        'multimer_assembly_defect': 0.20,
    },
    'C1-C2': {
        'collagen_binding_loss': 0.20,
        'gp1b_loss_of_function': 0.15,
    },
    'CT': {
        'multimer_assembly_defect': 0.28,
        'secretion_defect': 0.24,
    },
}

CLINICAL_RECOMMENDATIONS = {
    'type2A': [
        '建议优先做 VWF multimer analysis',
        '必要时补做 ADAMTS13-related phenotype / proteolysis 验证',
    ],
    'type2B': [
        '建议优先做 low-dose RIPA',
        '联合评估 platelet count / thrombocytopenia',
    ],
    'type2M': [
        '建议优先做 VWF activity/Ag ratio',
        '补做 collagen binding 或 GPIb binding 功能实验',
    ],
    'type2N': [
        '建议优先做 FVIII binding assay',
        '同时与轻型血友病A做鉴别',
    ],
    'uncertain': [
        '建议优先补充实验表型后再判断亚型',
        '如怀疑剪接/调控异常，可单独走 AlphaGenome / RNA 路径解释',
    ],
}


@dataclass
class VariantEvidence:
    position: int
    domain: str
    aa_change: str | None
    zygosity: str | None
    splice_score: float | None
    pae_delta: float | None
    notes: list[str]


def clamp(value: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(high, value))


def get_domain(position: int) -> str:
    for name, (start, end) in DOMAINS.items():
        if start <= position <= end:
            return name
    return 'Unknown'


def normalize_domain(value: Any) -> str | None:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    text = str(value).strip()
    if not text:
        return None

    aliases = {
        "D'D3": "D'D3",
        "D'D3 (FVIII binding)": "D'D3",
        'A1': 'A1',
        'A1 (Platelet binding)': 'A1',
        'A2': 'A2',
        'A2 (ADAMTS13 cleavage)': 'A2',
        'A3': 'A3',
        'A3 (Collagen binding)': 'A3',
        'D4': 'D4',
        'D4 (Multimerization)': 'D4',
        'CT': 'CT',
        'CT (C-terminal)': 'CT',
        'C1-C2': 'C1-C2',
        'C1-C2 (Integrin binding)': 'C1-C2',
        'D1-D2': 'D1-D2',
        'D1-D2 (Propeptide)': 'D1-D2',
    }
    return aliases.get(text, text)


def normalize_zygosity(value: Any) -> str | None:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    text = str(value).strip().lower()
    if not text:
        return None
    if text in {'het', 'heterozygous', '0/1', '1/0'}:
        return 'heterozygous'
    if text in {'hom', 'homozygous', '1/1'}:
        return 'homozygous'
    if 'compound' in text:
        return 'compound_heterozygous'
    return text


def parse_aa_change(aa_change: str | None) -> dict[str, Any]:
    if not aa_change:
        return {'raw': None, 'wt': None, 'mut': None, 'is_cysteine_disruption': False}

    text = str(aa_change).strip()
    match = re.search(r'([A-Z])(?:[a-z]{2})?(\d+)([A-Z*])', text)
    if not match:
        match = re.search(r'([A-Z])(\d+)([A-Z*])', text)

    if not match:
        return {'raw': text, 'wt': None, 'mut': None, 'is_cysteine_disruption': False}

    wt, _, mut = match.groups()
    return {
        'raw': text,
        'wt': wt,
        'mut': mut,
        'is_cysteine_disruption': wt == 'C' or mut == 'C',
        'introduces_stop': mut == '*',
    }


def extract_evidence(
    position: int,
    domain_name: str | None = None,
    aa_change: str | None = None,
    zygosity: str | None = None,
    splice_score: float | None = None,
    pae_delta: float | None = None,
) -> VariantEvidence:
    domain = normalize_domain(domain_name) or get_domain(position)
    notes = [f'变异位于 {domain} 功能域']

    if splice_score is not None and splice_score >= 0.2:
        notes.append('存在较高剪接/调控风险，需警惕并行的非经典 Type-2 解释路径')
    if pae_delta is not None and pae_delta >= 1.0:
        notes.append('局部结构扰动较明显，可作为蛋白功能异常的辅助证据')
    if normalize_zygosity(zygosity) in {'homozygous', 'compound_heterozygous'}:
        notes.append('双等位/复合杂合背景增强 2N 解释的可信度')

    return VariantEvidence(
        position=position,
        domain=domain,
        aa_change=aa_change,
        zygosity=normalize_zygosity(zygosity),
        splice_score=splice_score,
        pae_delta=pae_delta,
        notes=notes,
    )


def infer_mechanism_scores(evidence: VariantEvidence) -> dict[str, float]:
    scores = {mechanism: 0.0 for mechanism in MECHANISM_TO_SUBTYPE}
    priors = DOMAIN_PRIORS.get(evidence.domain, {})
    scores.update(priors)

    aa_info = parse_aa_change(evidence.aa_change)

    if aa_info.get('is_cysteine_disruption'):
        scores['multimer_assembly_defect'] += 0.20
        scores['secretion_defect'] += 0.15
        evidence.notes.append('涉及 cysteine/disulfide 改变，支持多聚体装配/分泌异常')

    if aa_info.get('introduces_stop'):
        scores['secretion_defect'] += 0.25
        scores['splice_regulatory_risk'] += 0.10
        evidence.notes.append('存在终止密码子改变，支持更广泛的功能破坏')

    if evidence.splice_score is not None:
        if evidence.splice_score >= 0.5:
            scores['splice_regulatory_risk'] += 0.65
        elif evidence.splice_score >= 0.2:
            scores['splice_regulatory_risk'] += 0.35

    if evidence.pae_delta is not None:
        if evidence.pae_delta >= 1.5:
            scores['multimer_assembly_defect'] += 0.15
            scores['adamts13_susceptibility'] += 0.10
        elif evidence.pae_delta >= 0.8:
            scores['multimer_assembly_defect'] += 0.08

    if evidence.domain == 'A1':
        if evidence.zygosity in {'heterozygous', None}:
            scores['gp1b_gain_of_function'] += 0.08
        if evidence.zygosity in {'homozygous', 'compound_heterozygous'}:
            scores['fviii_binding_loss'] += 0.10
    elif evidence.domain == "D'D3" and evidence.zygosity in {'homozygous', 'compound_heterozygous'}:
        scores['fviii_binding_loss'] += 0.20
    elif evidence.domain == 'A3':
        scores['collagen_binding_loss'] += 0.08
    elif evidence.domain == 'A2':
        scores['adamts13_susceptibility'] += 0.10

    return {k: round(clamp(v), 3) for k, v in scores.items()}


def infer_subtype_probabilities(mechanism_scores: dict[str, float]) -> dict[str, float]:
    subtype_scores = {subtype: 0.0 for subtype in SUBTYPES}

    for mechanism, score in mechanism_scores.items():
        for subtype, weight in MECHANISM_TO_SUBTYPE[mechanism].items():
            subtype_scores[subtype] += score * weight

    if subtype_scores['uncertain'] < 0.15:
        subtype_scores['uncertain'] = 0.15

    total = sum(subtype_scores.values()) or 1.0
    probabilities = {k: round(v / total, 3) for k, v in subtype_scores.items()}
    return dict(sorted(probabilities.items(), key=lambda item: item[1], reverse=True))


def mechanism_interpretation(mechanism_scores: dict[str, float]) -> list[dict[str, Any]]:
    ranked = sorted(mechanism_scores.items(), key=lambda item: item[1], reverse=True)
    top = []
    mechanism_labels = {
        'multimer_assembly_defect': '多聚体装配异常',
        'adamts13_susceptibility': 'ADAMTS13 易裂解',
        'gp1b_gain_of_function': 'GPIb gain-of-function',
        'gp1b_loss_of_function': 'GPIb loss-of-function',
        'collagen_binding_loss': '胶原结合减弱',
        'fviii_binding_loss': 'FVIII 结合减弱',
        'secretion_defect': '分泌/胞内转运异常',
        'splice_regulatory_risk': '剪接/调控异常风险',
    }
    for name, score in ranked[:3]:
        if score <= 0:
            continue
        top.append({'mechanism': name, 'label': mechanism_labels[name], 'score': round(score, 3)})
    return top


def subtype_explanation(predicted_subtype: str, evidence: VariantEvidence, top_mechanisms: list[dict[str, Any]]) -> str:
    domain_text = f'该变异位于 {evidence.domain} 域'
    mechanism_text = '；'.join([f"{m['label']}({m['score']:.2f})" for m in top_mechanisms]) or '机制信号不足'

    subtype_templates = {
        'type2A': f'{domain_text}，当前更支持多聚体缺失/易裂解相关的 2A-like 机制：{mechanism_text}。',
        'type2B': f'{domain_text}，当前更支持 GPIb 亲和力异常升高的 2B-like 机制：{mechanism_text}。',
        'type2M': f'{domain_text}，当前更支持功能下降但非典型多聚体缺失的 2M-like 机制：{mechanism_text}。',
        'type2N': f'{domain_text}，当前更支持 FVIII 结合受损的 2N-like 机制：{mechanism_text}。',
        'uncertain': f'{domain_text}，目前机制信号混合或偏弱：{mechanism_text}，应先补充分子/功能实验。',
    }
    return subtype_templates[predicted_subtype]


def predict_single(
    position: int,
    domain_name: str | None = None,
    aa_change: str | None = None,
    zygosity: str | None = None,
    splice_score: float | None = None,
    pae_delta: float | None = None,
) -> dict[str, Any]:
    evidence = extract_evidence(
        position=position,
        domain_name=domain_name,
        aa_change=aa_change,
        zygosity=zygosity,
        splice_score=splice_score,
        pae_delta=pae_delta,
    )
    mechanism_scores = infer_mechanism_scores(evidence)
    subtype_probabilities = infer_subtype_probabilities(mechanism_scores)
    predicted_subtype = next(iter(subtype_probabilities))
    top_mechanisms = mechanism_interpretation(mechanism_scores)

    return {
        'position': position,
        'domain': evidence.domain,
        'aa_change': aa_change,
        'zygosity': evidence.zygosity,
        'splice_score': splice_score,
        'pae_delta': pae_delta,
        'predicted_subtype': predicted_subtype,
        'confidence': subtype_probabilities[predicted_subtype],
        'subtype_probabilities': subtype_probabilities,
        'top_mechanisms': top_mechanisms,
        'explanation': subtype_explanation(predicted_subtype, evidence, top_mechanisms),
        'recommended_experiments': CLINICAL_RECOMMENDATIONS[predicted_subtype],
        'notes': evidence.notes,
        'method': 'mechanism_first_rule_based_triage',
    }


def row_to_prediction(row: pd.Series) -> dict[str, Any]:
    def pick(*names: str):
        for name in names:
            if name in row and pd.notna(row[name]):
                return row[name]
        return None

    position = pick('Position', 'position', 'AA_Position', 'Residue_Position')
    if position is None:
        raise ValueError("输入文件缺少 Position/position/AA_Position/Residue_Position 列")

    prediction = predict_single(
        position=int(position),
        domain_name=pick('VWF_Domain', 'domain', 'Domain'),
        aa_change=pick('AA_Change', 'aa_change', 'Protein_Change', 'HGVSp'),
        zygosity=pick('Zygosity', 'zygosity', 'Genotype'),
        splice_score=pick('Splice_Score', 'splice_score', 'AlphaGenome_Splice_Score'),
        pae_delta=pick('PAE_Delta_Local', 'pae_delta', 'PAE_Delta'),
    )
    return prediction


def batch_predict(input_file: str, output_file: str = 'predictions.csv') -> pd.DataFrame:
    df = pd.read_csv(input_file)
    predictions = [row_to_prediction(row) for _, row in df.iterrows()]
    pred_df = pd.DataFrame(predictions)

    for subtype in SUBTYPES:
        pred_df[f'prob_{subtype}'] = pred_df['subtype_probabilities'].apply(lambda d: d.get(subtype, 0.0))
    pred_df['top_mechanism_1'] = pred_df['top_mechanisms'].apply(lambda x: x[0]['label'] if x else None)
    pred_df['top_mechanism_2'] = pred_df['top_mechanisms'].apply(lambda x: x[1]['label'] if len(x) > 1 else None)
    pred_df['recommended_experiments_text'] = pred_df['recommended_experiments'].apply(lambda x: ' | '.join(x))
    pred_df['notes_text'] = pred_df['notes'].apply(lambda x: ' | '.join(x))

    output_df = pd.concat([df.reset_index(drop=True), pred_df], axis=1)
    output_df.to_csv(output_file, index=False)
    print(f'预测完成，结果保存到: {output_file}')
    return output_df


def print_single_result(result: dict[str, Any]):
    print('\n' + '=' * 68)
    print('VWF Type-2 可解释分诊结果')
    print('=' * 68)
    print(f"变异位置: {result['position']}")
    print(f"功能域: {result['domain']}")
    if result['aa_change']:
        print(f"氨基酸变化: {result['aa_change']}")
    if result['zygosity']:
        print(f"合子状态: {result['zygosity']}")
    print(f"预测分型: {result['predicted_subtype']}")
    print(f"置信度: {result['confidence']:.1%}")
    print(f"方法: {result['method']}")

    print('\n亚型概率:')
    for subtype, prob in result['subtype_probabilities'].items():
        print(f'  - {subtype}: {prob:.1%}')

    print('\n主要机制:')
    for item in result['top_mechanisms']:
        print(f"  - {item['label']}: {item['score']:.2f}")

    print(f"\n解释: {result['explanation']}")

    print('\n推荐实验:')
    for rec in result['recommended_experiments']:
        print(f'  - {rec}')

    if result['notes']:
        print('\n补充说明:')
        for note in result['notes']:
            print(f'  - {note}')
    print('=' * 68)


def main():
    parser = argparse.ArgumentParser(description='VWF Type-2 可解释分诊工具')
    parser.add_argument('--position', type=int, help='氨基酸位置 (1-2813)')
    parser.add_argument('--domain', type=str, help='VWF功能域 (可选)')
    parser.add_argument('--aa-change', type=str, help='氨基酸变化，如 p.Arg1308Cys')
    parser.add_argument('--zygosity', type=str, help='合子状态，如 heterozygous / homozygous / compound_heterozygous')
    parser.add_argument('--splice-score', type=float, help='剪接/调控风险分数')
    parser.add_argument('--pae-delta', type=float, help='局部 PAE 扰动值')
    parser.add_argument('--batch-predict', type=str, help='批量预测 CSV 文件路径')
    parser.add_argument('--output', type=str, default='predictions.csv', help='输出 CSV 文件路径')
    parser.add_argument('--json', action='store_true', help='单变异结果以 JSON 打印')
    args = parser.parse_args()

    if args.batch_predict:
        batch_predict(args.batch_predict, args.output)
        return

    if args.position:
        result = predict_single(
            position=args.position,
            domain_name=args.domain,
            aa_change=args.aa_change,
            zygosity=args.zygosity,
            splice_score=args.splice_score,
            pae_delta=args.pae_delta,
        )
        if args.json:
            print(json.dumps(result, ensure_ascii=False, indent=2))
        else:
            print_single_result(result)
        return

    parser.print_help()


if __name__ == '__main__':
    main()
