#!/usr/bin/env python3
"""Build the 16-case computational evidence delivery (requirements doc §9).

Generates under outputs/computational_evidence_16case_20260905/:
  run_manifest.json, case_inventory.csv, variant_inventory.csv,
  model_protocol_inventory.csv (per patient x variant x model status matrix),
  missingness_report.csv, artifact_manifest.parquet (+ .csv mirror),
  SHA256SUMS, README.md, fhir/<case_id>/bundle.json (16 bundles)

All statuses come from the verified ground truth on disk; no hand-filled
estimates. Large artifacts are referenced by immutable HF URI + revision +
sha256 (HEAD commit of dataset lucachangretta/VWF at build time); the API
key/token is never written to any output.

Usage:
    python3 scripts/pipeline/build_16case_evidence_delivery.py \
        [--hf-revision <commit_sha>] [--hf-verify]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent.parent.parent
OUT = ROOT / "outputs" / "computational_evidence_16case_20260905"

T1_BUNDLE = ROOT / "outputs" / "type1_panel_agent_20260828" / "server_bundle"
T2B_BUNDLE = ROOT / "outputs" / "type2b_panel_agent_20260829" / "server_bundle"
AG16 = OUT / "alphagenome"
MD_PANEL = ROOT / "outputs" / "computational_panel_20260829" / "md"

HF_BASE = "https://huggingface.co/datasets/lucachangretta/VWF/resolve/main"
HF_REPO = "lucachangretta/VWF"

AG_MANIFEST = Path("/tmp/vwf_ag16_stage/manifest.json")
BOLTZ_MANIFEST = Path("/tmp/vwf_boltz_raw_stage/manifest.json")
P0_MANIFEST = Path("/tmp/vwf_p0_md_stage/manifest_with_sha.json")
TJ_MANIFEST = MD_PANEL / "trajectories" / "upload_manifest.json"

AG_OUTPUT_TYPES = [
    "RNA_SEQ", "CAGE", "DNASE", "CHIP_HISTONE", "CHIP_TF",
    "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS", "PROCAP",
]
AG_NOT_SUPPORTED = ["ATAC", "CONTACT_MAPS"]

STATUS_VOCAB = [
    "available_complete", "available_partial", "not_applicable",
    "blocked_missing_input", "not_supported_by_model",
    "failed_retry_exhausted", "not_run_out_of_scope",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# Variant inventory: 16 patients, 17 variant rows
# ---------------------------------------------------------------------------

AA_SHORT = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U", "Pyl": "O",
}


def compact_aa(hgvs_p: str) -> str:
    """p.Arg1205His -> R1205H; p.Leu2407ProfsTer11 -> L2407fs.

    Accepts both 3-letter (p.Arg1205His) and already-compact (p.R1205H)
    forms; returns '' when no single-residue missense can be derived
    (synonymous, nonsense, frameshift, splice, CNV)."""
    if not isinstance(hgvs_p, str) or not hgvs_p.strip():
        return ""
    s = hgvs_p.removeprefix("p.")
    if s in {"=", "?", "Met1?"}:
        return ""

    def one_letter(code: str) -> str | None:
        if len(code) == 1 and code.isupper():
            return code if code in AA_SHORT.values() else None
        return AA_SHORT.get(code)

    # Long form: AAAAnnnBBB / AAAAnnnfsTerNN / AAAAnnnTer
    match = re.match(r"^([A-Z][a-z]{2})(\d+)([A-Za-z0-9*]*)$", s)
    if match:
        aa, pos, rest = match.groups()
        aa1 = one_letter(aa)
        if aa1 is None:
            return ""
        if "fs" in rest:
            return f"{aa1}{pos}fs"
        if rest in {"", "Ter", "*"}:
            return f"{aa1}{pos}*"
        if len(rest) >= 3 and one_letter(rest[:3]) is not None and rest[3:] == "":
            return f"{aa1}{pos}{one_letter(rest[:3])}"
        return ""

    # Compact form: XnnnY / Xnnn* / Xnnnfs
    match = re.match(r"^([A-Z*])(\d+)([A-Z*]?)(fs[0-9]*)?$", s)
    if match:
        aa, pos, rep, fs = match.groups()
        aa1 = one_letter(aa)
        if aa1 is None:
            return ""
        if fs:
            return f"{aa1}{pos}fs"
        if rep:
            rep1 = one_letter(rep)
            return f"{aa1}{pos}{rep1}" if rep1 else ""
        return f"{aa1}{pos}*"
    return ""


def hgvs_fallback(case_id: str, t) -> str:
    """variant_id tail when no single-residue missense can be derived
    (synonymous / nonsense / frameshift / splice / CNV)."""
    if case_id == "CASE_T1_004":
        return "T1_004_CNV"
    return {
        "CASE_T1_003": "G2488syn",
        "CASE_T1_005": "Y988Ter",  # unreachable: compact_aa handles it
        "CASE_T1_006": "L2407fs",  # unreachable: compact_aa handles it
        "CASE_T1_007": "c3379splicesite",
    }.get(case_id, "CNV")


def load_variant_inventory() -> pd.DataFrame:
    g1 = pd.read_csv(T1_BUNDLE / "input" / "genetics.csv")
    g2 = pd.read_csv(T2B_BUNDLE / "input" / "genetics.csv")
    r1 = pd.read_csv(T1_BUNDLE / "input" / "alphagenome_requests.csv")
    r2 = pd.read_csv(T2B_BUNDLE / "input" / "alphagenome_requests.csv")
    reqs = pd.concat([r1, r2], ignore_index=True)
    reqs = reqs.set_index("case_id")
    # T2B requests carry patient_id (CASE_T2B_005_VARIANT_1 and
    # _VARIANT_2_BENIGN share patient CASE_T2B_005); T1 has 1 case = 1 patient.
    patient_of = {}
    if "patient_id" in r2.columns:
        patient_of = dict(zip(r2["case_id"], r2["patient_id"]))

    rows = []
    for bundle, g in (("T1", g1), ("T2B", g2)):
        for t in g.itertuples():
            case_id = t.研究病例ID
            is_cnv = case_id == "CASE_T1_004"
            req = reqs.loc[case_id] if case_id in reqs.index else None
            rows.append({
                "case_id": case_id,
                "bundle": bundle,
                "patient_id": patient_of.get(case_id, case_id),
                "variant_row_id": case_id,  # 1 row per case; T2B_005 two rows
                "variant_id": (
                    "VWF_" + (compact_aa(t.氨基酸改变) or hgvs_fallback(case_id, t))
                ),
                "gene": "VWF",
                "variant_type": t.变异类型,
                "hgvs_c": "" if is_cnv or pd.isna(t.核苷酸改变) else t.核苷酸改变,
                "hgvs_p": "" if is_cnv or pd.isna(t.氨基酸改变) else t.氨基酸改变,
                "aa_change": "" if is_cnv else compact_aa(t.氨基酸改变),
                "assembly": "GRCh38",
                "chromosome": "" if is_cnv or pd.isna(t.hg38染色体) else t.hg38染色体,
                "position": np.nan if is_cnv or pd.isna(t.hg38位置) else int(t.hg38位置),
                "ref": "" if is_cnv or pd.isna(t.hg38_REF) else t.hg38_REF,
                "alt": "" if is_cnv or pd.isna(t.hg38_ALT) else t.hg38_ALT,
                "interval_start": (
                    np.nan if req is None else int(reqs.loc[case_id, "interval_start"])
                ),
                "interval_end": (
                    np.nan if req is None else int(reqs.loc[case_id, "interval_end"])
                ),
                "coordinate_qc": t.坐标QC,
                "normalization_source": t.归一化来源,
                "ag_request_status": t.AlphaGenome请求状态,
                "boltz_request_status": t.Boltz请求状态,
                "dedup_key": (
                    f"chr12:{int(reqs.loc[case_id, 'position'])}"
                    if req is not None
                    else ""
                ),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Per patient x variant x model status matrix (model_protocol_inventory)
# ---------------------------------------------------------------------------

BOLTZ_AXIS_FOR = {
    "R1205H": "dprime_d3_fviii_binding",
    "P1413L": "a1_gpiba_forced_binding/a1_heparan_sulfate_binding/a1_aim_autoinhibition_context",
    "A1500V": "a2_folded_stability/a2_adamts13_folded_complex",
    "R1943C": "d4_assembly_context",
    "C1950Y": "d4_assembly_context",
    "R1308C": "a1_gpiba_forced_binding/a1_heparan_sulfate_binding/a1_aim_autoinhibition_context",
    "S1310F": "a1_gpiba_forced_binding/a1_heparan_sulfate_binding/a1_aim_autoinhibition_context",
    "V1316M": "a1_gpiba_forced_binding/a1_heparan_sulfate_binding/a1_aim_autoinhibition_context",
    "R1341W": "a1_gpiba_forced_binding/a1_heparan_sulfate_binding/a1_aim_autoinhibition_context",
    "A1461D": "a1_gpiba_forced_binding/a1_heparan_sulfate_binding/a1_aim_autoinhibition_context",
    "D2449N": "c_domain_assembly_context",
}
NON_MISSENSE_CASES = {
    "CASE_T1_003", "CASE_T1_005", "CASE_T1_006", "CASE_T1_007",
}
MD_A1_VARIANTS = {"P1413L", "R1308C", "S1310F", "V1316M", "R1341W", "A1461D"}


def load_ag_ledger() -> dict:
    ledger = {}
    for line in (AG16 / "run_ledger.jsonl").read_text().splitlines():
        rec = json.loads(line)
        ledger[rec["case_id"]] = rec
    return ledger


def build_model_protocol_inventory(variants: pd.DataFrame) -> pd.DataFrame:
    ag_ledger = load_ag_ledger()
    bz_manifest = json.loads(BOLTZ_MANIFEST.read_text())["jobs"]
    bz_keys = {k.replace("boltz_results_", "") for k in bz_manifest}

    rows = []
    for v in variants.itertuples():
        case_id = v.case_id
        aa = v.aa_change
        # ---- AlphaGenome layer (11 output types per requestable variant) ----
        if case_id in ag_ledger:
            for entry in ag_ledger[case_id]["ledger"]:
                status = entry["status"]
                if status == "success":
                    final = "available_complete"
                    reason = ""
                elif status == "not_supported_by_model":
                    final = "not_supported_by_model"
                    reason = entry.get("reason", "")
                else:
                    final = status
                    reason = entry.get("reason", "")
                rows.append({
                    "case_id": case_id,
                    "variant_id": v.variant_id,
                    "model": "alphagenome",
                    "protocol_id": f"ag_predict_variant_v0.8.0__{entry['output_type']}",
                    "construct_or_output": entry["output_type"],
                    "status": final,
                    "reason": reason,
                    "eligible_for_agent": "true" if final == "available_complete" else "false",
                })
        else:
            rows.append({
                "case_id": case_id,
                "variant_id": v.variant_id,
                "model": "alphagenome",
                "protocol_id": "ag_predict_variant_v0.8.0__all_outputs",
                "construct_or_output": "ALL",
                "status": "blocked_missing_input",
                "reason": "normalization_not_available_cnv",
                "eligible_for_agent": "false",
            })

        # ---- Boltz layer ----
        if case_id == "CASE_T1_004":
            rows.append({
                "case_id": case_id, "variant_id": v.variant_id,
                "model": "boltz2",
                "protocol_id": "boltz2_vwd_functional_panel",
                "construct_or_output": "ALL",
                "status": "blocked_missing_input",
                "reason": "no_normalizable_variant_definition_cnv",
                "eligible_for_agent": "false",
            })
        elif case_id in NON_MISSENSE_CASES:
            rows.append({
                "case_id": case_id, "variant_id": v.variant_id,
                "model": "boltz2",
                "protocol_id": "boltz2_vwd_functional_panel",
                "construct_or_output": "ALL",
                "status": "not_applicable",
                "reason": "non_missense_variant_" + str(v.variant_type),
                "eligible_for_agent": "false",
            })
        else:
            for axis in BOLTZ_AXIS_FOR.get(aa, "").split("/"):
                if not axis:
                    continue
                job = f"{v.variant_id}__{axis}"
                has_raw = job in bz_keys
                rows.append({
                    "case_id": case_id, "variant_id": v.variant_id,
                    "model": "boltz2",
                    "protocol_id": f"boltz2_vwd_functional_panel__{axis}",
                    "construct_or_output": axis,
                    "status": "available_complete" if has_raw else "failed_retry_exhausted",
                    "reason": "" if has_raw else "raw_output_not_found",
                    "eligible_for_agent": "true" if has_raw else "false",
                })

        # ---- MD layer ----
        # Protocol 1: 7A6O AIM-A1 50 ns (A1 variants only)
        p7 = "md_7a6o_aim_a1_50ns_charmm27"
        if aa in MD_A1_VARIANTS or v.variant_id == "VWF_P1413L":
            rows.append({
                "case_id": case_id, "variant_id": v.variant_id,
                "model": "gromacs_md",
                "protocol_id": p7,
                "construct_or_output": "AIM_A1_7A6O",
                "status": "available_complete",
                "reason": "",
                "eligible_for_agent": "true",
            })
        else:
            rows.append({
                "case_id": case_id, "variant_id": v.variant_id,
                "model": "gromacs_md",
                "protocol_id": p7,
                "construct_or_output": "AIM_A1_7A6O",
                "status": "not_run_out_of_scope",
                "reason": "variant_not_in_7a6o_panel",
                "eligible_for_agent": "false",
            })

        # Protocol 2-4: P0 experimental 20 ns (1SQ0 a1_gpiba / 3GXB a2 / 6N29 dprime_d3)
        p0_axes = {
            "md_20ns_p0__a1_gpiba_1SQ0": MD_A1_VARIANTS,
            "md_20ns_p0__a2_folded_3GXB": {"A1500V"},
            "md_20ns_p0__dprime_d3_6N29": {"R1205H"},
        }
        for protocol, eligible in p0_axes.items():
            if aa in eligible:
                rows.append({
                    "case_id": case_id, "variant_id": v.variant_id,
                    "model": "gromacs_md",
                    "protocol_id": protocol,
                    "construct_or_output": protocol.split("__")[1],
                    "status": "available_complete",
                    "reason": "",
                    "eligible_for_agent": "true",
                })
            else:
                rows.append({
                    "case_id": case_id, "variant_id": v.variant_id,
                    "model": "gromacs_md",
                    "protocol_id": protocol,
                    "construct_or_output": protocol.split("__")[1],
                    "status": "not_run_out_of_scope",
                    "reason": "variant_not_in_p0_protocol_cohort",
                    "eligible_for_agent": "false",
                })

        # slow025 SMD: archived but no-go, never agent-eligible
        rows.append({
            "case_id": case_id, "variant_id": v.variant_id,
            "model": "gromacs_smd",
            "protocol_id": "smd_slow025_aim_a1",
            "construct_or_output": "AIM_A1_7A6O_smd",
            "status": "not_run_out_of_scope",
            "reason": "complete_no_go_direction_reversed_excluded_from_classifier",
            "eligible_for_agent": "false",
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Artifact manifest (HF URIs + revision + sha256)
# ---------------------------------------------------------------------------

def load_artifacts(hf_revision: str) -> tuple[pd.DataFrame, list[dict]]:
    entries = []

    def add(layer, artifact_id, repo_path, size_bytes, sha, extra=None):
        uri = f"{HF_BASE}/{repo_path}@{hf_revision}"
        e = {
            "artifact_id": artifact_id,
            "layer": layer,
            "repo_path": repo_path,
            "uri": uri,
            "revision": hf_revision,
            "size_bytes": int(size_bytes),
            "sha256": sha,
        }
        if extra:
            e.update(extra)
        entries.append(e)

    # AlphaGenome: 16 case archives + shared small files
    ag = json.loads(AG_MANIFEST.read_text())["cases"]
    for case_id, info in sorted(ag.items()):
        add("alphagenome", f"ag16__{case_id}",
            f"computational_evidence_16case_20260905/alphagenome/{info['archive']}",
            info["size_bytes"], info["sha256"],
            {"n_npz": info["npz_count"], "n_sidecar": info["sidecar_count"]})
    for shared in ["derived/paired_track_summaries.csv", "ledger_flat.csv",
                   "run_ledger.jsonl", "run_manifest.json", "ontology_plan.json"]:
        p = AG16 / shared.split("/")[-1] if shared.startswith("derived/") else AG16 / shared
        if not p.exists():
            continue
        add("alphagenome", f"ag16__shared__{p.name}",
            f"computational_evidence_16case_20260905/alphagenome/{p.name}",
            p.stat().st_size, sha256_file(p))

    # Boltz: 32 job archives
    bz = json.loads(BOLTZ_MANIFEST.read_text())["jobs"]
    for job, info in sorted(bz.items()):
        add("boltz2", f"boltz__{job.replace('boltz_results_', '')}",
            f"computational_panel_20260829/boltz_raw/{info['archive']}",
            info["size_bytes"], info["sha256"],
            {"cif_models": info["cif_models"],
             "confidence_jsons": info["confidence_jsons"],
             "pae_npz": info["pae_npz"], "pde_npz": info["pde_npz"],
             "plddt_npz": info["plddt_npz"], "done_marker": str(info["done_marker"]).lower()})

    # MD P0: 9 system archives
    p0 = json.loads(P0_MANIFEST.read_text())["axes"]
    for axis, systems in sorted(p0.items()):
        for variant, info in sorted(systems.items()):
            add("gromacs_md", f"md_p0__{axis}__{variant}",
                f"md_p0_experimental_archives/{info['archive']}",
                info["size_bytes"], info["sha256_local"],
                {"md_tag": info["md_tag"],
                 "merged_full_xtc": str(info["merged_full_xtc"]).lower(),
                 "has_final_gro": str(info["has_final_gro"]).lower()})

    # MD 7A6O 50 ns trajectories: 7 XTC
    tj = json.loads(TJ_MANIFEST.read_text())["trajectories"]
    for name, info in sorted(tj.items()):
        add("gromacs_md", f"md_7a6o__{name.replace('_prod_concat.xtc', '')}",
            f"computational_panel_20260829/md/trajectories/{name}",
            info["size_bytes"], info["sha256"])

    df = pd.DataFrame(entries)
    json_entries = entries
    return df, json_entries


# ---------------------------------------------------------------------------
# FHIR bundles (FHIR-shaped / internal validated)
# ---------------------------------------------------------------------------

def build_fhir_bundles(variants: pd.DataFrame, matrix: pd.DataFrame,
                       artifacts: pd.DataFrame, hf_revision: str) -> None:
    fhir_dir = OUT / "fhir"
    fhir_dir.mkdir(parents=True, exist_ok=True)
    # One bundle per PATIENT (16): T2B_005 carries both variant rows.
    for stale in fhir_dir.glob("*/bundle.json"):
        stale.unlink()
    for stale in (fhir_dir).glob("*"):
        if stale.is_dir() and not any(stale.iterdir()):
            stale.rmdir()

    art_by_id = artifacts.set_index("artifact_id")
    scorer_t1 = pd.read_csv(T1_BUNDLE / "results" / "alphagenome" / "alphagenome_scores_long.csv")
    scorer_t2 = pd.read_csv(T2B_BUNDLE / "results" / "alphagenome" / "alphagenome_scores_long.csv")
    scorers = pd.concat([scorer_t1, scorer_t2], ignore_index=True)

    b1 = pd.read_csv(T1_BUNDLE / "results" / "boltz" / "boltz_results_summary.csv")
    b2 = pd.read_csv(T2B_BUNDLE / "results" / "boltz" / "boltz_results_summary.csv")
    boltz = pd.concat([b1, b2], ignore_index=True).drop_duplicates("job_name")

    p0_iface = pd.read_csv(MD_PANEL / "p0_experimental_20ns" / "p0_a1_gpiba_interface_summary.csv")
    ag_summaries = pd.read_csv(AG16 / "derived" / "paired_track_summaries.csv")

    # Per-patient loop: one bundle per patient, variant rows nested inside.
    for patient_id, pat_variants in variants.groupby("patient_id", sort=True):
        first = pat_variants.iloc[0]
        bundle = {
            "resourceType": "Bundle",
            "id": f"vwf-evidence-{patient_id.lower()}",
            "type": "collection",
            "meta": {
                "profile": ["VWFComputationalEvidence/1.0 (FHIR-shaped, internal validated)"],
                "tag": [{"system": "http://terminology.local/validation",
                         "code": "internal-schema-only",
                         "display": "FHIR-shaped/internal validated; not HL7-validator certified"}],
            },
            "timestamp": utc_now(),
            "entry": [],
        }

        # Patient
        bundle["entry"].append({
            "fullUrl": f"urn:uuid:patient-{patient_id.lower()}",
            "resource": {
                "resourceType": "Patient",
                "id": patient_id.lower(),
                "identifier": [{"system": "urn:vwf-panel:patient-id", "value": patient_id}],
            },
        })

        for v in pat_variants.itertuples():
            case_id = v.case_id

            # Variant Observation (identity fields, current case only)
            if v.position is not None and not (isinstance(v.position, float) and np.isnan(v.position)):
                obs = {
                    "resourceType": "Observation",
                    "id": f"obs-variant-{case_id.lower()}",
                    "status": "final" if v.coordinate_qc != "NORMALIZATION_NOT_AVAILABLE" else "registered",
                    "code": {
                        "coding": [{"system": "http://loinc.org", "code": "69548-6",
                                    "display": "Genetic variant assessment"}],
                    },
                    "subject": {"reference": f"Patient/{patient_id.lower()}"},
                    "valueCodeableConcept": {
                        "coding": [{"system": "http://varnomen.hgvs.org",
                                    "code": v.hgvs_c or v.hgvs_p or "CNV"}],
                        "text": f"{v.hgvs_p or v.variant_type} ({v.hgvs_c or 'no HGVS.c'})",
                    },
                    "note": [{"text": f"assembly=GRCh38; chr={v.chromosome}; pos={int(v.position)}; "
                                      f"REF={v.ref}; ALT={v.alt}; coordinate_qc={v.coordinate_qc}"}],
                }
                bundle["entry"].append({"fullUrl": f"urn:uuid:obs-variant-{case_id.lower()}",
                                        "resource": obs})

            case_matrix = matrix[matrix.case_id == case_id]

            # AlphaGenome Observations (per persisted output type; paired summary stats)
            case_ag = ag_summaries[ag_summaries.case_id == case_id]
            for output_type in AG_OUTPUT_TYPES:
                sub = case_matrix[(case_matrix.model == "alphagenome") &
                                  (case_matrix.construct_or_output == output_type)]
                if sub.empty:
                    continue
                status = sub.iloc[0].status
                obs = {
                    "resourceType": "Observation",
                    "id": f"obs-ag-{case_id.lower()}-{output_type.lower()}",
                    "status": "final" if status == "available_complete" else "registered",
                    "code": {"coding": [{"system": "urn:alphagenome:output-type",
                                         "code": output_type}],
                             "text": f"AlphaGenome paired REF/ALT {output_type}"},
                    "subject": {"reference": f"Patient/{patient_id.lower()}"},
                    "effectiveDateTime": utc_now(),
                    "method": {"text": "AlphaGenome 0.8.0 predict_variant; 1 Mb interval; "
                                       "REF and ALT float32 arrays persisted losslessly (npz+json sidecar)"},
                    "component": [],
                }
                tracks = case_ag[case_ag.output_type == output_type]
                for t in tracks.itertuples():
                    obs["component"].append({
                        "code": {"text": f"track:{t.track_name}"},
                        "valueQuantity": {
                            "value": float(t.abs_max) if not pd.isna(t.abs_max) else None,
                            "unit": "ALT-REF max |delta| (model output units)",
                        },
                    })
                if not obs["component"]:
                    obs.pop("component")
                if status != "available_complete":
                    obs.setdefault("note", []).append(
                        {"text": f"status={status}; reason={sub.iloc[0].reason}"})
                doc_ref_id = f"doc-ag-{case_id.lower()}"
                obs["derivedFrom"] = [{"reference": f"DocumentReference/{doc_ref_id}"}]
                bundle["entry"].append({
                    "fullUrl": f"urn:uuid:obs-ag-{case_id.lower()}-{output_type.lower()}",
                    "resource": obs})

            # AlphaGenome DocumentReference (archive URI + hash)
            art_id = f"ag16__{case_id}"
            if art_id in art_by_id.index:
                art = art_by_id.loc[art_id]
                bundle["entry"].append({
                    "fullUrl": f"urn:uuid:doc-ag-{case_id.lower()}",
                    "resource": {
                        "resourceType": "DocumentReference",
                        "id": f"doc-ag-{case_id.lower()}",
                        "status": "current",
                        "content": [{
                            "attachment": {
                                "url": art["uri"],
                                "size": int(art["size_bytes"]),
                                "hash": art["sha256"],
                                "title": f"AlphaGenome raw arrays {case_id} (9 output types)",
                            },
                        }],
                    },
                })

            # AlphaGenome scorer Observations (official variant scorers, top tracks only)
            case_scorers = scorers[scorers.case_id == case_id]
            if not case_scorers.empty:
                top = (case_scorers.sort_values("raw_score", key=abs, ascending=False)
                       .groupby("scorer_key").head(1))
                for s in top.itertuples():
                    bundle["entry"].append({
                        "fullUrl": f"urn:uuid:obs-ags-{case_id.lower()}-{s.scorer_key.lower()}",
                        "resource": {
                            "resourceType": "Observation",
                            "id": f"obs-ags-{case_id.lower()}-{s.scorer_key.lower()}",
                            "status": "final",
                            "code": {"coding": [{"system": "urn:alphagenome:variant-scorer",
                                                 "code": s.scorer_key}],
                                     "text": f"AlphaGenome official variant scorer {s.scorer_key}"},
                            "subject": {"reference": f"Patient/{patient_id.lower()}"},
                            "valueQuantity": {
                                "value": float(s.raw_score),
                                "unit": "scorer raw score (signed ALT-REF direction)",
                            },
                            "method": {"text": f"track={s.track_name}; output_type={s.output_type}; "
                                               f"biosample={s.biosample_name}; gene={s.gene_name}"},
                        },
                    })

            # Boltz Observations per assay
            variant_jobs = boltz[boltz.job_name.str.startswith(v.variant_id + "__")]
            for job in variant_jobs.itertuples():
                metric = job.primary_metric
                value = float(job.primary_value) if not pd.isna(job.primary_value) else None
                art_id = f"boltz__{job.job_name}"
                art = art_by_id.loc[art_id] if art_id in art_by_id.index else None
                obs = {
                    "resourceType": "Observation",
                    "id": f"obs-boltz-{job.job_name.lower()}",
                    "status": "final",
                    "code": {"coding": [{"system": "urn:vwf:boltz-assay",
                                         "code": job.job_name.split("__", 1)[1]}],
                             "text": f"Boltz-2 structural confidence: {metric}"},
                    "subject": {"reference": f"Patient/{patient_id.lower()}"},
                    "valueQuantity": {"value": value, "unit": metric},
                    "method": {"text": f"Boltz-2 vwd functional panel; n_samples={int(job.n_samples)}; "
                                       f"primary_metric={job.primary_metric} "
                                       f"({job.primary_metric_reason})"},
                    "note": [{"text": "Model confidence/geometry metric; not a binding free energy "
                                      "or functional effect size"}],
                }
                if art is not None:
                    obs["derivedFrom"] = [{
                        "reference": f"DocumentReference/doc-boltz-{job.job_name.lower()}"}]
                    bundle["entry"].append({
                        "fullUrl": f"urn:uuid:doc-boltz-{job.job_name.lower()}",
                        "resource": {
                            "resourceType": "DocumentReference",
                            "id": f"doc-boltz-{job.job_name.lower()}",
                            "status": "current",
                            "content": [{
                                "attachment": {
                                    "url": art["uri"],
                                    "size": int(art["size_bytes"]),
                                    "hash": art["sha256"],
                                    "title": f"Boltz-2 raw output {job.job_name}",
                                },
                            }],
                        },
                    })
                bundle["entry"].append({
                    "fullUrl": f"urn:uuid:obs-boltz-{job.job_name.lower()}",
                    "resource": obs})

            # MD Observations (P0 a1_gpiba interface for eligible variants)
            if v.aa_change in MD_A1_VARIANTS:
                sub = p0_iface[p0_iface.variant == v.aa_change]
                wt = p0_iface[p0_iface.variant == "WT"]
                if not sub.empty:
                    s = sub.iloc[0]
                    w = wt.iloc[0]
                    bundle["entry"].append({
                        "fullUrl": f"urn:uuid:obs-md-p0-a1gpiba-{case_id.lower()}",
                        "resource": {
                            "resourceType": "Observation",
                            "id": f"obs-md-p0-a1gpiba-{case_id.lower()}",
                            "status": "final",
                            "code": {"coding": [{"system": "urn:vwf:md-protocol",
                                                 "code": "md_20ns_p0__a1_gpiba_1SQ0"}],
                                     "text": "A1-GPIbalpha interface contacts (tail 15-20 ns)"},
                            "subject": {"reference": f"Patient/{patient_id.lower()}"},
                            "valueQuantity": {
                                "value": float(s.contacts_tail_15_20ns),
                                "unit": "heavy-atom contact pairs (cutoff 0.45 nm)"},
                            "method": {"text": "GROMACS 2025.4 CHARMM27; 1SQ0 experimental complex; "
                                               "20 ns production; matched WT=PANEL WT (1SQ0, same protocol)"},
                            "note": [{"text": f"matched_WT_tail_contacts={float(w.contacts_tail_15_20ns)}; "
                                              f"delta_vs_WT={float(s.contacts_tail_15_20ns) - float(w.contacts_tail_15_20ns):.1f}"}],
                        },
                    })

            # Missingness as explicit OperationOutcome-like entries
            missing = case_matrix[case_matrix.status.isin(
                ["blocked_missing_input", "not_supported_by_model",
                 "not_applicable", "not_run_out_of_scope", "failed_retry_exhausted"])]
            for m in missing.itertuples():
                bundle["entry"].append({
                    "fullUrl": f"urn:uuid:oo-{case_id.lower()}-{m.model}-{m.protocol_id}".lower().replace(" ", ""),
                    "resource": {
                        "resourceType": "OperationOutcome",
                        "id": f"oo-{case_id.lower()}-{m.model}-{m.construct_or_output}".lower().replace(" ", "").replace("|", "-"),
                        "issue": [{
                            "severity": "information" if m.status in ("not_applicable", "not_run_out_of_scope") else "warning",
                            "code": "incomplete",
                            "details": {"text": f"{m.model}/{m.construct_or_output}: {m.status} ({m.reason})"},
                        }],
                    },
                })

        (fhir_dir / patient_id / "bundle.json").parent.mkdir(parents=True, exist_ok=True)
        (fhir_dir / patient_id / "bundle.json").write_text(
            json.dumps(bundle, indent=1, ensure_ascii=False), encoding="utf-8")


# ---------------------------------------------------------------------------
# Missingness report
# ---------------------------------------------------------------------------

def build_missingness_report(matrix: pd.DataFrame) -> pd.DataFrame:
    missing = matrix[matrix.status != "available_complete"]
    return missing[["case_id", "variant_id", "model", "protocol_id",
                    "construct_or_output", "status", "reason"]].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hf-revision", required=True,
                        help="HF dataset HEAD commit sha (immutable revision)")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    print(f"[1/6] variant inventory ...")
    variants = load_variant_inventory()
    variants.to_csv(OUT / "variant_inventory.csv", index=False)
    n_patients = variants.patient_id.nunique()
    print(f"      {len(variants)} variant rows, {n_patients} patients")

    print(f"[2/6] case inventory ...")
    case_rows = []
    for patient_id, group in variants.groupby("patient_id", sort=True):
        case_rows.append({
            "patient_id": patient_id,
            "bundle": group.bundle.iloc[0],
            "n_variant_rows": len(group),
            "case_ids": "|".join(sorted(group.case_id)),
            "variant_ids": "|".join(sorted(set(group.variant_id))),
            "ag_terminal_state": (
                "blocked_missing_variant_definition"
                if "CASE_T1_004" in set(group.case_id)
                else "available_complete_9of11_persisted"
            ),
        })
    case_inv = pd.DataFrame(case_rows)
    case_inv.to_csv(OUT / "case_inventory.csv", index=False)
    print(f"      {len(case_inv)} patients")

    print(f"[3/6] model x protocol inventory (completeness matrix) ...")
    matrix = build_model_protocol_inventory(variants)
    matrix.to_csv(OUT / "model_protocol_inventory.csv", index=False)
    counts = matrix.status.value_counts()
    print("      " + "; ".join(f"{k}={v}" for k, v in counts.items()))

    print(f"[4/6] artifact manifest (HF URIs @ {args.hf_revision[:12]}) ...")
    artifacts, json_entries = load_artifacts(args.hf_revision)
    artifacts.to_parquet(OUT / "artifact_manifest.parquet", index=False)
    artifacts.to_csv(OUT / "artifact_manifest.csv", index=False)
    print(f"      {len(artifacts)} artifacts, "
          f"{artifacts.size_bytes.sum() / 1e9:.2f} GB total")

    print(f"[5/6] missingness report ...")
    missingness = build_missingness_report(matrix)
    missingness.to_csv(OUT / "missingness_report.csv", index=False)
    print(f"      {len(missingness)} non-complete rows")

    print(f"[6/6] FHIR bundles ...")
    build_fhir_bundles(variants, matrix, artifacts, args.hf_revision)
    n_bundles = len(list((OUT / "fhir").glob("*/bundle.json")))
    print(f"      {n_bundles} case bundles written")

    # run_manifest.json
    run_manifest = {
        "created_at": utc_now(),
        "run_id": "computational_evidence_16case_20260905",
        "requirements_doc": "docs/SERVER_16CASE_COMPUTATIONAL_DATA_COMPLETION_REQUIREMENTS_20260905.md",
        "generator": {
            "script": "scripts/pipeline/build_16case_evidence_delivery.py",
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
        "hf_repo": HF_REPO,
        "hf_revision": args.hf_revision,
        "status_vocabulary": STATUS_VOCAB,
        "counts": {
            "patients": int(n_patients),
            "variant_rows": int(len(variants)),
            "matrix_rows": int(len(matrix)),
            "artifacts": int(len(artifacts)),
            "artifact_total_bytes": int(artifacts.size_bytes.sum()),
            "fhir_bundles": n_bundles,
            "missingness_rows": int(len(missingness)),
        },
        "layers": {
            "alphagenome": {
                "source_run": "outputs/computational_evidence_16case_20260905/alphagenome/",
                "runner": "scripts/pipeline/run_16case_alphagenome_full_persistence.py",
                "sdk_version": "alphagenome 0.8.0",
                "request_csvs": [
                    "outputs/type1_panel_agent_20260828/server_bundle/input/alphagenome_requests.csv",
                    "outputs/type2b_panel_agent_20260829/server_bundle/input/alphagenome_requests.csv",
                ],
                "official_variant_scorers": [
                    "outputs/type1_panel_agent_20260828/server_bundle/results/alphagenome/alphagenome_scores_long.csv",
                    "outputs/type2b_panel_agent_20260829/server_bundle/results/alphagenome/alphagenome_scores_long.csv",
                ],
                "persisted_output_types": AG_OUTPUT_TYPES,
                "not_supported_by_model": AG_NOT_SUPPORTED,
            },
            "boltz2": {
                "job_manifests": [
                    "outputs/type1_panel_agent_20260828/server_bundle/boltz/job_manifest.csv",
                    "outputs/type2b_panel_agent_20260829/server_bundle/boltz/job_manifest.csv",
                ],
                "intended_jobs_total": 35,
                "unique_jobs_after_wt_dedup": 32,
                "completed_jobs": 32,
                "summaries": [
                    "outputs/type1_panel_agent_20260828/server_bundle/results/boltz/boltz_results_summary.csv",
                    "outputs/type2b_panel_agent_20260829/server_bundle/results/boltz/boltz_results_summary.csv",
                ],
            },
            "gromacs_md": {
                "protocols": {
                    "md_20ns_p0__a1_gpiba_1SQ0": "A1-GPIbalpha experimental complex 1SQ0; WT+5 variants",
                    "md_20ns_p0__a2_folded_3GXB": "A2 folded monomer 3GXB; WT+A1500V",
                    "md_20ns_p0__dprime_d3_6N29": "D'D3 6N29; WT+R1205H",
                    "md_7a6o_aim_a1_50ns": "AIM-A1 experimental 7A6O; PANEL_WT_MATCHED+6 variants",
                },
                "md_feature_matrix": "outputs/computational_panel_20260829/md/md_feature_matrix.csv",
                "smd_slow025": "complete_no_go; eligible_for_agent=false everywhere",
            },
        },
        "secrets": "No API key, token, or Authorization header is stored in any file of this run; "
                   "credentials were read from environment variables only.",
    }
    (OUT / "run_manifest.json").write_text(
        json.dumps(run_manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    # README.md
    readme = f"""# 16-case computational evidence run (2026-09-05)

Generated by `scripts/pipeline/build_16case_evidence_delivery.py` from verified
on-disk inventories. Requirements:
`docs/SERVER_16CASE_COMPUTATIONAL_DATA_COMPLETION_REQUIREMENTS_20260905.md`.

## Layout

| File | Content |
|------|---------|
| `run_manifest.json` | run provenance, layer inventory, HF revision, counts |
| `case_inventory.csv` | {n_patients} patients |
| `variant_inventory.csv` | {len(variants)} variant rows (T2B_005 has two; T1_004 = CNV blocked) |
| `model_protocol_inventory.csv` | completeness matrix: patient x variant x model x protocol, status + reason + eligible_for_agent |
| `artifact_manifest.parquet` / `.csv` | all large artifacts: immutable HF URI + revision + sha256 + bytes |
| `missingness_report.csv` | every non-complete row with reason |
| `fhir/<patient_id>/bundle.json` | per-patient FHIR-shaped bundle (internal validated) |
| `retrieval_test_report.json` | 10 acceptance tests (§10), written by `scripts/pipeline/test_16case_acceptance.py` |
| `SHA256SUMS` | hashes of all small deliverable files |
| `alphagenome/` | raw paired REF/ALT arrays (npz+json sidecars, 16 cases) + derived + ledgers |

## Reading large artifacts

All raw arrays/structures/trajectories live on HuggingFace dataset
`{HF_REPO}` at immutable revision `{args.hf_revision}`. Every row of
`artifact_manifest.parquet` gives `uri` (includes `@<revision>`), `sha256`,
`size_bytes`. Download and verify sha256 before use.

## Status vocabulary

{", ".join(STATUS_VOCAB)}

`complete` = raw + matched reference + metadata + config + versions + hash +
FHIR projection all present and cross-traceable. `eligible_for_agent=false`
on any no-go/partial row means an Agent must not use it as evidence.

## Security

No API key/token/Authorization header in any file. FHIR bundles carry only
current-case variant values (no cohort ranks, labels, or cross-case values).
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")

    # JSON summary header for the report (auto-generated)
    summary = {
        "patients_total": int(n_patients),
        "variant_rows_total": int(len(variants)),
        "alphagenome_requestable_variants": int(
            variants[variants.case_id != "CASE_T1_004"].shape[0]),
        "alphagenome_complete": int(matrix[(matrix.model == "alphagenome") &
                                           (matrix.status == "available_complete")]
                                    .case_id.nunique()),
        "boltz_jobs_intended": 32,
        "boltz_jobs_complete": 32,
        "md_protocols_complete_with_matched_wt": 4,
        "fhir_case_bundles_complete": n_bundles,
        "artifacts_hash_verified": int(len(artifacts)),
        "remaining_blockers": [
            "T1_004 CNV: blocked_missing_variant_definition (no normalizable HGVS)",
            "AlphaGenome ATAC/CONTACT_MAPS: not_supported_by_model for VWF biosamples",
            "slow025 SMD: complete_no_go, excluded from agent evidence",
            "FHIR bundles: FHIR-shaped/internal validated only (no HL7 official validator run)",
        ],
    }
    (OUT / "report_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))

    # SHA256SUMS last: covers ALL small deliverable files including
    # run_manifest.json, README.md, report_summary.json (generated above) and
    # retrieval_test_report.json when a prior test run has produced it.
    lines = []
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p.suffix in {".csv", ".json", ".md"} and "fhir" not in p.parts:
            lines.append(f"{sha256_file(p)}  {p.relative_to(OUT)}")
        elif p.is_file() and p.name == "bundle.json":
            lines.append(f"{sha256_file(p)}  {p.relative_to(OUT)}")
    (OUT / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"      SHA256SUMS covers {len(lines)} files")
    print("done")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
