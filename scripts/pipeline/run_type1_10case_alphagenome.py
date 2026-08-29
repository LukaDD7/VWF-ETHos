#!/usr/bin/env python3
"""Run the model-ready Type-1 cases through AlphaGenome's complete profile.

The official API exposes 11 raw OutputType values.  This runner queries model
metadata first, chooses VWF-relevant endothelial/vascular tracks independently
for every output type, keeps non-biosample outputs as global, and records an
explicit missing reason instead of substituting an unrelated cell type.

It also requests all official recommended variant scorers (currently 19 score
views over the 11 raw outputs, including active-track and polyadenylation
derivatives).  API credentials are read only from ALPHAGENOME_API_KEY.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
from typing import Any, Iterable

import numpy as np
import pandas as pd


RAW_OUTPUTS: tuple[str, ...] = (
    "ATAC", "CAGE", "DNASE", "RNA_SEQ", "CHIP_HISTONE", "CHIP_TF",
    "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS",
    "CONTACT_MAPS", "PROCAP",
)
OUTPUT_ATTR = {name: name.casefold() for name in RAW_OUTPUTS}
DEFAULT_BIOSAMPLE_REGEX = (
    r"endothelial|vascular|blood vessel|umbilical vein|huvec|"
    r"megakaryocyte|platelet"
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_requests(path: Path, expected_requests: int | None = None) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {
        "case_id", "chromosome", "position", "ref", "alt",
        "interval_start", "interval_end",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Missing AlphaGenome request columns: {missing}")
    if frame["case_id"].duplicated().any():
        raise ValueError("Duplicate case_id in AlphaGenome requests")
    if expected_requests is not None and len(frame) != expected_requests:
        raise ValueError(f"Expected {expected_requests} model-ready requests; found {len(frame)}")
    requested = (
        frame["requested_outputs"].astype(str)
        if "requested_outputs" in frame
        else pd.Series("|".join(RAW_OUTPUTS), index=frame.index)
    )
    if set(requested.str.split("|").explode()) != set(RAW_OUTPUTS):
        raise ValueError("Request profile must contain all 11 official raw OutputTypes")
    return frame


def _metadata_frame(metadata: Any, output_name: str) -> pd.DataFrame:
    value = getattr(metadata, OUTPUT_ATTR[output_name], None)
    if isinstance(value, pd.DataFrame):
        return value.copy()
    return pd.DataFrame()


def build_ontology_plan(
    metadata: Any,
    preferred_terms: Iterable[str],
    biosample_regex: str = DEFAULT_BIOSAMPLE_REGEX,
    max_terms_per_output: int = 12,
) -> tuple[dict[str, dict[str, Any]], pd.DataFrame]:
    preferred = list(dict.fromkeys(preferred_terms))
    pattern = re.compile(biosample_regex, flags=re.IGNORECASE)
    plan: dict[str, dict[str, Any]] = {}
    inventory: list[pd.DataFrame] = []

    for output_name in RAW_OUTPUTS:
        frame = _metadata_frame(metadata, output_name)
        frame.insert(0, "output_type", output_name)
        if frame.empty:
            plan[output_name] = {
                "strategy": "unavailable",
                "ontology_terms": [],
                "track_count": 0,
                "reason": "metadata_empty",
            }
            inventory.append(frame)
            continue

        ontology = frame.get("ontology_curie", pd.Series(index=frame.index, dtype="object"))
        named = pd.concat(
            [
                frame.get("biosample_name", pd.Series("", index=frame.index)).fillna("").astype(str),
                frame.get("name", pd.Series("", index=frame.index)).fillna("").astype(str),
            ],
            axis=1,
        ).agg(" ".join, axis=1)
        exact_mask = ontology.astype(str).isin(preferred)
        relevant_mask = named.str.contains(pattern, na=False)
        selected_mask = exact_mask | relevant_mask

        if ontology.isna().all() or ontology.astype(str).str.strip().isin({"", "nan", "None"}).all():
            strategy = "global_non_biosample"
            selected_terms: list[str] = []
            selected_mask = pd.Series(True, index=frame.index)
            reason = "output_has_no_biosample_ontology"
        else:
            terms = [
                value for value in ontology[selected_mask].dropna().astype(str).tolist()
                if value and value not in {"nan", "None"}
            ]
            selected_terms = list(dict.fromkeys(terms))[:max_terms_per_output]
            selected_mask &= ontology.astype(str).isin(selected_terms)
            if selected_terms:
                strategy = "vwf_relevant_biosample"
                reason = "preferred_or_endothelial_vascular_megakaryocyte_track"
            else:
                strategy = "unavailable"
                reason = "no_vwf_relevant_biosample_track"

        frame["selected_for_type1_panel"] = selected_mask
        frame["selection_reason"] = np.where(selected_mask, reason, "not_selected")
        inventory.append(frame)
        plan[output_name] = {
            "strategy": strategy,
            "ontology_terms": selected_terms,
            "track_count": int(selected_mask.sum()),
            "reason": reason,
        }

    return plan, pd.concat(inventory, ignore_index=True, sort=False)


def _values(value: Any) -> np.ndarray | None:
    if value is None:
        return None
    raw = getattr(value, "values", value)
    try:
        array = np.asarray(raw)
    except Exception:
        return None
    return array if array.size else None


def summarize_raw_delta(case_id: str, output_name: str, ref: Any, alt: Any) -> dict[str, Any]:
    ref_values = _values(ref)
    alt_values = _values(alt)
    base = {"case_id": case_id, "output_type": output_name}
    if ref_values is None or alt_values is None:
        return {**base, "status": "not_returned", "reason": "empty_or_missing_output"}
    if ref_values.shape != alt_values.shape:
        return {
            **base,
            "status": "failed_qc",
            "reason": f"ref_alt_shape_mismatch:{ref_values.shape}!={alt_values.shape}",
        }
    delta = alt_values.astype(np.float64) - ref_values.astype(np.float64)
    abs_flat_index = int(np.nanargmax(np.abs(delta)))
    peak_index = np.unravel_index(abs_flat_index, delta.shape)
    return {
        **base,
        "status": "success",
        "reason": "",
        "shape": "x".join(map(str, delta.shape)),
        "signed_min": float(np.nanmin(delta)),
        "signed_max": float(np.nanmax(delta)),
        "abs_max": float(np.nanmax(np.abs(delta))),
        "signed_mean": float(np.nanmean(delta)),
        "positive_sum": float(delta[delta > 0].sum(initial=0.0)),
        "negative_sum": float(delta[delta < 0].sum(initial=0.0)),
        "peak_index": "|".join(map(str, peak_index)),
    }


def _filter_tidy_scores(
    tidy: pd.DataFrame,
    case_id: str,
    scorer_key: str,
    ontology_plan: dict[str, dict[str, Any]],
) -> pd.DataFrame:
    if tidy.empty:
        return tidy
    tidy = tidy.copy()
    tidy.insert(0, "case_id", case_id)
    tidy.insert(1, "scorer_key", scorer_key)
    if "gene_name" in tidy.columns and tidy["gene_name"].notna().any():
        tidy = tidy[tidy["gene_name"].astype(str).str.upper().eq("VWF")]
    output_name = str(tidy["output_type"].iloc[0]) if not tidy.empty else ""
    selected_terms = ontology_plan.get(output_name, {}).get("ontology_terms", [])
    strategy = ontology_plan.get(output_name, {}).get("strategy")
    if selected_terms and "ontology_curie" in tidy.columns:
        tidy = tidy[tidy["ontology_curie"].astype(str).isin(selected_terms)]
    elif strategy == "unavailable":
        tidy = tidy.iloc[0:0]
    return tidy.reset_index(drop=True)


def _write_case_outputs(case_dir: Path, raw_rows: list[dict[str, Any]], score_rows: list[pd.DataFrame]) -> None:
    case_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(raw_rows).to_csv(case_dir / "raw_summary.csv", index=False)
    long = pd.concat(score_rows, ignore_index=True, sort=False) if score_rows else pd.DataFrame()
    long.to_csv(case_dir / "scores_long.csv", index=False)


def _finalize(output_dir: Path, request_frame: pd.DataFrame, scorer_keys: list[str]) -> None:
    raw_files = sorted((output_dir / "cases").glob("*/raw_summary.csv"))
    score_files = sorted((output_dir / "cases").glob("*/scores_long.csv"))
    raw = pd.concat([pd.read_csv(path) for path in raw_files], ignore_index=True, sort=False)
    nonempty_scores = [frame for path in score_files if not (frame := pd.read_csv(path)).empty]
    scores = pd.concat(nonempty_scores, ignore_index=True, sort=False) if nonempty_scores else pd.DataFrame()
    raw.to_csv(output_dir / "alphagenome_raw_summary.csv", index=False)
    scores.to_csv(output_dir / "alphagenome_scores_long.csv", index=False)
    try:
        scores.to_parquet(output_dir / "alphagenome_scores_long.parquet", index=False)
    except Exception as exc:
        (output_dir / "parquet_warning.txt").write_text(str(exc), encoding="utf-8")

    base_columns = ["case_id", "chromosome", "position", "ref", "alt", "hgvs_c", "hgvs_p"]
    if "patient_id" in request_frame.columns:
        base_columns.insert(1, "patient_id")
    base = request_frame[base_columns].copy()
    feature_rows: list[dict[str, Any]] = []
    for case_id in base["case_id"]:
        row: dict[str, Any] = {"case_id": case_id}
        subset = scores[scores["case_id"].eq(case_id)] if not scores.empty else pd.DataFrame()
        for key in scorer_keys:
            values = pd.to_numeric(subset.loc[subset["scorer_key"].eq(key), "raw_score"], errors="coerce").dropna()
            prefix = f"ag_{key.casefold()}"
            row[f"{prefix}_signed_min"] = float(values.min()) if not values.empty else np.nan
            row[f"{prefix}_signed_max"] = float(values.max()) if not values.empty else np.nan
            row[f"{prefix}_abs_max"] = float(values.abs().max()) if not values.empty else np.nan
        feature_rows.append(row)
    features = base.merge(pd.DataFrame(feature_rows), on="case_id", how="left")
    features.to_csv(output_dir / "alphagenome_agent_features.csv", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--requests", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--preferred-ontology", action="append", default=["CL:0000115"])
    parser.add_argument("--biosample-regex", default=DEFAULT_BIOSAMPLE_REGEX)
    parser.add_argument("--expected-requests", type=int)
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--plan-only", action="store_true", help="Validate the request/11-output inventory without API access.")
    parser.add_argument("--no-raw", action="store_true")
    parser.add_argument("--no-scores", action="store_true")
    args = parser.parse_args()

    requests = validate_requests(args.requests, args.expected_requests)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    static_plan = {
        "created_at": utc_now(),
        "request_sha256": sha256(args.requests),
        "case_count": len(requests),
        "raw_output_count": len(RAW_OUTPUTS),
        "raw_outputs": list(RAW_OUTPUTS),
        "expected_raw_case_output_pairs": len(requests) * len(RAW_OUTPUTS),
        "credential_source": "ALPHAGENOME_API_KEY environment variable",
    }
    (args.output_dir / "static_run_plan.json").write_text(json.dumps(static_plan, indent=2), encoding="utf-8")
    if args.plan_only:
        print(json.dumps(static_plan, indent=2))
        return 0

    api_key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not api_key:
        raise SystemExit("ALPHAGENOME_API_KEY is required")
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    model = dna_client.create(api_key)
    metadata = model.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
    ontology_plan, metadata_inventory = build_ontology_plan(
        metadata,
        preferred_terms=args.preferred_ontology,
        biosample_regex=args.biosample_regex,
    )
    metadata_inventory.to_csv(args.output_dir / "output_metadata_inventory.csv", index=False)
    (args.output_dir / "ontology_plan.json").write_text(
        json.dumps(ontology_plan, indent=2), encoding="utf-8"
    )
    if args.metadata_only:
        print(json.dumps(ontology_plan, indent=2))
        return 0

    scorer_map = dict(variant_scorers.RECOMMENDED_VARIANT_SCORERS)
    scorer_keys = list(scorer_map)
    completed: set[str] = set()
    checkpoint = args.output_dir / "checkpoint.jsonl"
    if checkpoint.exists():
        for line in checkpoint.read_text(encoding="utf-8").splitlines():
            record = json.loads(line)
            if record.get("status") == "success":
                completed.add(record["case_id"])

    for _, request in requests.iterrows():
        case_id = str(request["case_id"])
        if case_id in completed:
            continue
        variant = genome.Variant(
            chromosome=str(request["chromosome"]),
            position=int(request["position"]),
            reference_bases=str(request["ref"]),
            alternate_bases=str(request["alt"]),
        )
        interval = genome.Interval(
            chromosome=str(request["chromosome"]),
            start=int(request["interval_start"]),
            end=int(request["interval_end"]),
        )
        raw_rows: list[dict[str, Any]] = []
        score_rows: list[pd.DataFrame] = []
        try:
            if not args.no_raw:
                grouped: dict[tuple[str, ...] | None, list[str]] = defaultdict(list)
                for output_name, profile in ontology_plan.items():
                    if profile["strategy"] == "unavailable":
                        raw_rows.append({
                            "case_id": case_id,
                            "output_type": output_name,
                            "status": "not_available_for_selected_cells",
                            "reason": profile["reason"],
                        })
                        continue
                    terms = tuple(profile["ontology_terms"]) or None
                    grouped[terms].append(output_name)
                for terms, output_names in grouped.items():
                    result = model.predict_variant(
                        interval=interval,
                        variant=variant,
                        requested_outputs=[getattr(dna_client.OutputType, name) for name in output_names],
                        ontology_terms=list(terms) if terms is not None else None,
                    )
                    for output_name in output_names:
                        attr = OUTPUT_ATTR[output_name]
                        raw_rows.append(summarize_raw_delta(
                            case_id,
                            output_name,
                            getattr(result.reference, attr, None),
                            getattr(result.alternate, attr, None),
                        ))
            if not args.no_scores:
                score_outputs = model.score_variant(
                    interval=interval,
                    variant=variant,
                    variant_scorers=list(scorer_map.values()),
                )
                for scorer_key, adata in zip(scorer_keys, score_outputs, strict=True):
                    tidy = variant_scorers.tidy_anndata(adata)
                    filtered = _filter_tidy_scores(tidy, case_id, scorer_key, ontology_plan)
                    if not filtered.empty:
                        score_rows.append(filtered)
            _write_case_outputs(args.output_dir / "cases" / case_id, raw_rows, score_rows)
            record = {"case_id": case_id, "status": "success", "finished_at": utc_now()}
        except Exception as exc:
            record = {"case_id": case_id, "status": "failed", "finished_at": utc_now(), "error": str(exc)}
        with checkpoint.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, ensure_ascii=False) + "\n")
        print(json.dumps(record, ensure_ascii=False), flush=True)

    successful_cases: set[str] = set()
    failed_cases: list[dict[str, Any]] = []
    if checkpoint.exists():
        for line in checkpoint.read_text(encoding="utf-8").splitlines():
            record = json.loads(line)
            if record.get("status") == "success":
                successful_cases.add(str(record.get("case_id")))
            else:
                failed_cases.append(record)
    missing_cases = sorted(set(requests["case_id"].astype(str)) - successful_cases)
    if missing_cases or failed_cases:
        raise RuntimeError(
            json.dumps(
                {"missing_cases": missing_cases, "failed_cases": failed_cases},
                ensure_ascii=False,
            )
        )

    _finalize(args.output_dir, requests, scorer_keys)
    manifest = {
        **static_plan,
        "finished_at": utc_now(),
        "recommended_scorer_count": len(scorer_keys),
        "recommended_scorers": scorer_keys,
        "sign_convention": "raw deltas and signed scorer values retain ALT minus REF direction where the official scorer is signed",
        "ontology_plan": ontology_plan,
        "outputs": {
            "agent_features": "alphagenome_agent_features.csv",
            "raw_summary": "alphagenome_raw_summary.csv",
            "scores_long": "alphagenome_scores_long.csv",
        },
    }
    (args.output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
