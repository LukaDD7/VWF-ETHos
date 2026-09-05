#!/usr/bin/env python3
"""16-case AlphaGenome evidence run with full paired REF/ALT persistence.

Implements the requirements in
docs/SERVER_16CASE_COMPUTATIONAL_DATA_COMPLETION_REQUIREMENTS_20260905.md §4:

- Request source: the two frozen server_bundle input CSVs (9 Type 1 + 7 Type 2B
  variant rows = 16 normalizable variants, 16 patients; T1_004 is a CNV and is
  NOT requested here - it gets blocked_missing_variant_definition upstream).
- For every requested output type, persists BOTH complete reference and
  alternate value arrays plus track metadata, interval, resolution, dtype,
  and an ALT-REF alignment check (persist_paired_tracks).
- Groups outputs by ontology-term tuple so one API call serves several
  modalities sharing the same cell-type selection, and records per-call
  request/response metadata (no secrets).
- ATAC/CONTACT_MAPS have no VWF-relevant biosample in official metadata;
  they are recorded as not_supported_by_model with the API-observed reason
  (empty track arrays), never zero-filled.
- Raw (npz + JSON sidecar) and derived (per-track paired summaries) layers
  are written separately under output_dir/raw/ and output_dir/derived/.

API key is read ONLY from the ALPHAGENOME_API_KEY environment variable and is
never written to any output file.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

RAW_OUTPUTS = (
    "ATAC", "CAGE", "DNASE", "RNA_SEQ", "CHIP_HISTONE", "CHIP_TF",
    "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS",
    "CONTACT_MAPS", "PROCAP",
)
OUTPUT_ATTR = {name: name.casefold() for name in RAW_OUTPUTS}
PREFERRED_TERMS = ["CL:0000115"]  # endothelial cell - VWF-expressing
BIOSAMPLE_REGEX = (
    r"endothelial|vascular|blood vessel|umbilical vein|huvec|"
    r"megakaryocyte|platelet"
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _values(track_data) -> np.ndarray | None:
    if track_data is None:
        return None
    raw = getattr(track_data, "values", None)
    if raw is None:
        return None
    try:
        array = np.asarray(raw)
    except Exception:
        return None
    return array if array.size else None


def _metadata_records(track_data) -> list[dict]:
    meta = getattr(track_data, "metadata", None)
    if meta is None or not isinstance(meta, pd.DataFrame) or meta.empty:
        return []
    return json.loads(meta.to_json(orient="records"))


def _safe_uns(track_data) -> dict:
    uns = getattr(track_data, "uns", None)
    if isinstance(uns, dict):
        return {str(k): str(v) for k, v in uns.items()}
    return {}


def build_ontology_plan(metadata, preferred=PREFERRED_TERMS,
                        regex=BIOSAMPLE_REGEX, max_terms=12):
    """Select VWF-relevant biosample tracks per output type from metadata."""
    import re
    pattern = re.compile(regex, flags=re.IGNORECASE)
    plan = {}
    for output_name in RAW_OUTPUTS:
        frame = getattr(metadata, OUTPUT_ATTR[output_name], None)
        if frame is None or (isinstance(frame, pd.DataFrame) and frame.empty):
            plan[output_name] = {
                "strategy": "unavailable",
                "ontology_terms": [],
                "track_count": 0,
                "reason": "metadata_empty",
            }
            continue
        ontology = frame.get("ontology_curie")
        if ontology is None or ontology.isna().all() or ontology.astype(str).str.strip().isin(["", "nan", "None"]).all():
            plan[output_name] = {
                "strategy": "global_non_biosample",
                "ontology_terms": [],
                "track_count": int(len(frame)),
                "reason": "output_has_no_biosample_ontology",
            }
            continue
        named = pd.concat(
            [
                frame.get("biosample_name", pd.Series("", index=frame.index)).fillna("").astype(str),
                frame.get("name", pd.Series("", index=frame.index)).fillna("").astype(str),
            ],
            axis=1,
        ).agg(" ".join, axis=1)
        selected = ontology.astype(str).isin(preferred) | named.str.contains(pattern, na=False)
        terms = [t for t in ontology[selected].dropna().astype(str).tolist() if t and t not in {"nan", "None"}]
        terms = list(dict.fromkeys(terms))[:max_terms]
        if terms:
            plan[output_name] = {
                "strategy": "vwf_relevant_biosample",
                "ontology_terms": terms,
                "track_count": int(selected.sum()),
                "reason": "preferred_or_endothelial_vascular_megakaryocyte_track",
            }
        else:
            plan[output_name] = {
                "strategy": "unavailable",
                "ontology_terms": [],
                "track_count": 0,
                "reason": "no_vwf_relevant_biosample_track",
            }
    return plan


def persist_paired_tracks(case_id: str, output_name: str, result,
                           raw_dir: Path) -> dict:
    """Losslessly persist REF and ALT arrays + metadata for one output type.

    Writes raw/<case_id>__<output_name>.npz (reference/alternate arrays) and a
    JSON sidecar with track metadata, interval, resolution, dtype, shape, and
    the alignment check.  Returns the row for the run ledger.
    """
    attr = OUTPUT_ATTR[output_name]
    ref_td = getattr(result.reference, attr, None)
    alt_td = getattr(result.alternate, attr, None)
    ref_values = _values(ref_td)
    alt_values = _values(alt_td)

    row = {
        "case_id": case_id,
        "output_type": output_name,
        "requested_at": utc_now(),
    }
    if ref_values is None and alt_values is None:
        row.update({
            "status": "not_returned",
            "reason": "empty_or_missing_output_for_selected_tracks",
            "arrays_persisted": False,
        })
        return row
    if ref_values is not None and alt_values is not None and ref_values.shape != alt_values.shape:
        row.update({
            "status": "failed_qc",
            "reason": f"ref_alt_shape_mismatch:{ref_values.shape}!={alt_values.shape}",
            "arrays_persisted": False,
        })
        return row

    interval = getattr(ref_td, "interval", None) if ref_td is not None else None
    interval_str = str(interval) if interval is not None else None
    resolution = getattr(ref_td, "resolution", None) if ref_td is not None else None
    metadata_records = _metadata_records(ref_td)
    uns = _safe_uns(ref_td)

    # Track identity alignment: metadata row order must match array columns.
    n_tracks = ref_values.shape[-1] if ref_values is not None else 0
    identity_aligned = bool(metadata_records and len(metadata_records) == n_tracks)

    npz_path = raw_dir / f"{case_id}__{output_name}.npz"
    arrays = {}
    if ref_values is not None:
        arrays["reference"] = ref_values.astype(np.float32, copy=False)
    if alt_values is not None:
        arrays["alternate"] = alt_values.astype(np.float32, copy=False)
    np.savez_compressed(npz_path, **arrays)

    sidecar = {
        "case_id": case_id,
        "output_type": output_name,
        "npz_file": npz_path.name,
        "npz_bytes": npz_path.stat().st_size,
        "interval": interval_str,
        "resolution": resolution,
        "shape": list(ref_values.shape) if ref_values is not None else None,
        "dtype": "float32",
        "track_metadata": metadata_records,
        "uns": uns,
        "track_identity_aligned": identity_aligned,
        "saved_keys": sorted(arrays.keys()),
    }
    sidecar_path = raw_dir / f"{case_id}__{output_name}.json"
    sidecar_path.write_text(json.dumps(sidecar, indent=2), encoding="utf-8")

    row.update({
        "status": "success",
        "reason": "",
        "shape": "x".join(map(str, ref_values.shape)) if ref_values is not None else "x".join(map(str, alt_values.shape)),
        "n_tracks": int(n_tracks),
        "interval": interval_str,
        "resolution": resolution,
        "track_identity_aligned": identity_aligned,
        "npz_bytes": npz_path.stat().st_size,
        "arrays_persisted": True,
    })
    return row


def summarize_paired_tracks(case_id: str, output_name: str, npz_path: Path,
                            track_metadata: list[dict]) -> list[dict]:
    """Derive per-track paired summaries from the persisted raw arrays."""
    with np.load(npz_path) as data:
        ref = data["reference"] if "reference" in data.files else None
        alt = data["alternate"] if "alternate" in data.files else None
    rows = []
    n_tracks = (ref.shape[-1] if ref is not None else alt.shape[-1])
    for track_index in range(n_tracks):
        track_name = (
            track_metadata[track_index].get("name", f"track_{track_index}")
            if track_metadata and track_index < len(track_metadata)
            else f"track_{track_index}"
        )
        ref_col = ref[:, track_index] if ref is not None else None
        alt_col = alt[:, track_index] if alt is not None else None
        row = {
            "case_id": case_id,
            "output_type": output_name,
            "track_index": track_index,
            "track_name": track_name,
        }
        if ref_col is not None and alt_col is not None:
            delta = alt_col.astype(np.float64) - ref_col.astype(np.float64)
            peak = int(np.nanargmax(np.abs(delta)))
            row.update({
                "status": "paired",
                "signed_min": float(np.nanmin(delta)),
                "signed_max": float(np.nanmax(delta)),
                "abs_max": float(np.nanmax(np.abs(delta))),
                "signed_mean": float(np.nanmean(delta)),
                "peak_position": peak,
            })
        else:
            row.update({"status": "single_side_only"})
        rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--requests", nargs=2, required=True,
                        help="Two request CSVs (type1 bundle, type2b bundle)")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--case", help="Run a single case_id only (debug)")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = args.output_dir / "raw"
    derived_dir = args.output_dir / "derived"
    raw_dir.mkdir(parents=True, exist_ok=True)
    derived_dir.mkdir(parents=True, exist_ok=True)

    requests = pd.concat(
        [pd.read_csv(p) for p in args.requests], ignore_index=True
    )
    requests = requests.astype({"case_id": str, "chromosome": str,
                                "position": int, "ref": str, "alt": str,
                                "interval_start": int, "interval_end": int})
    if args.case:
        requests = requests[requests["case_id"].eq(args.case)]
    print(f"cases to run: {len(requests)}")

    # API key: environment only, never persisted.
    api_key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not api_key:
        raise SystemExit("ALPHAGENOME_API_KEY is required")

    sys.path.insert(0, "/tmp/alphagenome_env")
    import alphagenome
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    model = dna_client.create(api_key)

    # --- metadata + ontology plan (shared across cases) ---
    metadata = model.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
    plan = build_ontology_plan(metadata)
    (args.output_dir / "ontology_plan.json").write_text(
        json.dumps(plan, indent=2), encoding="utf-8")
    print("ontology plan written")
    if args.metadata_only:
        print(json.dumps(plan, indent=2))
        return 0

    environment = {
        "alphagenome_version": alphagenome.__version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "requested_outputs": list(RAW_OUTPUTS),
        "preferred_ontology_terms": PREFERRED_TERMS,
    }

    ledger_path = args.output_dir / "run_ledger.jsonl"
    derived_rows_all: list[dict] = []

    for _, request in requests.iterrows():
        case_id = str(request["case_id"])
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
        request_record = {
            "case_id": case_id,
            "request": {
                "assembly": str(request.get("genome_build", "GRCh38")),
                "chromosome": str(request["chromosome"]),
                "position": int(request["position"]),
                "reference_bases": str(request["ref"]),
                "alternate_bases": str(request["alt"]),
                "interval_start": int(request["interval_start"]),
                "interval_end": int(request["interval_end"]),
                "requested_outputs": list(RAW_OUTPUTS),
                "hgvs_c": request.get("hgvs_c", ""),
                "hgvs_p": request.get("hgvs_p", ""),
            },
            "requested_at": utc_now(),
        }

        # Group output types by ontology-term tuple → one API call per group.
        grouped: dict[tuple[str, ...] | None, list[str]] = {}
        skipped: list[dict] = []
        for output_name in RAW_OUTPUTS:
            profile = plan[output_name]
            if profile["strategy"] == "unavailable":
                skipped.append({
                    "case_id": case_id,
                    "output_type": output_name,
                    "status": "not_supported_by_model",
                    "reason": profile["reason"],
                    "arrays_persisted": False,
                })
                continue
            terms = tuple(profile["ontology_terms"]) or None
            grouped.setdefault(terms, []).append(output_name)

        ledger_rows: list[dict] = list(skipped)
        api_calls: list[dict] = []
        t0 = time.time()
        for terms, output_names in sorted(grouped.items(), key=lambda kv: (kv[0] is None, kv[0] or ())):
            call_record = {
                "ontology_terms": list(terms) if terms else None,
                "output_types": output_names,
            }
            attempts = []
            result = None
            for attempt in range(1, 4):
                try:
                    result = model.predict_variant(
                        interval=interval,
                        variant=variant,
                        requested_outputs=[
                            getattr(dna_client.OutputType, name) for name in output_names
                        ],
                        ontology_terms=list(terms) if terms is not None else None,
                    )
                    attempts.append({"attempt": attempt, "status": "success"})
                    break
                except Exception as exc:
                    attempts.append({
                        "attempt": attempt,
                        "status": "error",
                        "error_type": type(exc).__name__,
                        "error": str(exc)[:500],
                    })
                    time.sleep(min(2 ** attempt, 30))
            call_record["attempts"] = attempts
            if result is None:
                for output_name in output_names:
                    ledger_rows.append({
                        "case_id": case_id,
                        "output_type": output_name,
                        "status": "failed_retry_exhausted",
                        "reason": "api_error_after_3_attempts",
                        "arrays_persisted": False,
                    })
                api_calls.append(call_record)
                continue

            for output_name in output_names:
                row = persist_paired_tracks(case_id, output_name, result, raw_dir)
                ledger_rows.append(row)
                if row.get("arrays_persisted"):
                    sidecar = json.loads(
                        (raw_dir / f"{case_id}__{output_name}.json").read_text()
                    )
                    derived_rows_all.extend(
                        summarize_paired_tracks(
                            case_id, output_name,
                            raw_dir / f"{case_id}__{output_name}.npz",
                            sidecar["track_metadata"],
                        )
                    )
            api_calls.append(call_record)

        elapsed = time.time() - t0
        case_record = {
            **request_record,
            "api_calls": api_calls,
            "ledger": ledger_rows,
            "elapsed_seconds": round(elapsed, 1),
            "finished_at": utc_now(),
        }
        with ledger_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(case_record, ensure_ascii=False) + "\n")
        statuses = [r["status"] for r in ledger_rows]
        print(f"  {case_id}: {statuses.count('success')}/11 persisted "
              f"({elapsed:.0f}s) " +
              ("; ".join(f"{s}x{statuses.count(s)}" for s in set(statuses) if s != "success")))

    pd.DataFrame(derived_rows_all).to_csv(
        derived_dir / "paired_track_summaries.csv", index=False)

    # Flat ledger CSV for the completeness matrix
    flat = []
    if ledger_path.exists():
        for line in ledger_path.read_text().splitlines():
            record = json.loads(line)
            for row in record["ledger"]:
                flat.append({
                    "case_id": row["case_id"],
                    "output_type": row["output_type"],
                    "status": row.get("status", ""),
                    "reason": row.get("reason", ""),
                    "n_tracks": row.get("n_tracks", ""),
                    "npz_bytes": row.get("npz_bytes", ""),
                    "track_identity_aligned": row.get("track_identity_aligned", ""),
                })
    pd.DataFrame(flat).to_csv(args.output_dir / "ledger_flat.csv", index=False)

    (args.output_dir / "run_manifest.json").write_text(
        json.dumps({
            "created_at": utc_now(),
            "case_count": len(requests),
            "environment": environment,
            "ontology_plan": plan,
            "raw_layer": "raw/<case_id>__<output_type>.npz + .json sidecar",
            "derived_layer": "derived/paired_track_summaries.csv",
        }, indent=2), encoding="utf-8")
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
