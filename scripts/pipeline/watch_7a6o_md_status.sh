#!/bin/bash
set -u
ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
MD_ROOT="$ROOT_DIR/output/gromacs_md_autoinhib"
WT_DIR="$MD_ROOT/7A6O_WT/md_7a6o"
RELAX_LOG="$MD_ROOT/_mutant_relax_logs"

echo "=== 7A6O MD status: $(date '+%F %T') ==="
echo

echo "[Processes]"
pgrep -af 'run_7a6o_autoinhib_md_batch|run_7a6o_variant_direct|launch_7a6o_followup_queue|gmx mdrun|md_prod' || true

echo
echo "[WT production files]"
if [ -d "$WT_DIR" ]; then
  ls -lh --time-style=long-iso "$WT_DIR"/md_prod.* "$WT_DIR"/md_prod_run.* 2>/dev/null || true
else
  echo "missing: $WT_DIR"
fi

echo
echo "[WT production log tail]"
wt_log=""
for cand in "$WT_DIR"/md_prod_run.part*.log "$WT_DIR/md_prod.log"; do
  [ -f "$cand" ] && wt_log="$cand"
done
if [ -n "$wt_log" ]; then
  echo "log: $wt_log"
  tail -40 "$wt_log"
else
  echo "missing: $WT_DIR/md_prod.log or $WT_DIR/md_prod_run.part*.log"
fi

echo
echo "[Mutant relaxed structures]"
count=$(find "$MD_ROOT" -path '*/relax_pdb/solv_ions_em.gro' -type f 2>/dev/null | wc -l)
echo "solv_ions_em.gro count: $count"
find "$MD_ROOT" -path '*/relax_pdb/solv_ions_em.gro' -type f 2>/dev/null | sed "s#^$MD_ROOT/##" | sort

echo
echo "[Mutant production status]"
for variant_dir in "$MD_ROOT"/*; do
  [ -d "$variant_dir/md_7a6o" ] || continue
  variant="$(basename "$variant_dir")"
  [ "$variant" = "7A6O_WT" ] && continue
  work="$variant_dir/md_7a6o"
  if [ -f "$work/md_prod.gro" ] || compgen -G "$work/md_prod.part*.gro" >/dev/null; then
    state="complete"
  elif [ -f "$work/md_prod.cpt" ] || compgen -G "$work/md_prod.part*.log" >/dev/null; then
    state="running_or_resumable"
  else
    state="not_started"
  fi
  latest_log="$(ls -t "$work"/md_prod.part*.log "$work"/md_prod.log 2>/dev/null | head -1 || true)"
  latest_cpt=""
  for cand in "$work/md_prod.cpt" "$work"/md_prod.part*.cpt; do
    [ -f "$cand" ] && latest_cpt="$cand"
  done
  echo "=== $variant: $state ==="
  [ -n "$latest_cpt" ] && ls -lh --time-style=long-iso "$latest_cpt"
  [ -n "$latest_log" ] && {
    echo "log: $latest_log"
    grep -E 'Step|Writing checkpoint|Performance:|Finished mdrun|Fatal error|LINCS WARNING|cudaError' "$latest_log" | tail -8 || true
  }
done

echo
echo "[Relax loop tail]"
if [ -f /tmp/relax_loop.log ]; then
  tail -20 /tmp/relax_loop.log
else
  echo "missing: /tmp/relax_loop.log"
fi
