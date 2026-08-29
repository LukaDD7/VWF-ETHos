# Merged Type-1 + Type-2B offline GPU run

Date: 2026-08-29

## CPU-side work already completed

AlphaGenome was run on the network-enabled CPU instance:

| Panel | Requests | Successful requests | Raw output rows | Selected-track raw successes | Explicit unavailable outputs |
|---|---:|---:|---:|---:|---:|
| Type 1 | 9 | 9 | 99 | 81 | 18 |
| Type 2B | 7 | 7 | 77 | 63 | 14 |

The unavailable rows are ATAC and contact maps for biosamples relevant to VWF;
they are recorded as missing rather than replaced with unrelated cells.

Outputs:

```text
outputs/type1_panel_agent_20260828/server_bundle/results/alphagenome/
outputs/type2b_panel_agent_20260829/server_bundle/results/alphagenome/
```

The CPU instance also downloaded and cleaned experimental structure 7A6O and
prebuilt the shared WT plus six FoldX mutants:

```text
structures/7A6O_AIM_A1_clean.pdb
structures/7a6o_mutants/PANEL_WT_MATCHED.pdb
structures/7a6o_mutants/P1413L.pdb
structures/7a6o_mutants/R1308C.pdb
structures/7a6o_mutants/S1310F.pdb
structures/7a6o_mutants/V1316M.pdb
structures/7a6o_mutants/R1341W.pdb
structures/7a6o_mutants/A1461D.pdb
```

Therefore the GPU node does not need AlphaGenome, FoldX, 7A6O download, or any
network access.

## GPU-side one-command run

From the repository root on the shared GPU instance:

```bash
cd /lzy/projects/VWF-ETHos
GPUS=7 bash scripts/pipeline/run_type1_type2b_gpu_offline.sh
```

For four GPUs:

```bash
cd /lzy/projects/VWF-ETHos
GPUS=4 GPU_IDS=0,1,2,3 bash scripts/pipeline/run_type1_type2b_gpu_offline.sh
```

The script uses the shared environments:

```text
/lzy/envs/boltz2
/lzy/envs/gromacs
/lzy/envs/tools
```

## Workload

- Combined Boltz-2: 15 Type-1 jobs + 20 Type-2B jobs = 35 jobs.
- Targeted MD: one shared WT plus P1413L, R1308C, S1310F, V1316M, R1341W, and
  A1461D = seven 50-ns trajectories.
- Agent ingestion and rerun are executed automatically after Boltz and MD pass
  their inventory checks.

## Output locations

```text
outputs/computational_panel_20260829/boltz/
outputs/computational_panel_20260829/md/
outputs/type1_panel_agent_20260828/server_results/
outputs/type2b_panel_agent_20260829/server_results/
```

Expected final agent runs:

```text
outputs/type1_panel_agent_20260828/server_results/agent_run/
outputs/type2b_panel_agent_20260829/server_results/agent_run/
```
