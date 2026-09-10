# Offline GPU task-agent protocol

These rules apply to every file below `mechanism_tasks/`.

1. Treat `requests/<task_id>/request.json` and `task.fhir.json` as immutable.
   Never edit the request, registry, protocol version, controls, ROIs, simulation
   parameters, or requested measurement set on a compute branch.
2. Before computation, run:

   ```bash
   python scripts/mechanism_task.py validate-request <task_id>
   python scripts/mechanism_task.py result-template <task_id>
   ```

3. Execute only the protocol and execution tier named by the request. If an
   input, structure, executable, resource, or QC prerequisite is missing, return
   `inconclusive` or `failed` with explicit QC failures and limitations. Do not
   improvise a substitute protocol. Remove the `template_only` marker before
   validation; its presence intentionally makes an untouched template invalid.
4. Compare case and matched WT with identical starting scaffold, preparation,
   tier, analysis selections, and units. Run the locked positive/negative or
   boundary controls required by the protocol. Record software/model versions,
   random seeds, runtime, ROI definitions, and replicate-level values.
5. Write the final contract only to `results/<task_id>/result.json`. Keep raw
   trajectories outside Git. Commit only the compact result, checksummed
   manifests, representative structures/plots when requested, and logs needed
   to audit QC.
6. Validate before committing:

   ```bash
   python scripts/mechanism_task.py validate-result \
     <task_id> mechanism_tasks/results/<task_id>/result.json
   ```

7. Report only `supports`, `contradicts`, or `indeterminate` mechanism evidence.
   Never state that MD/AI confirms a VWD subtype, pathogenicity, or ACMG PS3.
