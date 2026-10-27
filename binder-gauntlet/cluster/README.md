# Running the GPU steps as batch jobs on a cluster

> **Status: written down, not yet run on a customer's cluster.** The design loop this drives has been run many times on single data-centre GPUs, which is the guide's Route A. What is untested is pushing each step through a cluster scheduler from a Claude Code session. Expect to adjust node selection, storage and the image for your site.

This is the [guide](../README.md)'s Route B. Claude Code still holds the loop; the only change from Route A is that each GPU step in [`steps.md`](../steps.md) becomes one batch job that the session writes, submits with `kubectl` or `sbatch`, waits for, and reads back. There is deliberately no wrapper script here: your platform team already has conventions for GPU jobs, and a session that calls `kubectl` and `sbatch` directly, under the rules below, adapts to them faster than a tool of ours would. This page is the cluster-side checklist, one example manifest per object, and those rules. Do Route A once first, so that your `TOOLS.md` is proven before the cluster adds its own unknowns.

## What the cluster side needs

### Kubernetes

- **Where Claude Code runs:** any machine with `kubectl` and a kubeconfig context for the cluster. It needs no GPU itself.
- **One namespace, and an identity scoped to it** that can create, get, list, watch and delete Jobs, Pods and PersistentVolumeClaims, read `pods/log` and use `pods/exec` — nothing cluster-wide, no node access. A Role of roughly this shape, bound to your user or service account, is enough on a default cluster; check it against your cluster's admission policies, since some sites also require quotas, priority classes or pod-security labels on the namespace:

  ```yaml
  kind: Role
  apiVersion: rbac.authorization.k8s.io/v1
  metadata: {name: binder-design, namespace: <your-namespace>}
  rules:
    - {apiGroups: ["batch"], resources: ["jobs"], verbs: ["create", "get", "list", "watch", "patch", "delete"]}
    - {apiGroups: [""], resources: ["pods", "persistentvolumeclaims"], verbs: ["create", "get", "list", "watch", "patch", "delete"]}
    - {apiGroups: [""], resources: ["pods/log"], verbs: ["get"]}
    - {apiGroups: [""], resources: ["pods/exec"], verbs: ["create", "get"]}
  ```
- **A container image,** pullable from that namespace, with the tools from the guide's section 4 and their weights baked in, plus `bash` and `tar`. Building it is a task for whoever administers your compute: each tool's own repository is its install guide, the weights go in at build time so that no job downloads anything, and your `TOOLS.md` records what went in and how each tool is called inside it.
- **A work volume** that every job mounts at `/work`: one PersistentVolumeClaim, on a ReadWriteMany storage class if your cluster offers one (then pods may land on any GPU node); otherwise ReadWriteOnce with every pod pinned to one GPU node by a `kubernetes.io/hostname` selector, which is fine for a first run. Tens of gigabytes is enough when the weights are in the image. A small utility pod that mounts the same volume and sleeps gives `kubectl cp` and `kubectl exec` somewhere to land; [`workspace.example.yaml`](workspace.example.yaml) is the shape of both.
- **The GPU resource name** your device plugin advertises (`nvidia.com/gpu` unless your platform team says otherwise) and **the node labels and taints** your GPU pool carries, which become the Job's `nodeSelector` and `tolerations`; these vary by cloud and by site. If your role may read nodes, `kubectl --context <ctx> get nodes --show-labels` shows them; otherwise ask.
- **A GPU quota you know.** Jobs beyond it sit `Pending`, and a Job's deadline counts that waiting time.

Write the context, namespace, image, GPU resource name, node selector and tolerations, storage class and a short run id into the "Where the GPU work runs" section of your `TOOLS.md` (the guide's [`TOOLS.md.example`](../TOOLS.md.example) has the lines). [`gpu-job.example.yaml`](gpu-job.example.yaml) is the shape of one step: the session copies it per step, fills in the name, command and resources, applies it, waits for it with `kubectl wait` or by polling `kubectl get job`, saves `kubectl logs` under `jobs/` in your working folder, and reads the result. Small inputs go onto the volume and small outputs come back with `kubectl cp` against the utility pod; short CPU commands (listing a folder, a scoring script over confidence files) run there with `kubectl exec`.

### Slurm

Run Claude Code on a login node, from a working folder on the filesystem the compute nodes share, with the tools installed in an environment or container the compute nodes can activate. Each GPU step is a short batch script the session writes under `jobs/` and submits with `sbatch --wait`, so the command returns when the step ends; the script names the partition, `--gres=gpu:<n>`, `--cpus-per-task`, `--mem`, `--time` and whatever account or QoS options your site requires, activates the environment on its first line, and writes its log beside itself. `squeue -u $USER` and `sacct -j <id> --format=State,Elapsed,MaxRSS` give status and the reason for a failure; inputs and outputs are already on the shared filesystem, so nothing is copied. Put the partition, the extra `sbatch` options and the activation line in `TOOLS.md`. Check the `--gres` spelling and the time format your site expects.

## Rules the session follows on a cluster

Paste these rules into your first message.

1. **Stay in the box.** Every `kubectl` command carries `--context <ctx> -n <namespace>` explicitly, never whatever the current context happens to be; every `sbatch` names the partition. No cluster-scoped command, ever.
2. **Name, label, and remove narrowly.** Every object is named `<run-id>-<step>-<suffix>` and labelled `app=binder-design,run=<run-id>`. Objects are removed by exact name, after listing what carries both of this run's labels; never by one broad label on its own, never across namespaces. In a shared namespace a broad selector removes other people's work.
3. **Bounded jobs only.** Every Job has `backoffLimit: 0`, an `activeDeadlineSeconds`, CPU, memory and GPU limits, `automountServiceAccountToken: false`, no secrets mounted and no network use; its command tees its own output to `/work/logs/` so the log survives a deadline kill. Slurm scripts carry `--time` and write their log to the shared folder.
4. **One job per tool per batch of candidates,** never one per candidate: loading a model costs more than scoring one more complex. Several jobs at once only up to the GPU quota you state; long steps run in the background while the session reads earlier results and keeps the report current.
5. **On failure, change something.** Read the Job's status reason (`DeadlineExceeded`, `OOMKilled`, a non-zero exit) and the log tail, say in one line what went wrong and what will be different, then resubmit. Never the same command into the same failure.
6. **Pull back tables and selected models,** not whole prediction folders; keep `designs.csv` and `report.md` current in the working folder as the guide asks. Finish with nothing running, remove this run's utility pod by name, and leave the volume until you have copied what you need; remove it, by name, only when you say so.

## Security and data handling

- **Least privilege.** A namespace-scoped identity; no node access; no service-account token in any pod; no secrets mounted, because no step needs one.
- **No network from jobs.** With weights in the image and inputs on the volume, no step makes an outbound call. Where your network plugin enforces NetworkPolicy, a default-deny egress policy on the namespace turns that from a property of the tools into a guarantee.
- **What is on the volume** — target files, every candidate sequence and model, the scores — is, on your own targets, the confidential material the guide's section 10 describes. Put the namespace and its storage where such material may live, and remove the volume when the campaign is archived.
- **What goes to the model** is only what the session reads: job status, log tails, and the tables and files it opens. The bulk stays on the volume or the shared filesystem.

## Troubleshooting

| Symptom | Likely cause | What to do |
|---------|-------------|-----------|
| A Job's pod sits `Pending`, then the Job ends `DeadlineExceeded` | No free GPU of the kind you selected, a node selector or toleration that matches no node, or the namespace quota is used up | While it is `Pending`: `kubectl --context <ctx> -n <ns> describe pod <pod-name>`; the Events section names the reason. Fix the selector or tolerations in `TOOLS.md`, wait for quota, or ask for fewer GPUs at once. |
| `OOMKilled` | The container hit its memory limit (host memory, not GPU memory) | Raise the memory limit; alignment-based and diffusion-based predictors need more. GPU out-of-memory shows in the log instead: smaller batches, or one predictor per job. |
| `ErrImagePull` / `ImagePullBackOff` | Wrong image name, or the namespace cannot read your registry | Check the image line in `TOOLS.md`; add `imagePullSecrets` to both manifests. |
| The utility pod is `Pending` with a volume event | The storage class does not support the access mode asked for, or a ReadWriteOnce volume is attached to another node | See "A work volume" above; `describe pvc` and `describe pod` say which. |
| Permission denied writing under `/work` | The image runs as a non-root user and the fresh volume is root-owned | Set `securityContext.fsGroup` to the image user's group in both manifests. |
| A tool cannot find its weights, or tries to download them | The weights path inside the image differs from what `TOOLS.md` says, or the tool defaults to fetching | Fix `TOOLS.md`; point the tool's weights or cache variable at the baked-in path. |
| Slurm job ends `TIMEOUT` or `OUT_OF_MEMORY` in `sacct` | `--time` or `--mem` too small for the batch | Time one small batch first and size from that; split the batch. |
