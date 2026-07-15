# Vast Cloud Sync — offload without GPU on

**Docs:** [Data movement](https://docs.vast.ai/guides/instances/storage/data-movement),
[Cloud Sync](https://docs.vast.ai/guides/instances/storage/cloud-sync),
CLI `vastai cloud copy`.

## What exists

Vast **Cloud Sync / Cloud Copy** moves data between instance **local disk** and
cloud providers **even when the instance is stopped** (GPU off; storage still
billed until destroy).

Supported: **Backblaze B2, Amazon S3, Google Drive, Dropbox** (Docker
instances; not KVM/VM).

CLI examples (from docs):

```bash
vastai cloud copy ...          # after account-level cloud integration
vastai copy s3.101:/data/ C.INSTANCE:/workspace/
vastai copy C.INSTANCE:/path/ drive:/folder/
vastai execute INSTANCE_ID 'ls -l'   # browse stopped instance
```

GUI: Account → Connect Backblaze/S3/… → instance **Cloud Copy** button.

## Fit for SCP

| Path | Pros | Cons |
|------|------|------|
| **Current** `sim_download` rsync while running | Simple; in runner today | Needs instance **up**; slow multi-GB; burns GPU idle |
| **rclone on instance → scpsfa:B2** while running | Already used; same B2 remote | Needs instance **up** |
| **Vast Cloud Sync → B2** after **stop** | GPU $ off during multi-GB transfer; API/CLI | Account integration; Docker-only; still pay disk until destroy |
| Pull only **tracks + volview webp** | Tiny; science SoT for force/orbit | No re-slice later |

**Recommendation:** **(B) document + thin runner hook later**, not a big
rewrite now.

1. Connect **Backblaze** on Vast account to the same `scpsfa` bucket (dedicated
   app key; least privilege).
2. After campaign: `stop` instance → `vastai cloud copy` SFAs → `scpsfa:…/v75/` →
   verify → **destroy**.
3. Optional runner tools later:
   - `sim_cloud_sync(remote_path, cloud=b2, dest=…)` wrapping `vastai cloud copy`
   - `sim_remote_snapshot(view, frame)` → `volview -snapshot` on instance →
     download webp only

## Does the SFA have value?

| Need | Full SFA? |
|------|-----------|
| Force/orbit **D(t), pass/fail** | **No** — track+diag enough |
| Re-run `mf_pair_track` / new metrics | Yes, or re-sim cheaper sometimes |
| **Contact morphology** (G3 t~100) | **Snapshots** better than 23G pull |
| Paper figures / volume render | Remote volview frames |
| Checkpoint resume | Seed or mid-SFA if format supports |

**Verdict:** Keep SFAs only if you will re-analyze fields. For F11, **tracks are
SoT**. If archiving: Cloud Sync one copy to B2, then destroy; **do not** pull
23G to local disk by default. Prefer **in-place renders** of frames near
D_min for G3.

## Config steps (manual, first time)

1. B2 app key with R/W on `scpsfa` (or a vast-only prefix).
2. Vast account → Connect Backblaze (keyId + applicationKey). Connection name **`scp`** already exists.
3. Stop `v75mf` when GPU idle.
4. Cloud Copy `/root/mf_headon_D20.sfa` (and force if wanted) → B2 prefix `v75/`.
5. `rclone ls scpsfa:scpsfa/v75/` verify from home.
6. Destroy instance.

## 2026-07-13 test (instance 44719063)

| Method | Result |
|--------|--------|
| `vastai show connections` | **401** — stored API key is restricted (`run_instance`), lacks `api.user.cloud_integrations` |
| `vastai cloud copy … --connection scp` | Connection name accepted (not “No Cloud Operation”); API returns **500 server_error** |
| **On-instance `rclone → scpsfa:`** (same B2) | **Works** — diags/tracks/renders uploaded to `scpsfa:scpsfa/v75/mf/`; headon SFA copy started |

**Unblock Vast Cloud Copy:** regenerate or store a primary/full API key (or key with cloud_integrations + commands/rclone), then `vastai show connections` to get numeric ID, retest `cloud copy` after **stop**. Until then, **instance rclone → scpsfa** is the reliable path.
'''