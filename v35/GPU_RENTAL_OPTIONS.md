# GPU Rental Options for SCP Field Simulations

Comparison of two providers for running CUDA-based soliton field simulations.
Requirement: V100 32GB (or better), SSH access, Ubuntu + CUDA, pay-per-hour.

---

## Provider 1: OVHcloud

**Website**: https://us.ovhcloud.com/public-cloud/gpu/

### V100S Instance Lineup

| Instance   | GPUs         | GPU VRAM | vCPUs | RAM    | Storage   | Price/hr |
|------------|-------------|----------|-------|--------|-----------|----------|
| t2-le-45   | 1x V100S    | 32 GB    | 15    | 45 GB  | 300 GB SSD | $0.88   |
| t2-le-90   | 2x V100S    | 64 GB    | 30    | 90 GB  | 500 GB SSD | $1.76   |
| t2-le-180  | 4x V100S    | 128 GB   | 60    | 180 GB | 500 GB SSD | $3.53   |

Premium variants (t2-45, t2-90, t2-180) also exist at higher prices with enhanced storage.

The V100S has 5,120 CUDA cores, 640 Tensor cores, 32 GB HBM2, and 1,134 GB/s memory bandwidth.

### Other GPU Options on OVHcloud

| GPU     | VRAM  | Starting Price/hr |
|---------|-------|-------------------|
| L4      | 24 GB | $1.00             |
| L40S    | 48 GB | $1.80             |
| V100S   | 32 GB | $0.88             |

A100 and H100 are mentioned but listed only in GRA11 (EU region); pricing not publicly listed for US.

### Regions

- V100S: GRA7, GRA9, GRA11 (EU), BHS5 (Canada), plus Vint Hill VA and Hillsboro OR (US)
- Newer GPUs (A100, H100, L4, L40S): GRA11 only

### OS Images

- All standard OVHcloud Linux images work (Ubuntu, Debian, CentOS, etc.)
- Ubuntu 22.04 and 24.04 available
- CUDA is NOT pre-installed -- you install NVIDIA drivers + CUDA toolkit after boot
- Windows requires special UEFI images (not relevant for us)

### Billing

- Hourly billing, no commitment required
- Data transfer: unlimited (included)
- Can delete instance at any time to stop charges

### Pros

- Fixed, predictable pricing ($0.88/hr for 1x V100S)
- Unlimited bandwidth included -- good for streaming large SFA files
- Managed infrastructure, high reliability
- Standard OpenStack API for automation
- US data centers available (low latency)

### Cons

- CUDA not pre-installed (15-20 min setup on first boot)
- Must use OVHcloud web console or OpenStack CLI to manage
- More expensive than marketplace providers like Vast.ai
- V100S only in older regions; newer GPUs limited to EU

---

## Provider 2: Vast.ai

**Website**: https://vast.ai/
**Pricing**: https://vast.ai/pricing
**CLI docs**: https://docs.vast.ai/cli/get-started

### V100 32GB Pricing (Marketplace -- prices fluctuate)

| GPU        | On-Demand    | Interruptible      |
|------------|-------------|---------------------|
| V100 16GB  | $0.11-0.33  | ~50% less           |
| V100 32GB  | $0.17-0.23  | ~$0.08-0.12         |
| A100 40GB  | $0.29-0.79  | ~$0.15-0.40         |
| A100 80GB  | $0.94-1.40  | ~$0.50-0.70         |

### Pricing Tiers

| Tier           | Description                                    | Discount     |
|----------------|------------------------------------------------|-------------|
| On-demand      | Guaranteed uptime, best for production         | Base price   |
| Interruptible  | May be reclaimed, best for batch jobs          | ~50%+ off    |
| Reserved       | 1/3/6 month commitment, guaranteed capacity    | Up to 50% off|

### Cost Components

- **Compute**: Charged per second while running
- **Storage**: Charged continuously even when stopped (until instance destroyed)
- **Bandwidth**: Per-byte charges for upload/download

### How It Works

Vast.ai is a GPU marketplace: individual hosts list their machines and set prices.
You pick an offer, launch a Docker container on it, and get SSH access.

### Docker Images

- Use any Docker image (PyTorch, TensorFlow, CUDA base, or custom)
- Pre-built templates: `vastai/tensorflow`, `pytorch/pytorch`, `nvidia/cuda:12.2.0-devel-ubuntu22.04`
- For our use case: `nvidia/cuda:12.2.0-devel-ubuntu22.04` gives us CUDA toolkit + Ubuntu

### Pros

- 2-5x cheaper than OVHcloud ($0.17-0.23/hr vs $0.88/hr for V100 32GB)
- Per-second billing
- Docker-based -- can use pre-built CUDA images (zero setup time)
- Large inventory (10,000+ GPUs, 68+ GPU types)
- Full CLI and API for automation
- V100 32GB specifically available

### Cons

- Marketplace pricing is volatile
- Storage charges even when stopped
- Bandwidth costs extra (matters for large SFA file downloads)
- Variable reliability (host-dependent)
- Interruptible instances may be reclaimed mid-job
- Less predictable than traditional cloud

---

## Head-to-Head Comparison

| Feature                | OVHcloud t2-le-45        | Vast.ai V100 32GB       |
|------------------------|--------------------------|--------------------------|
| GPU                    | 1x V100S 32GB            | 1x V100 32GB            |
| Price/hr               | $0.88 (fixed)            | $0.17-0.23 (market)     |
| Monthly (24/7)         | ~$634                    | ~$122-$166              |
| CUDA pre-installed     | No (manual install)      | Yes (via Docker image)  |
| SSH access             | Yes (standard)           | Yes (port-forwarded)    |
| Bandwidth              | Unlimited (included)     | Metered (extra cost)    |
| Reliability            | High (managed cloud)     | Variable (marketplace)  |
| Setup time             | ~20 min (first boot)     | ~5 min (Docker pull)    |
| Billing granularity    | Hourly                   | Per-second              |
| API/CLI                | OpenStack CLI             | vastai CLI (pip)        |
| Commitment             | None                     | None (on-demand)        |

**Recommendation**: For development/iteration, use Vast.ai (4x cheaper, faster setup).
For long production runs where reliability matters, OVHcloud is safer.

---

## Setup Guide: OVHcloud

### 1. Prerequisites

```bash
# Install OpenStack CLI
pip install python-openstackclient

# Download your OpenStack RC file from OVHcloud control panel:
# Control Panel > Public Cloud > Project > Users & Roles > Download RC file
source openrc.sh
```

Create an SSH key in the OVHcloud control panel (Public Cloud > Key Pairs)
or upload your existing public key.

### 2. Launch an Instance

**Via Web Console:**
1. Go to https://us.ovhcloud.com/ > Public Cloud > your project
2. Compute > Instances > Create an instance
3. Select model: t2-le-45 (1x V100S, 45GB RAM)
4. Select region: US-EAST-VA-1 or US-WEST-OR-1
5. Select image: Ubuntu 22.04
6. Select your SSH key
7. Launch

**Via OpenStack CLI:**
```bash
# List available flavors
openstack flavor list | grep t2

# List available images
openstack image list | grep -i ubuntu

# Create the instance
openstack server create \
  --flavor t2-le-45 \
  --image "Ubuntu 22.04" \
  --key-name my-ssh-key \
  --network Ext-Net \
  scp-gpu-01

# Get the IP address
openstack server show scp-gpu-01 -f value -c addresses
```

### 3. SSH In and Install CUDA

```bash
# SSH into the instance
ssh ubuntu@<INSTANCE_IP>

# Install NVIDIA drivers + CUDA toolkit
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get install -y cuda-toolkit-12-2 nvidia-driver-535
sudo reboot

# After reboot, verify
ssh ubuntu@<INSTANCE_IP>
nvidia-smi
nvcc --version
```

### 4. Upload Simulation Binary and Run

```bash
# Upload binary
scp ./scp_simulation ubuntu@<INSTANCE_IP>:~/

# Upload any input files
scp ./input_config.dat ubuntu@<INSTANCE_IP>:~/

# SSH in and run
ssh ubuntu@<INSTANCE_IP>
chmod +x ~/scp_simulation
./scp_simulation --config input_config.dat --output /tmp/results/
```

### 5. Download Results

```bash
# Download result files
scp -r ubuntu@<INSTANCE_IP>:/tmp/results/ ./local_results/

# Or use rsync for large transfers (resumable)
rsync -avz --progress ubuntu@<INSTANCE_IP>:/tmp/results/ ./local_results/
```

### 6. Tear Down

```bash
# Via CLI
openstack server delete scp-gpu-01

# Or via web console: Instances > select > Delete
```

---

## Setup Guide: Vast.ai

### 1. Install CLI and Set API Key

```bash
# Install the CLI
pip install vastai

# Get your API key from https://cloud.vast.ai/cli/
vastai set api-key YOUR_API_KEY_HERE
```

### 2. Search for V100 32GB Instances

```bash
# Search for V100 32GB offers with good reliability
vastai search offers 'gpu_name=V100 gpu_ram>=32 reliability>0.95 num_gpus=1' \
  -o 'dph_total+'

# More specific: also require direct SSH and decent upload speed
vastai search offers \
  'gpu_name=V100 gpu_ram>=32 reliability>0.95 num_gpus=1 direct_port_count>0 inet_up>100' \
  -o 'dph_total+'

# Search for A100 40GB as alternative
vastai search offers 'gpu_name=A100 gpu_ram>=40 reliability>0.95 num_gpus=1' \
  -o 'dph_total+'
```

Output will show offer IDs, prices, specs, and location.

### 3. Launch an Instance

```bash
# Create instance from offer ID with CUDA+Ubuntu image and SSH
# Replace OFFER_ID with the ID from search results
vastai create instance OFFER_ID \
  --image nvidia/cuda:12.2.0-devel-ubuntu22.04 \
  --disk 50 \
  --ssh \
  --direct

# Check instance status
vastai show instances
```

### 4. Connect via SSH

```bash
# Get SSH connection details
vastai ssh-url INSTANCE_ID

# Output will be something like: ssh -p 12345 root@123.45.67.89
# Connect:
ssh -p 12345 root@123.45.67.89

# Make sure your SSH public key is uploaded to your Vast.ai account
# (Account > SSH Keys on the web console)
```

### 5. Upload Simulation Binary and Run

```bash
# Upload using SCP (note the -P flag for port)
scp -P 12345 ./scp_simulation root@123.45.67.89:~/

# Or use vastai copy
vastai copy ./scp_simulation INSTANCE_ID:/root/

# SSH in and run
ssh -p 12345 root@123.45.67.89
chmod +x ~/scp_simulation
nvidia-smi   # verify GPU is visible
nvcc --version  # verify CUDA toolkit
./scp_simulation --config input_config.dat --output /tmp/results/
```

### 6. Download Results

```bash
# Download via SCP
scp -P 12345 -r root@123.45.67.89:/tmp/results/ ./local_results/

# Or use vastai copy
vastai copy INSTANCE_ID:/tmp/results/ ./local_results/

# For large files, rsync with port specification
rsync -avz --progress -e "ssh -p 12345" root@123.45.67.89:/tmp/results/ ./local_results/
```

### 7. Stop or Destroy

```bash
# Stop (pauses compute charges, storage still billed)
vastai stop instance INSTANCE_ID

# Restart later
vastai start instance INSTANCE_ID

# Destroy (stops all charges, data lost)
vastai destroy instance INSTANCE_ID
```

---

## Cost Estimate for Typical Simulation Run

Assuming a 24-hour simulation producing ~20 GB of SFA output:

| Cost Component         | OVHcloud t2-le-45 | Vast.ai V100 32GB (on-demand) |
|------------------------|--------------------|-------------------------------|
| Compute (24 hrs)       | $21.12             | $4.08-5.52                    |
| Storage (50 GB)        | included           | ~$0.10-0.50                   |
| Bandwidth (20 GB out)  | included           | ~$0.40-2.00                   |
| **Total**              | **$21.12**         | **$4.58-8.02**                |

For a week-long campaign (168 hrs):

| Cost Component         | OVHcloud           | Vast.ai                       |
|------------------------|--------------------|-------------------------------|
| Compute                | $147.84            | $28.56-38.64                  |
| Storage + bandwidth    | included           | ~$2-5                         |
| **Total**              | **$147.84**        | **$30.56-43.64**              |

---

## Quick Decision Matrix

| If you need...                        | Use          |
|---------------------------------------|-------------|
| Cheapest option                       | Vast.ai     |
| Guaranteed uptime for multi-day run   | OVHcloud    |
| Fastest setup (zero CUDA install)     | Vast.ai     |
| Large data downloads without fees     | OVHcloud    |
| Per-second billing for short tests    | Vast.ai     |
| US-based with low latency             | Either      |
| A100 40GB upgrade path                | Vast.ai     |

---

## Bandwidth Consideration for SFA Streaming

Requirement: ~10 MB/s sustained for real-time frame streaming.

- **OVHcloud**: Up to 25 Gbps network, unlimited bandwidth included. No issue.
- **Vast.ai**: Host-dependent. Filter with `inet_up>100` (100 Mbps) in search.
  Many hosts have 1 Gbps+, but verify before renting. Bandwidth is metered.

For 1 frame/sec at 20 MB compressed per frame:
- Need ~160 Mbps sustained
- OVHcloud: easily handled
- Vast.ai: filter for `inet_up>200` to ensure headroom; expect ~$0.02/GB bandwidth cost
  (~$1.40/hr at 20 MB/s, which would ADD significantly to compute cost)

**If streaming is critical**, OVHcloud's unlimited bandwidth is a major advantage.
For batch jobs (run, then download), Vast.ai is cheaper overall.

---

*Last updated: 2026-03-20*
*Sources: us.ovhcloud.com, vast.ai, docs.vast.ai, getdeploying.com, vpsbenchmarks.com*
