# Fractal Hyperspace Simulation Sandbox

This repo implements a small, local simulator inspired by the `spec.md` outline:

- **Kuramoto synchronization** on a network whose couplings are weighted by an **ordered-factorization metric (OFM)**.
- **1D wave-packet tunneling** using a **split-step Fourier** method, with transmission \(T\) measured to the right of the barrier.

## Setup (Windows PowerShell)

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
py -m pip install --upgrade pip
py -m pip install -r .\requirements.txt
py .\app.py
```

## If the GUI doesn't start (Tcl/Tk / init.tcl error)

Some Windows Python installs are missing Tcl/Tk. In that case, run headless mode:

```powershell
py .\app.py --nogui --mode both --out outputs
```

This prints `final r` and `T` and saves plots to `.\outputs\kuramoto.png` and `.\outputs\tunneling.png`.

## Latest findings (Phase 1 baseline run)

Latest full run command (headless):

```powershell
py .\app.py --nogui --mode all --out outputs --unified-demo --chsh-noise-sweep --cantor-t-sweep --retro-switch --retro-vis-sweep --retro-sweep-step 5 --retro-phi-sweep --retro-phi-step 20
```

Outputs (from the latest run):

- **Kuramoto**: `r_final ≈ 0.997` (near-global synchronization)
- **Rectangular tunneling**: `T ≈ 0.083` (matches the ~0.082 benchmark used in `spec.md`)
- **CHSH (noise=0)**: `S ≈ 2.828` (Tsirelson bound)
- **CHSH noise sweep** (`outputs\chsh_noise_sweep.csv`): `S(noise=3.0) ≈ 2.052` (near-classical)
- **Retro (CPTP switch)**: `max_abs_delta ≈ 0.250`, `spec_dist ≈ 0.165`, `switch_vmax_max ≈ 0.500`
- **Retro visibility vs θ**: plots in `outputs\retro_vis_vs_theta.png`; raw numbers in `outputs\retro_vis_vs_theta.csv`
- **Cantor tunneling sweep + fit**: `d_eff ≈ -0.251`, `R2 ≈ 0.638` (see `outputs\cantor_t_sweep.png` and `outputs\cantor_t_sweep.csv`)

## CHSH / Entanglement (headless)

This repo also includes a simple **CHSH** computation on the Bell state \(|\Phi^+\rangle\) (pure NumPy; no QuTiP).

```powershell
py .\app.py --nogui --mode chsh --out outputs
```

It saves `.\outputs\chsh.png`. Use `--mode all` to run everything.

### Phase 1: CHSH → classical limit (noise sweep)

```powershell
py .\app.py --nogui --mode chsh --out outputs --chsh-noise-sweep --chsh-ofm-n 20
```

This writes `.\outputs\chsh_noise_sweep.png` for noise = `[0, 0.5, 1, 1.5, 2, 3]` and (with current calibration) drives \(S\) down to about ~2.05 at noise=3.

It also writes the raw numbers to `.\outputs\chsh_noise_sweep.csv`.

## Retrocausality / non-commuting projectors (headless)

This is a **toy model** showing that when projectors don’t commute (and/or don’t commute with the Hamiltonian), the **joint probability depends on the order**.

```powershell
py .\app.py --nogui --mode retro --out outputs
```

It saves:

- `.\outputs\retro.png` (order-dependent probabilities + delta + time-series)
- `.\outputs\retro_spectrum.png` (FFT magnitude)

and prints `max_abs_delta`, frequency estimates, and commutator norms.

You can tune parameters in headless mode too:

```powershell
py .\app.py --nogui --mode retro --out outputs --retro-theta 30 --retro-t2 8 --retro-dt 0.02 --retro-strobe 0
```

### Quantum-switch-style (indefinite causal order) variant

```powershell
py .\app.py --nogui --mode retro --out outputs --retro-switch --retro-switch-phase 45 --retro-theta 45
```

This enables a **fully CPTP quantum switch** (no postselection). It writes `.\outputs\retro_switch.png` and prints `switch_vmax_max`, where:

- In **CPTP mode**: `switch_vmax_max` is based on the control coherence, \(V(t_1)=2|\rho_{c,01}|\).
- A **classical mixture of orders** has \(V=0\) (baseline shown on the plot).

If you want the older **postselected overlap** model (legacy), use:

```powershell
py .\app.py --nogui --mode retro --out outputs --retro-switch-post --retro-switch-phase 45 --retro-theta 45
```

### Theta sweep (batch mode)

```powershell
py .\app.py --nogui --mode retro --out outputs --retro-sweep --retro-sweep-step 5
```

This writes `.\outputs\retro_theta_sweep.png`.

## CPTP quantum switch: visibility sweeps + Bloch animation

### Visibility vs theta (primary)

```powershell
py .\app.py --nogui --mode retro --out outputs --retro-switch --retro-vis-sweep --retro-sweep-step 5
```

Writes `.\outputs\retro_vis_vs_theta.png`.

It also writes the raw numbers to `.\outputs\retro_vis_vs_theta.csv`.

### Phase sweep (interference lobe)

```powershell
py .\app.py --nogui --mode retro --out outputs --retro-switch --retro-phi-sweep --retro-phi-step 10 --retro-theta 45
```

Writes `.\outputs\retro_vis_vs_phi.png`.

### Control Bloch animation

```powershell
py .\app.py --nogui --mode retro --out outputs --retro-switch --retro-control-bloch --retro-t2 4 --retro-dt 0.05
```

Writes `.\outputs\retro_control_bloch.gif`.

### Optional: OFM modulation

```powershell
py .\app.py --nogui --mode retro --out outputs --retro-switch --retro-vis-sweep --retro-ofm --retro-ofm-target phi --retro-ofm-n 20 --retro-ofm-strength-deg 30
```

## Kuramoto animation (headless GIF)

```powershell
py .\app.py --nogui --mode kuramoto_anim --out outputs
```

This writes `.\outputs\kuramoto.gif`.

## Sacred-geometry-inspired Kuramoto topologies

The Kuramoto sim supports alternative coupling topologies in code via `KuramotoConfig(topology=...)`:

- `ofm`: OFM-weighted dense coupling (default)
- `torus`: nearest-neighbor coupling on a periodic 2D torus grid
- `metatron`: a simple 13-node “Metatron-like” graph (center + ring + chords)
- `pyramid`: Great Pyramid-inspired internal topology (approximate chambers + shafts)

You can also set topology in headless mode:

```powershell
py .\app.py --nogui --mode kuramoto --out outputs --kuramoto-topology torus --torus-m 6
```

Pyramid example (approximate geometry; distances scaled in meters):

```powershell
py .\app.py --nogui --mode kuramoto --out outputs --kuramoto-topology pyramid --kuramoto-j 1.0 --pyramid-scale-m 35 --pyramid-shaft-weight 1.15
```

Custom pyramid geometry (CSV-driven):

- Nodes CSV columns: `name,x,y,z`
- Edges CSV columns: `u,v,medium,weight`
  - `medium` is optional and uses `--pyramid-medium-weights`
  - `weight` is optional; if omitted, distance-based weight is used

```powershell
py .\app.py --nogui --mode kuramoto --out outputs --kuramoto-topology pyramid `
  --pyramid-nodes-csv .\inputs\pyramid_nodes.csv `
  --pyramid-edges-csv .\inputs\pyramid_edges.csv `
  --pyramid-medium-weights "stone=1.0,air_shaft=1.2,water=1.1,granite=1.05"
```

Second-order (inertia) Kuramoto:

```powershell
py .\app.py --nogui --mode kuramoto --out outputs --kuramoto-topology pyramid `
  --kuramoto-model second --kuramoto-mass 2.0 --kuramoto-damping 0.8
```

Longer run + finer sampling (useful for breathing analysis):

```powershell
py .\app.py --nogui --mode kuramoto --out outputs --kuramoto-topology pyramid `
  --kuramoto-t-end 60 --kuramoto-dt 0.05
```

Shaft-weight sweep (writes CSV + PNG):

```powershell
py .\app.py --nogui --mode kuramoto --out outputs --kuramoto-topology pyramid `
  --pyramid-shaft-sweep "0.8,0.9,1.0,1.1,1.2,1.3,1.4"
```

## Notes

## Phase 1: Cantor tunneling sweep + log-log fit

```powershell
py .\app.py --nogui --mode tunneling --out outputs --barrier cantor --cantor-t-sweep
```

This writes `.\outputs\cantor_t_sweep.png` and prints a log-log fit exponent `d_eff` and `R2`.

It also writes the raw numbers to `.\outputs\cantor_t_sweep.csv`.

For a power-of-two sweep (cleaner spacing), you can run:

```powershell
py .\app.py --nogui --mode tunneling --out outputs --barrier cantor --cantor-t-sweep --cantor-levels "1,2,4,8,16"
```

Note: the CSV includes `min_feature_over_dx` (smallest Cantor feature size divided by grid step). Values < 1 mean that level is **not resolved** on the chosen grid, so any “fractal dimension” fit is not physically meaningful at that level.

## Optional: Unified demo (tie Kuramoto → retro)

```powershell
py .\app.py --nogui --mode all --out outputs --unified-demo
```

This maps `kuramoto r_final` into the retro switch phase via:

- `phi = 45 + 30*(1 - r_final)`

and applies control dephasing proportional to `(1 - r_final)` so the **CPTP switch visibility actually drops** when synchronization is poor.

Low-r example:

```powershell
py .\app.py --nogui --mode all --out outputs --unified-demo --retro-switch --kuramoto-omega-std 2.0 --kuramoto-j 0.2
```

## Phase 1 PDF bundle

After generating plots/CSVs, you can bundle them into one PDF:

```powershell
py .\phase1_report.py --outdir outputs
```

This writes `.\outputs\phase1_report.pdf`.

## Optional: FDTD (Meep) Cantor barrier spectrum (pre-hardware target)

This runs a **1D FDTD** simulation (Meep) and outputs a normalized transmission spectrum for a Cantor barrier with index contrast (default n=1.5 vs air).

### Important (Windows)

On Windows, `pip install meep` usually installs an **unrelated** package (name collision), and `conda-forge` often does **not** provide Meep FDTD builds for `win-64`.

The most reliable Windows path is to run Meep in a **Linux container** via Docker.

### Docker (recommended)

Build the Meep image once:

```powershell
docker build -f .\Dockerfile.meep -t qinvestigation-meep .
```

Sanity-check Meep imports:

```powershell
docker run --rm qinvestigation-meep
```

Run the Cantor spectrum job (writes into your local `.\outputs\` folder):

```powershell
docker run --rm -v ${PWD}:/work -w /work qinvestigation-meep micromamba run -n meep python .\meep_cantor.py --levels 4 --resolution 80 --nfreq 200 --out outputs\cantor_transmission.png
```

Run (wideband overlay + slab baseline):

```powershell
docker run --rm -v ${PWD}:/work -w /work qinvestigation-meep micromamba run -n meep python ./meep_cantor.py --wl-min 0.5 --wl-max 2.0 --levels-list 2,3,4,5 --baseline-slab --resolution 40 --nfreq 160 --cell-y 2 --out outputs/cantor_transmission_wideband.png
```

It writes:

- `outputs\cantor_transmission_wideband.png` (overlay plot)
- `outputs\cantor_transmission_wideband.npz` (combined arrays)
- `outputs\cantor_transmission_wideband_L*.npz` (per-level arrays)
- `outputs\cantor_transmission_wideband_slab_match_*.npz` (baseline slab spectrum)

### Text-only logging (no charts)

After running Meep, you can print dip targets and summary stats:

```powershell
.\.venv\Scripts\python.exe .\summarize_meep.py --npz outputs\cantor_transmission_wideband.npz --dip-prominence 0.02 --top-dips 5
```

### Animation (how the FDTD sim works)

This produces a GIF showing the **Ez pulse** propagating and interacting with the Cantor barrier (shaded region).

```powershell
docker run --rm -v ${PWD}:/work -w /work qinvestigation-meep micromamba run -n meep python ./meep_cantor.py --animate --levels 3 --wl-min 0.6 --wl-max 1.6 --no-dispersion --resolution 30 --cell-x 40 --cell-y 2 --pml 2 --barrier-width 4 --src-x -15 --anim-t-end 160 --anim-dt 2.0 --npts 4000 --anim-out outputs/cantor_field_L3.gif
```

### Build MVP artifacts (prediction + fab spec + measurement)

- **Prediction text (no charts)**: `outputs\prediction_L3_disp.txt`
- **Fab spec PDF**: `outputs\fab_spec_cantor_L3.pdf` (generate via `.\.venv\Scripts\python.exe .\fab_spec_cantor_L3.py`)
- **Measurement protocol**: `measurement_protocol.md`
- **Measurement analysis**: `analyze_measurement.py` (normalizes sample/ref, prints top dips)

- The GUI uses **Tkinter** (bundled with standard Python on Windows) and embeds **Matplotlib** plots.
- The OFM here is a practical approximation based on **ordered factorizations** \(F(n)\), used to build distances \(d(i,j)=|F(i)-F(j)|\) and couplings \(J_{ij}=J/(1+d)\), per the spec’s suggested approach.


## Hardware: Substrate Field Transceiver

A physical device for generating OFM-derived multi-frequency magnetic fields. Based on the theory that phase-coherent harmonic field patterns may couple to quantum/consciousness substrates.

### Files

- **`substrate_transceiver_build_guide.md`** - Full build documentation
- **`substrate_transceiver.ino`** - Arduino firmware
- **`SHOPPING_LIST.txt`** - Parts list with links/prices
- **`QUICKSTART.md`** - 1-page getting started guide

### Quick Overview

**Hardware:** Arduino Nano + AD9833 DDS modules + audio amplifier + Helmholtz coils

**Cost:** ~$60-150 depending on options

**Key Features:**
- Phase-coherent multi-frequency generation (all derived from same clock)
- OFM-derived harmonic ratios (prime factorization structure)
- Preset frequencies: OR timescale (6.25 Hz), Schumann (7.83 Hz), Theta (8 Hz), Gamma (40 Hz), Bandyopadhyay resonance (111 Hz)
- Sweep and pulse modes for frequency search
- Serial control interface

**Theory:** If nuclear spin affects consciousness (xenon isotope evidence), and nuclear spin couples to magnetic fields (NMR principle), then tuned magnetic fields should couple to consciousness-relevant processes. Multi-frequency patterns following number-theoretic structure (OFM) may create resonance with microtubule/substrate geometry.
