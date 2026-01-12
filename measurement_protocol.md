# Measurement Protocol (MVP): Optical Cantor Barrier (free-space)

## Goal
Measure the **normalized transmission spectrum** \(T(\lambda)\) of the Cantor sample and compare it to a **matched slab control**.

## Equipment
- **Light source**: tunable laser **or** broadband source + spectrometer
- **Detector**: power meter or spectrometer detector
- **Optics**: collimation + aperture/iris to ensure stable beam
- **Mount**: kinematic mount for repeatable normal incidence

## Data to capture
For each wavelength \(\lambda\) (or each spectrometer bin):
- `P_ref(λ)`: power without sample (or through blank substrate)
- `P_sample(λ)`: power through Cantor pattern
- `P_slab(λ)`: power through slab-control sample (same total high-index length)

## Steps
1. **Warmup & stability**
   - Warm the laser/source to steady-state.
2. **Reference sweep**
   - Record `P_ref(λ)` across your target band (e.g. 0.6–1.6 µm).
3. **Cantor sweep**
   - Insert Cantor sample at normal incidence. Record `P_sample(λ)` using the same alignment.
4. **Slab-control sweep**
   - Replace with slab-control. Record `P_slab(λ)`.
5. **Repeatability**
   - Repeat the Cantor sweep once (re-seat sample) to estimate alignment sensitivity.

## Compute
- `T_cantor(λ) = P_sample(λ) / P_ref(λ)`
- `T_slab(λ)   = P_slab(λ) / P_ref(λ)`

Use `analyze_measurement.py` to compute and print “dip targets” (top N deepest dips).

