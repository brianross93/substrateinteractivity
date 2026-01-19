# Substrate Field Transceiver - Quick Start

## 1. Wire It Up (30 min)

```
Arduino Nano          AD9833 #1        AD9833 #2
------------          ---------        ---------
5V  ───────────────── VCC              VCC
GND ───────────────── GND              GND
D13 ───────────────── SCLK             SCLK     (shared)
D11 ───────────────── SDATA            SDATA    (shared)
D10 ───────────────── FSYNC
D9  ──────────────────────────────────  FSYNC


AD9833 #1 OUT ──┐
                ├──→ TDA7297 LEFT IN
AD9833 #2 OUT ──┘


TDA7297             COILS                    POWER
-------             -----                    -----
L+ ─────────────→ Coil1+ ─→ Coil1- ─→ Coil2+ ─→ Coil2- ─→ L-
VCC ←─────────────────────────────────────────────────── 12V+
GND ←─────────────────────────────────────────────────── 12V-
```

## 2. Build Coils (1 hr)

1. Wind 100 turns of magnet wire around 6" PVC ring
2. Repeat for second coil
3. Mount coils 7.5 cm apart (spacing = radius), parallel
4. Wire in series (fields add)

## 3. Upload Code

1. Open Arduino IDE
2. Open `substrate_transceiver.ino`
3. Select Board: Arduino Nano
4. Select Port: (your COM port)
5. Upload

## 4. Test

1. Open Serial Monitor (115200 baud)
2. Type `h` for help
3. Type `2` for Theta preset (8 Hz)
4. Type `s` for frequency sweep

## 5. Commands

| Key | Action |
|-----|--------|
| 0-6 | Select preset |
| +/- | Adjust frequency |
| s | Sweep 1-100 Hz |
| c | OFM chord |
| r | Status |
| h | Help |

## Presets

| # | Name | Base Hz | Use |
|---|------|---------|-----|
| 0 | OR Timescale | 6.25 | Consciousness collapse cycle |
| 1 | Schumann | 7.83 | Earth resonance |
| 2 | Theta | 8.0 | Default, meditation |
| 3 | Alpha | 10.0 | Relaxed awareness |
| 4 | Gamma | 40.0 | Conscious binding |
| 5 | Bandyopadhyay | 111.0 | Microtubule resonance |

## First Experiment

1. Set preset 0 (OR timescale)
2. Place target in center of coils (plant, RNG, or yourself)
3. Run for 30+ minutes
4. Document observations
5. Compare to control (no field)

---

**Full documentation:** `substrate_transceiver_build_guide.md`
