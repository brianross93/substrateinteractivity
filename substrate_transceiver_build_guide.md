# Substrate Field Transceiver - Build Guide

## Overview

This device generates phase-coherent multi-frequency magnetic fields based on OFM (ordered-factorization metric) harmonic structures. The goal is to create tunable field patterns that may couple to quantum/consciousness substrates as described in the theoretical framework.

**Estimated Build Cost:** $150-400 depending on options
**Build Time:** 1-2 weekends
**Skill Level:** Intermediate electronics (soldering, basic circuits)

---

## Theory of Operation

The device generates multiple simultaneous frequencies derived from a master clock, ensuring phase coherence. Frequency ratios follow number-theoretic patterns (OFM structure) that may correspond to microtubule resonances and neural oscillation harmonics.

**Key frequencies:**
- ~6 Hz (OR collapse timescale, 1/160ms)
- ~8 Hz (theta rhythm, Schumann resonance range)
- ~40 Hz (gamma, consciousness correlate)
- Higher harmonics following prime-ratio structure

All frequencies are phase-locked (derived from same clock source).

---

## Bill of Materials

### Core Components (Required)

| Item | Description | Qty | Est. Cost | Source |
|------|-------------|-----|-----------|--------|
| Arduino Nano | Microcontroller for signal generation | 1 | $5-20 | Amazon, AliExpress |
| AD9833 DDS Module | Precision frequency synthesizer (0-12.5 MHz) | 2-4 | $3-8 each | AliExpress, eBay |
| TDA7297 Audio Amplifier | 15W dual channel amp module | 1 | $5-10 | Amazon, AliExpress |
| 12V 2A Power Supply | Powers amplifier | 1 | $8-15 | Amazon |
| 5V USB Power Supply | Powers Arduino | 1 | $5 | (any USB adapter) |
| Magnet Wire 22-26 AWG | For Helmholtz coils (~100m) | 1 roll | $10-20 | Amazon, hardware store |
| PVC Pipe 6" diameter | Coil forms | 2 | $5-10 | Hardware store |
| 3.5mm Audio Cables | Connections | 3 | $5 | Amazon |
| Breadboard | Prototyping | 1 | $5 | Amazon |
| Jumper Wires | Connections | 1 pack | $5 | Amazon |
| 10kΩ Potentiometers | Level control | 4 | $3 | Amazon |
| 0.1µF Capacitors | Filtering | 10 | $2 | Amazon |

**Subtotal Core: ~$70-120**

### Detection/Measurement (Choose One or More)

| Item | Description | Est. Cost | Notes |
|------|-------------|-----------|-------|
| HMC5883L Magnetometer Module | Measures field strength/direction | $5-10 | Verify field is as expected |
| ADS1115 16-bit ADC | High-precision analog input | $5-10 | For sensitive measurements |
| Hardware RNG (OneRNG or similar) | True random number generator | $50 | Anomaly detection |
| Muse 2 EEG Headband | Brain activity monitoring | $250 | If testing on humans |

### Optional Enhancements

| Item | Description | Est. Cost | Notes |
|------|-------------|-----------|-------|
| OLED Display (SSD1306) | Show current frequencies | $5 | Nice to have |
| Rotary Encoder | Frequency adjustment knob | $3 | Real-time tuning |
| SD Card Module | Log data | $5 | Record experiments |
| Aluminum Project Box | Shielding/housing | $15-30 | Cleaner build |
| BNC Connectors | Professional connections | $10 | Optional upgrade |

---

## Circuit Schematic

```
                                    ┌─────────────────┐
                                    │   HELMHOLTZ     │
                                    │   COIL PAIR     │
                                    │                 │
┌──────────────┐                    │    ┌─────┐      │
│              │                    │    │     │      │
│   ARDUINO    │   ┌──────────┐     │  ┌─┴─┐ ┌─┴─┐    │
│   NANO       │   │ AD9833   │     │  │ L │ │ R │    │
│              │   │ DDS #1   ├──┐  │  │   │ │   │    │
│         D13 ├───┤ SCLK     │  │  │  └─┬─┘ └─┬─┘    │
│         D11 ├───┤ SDATA    │  │  │    │     │      │
│         D10 ├───┤ FSYNC    │  │  └────┴─────┴──────┘
│              │   │     OUT ├──┤           │
│              │   └──────────┘  │           │
│              │                 │    ┌──────┴──────┐
│              │   ┌──────────┐  │    │   TDA7297   │
│              │   │ AD9833   │  ├────┤   AMP       │
│          D9 ├───┤ FSYNC    │  │    │             │
│              │   │ DDS #2   ├──┘    │        L+ ├───→ Coil L+
│              │   │     OUT ├────────┤        L- ├───→ Coil L-
│              │   └──────────┘       │        R+ ├───→ Coil R+
│              │                      │        R- ├───→ Coil R-
│              │   (Add more DDS      │             │
│              │    for more freqs)   │    12V ────┤ VCC
│              │                      │    GND ────┤ GND
│          A0 ├───← Pot (base freq)   └─────────────┘
│          A1 ├───← Pot (mix level)
│              │
│         GND ├───→ Common Ground
│         5V  ├───→ DDS VCC
│         VIN ├───← 7-12V (or USB)
│              │
└──────────────┘

         ┌─────────────────────────────────────┐
         │         MAGNETOMETER (optional)      │
         │                                      │
         │  HMC5883L                            │
         │  VCC ←─── 3.3V/5V                    │
         │  GND ←─── GND                        │
         │  SDA ←─── Arduino A4                 │
         │  SCL ←─── Arduino A5                 │
         └─────────────────────────────────────┘
```

---

## Helmholtz Coil Construction

Helmholtz coils produce a uniform magnetic field in the center region. Two identical coils separated by their radius.

### Specifications

| Parameter | Value | Notes |
|-----------|-------|-------|
| Coil diameter | 15 cm (6 inches) | PVC pipe form |
| Coil separation | 7.5 cm (= radius) | Critical for uniformity |
| Wire gauge | 24 AWG magnet wire | Good balance of resistance/turns |
| Turns per coil | 100 turns | Adjustable |
| Resistance per coil | ~5-10Ω | Depends on wire gauge |
| Max field strength | ~1-5 mT at center | With TDA7297 amp |

### Construction Steps

1. **Cut PVC forms**
   - Cut two 3cm wide rings from 6" PVC pipe
   - Sand edges smooth

2. **Wind coils**
   - Secure wire end with tape
   - Wind 100 turns evenly around PVC ring
   - Keep turns tight and even
   - Secure end, leave 30cm leads
   - Repeat for second coil

3. **Mount coils**
   - Build frame from wood/PVC to hold coils
   - Spacing between coils = 7.5cm (equal to radius)
   - Coils must be parallel and coaxial

4. **Wire coils in series (same polarity)**
   - Connect: Amp L+ → Coil1+ → Coil1- → Coil2+ → Coil2- → Amp L-
   - Field adds constructively in center

```
    SIDE VIEW:
    
    ════════════════  ← Coil 1
          │
          │ 7.5 cm
          │
    ════════════════  ← Coil 2
    
    ←───15 cm───→
    
    Uniform field region in center (~5cm diameter sphere)
```

---

## Software

### Arduino Code (substrate_transceiver.ino)

Save this file and upload to Arduino:

```cpp
/*
 * Substrate Field Transceiver
 * Generates OFM-derived multi-frequency magnetic field patterns
 * 
 * Hardware: Arduino Nano + AD9833 DDS modules + TDA7297 amp + Helmholtz coils
 */

#include <SPI.h>

// Pin definitions
#define DDS1_FSYNC 10
#define DDS2_FSYNC 9
#define DDS3_FSYNC 8  // Optional third DDS
#define DDS4_FSYNC 7  // Optional fourth DDS

#define POT_BASE_FREQ A0
#define POT_MIX A1

// AD9833 registers
#define AD9833_FREQ0  0x4000
#define AD9833_PHASE0 0xC000
#define AD9833_RESET  0x0100
#define AD9833_SINEWAVE 0x2000

// Master clock on AD9833
const float MCLK = 25000000.0;  // 25 MHz

// OFM-derived frequency ratios (relative to base)
// These follow prime-factorization harmonic structure
const float OFM_RATIOS[] = {
  1.0,      // f0: base frequency
  1.5,      // f1: 3/2 ratio (perfect fifth)
  2.0,      // f2: octave
  2.5,      // f3: 5/2
  3.0,      // f4: 3/1 (perfect twelfth)
  4.0,      // f5: 2 octaves
  5.0,      // f6: 5/1
  6.0,      // f7: 6/1 (=2*3)
  8.0       // f8: 3 octaves
};

// Preset frequency sets (base frequencies in Hz)
const float PRESETS[] = {
  6.25,     // OR timescale (160ms = 6.25 Hz)
  7.83,     // Schumann resonance
  8.0,      // Theta rhythm
  10.0,     // Alpha
  40.0,     // Gamma
  111.0,    // Bandyopadhyay's reported resonance
  432.0,    // "Natural" tuning pitch
};
const int NUM_PRESETS = 7;

int currentPreset = 0;
float baseFreq = 8.0;
bool running = true;

// Send 16-bit word to AD9833
void writeAD9833(int csPin, uint16_t data) {
  digitalWrite(csPin, LOW);
  SPI.transfer16(data);
  digitalWrite(csPin, HIGH);
}

// Set frequency on AD9833
void setFrequency(int csPin, float freq) {
  // Calculate frequency register value
  // FreqReg = (freq * 2^28) / MCLK
  uint32_t freqReg = (uint32_t)((freq * 268435456.0) / MCLK);
  
  // Split into two 14-bit words
  uint16_t lsb = (freqReg & 0x3FFF) | AD9833_FREQ0;
  uint16_t msb = ((freqReg >> 14) & 0x3FFF) | AD9833_FREQ0;
  
  // Write control word (reset)
  writeAD9833(csPin, AD9833_RESET);
  
  // Write frequency registers
  writeAD9833(csPin, lsb);
  writeAD9833(csPin, msb);
  
  // Write control word (sine wave output, no reset)
  writeAD9833(csPin, AD9833_SINEWAVE);
}

// Initialize AD9833
void initAD9833(int csPin) {
  pinMode(csPin, OUTPUT);
  digitalWrite(csPin, HIGH);
  writeAD9833(csPin, AD9833_RESET);
  delay(10);
}

void setup() {
  Serial.begin(115200);
  Serial.println("Substrate Field Transceiver v1.0");
  Serial.println("================================");
  
  // Initialize SPI
  SPI.begin();
  SPI.setDataMode(SPI_MODE2);  // AD9833 uses mode 2
  SPI.setBitOrder(MSBFIRST);
  
  // Initialize all DDS modules
  initAD9833(DDS1_FSYNC);
  initAD9833(DDS2_FSYNC);
  // initAD9833(DDS3_FSYNC);  // Uncomment if using
  // initAD9833(DDS4_FSYNC);  // Uncomment if using
  
  // Set initial frequencies
  updateFrequencies();
  
  printHelp();
}

void updateFrequencies() {
  // Read base frequency from pot (scale to 1-100 Hz range)
  int potVal = analogRead(POT_BASE_FREQ);
  float potFreq = map(potVal, 0, 1023, 1, 100);
  
  // Use pot value or preset
  if (potVal > 10) {
    baseFreq = potFreq;
  }
  
  // Set frequencies based on OFM ratios
  float f1 = baseFreq * OFM_RATIOS[0];  // Base
  float f2 = baseFreq * OFM_RATIOS[4];  // 3x (creates beat pattern)
  
  setFrequency(DDS1_FSYNC, f1);
  setFrequency(DDS2_FSYNC, f2);
  
  Serial.print("Frequencies: f1=");
  Serial.print(f1, 2);
  Serial.print(" Hz, f2=");
  Serial.print(f2, 2);
  Serial.println(" Hz");
}

void setPreset(int preset) {
  currentPreset = preset % NUM_PRESETS;
  baseFreq = PRESETS[currentPreset];
  
  Serial.print("Preset ");
  Serial.print(currentPreset);
  Serial.print(": ");
  Serial.print(baseFreq);
  Serial.println(" Hz base");
  
  updateFrequencies();
}

// Set full OFM chord (requires 4 DDS modules)
void setOFMChord(float base) {
  baseFreq = base;
  setFrequency(DDS1_FSYNC, base * 1.0);   // Fundamental
  setFrequency(DDS2_FSYNC, base * 2.0);   // Octave
  // setFrequency(DDS3_FSYNC, base * 3.0);   // Fifth above octave
  // setFrequency(DDS4_FSYNC, base * 5.0);   // Major third above 2 octaves
  
  Serial.print("OFM Chord: ");
  Serial.print(base);
  Serial.println(" Hz fundamental");
}

// Frequency sweep mode
void sweep(float startFreq, float endFreq, float duration_ms) {
  Serial.println("Starting frequency sweep...");
  
  unsigned long startTime = millis();
  unsigned long endTime = startTime + (unsigned long)duration_ms;
  
  while (millis() < endTime) {
    float progress = (float)(millis() - startTime) / duration_ms;
    float freq = startFreq + progress * (endFreq - startFreq);
    
    setFrequency(DDS1_FSYNC, freq);
    setFrequency(DDS2_FSYNC, freq * 3.0);  // Harmonic
    
    delay(50);  // Update rate
    
    if (Serial.available()) {
      char c = Serial.read();
      if (c == 'x' || c == 'X') {
        Serial.println("Sweep cancelled");
        return;
      }
    }
  }
  
  Serial.println("Sweep complete");
}

void printHelp() {
  Serial.println("\nCommands:");
  Serial.println("  0-6  : Select preset");
  Serial.println("  +/-  : Adjust base freq");
  Serial.println("  s    : Start sweep (1-100 Hz)");
  Serial.println("  c    : Set OFM chord");
  Serial.println("  r    : Read current settings");
  Serial.println("  h    : Help");
  Serial.println();
}

void loop() {
  if (Serial.available()) {
    char cmd = Serial.read();
    
    switch (cmd) {
      case '0': case '1': case '2': case '3':
      case '4': case '5': case '6':
        setPreset(cmd - '0');
        break;
      
      case '+':
      case '=':
        baseFreq *= 1.1;
        updateFrequencies();
        break;
      
      case '-':
        baseFreq /= 1.1;
        updateFrequencies();
        break;
      
      case 's':
      case 'S':
        sweep(1.0, 100.0, 60000);  // 1-100 Hz over 60 seconds
        break;
      
      case 'c':
      case 'C':
        setOFMChord(baseFreq);
        break;
      
      case 'r':
      case 'R':
        Serial.print("Current base: ");
        Serial.print(baseFreq);
        Serial.println(" Hz");
        break;
      
      case 'h':
      case 'H':
      case '?':
        printHelp();
        break;
    }
  }
  
  // Periodically check pot for manual tuning
  static unsigned long lastPotRead = 0;
  if (millis() - lastPotRead > 500) {
    int potVal = analogRead(POT_BASE_FREQ);
    if (potVal > 10) {  // Only if pot is turned
      float newFreq = map(potVal, 0, 1023, 1, 100);
      if (abs(newFreq - baseFreq) > 1.0) {
        baseFreq = newFreq;
        updateFrequencies();
      }
    }
    lastPotRead = millis();
  }
}
```

---

## Assembly Instructions

## Step-by-Step Assembly (Detailed)

### Step 0: Prep and Workspace

1. Lay out all parts and verify quantities.
2. Label modules with tape: `DDS1`, `DDS2`, `AMP`, `ARDUINO`.
3. Set up a clean work surface and common ground bus on the breadboard.

### Step 1: Prepare Arduino and Breadboard

1. Insert Arduino Nano into breadboard
2. Connect USB for power (or VIN to 7-12V)
3. Wire common ground rail

### Step 2: Wire AD9833 DDS Modules

For each AD9833 module:

| AD9833 Pin | Arduino Pin |
|------------|-------------|
| VCC | 5V |
| GND | GND |
| SCLK | D13 (SPI clock) |
| SDATA | D11 (SPI MOSI) |
| FSYNC | D10 (DDS #1), D9 (DDS #2), etc. |

**Note:** All DDS modules share SCLK and SDATA, only FSYNC is unique per module.

### Step 3: Wire Audio Amplifier

1. Connect DDS outputs to amp inputs:
   - DDS1 OUT → Amp LEFT IN (via 10kΩ pot for level control)
   - DDS2 OUT → Amp RIGHT IN (via 10kΩ pot)
   
2. Mix both DDS to single channel for single coil pair:
   - Use summing resistors (10kΩ from each DDS to amp input)

3. Power:
   - Amp VCC → 12V supply
   - Amp GND → Common ground

### Step 3.5: Confirm Signal Path (Quick Check)

1. Power only the Arduino and DDS modules (not the amp).
2. Upload `substrate_transceiver.ino`.
3. Open Serial Monitor and set preset `2` (Theta).
4. Use a multimeter (AC) or scope on DDS OUT pins to confirm signal.

### Step 4: Connect Helmholtz Coils

1. Amp LEFT+ → Coil 1 terminal A
2. Coil 1 terminal B → Coil 2 terminal A
3. Coil 2 terminal B → Amp LEFT-

Coils are in series, same winding direction (fields add).

### Step 5: Add Control Pots (Optional)

- 10kΩ pot between 5V and GND, wiper to A0 (base frequency)
- 10kΩ pot between 5V and GND, wiper to A1 (mix level)

### Step 6: Add Magnetometer (Optional)

| HMC5883L Pin | Arduino Pin |
|--------------|-------------|
| VCC | 3.3V or 5V |
| GND | GND |
| SDA | A4 |
| SCL | A5 |

---

## Final Assembly Diagram (Top View)

This is the **finished physical layout** you should aim for (single-station build):

```
                  [ COIL 1 ]           [ COIL 2 ]
               ┌─────────────┐       ┌─────────────┐
               │             │       │             │
               │   (field)   │  <--  │   (field)   │
               │             │       │             │
               └─────────────┘       └─────────────┘
                        │                   │
                        └───── series ──────┘
                              connection

    [Arduino + DDS modules] ----> [TDA7297 amp] ----> [Coils]
            (breadboard)                (module)        (frame)
```

**Field region:** The uniform field is in the **center gap** between the two coils.

---

## Final Assembly Diagram (Wiring Overview)

```
Arduino Nano        AD9833 #1         AD9833 #2
-----------         ---------         ---------
5V  -------------→ VCC              VCC
GND -------------→ GND              GND
D13 -------------→ SCLK             SCLK
D11 -------------→ SDATA            SDATA
D10 -------------→ FSYNC
D9  --------------------------------→ FSYNC

DDS1 OUT ----┐
             ├──→ AMP LEFT IN
DDS2 OUT ----┘

AMP LEFT OUT → Coil1 → Coil2 → AMP LEFT OUT (return)
AMP VCC → 12V+
AMP GND → Common GND
```

---

## Finished Build Checklist

- Arduino and DDS modules mounted on breadboard
- DDS outputs verified with multimeter/scope
- Amplifier powered from 12V supply
- Coils mounted with correct spacing (radius apart)
- Coils wired in series with same winding direction
- Field region clear for sample/target placement

---

## Testing Procedure

### Phase 1: Verify Electronics

1. **Upload code** to Arduino, open Serial Monitor (115200 baud)
2. **Check DDS output** with oscilloscope or multimeter (AC)
   - Should see sine wave at set frequency
   - Frequency should change with commands
3. **Check amplifier** - output should follow input
4. **Check coils** - feel slight warmth when driven (indicates current flow)

### Phase 2: Verify Magnetic Field

1. Place magnetometer in center of Helmholtz coils
2. Run sweep from 1-100 Hz
3. Magnetometer should show oscillating field
4. Field strength should be relatively uniform in center region

### Phase 3: Initial Experiments

**Experiment A: Find Resonance (with RNG)**
1. Place hardware RNG in center of coils
2. Log RNG output continuously
3. Sweep frequencies slowly (1 Hz per minute)
4. Look for statistical deviation from baseline at specific frequencies

**Experiment B: Biological Response**
1. Place plant seedling in center of coils
2. Run specific frequencies for hours/days
3. Compare growth to control (no field)
4. Document which frequencies affect growth

**Experiment C: Subjective (with EEG)**
1. Human subject with EEG headband
2. Head positioned near (not in) coil center
3. Run preset frequencies
4. Record subjective reports + EEG changes
5. Correlate frequency with experience

---

## Frequency Reference Table

| Frequency (Hz) | Ratio to 8 Hz Base | Significance |
|----------------|-------------------|--------------|
| 1.0 | 1/8 | Sub-delta |
| 2.0 | 1/4 | Deep delta |
| 4.0 | 1/2 | Delta |
| 6.25 | 0.78 | OR collapse timescale |
| 7.83 | 0.98 | Schumann resonance |
| 8.0 | 1 | Theta/base |
| 10.0 | 1.25 | Alpha |
| 12.0 | 1.5 | SMR |
| 16.0 | 2 | Beta low |
| 24.0 | 3 | Beta mid |
| 40.0 | 5 | Gamma |
| 64.0 | 8 | High gamma |
| 111.0 | 13.9 | Bandyopadhyay resonance |

**OFM Chord (based on 8 Hz):** 8, 16, 24, 40 Hz (ratios 1:2:3:5)

---

## Safety Notes

1. **Magnetic fields:** Field strengths in this device are ~1-5 mT, comparable to MRI (1.5-3 T) is 1000x stronger. Generally safe, but:
   - Keep away from pacemakers, implants
   - Remove magnetic media (credit cards, hard drives)
   - Pregnant women should not participate in experiments

2. **Electrical:** 
   - 12V is low voltage but can still cause burns with shorts
   - Use fused power supply
   - Don't operate near water

3. **Coils:**
   - Will get warm during extended operation
   - Don't cover or insulate (needs airflow)
   - Use appropriate wire gauge for expected current

4. **Human experiments:**
   - Always have second person present
   - Start with low amplitudes
   - Document any adverse effects
   - Stop immediately if discomfort occurs

---

## Expansion Ideas

1. **Add more DDS modules** for richer harmonic content
2. **Add second coil pair** oriented orthogonally (3D field)
3. **Implement modulation** (AM, FM, pulse) for information encoding
4. **Add feedback loop** - magnetometer → Arduino → adjust output
5. **Interface with PC** for data logging and analysis
6. **Synchronize with external signals** (heartbeat, breath, EEG)

---

## Troubleshooting

| Problem | Possible Cause | Solution |
|---------|---------------|----------|
| No output from DDS | Bad wiring, wrong SPI mode | Check connections, verify code |
| Weak magnetic field | Low amplifier gain, coil resistance | Increase amp volume, check coil continuity |
| Distorted waveform | Amp clipping, overdriven | Reduce input level |
| Overheating coils | Too much current, wire too thin | Use thicker wire, reduce amplitude |
| Arduino resets randomly | Power supply inadequate | Use separate supply for amp |
| Frequency inaccurate | Wrong MCLK value in code | Measure actual AD9833 clock |

---

## Resources

- AD9833 Datasheet: https://www.analog.com/media/en/technical-documentation/data-sheets/AD9833.pdf
- Helmholtz Coil Calculator: https://www.accelinstruments.com/Encyclopaedia/Helmholtz-Coils.html
- Bandyopadhyay Lab: https://sites.google.com/view/anirbanbandyopadhyay/
- Original Theory Document: See free_will_investigation_synthesis.md

---

*Version 1.0 - Built for substrate exploration research*
