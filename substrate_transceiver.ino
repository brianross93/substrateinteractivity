/*
 * Substrate Field Transceiver v1.0
 * 
 * Generates OFM-derived multi-frequency magnetic field patterns
 * for substrate coupling experiments.
 * 
 * Hardware: Arduino Nano + AD9833 DDS modules + TDA7297 amp + Helmholtz coils
 * 
 * Theory: Phase-coherent multi-frequency fields may couple to 
 * quantum/consciousness substrates via resonance patterns following
 * number-theoretic (OFM) harmonic structure.
 * 
 * See: substrate_transceiver_build_guide.md for full documentation
 */

#include <SPI.h>

// ============================================================================
// PIN DEFINITIONS
// ============================================================================

#define DDS1_FSYNC 10   // First DDS chip select
#define DDS2_FSYNC 9    // Second DDS chip select
#define DDS3_FSYNC 8    // Third DDS (optional)
#define DDS4_FSYNC 7    // Fourth DDS (optional)

#define POT_BASE_FREQ A0  // Potentiometer for base frequency
#define POT_MIX A1        // Potentiometer for mix level

// ============================================================================
// AD9833 CONSTANTS
// ============================================================================

#define AD9833_FREQ0     0x4000
#define AD9833_PHASE0    0xC000
#define AD9833_RESET     0x0100
#define AD9833_SINEWAVE  0x2000
#define AD9833_TRIANGLE  0x2002
#define AD9833_SQUARE    0x2028

const float MCLK = 25000000.0;  // 25 MHz master clock on AD9833

// ============================================================================
// OFM-DERIVED FREQUENCY RATIOS
// ============================================================================
// These ratios follow prime-factorization harmonic structure
// Based on ordered factorization metric principles

const float OFM_RATIOS[] = {
  1.0,      // f0: base frequency (fundamental)
  1.5,      // f1: 3/2 ratio (perfect fifth)
  2.0,      // f2: octave (2/1)
  2.5,      // f3: 5/2 
  3.0,      // f4: perfect twelfth (3/1)
  4.0,      // f5: two octaves (4/1)
  5.0,      // f6: 5/1
  6.0,      // f7: 6/1 = 2*3
  8.0,      // f8: three octaves (8/1)
  9.0,      // f9: 9/1 = 3^2
  12.0,     // f10: 12/1 = 2^2 * 3
  15.0      // f11: 15/1 = 3 * 5
};
const int NUM_RATIOS = 12;

// ============================================================================
// PRESET FREQUENCIES
// ============================================================================
// Key frequencies derived from theory and empirical findings

struct Preset {
  const char* name;
  float baseHz;
  float ratio1;  // Second frequency as ratio
  float ratio2;  // Third frequency as ratio (if using 3+ DDS)
};

const Preset PRESETS[] = {
  {"OR Timescale",    6.25,  3.0,  5.0},   // 160ms collapse cycle
  {"Schumann",        7.83,  2.0,  4.0},   // Earth-ionosphere resonance
  {"Theta",           8.0,   3.0,  5.0},   // Theta rhythm
  {"Alpha",           10.0,  2.0,  4.0},   // Alpha rhythm
  {"Gamma",           40.0,  1.5,  2.0},   // Consciousness correlate
  {"Bandyopadhyay",   111.0, 2.0,  3.0},   // Microtubule resonance
  {"Custom",          8.0,   3.0,  5.0}    // User-adjustable
};
const int NUM_PRESETS = 7;

// ============================================================================
// STATE VARIABLES
// ============================================================================

int currentPreset = 0;
float baseFreq = 8.0;
float ratio1 = 3.0;
float ratio2 = 5.0;
bool running = true;
int waveform = AD9833_SINEWAVE;

// ============================================================================
// AD9833 FUNCTIONS
// ============================================================================

void writeAD9833(int csPin, uint16_t data) {
  digitalWrite(csPin, LOW);
  delayMicroseconds(1);
  SPI.transfer16(data);
  delayMicroseconds(1);
  digitalWrite(csPin, HIGH);
}

void setFrequency(int csPin, float freq) {
  if (freq < 0.1) freq = 0.1;
  if (freq > 12500000) freq = 12500000;
  
  // Calculate frequency register: FreqReg = (freq * 2^28) / MCLK
  uint32_t freqReg = (uint32_t)((freq * 268435456.0) / MCLK);
  
  // Split into two 14-bit words
  uint16_t lsb = (freqReg & 0x3FFF) | AD9833_FREQ0;
  uint16_t msb = ((freqReg >> 14) & 0x3FFF) | AD9833_FREQ0;
  
  // Reset, write frequency, then set output mode
  writeAD9833(csPin, AD9833_RESET);
  writeAD9833(csPin, lsb);
  writeAD9833(csPin, msb);
  writeAD9833(csPin, waveform);
}

void initAD9833(int csPin) {
  pinMode(csPin, OUTPUT);
  digitalWrite(csPin, HIGH);
  writeAD9833(csPin, AD9833_RESET);
  delay(10);
}

void stopOutput(int csPin) {
  writeAD9833(csPin, AD9833_RESET);
}

// ============================================================================
// FREQUENCY CONTROL
// ============================================================================

void updateAllFrequencies() {
  float f1 = baseFreq;
  float f2 = baseFreq * ratio1;
  float f3 = baseFreq * ratio2;
  
  setFrequency(DDS1_FSYNC, f1);
  setFrequency(DDS2_FSYNC, f2);
  // Uncomment if using more DDS modules:
  // setFrequency(DDS3_FSYNC, f3);
  // setFrequency(DDS4_FSYNC, baseFreq * 8.0);
  
  Serial.print(F("Frequencies: "));
  Serial.print(f1, 2);
  Serial.print(F(" + "));
  Serial.print(f2, 2);
  Serial.print(F(" + "));
  Serial.print(f3, 2);
  Serial.println(F(" Hz"));
}

void setPreset(int index) {
  if (index < 0 || index >= NUM_PRESETS) return;
  
  currentPreset = index;
  baseFreq = PRESETS[index].baseHz;
  ratio1 = PRESETS[index].ratio1;
  ratio2 = PRESETS[index].ratio2;
  
  Serial.print(F("Preset "));
  Serial.print(index);
  Serial.print(F(": "));
  Serial.println(PRESETS[index].name);
  
  updateAllFrequencies();
}

// ============================================================================
// OFM CHORD GENERATOR
// ============================================================================

void setOFMChord(float base, int numHarmonics) {
  Serial.print(F("OFM Chord @ "));
  Serial.print(base);
  Serial.print(F(" Hz: "));
  
  // Set first DDS to fundamental
  setFrequency(DDS1_FSYNC, base);
  Serial.print(base);
  
  // Set subsequent DDS to OFM ratio harmonics
  if (numHarmonics >= 2) {
    setFrequency(DDS2_FSYNC, base * OFM_RATIOS[4]);  // 3x
    Serial.print(F(", "));
    Serial.print(base * OFM_RATIOS[4]);
  }
  
  Serial.println();
}

// ============================================================================
// SWEEP MODES
// ============================================================================

void linearSweep(float startHz, float endHz, unsigned long durationMs) {
  Serial.print(F("Sweeping "));
  Serial.print(startHz);
  Serial.print(F(" -> "));
  Serial.print(endHz);
  Serial.print(F(" Hz over "));
  Serial.print(durationMs / 1000);
  Serial.println(F(" sec (x to cancel)"));
  
  unsigned long startTime = millis();
  unsigned long endTime = startTime + durationMs;
  
  while (millis() < endTime) {
    float progress = (float)(millis() - startTime) / (float)durationMs;
    float currentFreq = startHz + progress * (endHz - startHz);
    
    setFrequency(DDS1_FSYNC, currentFreq);
    setFrequency(DDS2_FSYNC, currentFreq * ratio1);
    
    // Progress indicator every second
    static unsigned long lastPrint = 0;
    if (millis() - lastPrint > 1000) {
      Serial.print(F("  "));
      Serial.print(currentFreq, 2);
      Serial.println(F(" Hz"));
      lastPrint = millis();
    }
    
    delay(50);
    
    // Check for cancel
    if (Serial.available() && Serial.read() == 'x') {
      Serial.println(F("Sweep cancelled"));
      return;
    }
  }
  
  Serial.println(F("Sweep complete"));
  updateAllFrequencies();  // Restore preset
}

void logSweep(float startHz, float endHz, unsigned long durationMs) {
  Serial.print(F("Log sweep "));
  Serial.print(startHz);
  Serial.print(F(" -> "));
  Serial.print(endHz);
  Serial.println(F(" Hz (x to cancel)"));
  
  float logStart = log(startHz);
  float logEnd = log(endHz);
  
  unsigned long startTime = millis();
  unsigned long endTime = startTime + durationMs;
  
  while (millis() < endTime) {
    float progress = (float)(millis() - startTime) / (float)durationMs;
    float currentFreq = exp(logStart + progress * (logEnd - logStart));
    
    setFrequency(DDS1_FSYNC, currentFreq);
    setFrequency(DDS2_FSYNC, currentFreq * ratio1);
    
    delay(50);
    
    if (Serial.available() && Serial.read() == 'x') {
      Serial.println(F("Sweep cancelled"));
      return;
    }
  }
  
  Serial.println(F("Sweep complete"));
  updateAllFrequencies();
}

// ============================================================================
// MODULATION MODES
// ============================================================================

void amplitudeModulate(float carrierHz, float modHz, unsigned long durationMs) {
  Serial.print(F("AM: "));
  Serial.print(carrierHz);
  Serial.print(F(" Hz carrier, "));
  Serial.print(modHz);
  Serial.println(F(" Hz mod"));
  
  setFrequency(DDS1_FSYNC, carrierHz);
  setFrequency(DDS2_FSYNC, modHz);  // This creates beat frequency
  
  delay(durationMs);
}

void pulseMode(float freqHz, unsigned long onMs, unsigned long offMs, int pulses) {
  Serial.print(F("Pulse: "));
  Serial.print(freqHz);
  Serial.print(F(" Hz, "));
  Serial.print(onMs);
  Serial.print(F("/"));
  Serial.print(offMs);
  Serial.println(F(" ms on/off"));
  
  for (int i = 0; i < pulses; i++) {
    setFrequency(DDS1_FSYNC, freqHz);
    setFrequency(DDS2_FSYNC, freqHz * ratio1);
    delay(onMs);
    
    stopOutput(DDS1_FSYNC);
    stopOutput(DDS2_FSYNC);
    delay(offMs);
    
    if (Serial.available() && Serial.read() == 'x') {
      Serial.println(F("Pulse cancelled"));
      return;
    }
  }
  
  Serial.println(F("Pulse complete"));
}

// ============================================================================
// SERIAL INTERFACE
// ============================================================================

void printStatus() {
  Serial.println(F("\n=== STATUS ==="));
  Serial.print(F("Preset: "));
  Serial.print(currentPreset);
  Serial.print(F(" ("));
  Serial.print(PRESETS[currentPreset].name);
  Serial.println(F(")"));
  Serial.print(F("Base: "));
  Serial.print(baseFreq);
  Serial.println(F(" Hz"));
  Serial.print(F("Ratios: 1 : "));
  Serial.print(ratio1);
  Serial.print(F(" : "));
  Serial.println(ratio2);
  Serial.print(F("Waveform: "));
  Serial.println(waveform == AD9833_SINEWAVE ? "Sine" : 
                 waveform == AD9833_TRIANGLE ? "Triangle" : "Square");
  Serial.println();
}

void printHelp() {
  Serial.println(F("\n=== COMMANDS ==="));
  Serial.println(F("Presets:"));
  Serial.println(F("  0: OR Timescale (6.25 Hz)"));
  Serial.println(F("  1: Schumann (7.83 Hz)"));
  Serial.println(F("  2: Theta (8 Hz)"));
  Serial.println(F("  3: Alpha (10 Hz)"));
  Serial.println(F("  4: Gamma (40 Hz)"));
  Serial.println(F("  5: Bandyopadhyay (111 Hz)"));
  Serial.println(F("  6: Custom"));
  Serial.println(F("\nAdjust:"));
  Serial.println(F("  +/-: Base freq ±10%"));
  Serial.println(F("  [/]: Ratio1 ±0.5"));
  Serial.println(F("  f<num>: Set base freq (e.g., f40)"));
  Serial.println(F("\nModes:"));
  Serial.println(F("  s: Linear sweep 1-100 Hz"));
  Serial.println(F("  l: Log sweep 1-1000 Hz"));
  Serial.println(F("  p: Pulse mode (160ms on/off)"));
  Serial.println(F("  c: OFM chord"));
  Serial.println(F("\nWaveform:"));
  Serial.println(F("  w: Cycle waveform (sine/tri/square)"));
  Serial.println(F("\nInfo:"));
  Serial.println(F("  r: Status"));
  Serial.println(F("  h: Help"));
  Serial.println(F("  x: Cancel sweep/pulse"));
  Serial.println();
}

void processCommand(char cmd) {
  switch (cmd) {
    // Presets
    case '0': case '1': case '2': case '3':
    case '4': case '5': case '6':
      setPreset(cmd - '0');
      break;
    
    // Frequency adjustment
    case '+':
    case '=':
      baseFreq *= 1.1;
      updateAllFrequencies();
      break;
    
    case '-':
      baseFreq /= 1.1;
      updateAllFrequencies();
      break;
    
    // Ratio adjustment
    case '[':
      ratio1 = max(1.0, ratio1 - 0.5);
      updateAllFrequencies();
      break;
    
    case ']':
      ratio1 += 0.5;
      updateAllFrequencies();
      break;
    
    // Sweep modes
    case 's':
    case 'S':
      linearSweep(1.0, 100.0, 60000);
      break;
    
    case 'l':
    case 'L':
      logSweep(1.0, 1000.0, 120000);
      break;
    
    // Pulse mode
    case 'p':
    case 'P':
      pulseMode(baseFreq, 160, 160, 100);
      break;
    
    // OFM chord
    case 'c':
    case 'C':
      setOFMChord(baseFreq, 2);
      break;
    
    // Waveform cycle
    case 'w':
    case 'W':
      if (waveform == AD9833_SINEWAVE) waveform = AD9833_TRIANGLE;
      else if (waveform == AD9833_TRIANGLE) waveform = AD9833_SQUARE;
      else waveform = AD9833_SINEWAVE;
      Serial.print(F("Waveform: "));
      Serial.println(waveform == AD9833_SINEWAVE ? "Sine" : 
                     waveform == AD9833_TRIANGLE ? "Triangle" : "Square");
      updateAllFrequencies();
      break;
    
    // Status
    case 'r':
    case 'R':
      printStatus();
      break;
    
    // Help
    case 'h':
    case 'H':
    case '?':
      printHelp();
      break;
    
    // Direct frequency setting (f<num>)
    case 'f':
    case 'F':
      {
        delay(50);  // Wait for rest of number
        String numStr = "";
        while (Serial.available()) {
          char c = Serial.read();
          if (c >= '0' && c <= '9' || c == '.') {
            numStr += c;
          }
          delay(10);
        }
        if (numStr.length() > 0) {
          baseFreq = numStr.toFloat();
          updateAllFrequencies();
        }
      }
      break;
  }
}

// ============================================================================
// SETUP & LOOP
// ============================================================================

void setup() {
  Serial.begin(115200);
  while (!Serial) delay(10);
  
  Serial.println(F("\n================================"));
  Serial.println(F("  SUBSTRATE FIELD TRANSCEIVER"));
  Serial.println(F("  v1.0"));
  Serial.println(F("================================\n"));
  
  // Initialize SPI
  SPI.begin();
  SPI.setDataMode(SPI_MODE2);
  SPI.setBitOrder(MSBFIRST);
  
  // Initialize DDS modules
  Serial.println(F("Initializing DDS modules..."));
  initAD9833(DDS1_FSYNC);
  initAD9833(DDS2_FSYNC);
  // initAD9833(DDS3_FSYNC);
  // initAD9833(DDS4_FSYNC);
  
  // Set default preset
  setPreset(2);  // Theta (8 Hz)
  
  printHelp();
}

void loop() {
  // Process serial commands
  if (Serial.available()) {
    char cmd = Serial.read();
    if (cmd > 32) {  // Ignore whitespace
      processCommand(cmd);
    }
  }
  
  // Check analog pot for manual tuning
  static unsigned long lastPotRead = 0;
  if (millis() - lastPotRead > 250) {
    int potVal = analogRead(POT_BASE_FREQ);
    if (potVal > 20) {  // Only if pot is connected and turned
      float potFreq = 1.0 + (potVal / 1023.0) * 99.0;  // 1-100 Hz
      if (abs(potFreq - baseFreq) > 2.0) {
        baseFreq = potFreq;
        currentPreset = 6;  // Switch to custom
        updateAllFrequencies();
      }
    }
    lastPotRead = millis();
  }
}
