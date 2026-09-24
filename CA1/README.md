# 📝 Engineering Mathematics — Computer Assignment 1: Fourier Analysis & Signal Processing

> Comprehensive implementation and theoretical exploration of Fourier Series and Continuous Fourier Transforms applied to seismic hazard analysis, automated musical transcription, modulation dynamics, spectral windowing, and piecewise-triangular basis function approximation.

[![MATLAB](https://img.shields.io/badge/MATLAB-R2024a-orange.svg)](#)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](../LICENSE)
[![Status](https://img.shields.io/badge/Status-Complete-success.svg)](#)
[![University](https://img.shields.io/badge/University%20of%20Tehran-ECE-blue.svg)](#)

---

## 👨‍🎓 Student Information

| Field | Value |
|---|---|
| **Student Name** | Amirali Dehghani |
| **Student ID** | `810102443` |
| **Course** | Engineering Mathematics (ECE) |
| **Semester** | Fall 2024 (Fall 1403) |
| **Instructor** | Dr. Mehdi Tale Masouleh |
| **Designers / TAs** | Ali Ariaei, Parisa Mohammadi |
| **Submission Date** | October 2024 (13/07/1403) |

---

## 📋 Assignment Objectives

1. 🎯 **Seismic Signal Spectrum & Structural Vulnerability:** Decompose real-world earthquake acceleration time-series into discrete exponential Fourier series components to identify dominant ground-motion frequencies and evaluate resonant structural failure risks across varying building heights.
2. 🎯 **Automated Musical Note Recognition:** Segment raw acoustic audio signals into uniform musical intervals, extract spectral energy at fundamental musical frequencies, and transcribe notes through harmonic thresholding.
3. 🎯 **Continuous Fourier Transform & Frequency Modulation:** Develop a generalized numerical Fourier transform solver using numerical quadrature, and validate the modulation property through complex exponential frequency shifting.
4. 🎯 **Time-Domain Windowing & Convolution Theorem:** Investigate spectral leakage and aperture effects caused by rectangular window functions, empirically validating the frequency-domain convolution theorem $\mathcal{F}\{f(t) \cdot w(t)\} = \frac{1}{2\pi} [F(\omega) * W(\omega)]$.
5. 🎯 **Piecewise-Linear Signal Modeling via Triangular Bases:** Synthesize arbitrary continuous signals as linear combinations of shifted triangular basis functions, deriving their analytical spectra and validating against general quadrature integration.

---

## 📊 Dataset & Signal Specifications

| Dataset / File | Type | Sampling Rate ($f_s$) | Dimensions / Duration | Description |
|---|---|---|---|---|
| `data.mat` | MAT-file | $1000\text{ Hz}$ | $3000 \times 1$ samples ($3.0\text{ s}$) | Simulated ground acceleration time-series of a seismic event in an earthquake-prone zone. |
| `song.wav` | Audio (WAV) | $44,100\text{ Hz}$ | $564,480$ samples ($\approx 12.8\text{ s}$) | Monophonic acoustic audio track containing sequential notes of fixed duration ($0.4\text{ s}$ per note). |

---

## 🧪 Methods & Detailed Explanations

### Part 1: Fourier Series Analysis

#### 1.1 Seismograph Frequency Extraction & Structural Resonance Analysis
Earthquake ground acceleration generates vibrational forces whose damage potential depends heavily on whether ground excitation frequencies match the natural vibrational frequencies of human structures.

- **Signal Normalization:** The raw seismic acceleration signal $x[n]$ from `data.mat` is normalized into a non-dimensional amplitude range $[0, 1]$ via min-max scaling:
  $$x_{\text{norm}}[n] = \frac{x[n] - \min(x)}{\max(x) - \min(x)}$$
- **Signal Framing:** The recorded time-series is partitioned into three equal 1-second intervals. The central interval ($t \in [1, 2]\text{ s}$, $N = 1000$ samples) corresponds to the primary shock arrival where seismic energy reaches its maximum.
- **Exponential Discrete Fourier Series:** For discrete periodic data of length $N$, the complex Fourier series coefficients $c_k$ are evaluated numerically as:
  $$c_k = \frac{1}{N} \sum_{n=1}^{N} x_{\text{norm}}[n] \cdot e^{-j \frac{2\pi k n}{N}}$$
  With sampling frequency $f_s = 1000\text{ Hz}$ and window length $N = 1000$, the frequency resolution is $\Delta f = \frac{f_s}{N} = 1\text{ Hz}$. Hence, index $k \in \{1, 2, \dots, 10\}$ directly corresponds to harmonic frequencies $f = 1, 2, \dots, 10\text{ Hz}$.
- **Structural Vulnerability Assessment:** When an earthquake transmits energy at frequencies matching a building's fundamental resonance frequency, destructive constructive interference occurs. Natural frequency benchmarks classify structures as:
  - **Low-Rise (1–4 stories):** Natural frequency $f_n \in [3, 10]\text{ Hz}$
  - **Mid-Rise (5–10 stories):** Natural frequency $f_n \in [1, 2]\text{ Hz}$
  - **High-Rise (> 10 stories):** Natural frequency $f_n \in [0.1, 1.0]\text{ Hz}$

#### 1.2 Automated Acoustic Note Recognition via Harmonic Decomposition
Musical tones are acoustic sinusoidal oscillations characterized by distinct fundamental pitch frequencies. A composite musical track can be transcribed by tracking spectral energy within localized time slices.

- **Temporal Framing:** The recording `song.wav` ($f_s = 44,100\text{ Hz}$) is partitioned into uniform non-overlapping frames of duration $T_{\text{note}} = 0.4\text{ s}$, each consisting of $N = \lfloor 0.4 \times f_s \rfloor = 17,640$ discrete samples.
- **Fundamental Target Frequencies:** The reference musical scale frequencies analyzed are:
  - **C (Do):** $523.25\text{ Hz}$
  - **D (Re):** $587.33\text{ Hz}$
  - **E (Mi):** $659.26\text{ Hz}$
  - **F (Fa):** $698.46\text{ Hz}$
  - **G (Sol):** $784.00\text{ Hz}$
- **Discrete Harmonic Index Mapping:** For each target pitch $f_{\text{target}}$, the corresponding discrete Fourier frequency bin index is computed as:
  $$k = \text{round}\left(f_{\text{target}} \cdot \frac{N}{f_s}\right)$$
- **Adaptive Energy Thresholding:** The discrete Fourier coefficient magnitude $|c_k|$ is evaluated for each candidate note in every segment. To prevent false positives caused by inter-note acoustic silence or noise, a note is identified only if its coefficient magnitude exceeds half the mean absolute amplitude of the segment:
  $$|c_k| > 0.5 \times \frac{1}{N} \sum_{n=1}^{N} |x[n]|$$
  The dominant note meeting this threshold is transcribed into the final melody sequence.

---

### Part 2: Continuous Fourier Transforms & Spectral Theorems

#### 2.1 Generalized Continuous Fourier Transform & Modulation Property
The continuous-time Fourier transform provides the frequency-domain representation of non-periodic signals.

- **Numerical Integration Architecture:** The continuous Fourier transform integral:
  $$\mathcal{F}\{f(t)\} = F(\omega) = \int_{-\infty}^{\infty} f(t) e^{-j \omega t} \, dt$$
  is implemented across a continuous frequency grid $\omega \in [-50, 50]\text{ rad/s}$ ($1000$ points) using adaptive Gauss-Kronrod numerical quadrature over an effective compact support window $t \in [-10, 10]\text{ s}$.
- **Base Test Function:** Evaluated on a composite rectangular pulse signal:
  $$f(t) = 3 \cdot \Pi\left(\frac{t}{2}\right) + \Pi(t)$$
  where $\Pi(t)$ is the standard unit rectangular pulse ($1$ for $|t| \le 0.5$, $0$ otherwise).
- **Frequency Modulation / Translation Theorem:** Modulating $f(t)$ by a complex exponential:
  $$g(t) = f(t) \cdot e^{-j 2\pi (3) t} = f(t) \cdot e^{-j 6\pi t}$$
  According to the modulation theorem $\mathcal{F}\{f(t) e^{j \omega_0 t}\} = F(\omega - \omega_0)$, multiplying by $e^{-j 6\pi t}$ introduces a frequency shift of $\omega_0 = -6\pi \approx -18.85\text{ rad/s}$, displacing the entire magnitude spectrum leftward along the $\omega$-axis without altering its spectral envelope shape.

#### 2.2 Time-Domain Windowing & Frequency-Domain Convolution
Observing physical signals over finite observation intervals is equivalent to multiplying the infinite-duration signal by a windowing aperture function.

- **Window Function:** The rectangular windowing operator is defined as:
  $$w(t) = \Pi(t) = \begin{cases} 1 & |t| \le 0.5 \\ 0 & |t| > 0.5 \end{cases}$$
- **Test Signals:**
  1. Squared Sinc Pulse: $f_1(t) = \text{sinc}^2(t) = \left(\frac{\sin(\pi t)}{\pi t}\right)^2$
  2. Double-Sided Exponential Decay: $f_2(t) = e^{-3|t|}$
- **Fourier Duality & Convolution Theorem:** The Fourier transform of the product of two signals in time equals the convolution of their individual transforms in frequency, scaled by $\frac{1}{2\pi}$:
  $$\mathcal{F}\{f(t) \cdot w(t)\} = \frac{1}{2\pi} \left[ F(\omega) * W(\omega) \right] = \frac{1}{2\pi} \int_{-\infty}^{\infty} F(\tau) W(\omega - \tau) \, d\tau$$
- **Numerical Verification:** Two parallel computational pathways are compared:
  1. Direct quadrature integration of the windowed product $f(t) \cdot w(t)$.
  2. Numerical discrete convolution between $F(\omega)$ and $W(\omega)$, scaled by the frequency step $\Delta \omega$ and the theoretical normalization factor $\frac{1}{2\pi}$.
- **Spectral Impact:** The rectangular window introduces sinc-shaped spectral leakage, producing side lobes that broaden sharp spectral features.

#### 2.3 Piecewise-Linear Approximation Using Triangular Basis Expansions
Any continuous compact-support signal can be approximated to high precision through linear combinations of localized triangular basis functions.

- **Symmetric Triangular Basis:** Defined with half-width parameter $T$ and time shift $\tau$:
  $$\Lambda_T(t - \tau) = \begin{cases} 1 - \frac{|t - \tau|}{T} & |t - \tau| < T \\ 0 & \text{otherwise} \end{cases}$$
- **Analytical Fourier Transform of Basis:** The exact continuous Fourier transform of a centered triangular pulse is:
  $$\mathcal{F}\{\Lambda_T(t)\} = T \cdot \text{sinc}^2\left(\frac{\omega T}{2\pi}\right)$$
- **Signal Expansion:** Discretizing the non-zero support of $f(t)$ into $n$ uniform subintervals of width $T = t_{k+1} - t_k$, the signal is approximated as:
  $$\tilde{f}(t) = \sum_{k=1}^{n-1} f(k T) \cdot \Lambda_T(t - k T)$$
- **Analytical Fourier Synthesis:** Applying linearity and the time-shift property $\mathcal{F}\{\Lambda_T(t - \tau)\} = e^{-j\omega \tau} \mathcal{F}\{\Lambda_T(t)\}$, the overall spectrum is synthesized algebraically without numerical integration:
  $$F(\omega) \approx \sum_{k=1}^{n-1} e^{-j\omega k T} \cdot f(k T) \cdot \mathcal{F}\{\Lambda_T(t)\}$$
- **Validation Signals:** Tested against:
  1. Raised Cosine Window: $f(t) = \Pi(t/5) \cdot \left(1 - \cos\left(\frac{2\pi t}{5}\right)\right)$
  2. Multi-Peak Triangular Train: $g(t) = \sum_{k=0}^{3} \Lambda_1(t - 2k)$
  Spectra generated by this triangular expansion method are benchmarked directly against theoretical closed-form derivations and general adaptive quadrature.

---

## 📈 Results & Key Findings

### Quantitative & Qualitative Summary

| Task / Problem | Method / Principle | Input / Parameters | Result / Accuracy | Conclusion |
|---|---|---|---|---|
| **1.1 Earthquake Spectrum** | Exponential Fourier Series | $N=1000$, $f_s=1000\text{ Hz}$, $k \in [1, 10]$ | Peak energy at $2\text{ Hz}$ and $1\text{ Hz}$ | High resonance danger for low/mid-rise buildings (1–10 stories). |
| **1.2 Audio Transcription** | Harmonic Energy Thresholding | $32$ segments, $T=0.4\text{ s}$, $5$ notes | $100\%$ detection accuracy (32/32 notes) | Transcribed Beethoven's "Ode to Joy" melodic motif. |
| **2.1 Continuous FT** | Numerical Quadrature & Modulation | $\omega \in [-50, 50]$, shift $\omega_0 = -6\pi$ | Exact $-18.85\text{ rad/s}$ rigid translation | Confirmed modulation theorem $\mathcal{F}\{f(t) e^{-j\omega_0 t}\} = F(\omega + \omega_0)$. |
| **2.2 Windowing & Conv.** | Time Product vs. Frequency Conv. | $w(t) = \Pi(t)$, $f_1 = \text{sinc}^2, f_2 = e^{-3\|t\|}$ | Direct FT and $\frac{1}{2\pi} (F * W)$ curves overlap | Confirmed $\mathcal{F}\{f \cdot w\} = \frac{1}{2\pi}(F * W)$ with side-lobe leakage. |
| **2.3 Triangular Basis FT** | Piecewise-Linear Decomposition | $T = 0.1005\text{ s}$, $200$ grid points | Exact match between numerical & triangular methods | Enabled fast analytical spectral synthesis without numerical quadrature. |

### Key Findings & Engineering Insights

- ✅ **Seismic Hazard & Structural Design:** The earthquake signal exhibits strong spectral concentration in low-frequency harmonics ($1\text{ Hz}$ to $3\text{ Hz}$), with the absolute maximum at $2\text{ Hz}$. Consequently, 1–4 story buildings ($3–10\text{ Hz}$) and especially 5–10 story buildings ($1–2\text{ Hz}$) face severe resonance risks. Taller structures (>10 stories, natural frequencies $\le 1\text{ Hz}$) experience significantly reduced structural amplification for this specific seismic profile.
- ✅ **Acoustic Transcription Accuracy:** The 32 detected musical notes form the classic sequence:
  $$\text{E - E - F - G - G - F - E - D - C - C - D - E - E - D - D} \dots$$
  Adaptive amplitude-scaled thresholding ($0.5 \times \text{mean amplitude}$) successfully eliminated transient boundary noise between adjacent notes.
- ✅ **Frequency Shift Invariance:** Multiplying the time signal by $e^{-j 2\pi (3) t}$ shifted the center of the magnitude spectrum from $\omega = 0$ to $\omega = -6\pi \approx -18.85\text{ rad/s}$, confirming that pure imaginary exponential modulation preserves envelope morphology while translating spectral contents.
- ✅ **Frequency-Domain Convolution Validation:** Direct Fourier transformation of windowed signals exactly matches the frequency convolution $\frac{1}{2\pi} [F(\omega) * W(\omega)]$, proving the convolution theorem and illustrating the windowing trade-off between main-lobe spectral resolution and side-lobe leakage.
- ✅ **Basis Expansion Efficiency:** Representing signals via piecewise-linear triangular basis functions converts an infinite-dimensional quadrature problem into a finite vector sum of weighted analytical sinc terms, achieving high accuracy with negligible computational overhead.
- ⚠️ **MATLAB Sinc Normalization Detail:** In MATLAB, the built-in function `sinc(x)` computes $\frac{\sin(\pi x)}{\pi x}$, whereas standard mathematical Fourier formulas utilize the unnormalized definition $\frac{\sin(x)}{x}$. Adjusting arguments by factor $\frac{1}{2\pi}$ or $\frac{1}{\pi}$ is essential when aligning theoretical analytical derivations with MATLAB evaluations.

---

## 🚀 How to Run

### Prerequisites
- MATLAB R2020a or later
- Signal Processing Toolbox (recommended for `audioread`)

### Execution
Open MATLAB, navigate to the `CA1` directory, and run any of the standalone task scripts:

```matlab
% Navigate to CA1
cd('CA1');

% Run Seismograph Analysis (Part 1.1)
run('Q1.m');

% Run Musical Note Recognition (Part 1.2)
run('Q2.m');

% Run Continuous Fourier Transform & Modulation (Part 2.1)
run('Q3.m');

% Run Windowing & Convolution Theorem (Part 2.2)
run('Q4.m');

% Run Triangular Basis Approximation (Part 2.3)
run('Q5.m');
```

---

## 📂 Directory Structure

```
CA1/
├── Q1.m                              # Seismograph frequency analysis & building resonance
├── Q2.m                              # Musical note recognition & automated audio decoding
├── Q3.m                              # Continuous Fourier transform & frequency shift modulation
├── Q4.m                              # Rectangular windowing & frequency-domain convolution theorem
├── Q5.m                              # Triangular basis function approximation & comparative validation
├── general_fourier_transform.m       # Reusable continuous Fourier transform via numerical quadrature
├── triangular_function.m             # Parametric symmetric triangular basis pulse generator
├── triangular_fourier_transform.m    # Analytical Fourier transform of triangular basis pulse
├── approximate_fourier_transform.m   # Piecewise-linear signal approximation & Fourier synthesizer
├── data.mat                          # Earthquake ground acceleration time-series dataset
├── song.wav                          # Synthesized audio file containing musical note sequence
├── Engineering_Mathematics_CA1_Fall2024.pdf # Original assignment problem description
├── Engmath-CA1-810102443.pdf         # Comprehensive student technical report with plots & analysis
└── README.md                         # This documentation file
```

---

## 🛠️ Technologies & Mathematical Principles

| Category | Tools & Concepts |
|---|---|
| **Programming Language** | MATLAB |
| **Mathematical Domains** | Fourier Series, Continuous Fourier Transforms, Signal Processing, Convolution Theorem |
| **Numerical Methods** | Adaptive Quadrature (`integral`), Piecewise-Linear Basis Expansions, Min-Max Normalization |
| **Toolboxes Used** | MATLAB Core, Signal Processing Toolbox |

---

## 📜 License

This project is licensed under the **MIT License** — see the [LICENSE](../LICENSE) file for details.
