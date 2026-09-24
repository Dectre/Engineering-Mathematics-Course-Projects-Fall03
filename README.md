# 🎓 Engineering Mathematics — Fall 2024 (Fall 1403)

**Comprehensive collection of computational engineering mathematics projects at the University of Tehran, exploring Fourier analysis, seismic hazard evaluation, acoustic note recognition, and numerical solutions to partial differential equations (2D magnetostatics via FEM, 1D/2D heat diffusion via FDM and `pdepe`).**

[![University](https://img.shields.io/badge/University%20of%20Tehran-ECE-blue.svg)](#)
[![Semester](https://img.shields.io/badge/Semester-Fall%201403%20%2F%202024-orange.svg)](#)
[![Language](https://img.shields.io/badge/MATLAB-R2024a-green.svg)](#)
[![License](https://img.shields.io/badge/License-MIT-purple.svg)](./LICENSE)

---

### Department of Electrical and Computer Engineering — University of Tehran

> 👨‍💻 **Author:** Amirali Dehghani (`810102443`)  
> 🎓 **Course:** Engineering Mathematics — Dr. Mehdi Tale Masouleh  
> 🏢 **Faculty:** Faculty of Electrical and Computer Engineering, University of Tehran

---

## 🚀 Projects Overview

| Project | Topic | Key Highlights | Report | Docs |
|:-------:|-------|----------------|:------:|:----:|
| **CA1** | **Fourier Analysis & Signal Processing** | Seismic harmonic analysis & resonance, audio note recognition, continuous FT quadrature, windowing & convolution theorem, triangular basis expansion | [📄 Report](./CA1/Engmath-CA1-810102443.pdf) | [📖 README](./CA1/README.md) |
| **CA2** | **Numerical Partial Differential Equations (PDEs)** | 2D magnetostatic FEM in PDE Modeler, theoretical comparative survey (FDM vs. FVM vs. FEM), 1D heat equation via `pdepe`, 2D steady & transient heat diffusion via FDM | [📄 Report](./CA2/Engmath-CA2-810102443.pdf) | [📖 README](./CA2/README.md) |

---

## 🔹 CA1: Fourier Analysis, Earthquake Dynamics & Signal Processing

**Topic:** Harmonic decomposition, continuous Fourier analysis, spectral theorems, and basis function approximation with real-world applications in civil and acoustic engineering.

- **1.1 Seismograph Frequency Extraction & Structural Resonance:** Analyzed real-world ground acceleration time-series ($f_s = 1000\text{ Hz}$) from `data.mat`, calculated discrete exponential Fourier series coefficients for dominant frequencies ($1–10\text{ Hz}$), and evaluated resonance hazards across low-rise ($3–10\text{ Hz}$), mid-rise ($1–2\text{ Hz}$), and high-rise ($0.1–1\text{ Hz}$) structures.
- **1.2 Automated Acoustic Note Recognition:** Segmented a musical audio track (`song.wav`, $44,100\text{ Hz}$) into uniform $0.4\text{ s}$ intervals, evaluated discrete Fourier energy at fundamental notes (C, D, E, F, G), applied adaptive amplitude-based thresholding, and decoded the 32-note melody (Beethoven's "Ode to Joy").
- **2.1 Generalized Continuous Fourier Transform & Modulation:** Implemented continuous Fourier integration via adaptive numerical quadrature over compact support $[-10, 10]\text{ s}$, and experimentally demonstrated the modulation property $\mathcal{F}\{f(t) e^{j\omega_0 t}\} = F(\omega - \omega_0)$ by inducing a $-6\pi\text{ rad/s}$ rigid spectral translation.
- **2.2 Time-Domain Windowing & Convolution Theorem:** Evaluated rectangular aperture windowing on $\text{sinc}^2(t)$ and $e^{-3|t|}$, empirically validating the Fourier convolution theorem $\mathcal{F}\{f \cdot w\} = \frac{1}{2\pi} [F * W]$ and characterizing sinc-induced spectral leakage and main-lobe broadening.
- **2.3 Piecewise-Linear Signal Modeling via Triangular Bases:** Decomposed arbitrary continuous compact-support signals into shifted triangular basis pulses $\Lambda(t)$, analytically synthesized continuous spectra via linearity and time shifting, and benchmarked against theoretical solutions and general numerical integration.

**Key Concepts:** Exponential Fourier Series, Continuous Fourier Transform (CFT), Frequency Modulation, Rectangular Windowing, Convolution Theorem, Spectral Leakage, Structural Resonance, Triangular Basis Expansion.

**Notable Results:**
- Identified peak seismic energy at $2\text{ Hz}$ and $1\text{ Hz}$, showing that 5–10 story buildings are at highest risk of destructive resonance for the given earthquake profile.
- Successfully transcribed all 32 notes of the audio file with $100\%$ accuracy using an adaptive threshold of $0.5 \times \text{mean amplitude}$.
- Confirmed exact equivalence between time-domain window multiplication and normalized frequency-domain convolution $\frac{1}{2\pi} (F * W)$.
- Demonstrated that piecewise-linear triangular expansion provides analytical Fourier transform accuracy comparable to expensive adaptive quadrature integration.

👉 **[Read CA1 Full Documentation](./CA1/README.md)**

---

## 🔹 CA2: Numerical Partial Differential Equations (PDEs)

**Topic:** Numerical solutions to elliptic and parabolic partial differential equations using the Finite Element Method (FEM), Finite Difference Method (FDM), and MATLAB ODE/PDE solvers.

- **Part 1: 2D Magnetostatic Field Modeling (FEM via PDE Modeler):** Formulated an electromagnetic actuator system consisting of an iron core ($\mu_r = 2000$), copper excitation coil windings ($J_z = \pm 3 \times 10^7\text{ A/m}^2$), and a working air gap. Solved the vector Poisson equation $-\nabla \cdot \left(\frac{1}{\mu} \nabla A_z\right) = J_z$ with homogeneous Dirichlet boundary conditions ($A_z = 0$) on an unstructured triangular finite element mesh. Analyzed magnetic flux guiding, reluctance paths, and the impact of core degradation.
- **Part 2.1: Numerical PDE Methodology Benchmark:** Conducted a comprehensive theoretical comparison between the Finite Difference Method (FDM), Finite Volume Method (FVM), and Finite Element Method (FEM), comparing governing principles, mesh adaptability, local/global conservation laws, and computational overhead.
- **Part 2.2: 1D Parabolic Heat Equation via MATLAB `pdepe`:** Solved the transient 1D heat equation $\frac{\partial u}{\partial t} = \frac{\partial^2 u}{\partial x^2}$ with initial profile $u(x, 0) = \frac{2x}{1 + x^2}$ and fixed boundary temperatures ($u(0) = 0$, $u(1) = 1$), observing asymptotic relaxation toward the steady-state profile $u(x) = x$.
- **Part 2.3: 2D Steady-State Heat Conduction (Laplace) via FDM:** Discretized the elliptic equation $\nabla^2 T = 0$ on a $100 \times 100$ uniform Cartesian grid using a 5-point discrete Laplacian stencil with Jacobi iterative relaxation under non-uniform boundary conditions (Top: 75°C, Bottom: 130°C, Left: 110°C, Right: 45°C), converging to tolerance $\epsilon = 5 \times 10^{-4}$.
- **Part 2.4: 2D Transient Heat Diffusion via Explicit FTCS:** Discretized the 2D parabolic diffusion equation $\frac{\partial T}{\partial t} = \alpha \nabla^2 T$ with thermal diffusivity $\alpha = 0.1\text{ m}^2/\text{s}$ using explicit Forward-Time Central-Space finite differences. Enforced the CFL stability criterion ($r = \frac{\alpha \Delta t}{\Delta x^2} \approx 0.098 \le 0.25$) with $\Delta t = 0.0001\text{ s}$, and captured temporal thermal diffusion snapshots at $t = 0, 0.25, 0.5, 1.0\text{ s}$.

**Key Concepts:** Magnetostatics, Magnetic Vector Potential ($A_z$), Relative Permeability & Reluctance, Finite Element Method (FEM), Finite Difference Method (FDM), Finite Volume Method (FVM), 5-Point Laplacian Stencil, Jacobi Relaxation, Forward-Time Central-Space (FTCS), Von Neumann / CFL Stability.

**Notable Results:**
- Proved that the high-permeability iron core ($\mu_r = 2000$) channels almost all magnetic flux across the air gap, and demonstrated that reducing permeability leads to severe flux leakage and field dissipation.
- Validated that while magnetic flux density $\mathbf{B}$ concentrates inside the iron core, magnetic field intensity $\mathbf{H}$ peaks inside the air gap and coil due to the constitutive relation $\mathbf{H} = \mathbf{B}/\mu$.
- Verified numerical stability and smooth thermal evolution for the 2D explicit FTCS solver under the CFL threshold, matching the steady-state Laplace solution as $t \to 1.0\text{ s}$.

👉 **[Read CA2 Full Documentation](./CA2/README.md)**

---

## 🛠️ Tools & Technologies

| Category | Tools & Frameworks | Details |
|---|---|---|
| **Programming Language** | MATLAB (R2024a) | Primary computing and simulation platform |
| **Toolboxes** | PDE Toolbox, Signal Processing Toolbox | Mesh generation, finite element solver, ODE/PDE integrators, audio I/O |
| **Mathematical Domains** | Fourier Analysis, Numerical PDEs | Fourier Series, Continuous Fourier Transforms, Laplace & Heat Diffusion Equations, Maxwell's Equations |
| **Numerical Schemes** | Quadrature, FEM, FDM, FTCS | Gauss-Kronrod quadrature, triangular basis expansion, 5-point discrete Laplacian, explicit Euler time-stepping |
| **Visualizations** | 2D/3D MATLAB Plotting Engine | `pcolor`, `surf`, `shading interp`, `contour`, PDE Modeler GUI |

---

## 📂 Repository Structure

```
Engineering-Mathematics-Course-Projects-Fall03/
├── CA1/                                     # Computer Assignment 1: Fourier Analysis & Signal Processing
│   ├── README.md                            # Comprehensive CA1 documentation
│   ├── Q1.m                                 # Seismograph frequency analysis & building resonance
│   ├── Q2.m                                 # Automated musical note recognition & audio decoding
│   ├── Q3.m                                 # Continuous Fourier transform & frequency modulation
│   ├── Q4.m                                 # Rectangular windowing & convolution theorem verification
│   ├── Q5.m                                 # Piecewise-triangular basis expansion & spectral synthesis
│   ├── general_fourier_transform.m          # Numerical continuous Fourier transform function
│   ├── triangular_function.m                # Parametric symmetric triangular pulse generator
│   ├── triangular_fourier_transform.m       # Analytical Fourier transform of triangular basis
│   ├── approximate_fourier_transform.m      # Triangular basis expansion synthesizer
│   ├── data.mat                             # Earthquake ground acceleration time-series dataset
│   ├── song.wav                             # Acoustic musical note recording
│   ├── Engineering_Mathematics_CA1_Fall2024.pdf # Original CA1 problem description
│   └── Engmath-CA1-810102443.pdf            # Technical report with derivations and results
│
├── CA2/                                     # Computer Assignment 2: Numerical Partial Differential Equations
│   ├── README.md                            # Comprehensive CA2 documentation
│   ├── Part1.m                              # 2D Magnetostatic finite element model (PDE Modeler)
│   ├── Part2_2.m                            # 1D Parabolic Heat Equation solution via pdepe
│   ├── Part2_3.m                            # 2D Steady-State Laplace Heat Equation via iterative FDM
│   ├── Part2_4.m                            # 2D Transient Heat Diffusion via explicit FTCS FDM
│   ├── 1-2.png ~ 1-Q3-2.png                 # Simulation screenshots, meshes, and field distributions
│   ├── Engmath-CA2-810102443.pdf            # Technical report with theoretical derivations and figures
│   └── Engmath-CA2-810102443.docx           # Editable report document
│
├── LICENSE                                  # MIT License
└── README.md                                # Root repository documentation
```

---

## 🚀 How to Run

### Requirements
- MATLAB R2020a or later
- Partial Differential Equation Toolbox
- Signal Processing Toolbox

### Running CA1 (Fourier Analysis):
```matlab
cd('CA1');

% Run individual problems:
run('Q1.m');  % Earthquake analysis & structural resonance
run('Q2.m');  % Musical note recognition
run('Q3.m');  % Continuous FT & modulation property
run('Q4.m');  % Windowing & convolution theorem
run('Q5.m');  % Triangular basis approximation
```

### Running CA2 (Numerical PDEs):
```matlab
cd('CA2');

% Run individual models:
run('Part1.m');    % 2D Magnetostatic simulation in PDE Modeler
run('Part2_2.m');  % 1D Heat Equation using pdepe
run('Part2_3.m');  % 2D Steady-State Heat Equation (FDM Jacobi)
run('Part2_4.m');  % 2D Transient Heat Equation (FDM FTCS)
```

---

## 📊 Reports & Deliverables

Each assignment folder includes:
- 📄 **Technical PDF Report** — Detailed theoretical derivations, mathematical models, methodology descriptions, and labeled simulation figures.
- 💻 **Documented MATLAB Scripts** — Clean, modular code files implementing all numerical algorithms and visualization routines.
- 🔬 **Simulation Artifacts & Data** — Datasets (`data.mat`, `song.wav`) and exported PDE simulation figures.
- 📖 **Comprehensive Documentation** — Dedicated `README.md` files detailing assignment objectives, methodologies, and engineering conclusions.

---

## 👨‍💻 Author

**Amirali Dehghani**  
Department of Electrical and Computer Engineering  
University of Tehran  
Student ID: `810102443`  

---

## 📜 License

This repository is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
