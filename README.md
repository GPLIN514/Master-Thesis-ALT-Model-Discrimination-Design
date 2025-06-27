# Model Discrimination Design Generation for Accelerated Life Testing Experiments via Hybridized Optimization Algorithms

This repository contains the complete code and supplementary materials for my master's thesis titled:

**"Model Discrimination Design Generation for Accelerated Life Testing Experiments via Hybridized Optimization Algorithms"**
（中文題目：*以合成最佳化演算法生成加速壽命試驗之模型辨識設計*）

## 📄 Thesis Documents

- `Model Discrimination Design Generation for Accelerated Life Testing Experiments via Hybridized Optimization Algorithms.pdf` — English version of the full thesis  
- `以合成最佳化演算法生成加速壽命試驗之模型辨識設計.pdf` — 中文論文完整版

## 📁 Repository Structure

- `Thesis-code/`  
  Source code for all simulation experiments and Shiny applications. Contains:
  - `code/` — Core R scripts for model generation and divergence functions
  - `shiny-demo/` — Shiny UI & server code (Fidalgo, Meeker, Arrhenius modules)
  - `Appendix B example code.R` — Example code used in the appendix

- `Presentation/`  
  Presentation slides for thesis defense and midterm progress updates.

- `Record/`  
  Records of intermediate simulation results during the thesis study, including:
  - Search results
  - Design summary tables
  - Directional derivative plots

- `README.md` — You are here!

## 🧪 Features

### 🧮 Divergence-Based Design Optimization
- Support for four divergence criteria:  
  KL (CKL), Lin-Wong (CLW), Bhattacharyya (CB), and Chi-square

### 🧠 Hybrid Metaheuristic Algorithms
- Particle Swarm Optimization (PSO)  
- L-BFGS Gradient Search  
- Supports censored lognormal / Weibull models

### 💻 Interactive Shiny App
A fully interactive Shiny application allowing users to:
- Choose between Fidalgo, Arrhenius, and Meeker models
- Select divergence criteria (CKL, CLW, CB, Chi-square)
- Customize design space and view optimal designs and plots

---

## 🔧 Environment & Dependencies

- R (>= 4.2.0)
- Required R packages:
  - `shiny`, `shinyjs`, `ggplot2`, `DiscrimOD`, `DT`, `Rcpp`, `RcppArmadillo`

---

## 📚 Citation

If you use this code, model structure, or methodology, please cite:

> Lin Kuan-Yuan (2025). *Model Discrimination Design Generation for Accelerated Life Testing Experiments via Hybridized Optimization Algorithms*. Master's Thesis, National Taipei University.

---

## 📬 Contact

If you have any questions, feel free to contact me at:  
Email: [a0921129003@gmail.com](mailto:a0921129003@gmail.com)
