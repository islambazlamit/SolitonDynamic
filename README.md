# SolitonDynamic

[![C++ Standard](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![IDE](https://img.shields.io/badge/IDE-Xcode-informational)](https://developer.apple.com/xcode/)

`SolitonDynamic` is a high-precision C++ simulation framework designed to study soliton dynamics by solving the **1D Nonlinear Schrödinger Equation (NLSE) with dual nonlinearity and an external potential** using the advanced **Iterative Power Series (IPS)** method.

This project is a modern, ultra-high-precision evolution of the original FORTRAN implementation developed by **Prof. Usama Al Khawaja**. While the original FORTRAN version provided standard 16-digit (double) precision, this C++ version extends the computational limits significantly—offering **arbitrary high precision up to 50 to 100 digits of accuracy**.

---

## 🔬 Physics & Mathematical Model

The core solver simulates the generalized 1D NLSE, incorporating both cubic and quintic nonlinearities along with a spatially dependent external potential:

$$i \psi_t(x,t) + a_1 \psi_{xx}(x,t) + a_2 |\psi(x,t)|^2 \psi(x,t) + a_3 |\psi(x,t)|^4 \psi(x,t) + V(x) \psi(x,t) = 0$$

Where:
*   $\psi(x,t)$ is the complex wave function.
*   $a_1$ is the dispersion coefficient.
*   $a_2$ and $a_3$ represent the coefficients for the dual (cubic-quintic) nonlinearities.
*   $V(x)$ is the external potential.

### Key Features:
*   **Flexible Framework:** Fully capable of handling both **single** and **two-soliton** solutions.
*   **Customizable Parameters:** Easily choose and configure your own initial soliton solutions and external potentials $V(x)$.
*   **Scientific Validation:** The original codebase has been widely utilized in various peer-reviewed research papers. The specific solution configured in this release is referenced in:
    > **Reference Paper:** *Phys. Rev. E 107, 014202*  
    > 🔗 [DOI: 10.1103/PhysRevE.107.014202](https://doi.org/10.1103/PhysRevE.107.014202)

---

## ⚡ Key Enhancements in the C++ Version

*   **Ultra-High Precision:** Upgraded from 16-digit floating-point precision to custom multi-precision layouts capable of sustained **50–100 digit accuracy**, critical for analyzing chaotic dynamics or long-term stability in soliton interactions.
*   **Modern Software Stack:** Rebuilt from scratch using modern **C++20** standards.
*   **Iterative Power Series (IPS) Solver:** Uses an optimized recurrence routine to compute the time-evolution operator with extreme fidelity.

---

## 🛠️ Environment & Prerequisites

*   **IDE:** Apple Xcode
*   **Compiler:** C++20 compliant compiler (e.g., Clang/LLVM via Xcode)

---

## 🚧 Status Note
> ⚠️ **Under Construction:** The current version of the code is actively being refactored to optimize performance, enhance code clarity, and improve overall ease of use for the scientific community.

---

## 👤 Author & Contact

**Islam M. Albazlamit**  
*   **Email:** [islam.albazlamit@outlook.com](mailto:islam.albazlamit@outlook.com)  
*   **LinkedIn:** [linkedin.com/in/islambzl](https://www.linkedin.com/in/islambzl/)  

*Original Method & FORTRAN Implementation by:* **Prof. Usama Al Khawaja**

---

## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
