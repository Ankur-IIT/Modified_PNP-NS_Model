# 🔬 Modified_PNP-NS_Model

This repository implements the **Finite Element Method (FEM)** for the **Modified Poisson–Nernst–Planck/Navier-Stokes (PNP/NS) Model** from the paper:

> ["Finite Element Method for the Numerical Simulation of Modified PNP/NS Model"](https://www.arxiv.org/abs/2409.08746)

The model solves for parameters like **temperature** and **bulk modulus** to study system stability.

---

## 🚀 Installation

### 📦 Using Docker (Recommended)
Install Docker CE for your OS (Windows, Mac, Linux), then run:
```bash
docker run -ti -p 127.0.0.1:8000:8000 \
-v $(pwd):/root/shared \
-w /root/shared \
ghcr.io/scientificcomputing/fenics-gmsh:2024-05-30