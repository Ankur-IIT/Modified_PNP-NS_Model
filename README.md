# 🔬 Modified_PNP-NS_Model

This repository implements the **Finite Element Method (FEM)** for the **Modified Poisson–Nernst–Planck/Navier-Stokes (PNP/NS) Model** from the paper:

> <a href="https://www.arxiv.org/abs/2409.08746" target="_blank">Finite Element Method for the Numerical Simulation of Modified PNP/NS Model</a>

Twe simulate the present model using the Finite element scheme discussed in the previous section. We simulate the model for a wide range of parameter values, such as temperature and bulk modulus,  as they affect the numerical stability and qualitative behavior of the model.

---

## 🚀 Installation

### 📦 Using Docker (Recommended)
Install Docker CE for your OS (Windows, Mac, Linux), then run:
```bash
docker run -ti -p 127.0.0.1:8000:8000 \
-v $(pwd):/root/shared \
-w /root/shared \
ghcr.io/scientificcomputing/fenics-gmsh:2024-05-30
```

### 🪟 Windows 10
Enable Windows Subsystem for Linux (WSL) and install Ubuntu from the Microsoft Store. Then follow the Ubuntu installation instructions below.

### 🐧 Ubuntu
To install FEniCS on Ubuntu, run:
```bash
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install fenics
```