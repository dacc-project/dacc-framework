# Derived Adelic Cohomology Conjecture (DACC)

This repository contains the implementation of the Derived Adelic Cohomology Conjecture (DACC), a novel cohomological framework for understanding the Birch and Swinnerton-Dyer (BSD) conjecture for elliptic curves.

## Overview

The DACC framework provides a cohomological interpretation that explains both aspects of the BSD conjecture:

1. The equality between the order of vanishing of the L-function and the rank: 
   $$\text{ASI}(E) = \text{rank}(E) = \text{ord}_{s=1}L(s, E)$$

2. The precise formula for the leading coefficient: 
   $$\frac{L^{(r)}(E,1)}{r!} = \frac{\Omega_E \cdot R_E \cdot \prod_{p} c_p}{\text{Sha}(E)}$$

Our approach constructs a derived sheaf by gluing local arithmetic data at each place of $\mathbb{Q}$. The resulting adelic complex, equipped with a natural filtration, gives rise to a spectral sequence whose behavior directly encodes the BSD conjecture.

## Key Features

- **Unified Framework**: Explains both the rank equality and the special value formula simultaneously
- **Universal Applicability**: Works for curves of all ranks, not just ranks 0 and 1
- **Cohomological Interpretation**: Provides a natural interpretation of the Tate-Shafarevich group
- **Extensive Verification**: Numerically verified across hundreds of elliptic curves

## Installation Requirements

⚠️ **Important**: Building SageMath from source is strongly recommended for the DACC framework. Pre-built binaries may not include all necessary components or have version incompatibilities.

The `build.sh` script will help you build SageMath from source and set up the environment:

```bash
# Make the build script executable
chmod +x build.sh

# Run the build script
./build.sh