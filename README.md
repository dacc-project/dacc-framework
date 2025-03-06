# Derived Adelic Cohomology Conjecture (DACC)

This repository contains the implementation of the Derived Adelic Cohomology Conjecture (DACC), a novel cohomological framework for understanding the Birch and Swinnerton-Dyer (BSD) conjecture for elliptic curves.

## Overview

The DACC framework provides a cohomological interpretation that explains both aspects of the BSD conjecture:

1. The equality between the order of vanishing of the L-function and the rank: `ASI(E) = rank(E) = ord_{s=1}L(s, E)`
2. The precise formula for the leading coefficient: `L^(r)(E,1)/r! = (Ω_E·R_E·∏c_p)/#Sha(E)`

Our approach constructs a derived sheaf by gluing local arithmetic data at each place of Q. The resulting adelic complex, equipped with a natural filtration, gives rise to a spectral sequence whose behavior directly encodes the BSD conjecture.

## Installation Requirements

⚠️ **Important**: Building SageMath from source is strongly recommended for the DACC framework. Pre-built binaries may not include all necessary components or have version incompatibilities.

The `build.sh` script will help you build SageMath from source and set up the environment:

```bash
# Make the build script executable
chmod +x build.sh

# Run the build script
./build.sh