# Derived Adelic Cohomology Conjecture (DACC)

This repository contains the implementation of the Derived Adelic Cohomology Conjecture (DACC), a novel cohomological framework for understanding the Birch and Swinnerton-Dyer (BSD) conjecture for elliptic curves.

## Overview

The DACC framework provides a cohomological interpretation that explains both aspects of the BSD conjecture:

1. The equality between the order of vanishing of the L-function and the rank: 
   ```math\text{ASI}(E) = \text{rank}(E) = \text{ord}_{s=1}L(s, E)```

2. The precise formula for the leading coefficient: 
   ```math\frac{L^{(r)}(E,1)}{r!} = \frac{\Omega_E \cdot R_E \cdot \prod_{p} c_p}{\#\text{Sha}(E)}```

Our approach constructs a derived sheaf by gluing local arithmetic data at each place of $\mathbb{Q}$. The resulting adelic complex, equipped with a natural filtration, gives rise to a spectral sequence whose behavior directly encodes the BSD conjecture.

## Key Features

- **Unified Framework**: Explains both the rank equality and the special value formula simultaneously
- **Universal Applicability**: Works for curves of all ranks, not just ranks 0 and 1
- **Cohomological Interpretation**: Provides a natural interpretation of the Tate-Shafarevich group
- **Extensive Verification**: Numerically verified across hundreds of elliptic curves

## Repository Structure

- `dacc_master.sage`: Main orchestrator script for running all components  
- `dacc_curve.sage`: Analysis of individual elliptic curves  
- `dacc_utils.sage`: Utility functions for curve validation and LMFDB access  
- `dacc_config.json`: Configuration file with test curve families  
- `dacc_comprehensive_proof.sage`: Integration of all proof components  
- `dacc_derived_determinant.sage`: Implementation of determinant theory  
- `dacc_spectral_vanishing.sage`: Implementation of differential vanishing theorems  
- `dacc_theoretical_proof.sage`: Generation of theoretical proof documents  
- `dacc_comprehensive_test.sage`: Comprehensive testing framework  
- `dacc_simple_summary.sage`: Generation of simplified summary reports  

## Installation Requirements

### Core Requirements

- **SageMath**: Version 10.4 or later  
- **Python**: 3.8 or later (included with SageMath)  
- **NumPy**: For numerical computations (included with SageMath)  
- **Matplotlib**: For visualizations (included with SageMath)  

### Optional Dependencies

```bash
# Install within SageMath Python environment
sage -pip install plotly pandas scipy adjustText
```

- Plotly & Pandas: For interactive visualizations  
- SciPy: For advanced statistical analysis  
- adjustText: For improved label positioning in plots  

## Usage

### Running the Complete Framework

To run the entire DACC framework, including data retrieval, analysis, proof generation, and visualization:

```
sage dacc_master.sage
```

### Analyzing a Specific Curve

To analyze a single curve using the DACC framework:

```
sage dacc_curve.sage --curve=11.a1
sage dacc_curve.sage --curve=37.a1 --comprehensive  # With comprehensive tests
```

### Running Individual Components

# Run the determinant theory tests  
sage dacc_derived_determinant.sage  

# Generate theoretical proof document  
sage dacc_theoretical_proof.sage  

# Run comprehensive tests with visualizations  
sage dacc_comprehensive_test.sage  

## Examples

### Basic Analysis

```
# Analyze a rank 0 curve with non-trivial Sha
sage dacc_curve.sage --curve=571.a1

# Analyze a rank 2 curve
sage dacc_curve.sage --curve=389.a1 --comprehensive
```

### Customized Testing

You can modify \`dacc_config.json\` to specify custom test families of curves.

## Citing this Work

If you use the DACC framework in your research, please cite:

@article{Wachs2025,
title={The Derived Adelic Cohomology Conjecture for Elliptic Curves},
author={Wachs, Dane},
journal={Preprint},
year={2025},
month={March}
}

## License

This project is licensed under the MIT License - see the LICENSE file for details.
