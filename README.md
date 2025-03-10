```markdown
# Derived Adelic Cohomology Conjecture (DACC)

Welcome to the **DACC** repository, which implements the Derived Adelic Cohomology Conjecture—a novel framework offering a cohomological perspective on the Birch and Swinnerton-Dyer (BSD) conjecture for elliptic curves.

## Overview

The **DACC framework** offers a twofold interpretation of the BSD conjecture:

1. **Rank Equality:**
   The order of vanishing of the L-function equals the elliptic curve's rank:
   $$\text{ASI}(E) = \text{rank}(E) = \text{ord}_{s=1}L(s, E)$$

2. **Leading Coefficient Formula:**
   The leading coefficient of the L-series is given by:
   $$\frac{L^{(r)}(E,1)}{r!} = \frac{\Omega_E \cdot R_E \cdot \prod_{p} c_p}{\text{\#Sha}(E)}$$

The method constructs a derived sheaf by combining local arithmetic data from each place of \(\mathbb{Q}\). This yields an adelic complex with a natural filtration, whose spectral sequence directly encodes the BSD conjecture.

## Key Features

- **Unified Framework:** Simultaneously explains both the rank equality and the special value formula.
- **Universal Applicability:** Suitable for elliptic curves of all ranks.
- **Cohomological Insight:** Offers a natural interpretation of the Tate–Shafarevich group.
- **Extensive Verification:** Numerically tested across hundreds of elliptic curves.

## Repository Structure

- dacc_master.sage: Main script to orchestrate the framework.
- dacc_curve.sage: Analyzes individual elliptic curves.
- dacc_utils.sage: Utility functions for curve validation and LMFDB access.
- dacc_config.json: Configuration file for test curve families.
- dacc_comprehensive_proof.sage: Integrates all proof components.
- dacc_derived_determinant.sage: Implements determinant theory.
- dacc_spectral_vanishing.sage: Handles differential vanishing theorems.
- dacc_theoretical_proof.sage: Generates theoretical proof documents.
- dacc_comprehensive_test.sage: Runs comprehensive tests.
- dacc_simple_summary.sage: Generates simplified summary reports.

## Installation Requirements

### Core Requirements

- **SageMath**: Version 10.4 or later.
- **Python**: Version 3.8 or later (bundled with SageMath).
- **NumPy**: For numerical computations (bundled with SageMath).
- **Matplotlib**: For visualizations (bundled with SageMath).

### Optional Dependencies

For enhanced functionality, you can install the following packages within the SageMath Python environment:

```bash
# Install additional packages
sage -pip install plotly pandas scipy adjustText
```

- **Plotly & Pandas:** For interactive visualizations.
- **SciPy:** For advanced statistical analysis.
- **adjustText:** To improve label positioning in plots.

## Usage

### Running the Complete Framework

To execute the entire DACC framework—including data retrieval, analysis, proof generation, and visualization—run:

```
sage dacc_master.sage
```

### Analyzing a Specific Curve

To analyze an individual elliptic curve:

```
sage dacc_curve.sage --curve=11.a1
sage dacc_curve.sage --curve=37.a1 --comprehensive  # Includes comprehensive tests
```

### Running Individual Components

- **Determinant Theory Tests:**  
  `sage dacc_derived_determinant.sage`

- **Generate Theoretical Proof Document:**  
  `sage dacc_theoretical_proof.sage`

- **Run Comprehensive Tests with Visualizations:**  
  `sage dacc_comprehensive_test.sage`

## Examples

### Basic Analysis

```
# Analyze a rank 0 curve with a non-trivial Sha
sage dacc_curve.sage --curve=571.a1

# Analyze a rank 2 curve
sage dacc_curve.sage --curve=389.a1 --comprehensive
```

### Customized Testing

Customize the test families by modifying the `dacc_config.json` file.

## Citing this Work

If you use the DACC framework in your research, please cite:

```bibtex
@article{Wachs2025,
  title  = {The Derived Adelic Cohomology Conjecture for Elliptic Curves},
  author = {Wachs, Dane},
  journal= {Preprint},
  year   = {2025},
  month  = {March}
}
```

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

