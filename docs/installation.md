# Installation Guide for DACC Framework

## Prerequisites

- Linux or macOS environment (recommended; Windows users should use WSL2)
- Git
- Build tools: gcc, g++, make
- Required libraries: ssl, readline, etc.
- At least 8GB of RAM and 20GB of disk space for building SageMath

## Building SageMath from Source (Recommended)

Building SageMath from source is **strongly recommended** for the DACC framework to ensure all components function correctly.

### Step 1: Download SageMath Source

```bash
# Create a directory for SageMath
mkdir -p ~/sage-build
cd ~/sage-build

# Download SageMath 10.0 source (or later version)
git clone https://github.com/sagemath/sage.git
cd sage
git checkout 10.0