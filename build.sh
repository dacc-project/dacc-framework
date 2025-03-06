### build.sh

```bash
#!/bin/bash
# DACC Project Build Script
# This script helps set up the environment for the DACC framework

echo "DACC Framework Setup"
echo "===================="
echo "This script will help you set up the environment for running the DACC framework."
echo ""

# Check if SageMath is installed
if command -v sage &> /dev/null; then
    echo "✓ SageMath is already installed in your PATH"
    sage_path=$(which sage)
    echo "  Path: $sage_path"
    sage_version=$(sage --version)
    echo "  Version: $sage_version"
    
    read -p "Do you want to use this existing SageMath installation? (y/n) " use_existing
    if [[ $use_existing == "y" ]]; then
        echo "Using existing SageMath installation"
    else
        echo "Building SageMath from source is recommended for DACC framework"
        read -p "Build SageMath from source now? (y/n) " build_sage
        
        if [[ $build_sage == "y" ]]; then
            # Build SageMath from source
            echo "Building SageMath from source..."
            
            # Create directory for SageMath
            mkdir -p ~/sage-build
            cd ~/sage-build
            
            # Download SageMath source
            echo "Downloading SageMath source..."
            git clone https://github.com/sagemath/sage.git
            cd sage
            git checkout 10.0
            
            # Configure and build
            echo "Configuring and building SageMath (this will take some time)..."
            ./configure --with-python=3
            make -j4
            
            # Create configuration file
            echo "Creating configuration file..."
            cd ~/dacc-project
            echo "export SAGE_ROOT=~/sage-build/sage" > config.sh
            echo "export PATH=\$SAGE_ROOT:\$PATH" >> config.sh
            chmod +x config.sh
            
            echo "SageMath built successfully at ~/sage-build/sage"
            echo "Before running DACC scripts, source the configuration file:"
            echo "  source config.sh"
        else
            echo "Please install SageMath before running the DACC framework"
            exit 1
        fi
    fi
else
    echo "SageMath not found in your PATH"
    read -p "Build SageMath from source now? (y/n) " build_sage
    
    if [[ $build_sage == "y" ]]; then
        # Build SageMath from source (same as above)
        echo "Building SageMath from source..."
        
        # Create directory for SageMath
        mkdir -p ~/sage-build
        cd ~/sage-build
        
        # Download SageMath source
        echo "Downloading SageMath source..."
        git clone https://github.com/sagemath/sage.git
        cd sage
        git checkout 10.0
        
        # Configure and build
        echo "Configuring and building SageMath (this will take some time)..."
        ./configure --with-python=3
        make -j4
        
        # Create configuration file
        echo "Creating configuration file..."
        cd ~/dacc-project
        echo "export SAGE_ROOT=~/sage-build/sage" > config.sh
        echo "export PATH=\$SAGE_ROOT:\$PATH" >> config.sh
        chmod +x config.sh
        
        echo "SageMath built successfully at ~/sage-build/sage"
        echo "Before running DACC scripts, source the configuration file:"
        echo "  source config.sh"
    else
        echo "Please install SageMath before running the DACC framework"
        exit 1
    fi
fi

# Create output directory
echo "Creating output directory..."
mkdir -p dacc_output/dacc_plots

echo ""
echo "DACC framework is ready to use."
echo "To run the full analysis:"
echo "  sage dacc_master.sage"
echo ""
echo "Or try one of the examples:"
echo "  sage examples/rank0_curves.sage"