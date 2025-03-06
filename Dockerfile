FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    python3-dev \
    libssl-dev \
    libreadline-dev \
    libsqlite3-dev \
    libgmp-dev \
    libmpfr-dev \
    libflint-dev \
    git \
    wget \
    m4 \
    cmake \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create directory for SageMath
WORKDIR /sage-build

# Clone SageMath repository
RUN git clone https://github.com/sagemath/sage.git && \
    cd sage && \
    git checkout 10.0

# Build SageMath (this step takes a long time)
WORKDIR /sage-build/sage
RUN ./configure --with-python=3 && \
    make -j4

# Add SageMath to PATH
ENV PATH="/sage-build/sage:${PATH}"

# Create DACC project directory
WORKDIR /dacc-project

# Copy DACC project files
COPY . /dacc-project/

# Create output directory
RUN mkdir -p dacc_output/dacc_plots

# Set entrypoint
ENTRYPOINT ["/sage-build/sage/sage", "dacc_master.sage"]