# Dockerfile
FROM crazymax/osxcross:latest

# Install Miniconda
RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/conda/bin:$PATH"

# Create and activate a conda environment
RUN conda create -n build-wheels python=3.11 numpy -y && \
    echo "source activate build-wheels" > ~/.bashrc

# Install cibuildwheel
RUN /opt/conda/bin/conda run -n build-wheels pip install cibuildwheel==2.18.1

# Set environment variables for cibuildwheel
ENV CIBW_ARCHS_MACOS="x86_64"
ENV CIBW_BUILD_VERBOSITY=1
ENV CIBW_ENVIRONMENT="FCFLAGS='-arch x86_64' FFLAGS='-arch x86_64' LDFLAGS='-arch x86_64' MACOSX_DEPLOYMENT_TARGET=11.0"

# Copy the project files
COPY . /project

# Set the working directory
WORKDIR /project

# Build the wheels
CMD ["/opt/conda/bin/conda", "run", "-n", "build-wheels", "cibuildwheel", "--output-dir", "wheelhouse"]
