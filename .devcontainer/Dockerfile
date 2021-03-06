# To fully customize the contents of this image, use the following Dockerfile as a base and add the RUN statement from this file:
# https://github.com/microsoft/vscode-dev-containers/blob/v0.112.0/containers/debian-10-git/.devcontainer/Dockerfile
#FROM mcr.microsoft.com/vscode/devcontainers/base:0-debian-10
#FROM mcr.microsoft.com/vscode/devcontainers/base:ubuntu-18.04
FROM mcr.microsoft.com/vscode/devcontainers/miniconda:0-3

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# Configure apt and install packages
RUN apt-get update \
    # Install C++ tools
    && apt-get -y install build-essential cmake cppcheck valgrind \
            libglpk-dev libgsl-dev libgmp-dev libltdl-dev libmpfr-dev pkg-config \
    # Clean up
    && rm -f *.deb \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# RUN wget -q https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.17.tar.gz \
#     && tar zxvf ViennaRNA-2.4.17.tar.gz \
#     && cd ViennaRNA-2.4.17 \
#     && ./configure --without-perl --without-python --without-python3 --without-forester --without-rnalocmin \
#     && make && make install \
#     && cd .. && rm -rf ViennaRNA-2.4.17 ViennaRNA-2.4.17.tar.gz 

# install ViennaRNA package from bioconda
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda update --all && \
    conda install viennarna && \
    conda clean -afy
ENV PKG_CONFIG_PATH /opt/conda/lib/pkgconfig:$PKG_CONFIG_PATH

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog
