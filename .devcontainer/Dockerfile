#-------------------------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See https://go.microsoft.com/fwlink/?linkid=2090316 for license information.
#-------------------------------------------------------------------------------------------------------------

# To fully customize the contents of this image, use the following Dockerfile as a base and add the RUN statement from this file:
# https://github.com/microsoft/vscode-dev-containers/blob/v0.112.0/containers/debian-10-git/.devcontainer/Dockerfile
FROM mcr.microsoft.com/vscode/devcontainers/base:0-debian-10
#FROM mcr.microsoft.com/vscode/devcontainers/base:ubuntu-18.04

# This Dockerfile's base image has a non-root user with sudo access. Use the "remoteUser"
# property in devcontainer.json to use it. On Linux, the container user's GID/UIDs
# will be updated to match your local UID/GID (when using the dockerFile property).
# See https://aka.ms/vscode-remote/containers/non-root-user for details.
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# Configure apt and install packages
RUN apt-get update \
    #
    # Install C++ tools
    && apt-get -y install build-essential cmake cppcheck valgrind \
            libglpk-dev libgsl-dev libgmp-dev libltdl-dev libmpfr-dev pkg-config \
    #
    # [Optional] Update UID/GID if needed
    && if [ "$USER_GID" != "1000" ] || [ "$USER_UID" != "1000" ]; then \
        groupmod --gid $USER_GID $USERNAME \
        && usermod --uid $USER_UID --gid $USER_GID $USERNAME \
        && chown -R $USER_UID:$USER_GID /home/$USERNAME; \
    fi \
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
ENV PATH /opt/conda/bin:$PATH
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    mkdir -p /opt && \
    sh miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda update --all && \
    /opt/conda/bin/conda install viennarna && \
    /opt/conda/bin/conda clean -afy
ENV PKG_CONFIG_PATH /opt/conda/lib/pkgconfig:$PKG_CONFIG_PATH

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog
