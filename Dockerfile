# From ubuntu:18.04
#FROM debian:10
FROM continuumio/miniconda3

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /workspaces

RUN apt-get update \
    # Install C++ tools
    && apt-get -y install build-essential cmake cppcheck valgrind \
            libglpk-dev libgsl-dev libgmp-dev libltdl-dev libmpfr-dev pkg-config \
    # Clean up
    && rm -f *.deb \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog

# install ViennaRNA package from bioconda
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda update --all && \
    conda install viennarna && \
    conda clean -afy
ENV PKG_CONFIG_PATH /opt/conda/lib/pkgconfig:$PKG_CONFIG_PATH

# build from the source
# RUN wget -q https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.17.tar.gz \
#     && tar zxvf ViennaRNA-2.4.17.tar.gz \
#     && cd ViennaRNA-2.4.17 \
#     && ./configure --without-perl --without-python --without-python3 --without-forester --without-rnalocmin \
#     && make && make install \
#     && cd .. && rm -rf ViennaRNA-2.4.17 ViennaRNA-2.4.17.tar.gz 

COPY . .

RUN rm -rf build && mkdir build \
    && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release .. \
    #&& cmake -DCMAKE_BUILD_TYPE=Release -DSTATIC_BULD=ON .. \
    && make && make install