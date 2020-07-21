# From ubuntu:18.04

# ENV DEBIAN_FRONTEND=noninteractive

# WORKDIR /workspaces

# RUN apt-get update \
#     && apt-get -y install build-essential wget cmake \
#             libglpk-dev libgsl-dev libgmp-dev libltdl-dev libmpfr-dev pkg-config \
#     && apt-get autoremove -y \
#     && apt-get clean -y \
#     && rm -rf /var/lib/apt/lists/*

# # build from the source
# RUN wget -q https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz \
#     && tar zxvf ViennaRNA-2.4.14.tar.gz \
#     && cd ViennaRNA-2.4.14 \
#     && ./configure --without-perl --without-python --without-python3 --without-forester --without-rnalocmin \
#     && make && make install \
#     && cd .. && rm -rf ViennaRNA-2.4.14 ViennaRNA-2.4.14.tar.gz 

# use the official package
#RUN wget -q https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_18_04/viennarna_2.4.14-1_amd64.deb \
#        https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_18_04/viennarna-dev_2.4.14-1_amd64.deb \
#    && apt-get install ./*.deb \
#    && rm *.deb

From satoken/vsc-viennarna-dev:latest

WORKDIR /workspaces

COPY . .

RUN rm -rf build && mkdir build \
    && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release .. \
    #&& cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXE_LINKER_FLAGS='-static' -DCMAKE_FIND_LIBRARY_SUFFIXES='.a' .. \ # static link
    && make && make install
