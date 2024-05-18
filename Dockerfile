FROM python:3.12

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /workspaces
RUN apt update \
    && apt install -y gcc g++ cmake make ninja-build git pkg-config \
    && apt autoremove \
    && apt clean 

# HiGHS
RUN git clone https://github.com/ERGO-Code/HiGHS \
    && cd HiGHS \
    && cmake -DFAST_BUILD=ON -DCMAKE_BUILD_TYPE=Release -G Ninja -B build \
    && cmake --build build \
    && cmake --install build --strip \
    && cd ..

# ViennaRNA
COPY --from=satoken/viennarna:latest /usr/local/ /usr/local/
# build from the source
# RUN wget -q https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.17.tar.gz \
#     && tar zxvf ViennaRNA-2.4.17.tar.gz \
#     && cd ViennaRNA-2.4.17 \
#     && ./configure --without-perl --without-python --without-python2 --without-forester --without-rnalocmin \
#     && make && make install \
#     && cd .. && rm -rf ViennaRNA-2.4.17 ViennaRNA-2.4.17.tar.gz 

# MXfold2
RUN pip install --no-cache-dir "pybind11[global]"
RUN git clone https://github.com/mxfold/mxfold2.git \
    && cd mxfold2 \
    && git fetch origin with_ipknot \
    && git checkout with_ipknot \
    && pip install --no-cache-dir . \
    && cd .. \
    && rm -rf mxfold2

# IPknot
COPY . . 
RUN cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_HIGHS=ON -G Ninja -B build \
    && cmake --build build \
    && cmake --install build --strip \
    && pip install --no-cache-dir matplotlib seaborn \
    && cp utils/rna2dheatmap.py /usr/local/bin/ \
    && chmod +x /usr/local/bin/rna2dheatmap.py \
    && rm -rf /workspaces/*
