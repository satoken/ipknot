FROM alpine:latest AS builder

WORKDIR /workspaces
RUN apk add --no-cache musl-dev gcc g++ cmake make ninja git pkgconfig zlib-static

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

# IPknot
COPY . . 
RUN cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_HIGHS=ON -DSTATIC_BUILD=ON -G Ninja -B build \
    && cmake --build build \
    && cmake --install build --strip 

FROM alpine:latest
COPY --from=builder /usr/local/bin/ipknot /usr/local/bin/ipknot