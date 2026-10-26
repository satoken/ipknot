FROM mcr.microsoft.com/devcontainers/cpp:1-debian-12

ARG REINSTALL_CMAKE_VERSION_FROM_SOURCE="none"

# Optionally install the cmake for vcpkg
COPY ./reinstall-cmake.sh /tmp/

RUN if [ "${REINSTALL_CMAKE_VERSION_FROM_SOURCE}" != "none" ]; then \
        chmod +x /tmp/reinstall-cmake.sh && /tmp/reinstall-cmake.sh ${REINSTALL_CMAKE_VERSION_FROM_SOURCE}; \
    fi \
    && rm -f /tmp/reinstall-cmake.sh

# [Optional] Uncomment this section to install additional vcpkg ports.
# RUN su vscode -c "${VCPKG_ROOT}/vcpkg install <your-port-name-here>"

# [Optional] Uncomment this section to install additional packages.
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends libglpk-dev 

COPY --from=satoken/viennarna:latest /usr/local/ /usr/local/
COPY --from=python:3.12 /usr/local/ /usr/local/
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# MXfold2
RUN pip install --no-cache-dir "pybind11[global]"
RUN git clone https://github.com/mxfold/mxfold2.git \
    && cd mxfold2 \
    && git fetch origin with_ipknot \
    && git checkout with_ipknot \
    && pip install --no-cache-dir . \
    && cd .. \
    && rm -rf mxfold2