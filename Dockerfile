FROM ninai/pipeline:base
LABEL maintainer="Edgar Y. Walker, Fabian Sinz, Erick Cobos, Donnie Kim"

WORKDIR /data

# Upgrade datajoint
RUN pip3 install --upgrade datajoint~=0.12.7

# Install pipeline
COPY . /data/pipeline
RUN pip3 install -e pipeline/python/

ENTRYPOINT ["/bin/bash"]
