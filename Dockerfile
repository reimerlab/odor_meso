FROM ninai/pipeline:base
LABEL maintainer="Edgar Y. Walker, Fabian Sinz, Erick Cobos, Donnie Kim"

WORKDIR /data

# Upgrade datajoint
RUN pip3 install --upgrade datajoint~=0.12.7

# Upgrade pip
RUN pip3 install --upgrade pip

# Install matplotlib-scalebar
RUN pip3 install matplotlib-scalebar

# Install odor_meso
COPY ./odor_meso /data/odor_meso
RUN pip3 install -e /data/odor_meso/python/

# Install DataPlot
COPY ./DataPlot /data/DataPlot 
RUN pip3 install -e /data/DataPlot

ENTRYPOINT ["/bin/bash"]
