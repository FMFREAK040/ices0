FROM debian:bookworm-slim AS build

LABEL name="docker-ices" \
      maintainer="Jochem van den Anker" \
      description="Ices is free server software for streaming multimedia." \
      url="§https://icecast.org"

COPY rootfs /

RUN apt update && \
	apt install -y \
        tar \
        automake \
        autoconf \
        gcc \
        git \
        build-essential \
        checkinstall \
        libncursesw5-dev \
        libssl-dev \
        libsqlite3-dev \
        tk-dev \
        libgdbm-dev \
        libc6-dev \
        libbz2-dev \
        libffi-dev \
        wget \
        libtool \
        libshout-dev \
        libmp3lame-dev \
        libxml2-dev \
        pip \
        python3 \
        python3-dev
        
WORKDIR /
RUN git clone https://github.com/FMFREAK040/ices0.git
WORKDIR /ices0
RUN aclocal && \
    autoreconf -fi && \
    autoconf && \
    automake --add-missing
RUN ./configure --with-python=/usr/bin/python3 --with-perl=no && \
    make && \
    make install

FROM debian:bookworm-slim AS ices0
RUN apt update && \
    apt install -y\
    libshout-dev \
    python3 \
    python3-dev \
    libxml2-dev \
    libmp3lame-dev
COPY rootfs/config/ices.py /usr/local/etc/modules/ices.py
COPY --from=build /usr/local/bin/ices /usr/local/bin/ices
COPY rootfs /
ENTRYPOINT ["ices", "-c", "/config/ices.conf"]

