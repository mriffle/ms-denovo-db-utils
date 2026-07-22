FROM python:3.10-slim

ARG VERSION=dev
ARG VCS_REF=unknown

LABEL org.opencontainers.image.title="ms-denovo-db-utils" \
      org.opencontainers.image.description="Utilities for the ms-denovo-db pipeline" \
      org.opencontainers.image.source="https://github.com/mriffle/ms-denovo-db-utils" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}"

# procps supplies ps, which Nextflow shells out to for task resource metrics.
# hadolint ignore=DL3008
RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean

# Copied before the sources so the install layer stays cache-friendly.
COPY pyproject.toml README.md LICENSE /build/
COPY src /build/src
RUN pip install --no-cache-dir /build && rm -rf /build

# nf-ms-denovo-db invokes these absolute paths directly; the shims forward to
# ms_denovo_db_utils.cli.
COPY docker_bin/ /usr/local/bin/
COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod 755 /usr/local/bin/entrypoint.sh /usr/local/bin/*.py

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD []
