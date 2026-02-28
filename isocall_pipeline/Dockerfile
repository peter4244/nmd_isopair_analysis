FROM --platform=linux/amd64 ubuntu:22.04

# Install samtools dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    tabix \
    && rm -rf /var/lib/apt/lists/*

# Copy isocall binary
COPY bin/isocall /usr/local/bin/isocall
RUN chmod +x /usr/local/bin/isocall

# Verify binaries
RUN samtools --version | head -1 && isocall --version
