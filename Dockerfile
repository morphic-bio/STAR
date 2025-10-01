# Option 1: Use pre-built binaries (fastest)
# Build STAR first: ./build_star.sh --static --long
# Then: docker build -f Dockerfile.prebuilt -t star:latest .

# Option 2: Multi-stage build (slower but self-contained)
FROM debian:bookworm-slim as builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    make \
    xxd \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /build

# Copy source code
COPY . .

# Build STAR (static binaries for smaller runtime image)
RUN ./build_star.sh --static --long --clean

# Runtime stage
FROM debian:bookworm-slim

# Install minimal runtime dependencies
RUN apt-get update && apt-get install -y \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Create non-root user
RUN groupadd -r star && useradd -r -g star -s /bin/false star

# Copy binaries from builder
COPY --from=builder /build/bin/Linux_x86_64_static/STAR /usr/local/bin/STAR
COPY --from=builder /build/bin/Linux_x86_64_static/STARlong /usr/local/bin/STARlong

# Set permissions
RUN chmod +x /usr/local/bin/STAR /usr/local/bin/STARlong

# Create working directory for data
WORKDIR /data
RUN chown star:star /data

# Switch to non-root user
USER star

# Set default command
CMD ["STAR", "--help"]

# Labels
LABEL maintainer="STAR Docker Image" \
      description="STAR RNA-seq aligner" \
      version="2.7.11b"
