# ── Stage 1: build the C extension ───────────────────────────────────────────
FROM python:3.11-slim AS builder

WORKDIR /app

# Build dependencies for the C extension
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Install hseeker from source (compiles _hdna.c)
COPY pyproject.toml setup.py MANIFEST.in ./
COPY src/ ./src/
RUN pip install --no-cache-dir ".[app]"

# Install webapp dependencies
COPY webapp/requirements.txt ./webapp/requirements.txt
RUN pip install --no-cache-dir -r webapp/requirements.txt


# ── Stage 2: runtime image ────────────────────────────────────────────────────
FROM python:3.11-slim

WORKDIR /app

# Copy installed packages from builder
COPY --from=builder /usr/local/lib/python3.11 /usr/local/lib/python3.11
COPY --from=builder /usr/local/bin /usr/local/bin

# Copy webapp source only (hseeker package is already installed via wheel)
COPY webapp/ ./webapp/

# Non-root user for security
RUN useradd --create-home --shell /bin/bash appuser \
    && chown -R appuser:appuser /app
USER appuser

EXPOSE 8000

# Railway sets $PORT; fall back to 8000 for local docker run
CMD uvicorn webapp.main:app --host 0.0.0.0 --port "${PORT:-8000}"
