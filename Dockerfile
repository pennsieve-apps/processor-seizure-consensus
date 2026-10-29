FROM --platform=linux/amd64 python:3.12.15-slim@sha256:a0f57064cbdcd6285fa45f51237ab0b35ea8ebaed06a1f5973e89e19dfb86704

WORKDIR /processor

# Pinned Python dependencies (binary wheels, no apt packages), plus the Lambda runtime
# interface client for running as a Lambda function.
COPY requirements.txt /processor/requirements.txt
RUN pip install --no-cache-dir --only-binary=:all: -r /processor/requirements.txt \
    && find /usr/local/lib/python3.12 -depth -type d \( -name __pycache__ -o -name tests \) -path "*site-packages*" -exec rm -rf {} +

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    HOME=/tmp

# Runs as uid 1000 on the platform; nothing here needs root.
USER 1000

# Runs main() once as a container, or serves processor.handler as a Lambda
# function (see processor/entry.py).
ENTRYPOINT ["python3", "-m", "processor.entry"]
