FROM python:3.9-slim

# Install RDKit dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libboost-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
# Gunicorn handles multiple requests simultaneously (unlike Flask’s dev server).
CMD ["gunicorn", "--bind", "0.0.0.0:8080", "app:app"]