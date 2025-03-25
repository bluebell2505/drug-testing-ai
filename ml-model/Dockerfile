# ml-model/Dockerfile
FROM python:3.9-slim

# Install system dependencies for RDKit
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libboost-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
CMD ["python", "app.py"]