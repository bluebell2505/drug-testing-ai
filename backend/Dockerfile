FROM node:20-alpine
WORKDIR /app
COPY package*.json ./
RUN npm install
COPY . .
# Explicitly expose port 8080
EXPOSE 8080
CMD ["node", "server.js"]