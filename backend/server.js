// backend/server.js
const express = require('express');
const app = express();
const port = 5000;

app.get('/', (req, res) => {
  res.send('Backend is running!');
});

app.listen(port, () => {
  console.log(`Backend server running on port ${port}`);
});