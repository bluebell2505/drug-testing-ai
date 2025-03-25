const express = require('express');
const app = express();
const port = 8080;

app.get('/', (req, res) => {
  res.send('Backend is running!');
});

app.listen(port, () => {
  console.log(`Backend server running on port ${port}`);
});