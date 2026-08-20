import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

// The API has no CORS headers yet, so in dev we route calls through this proxy:
// the browser talks to the same-origin dev server, and Vite forwards to the API
// server-side (where CORS doesn't apply). Remove once the API sends CORS headers.
export default defineConfig({
  plugins: [react()],
  server: {
    proxy: {
      '/api': {
        target: 'https://9368clqa34.execute-api.eu-north-1.amazonaws.com',
        changeOrigin: true,
        secure: true,
        rewrite: (path) => path.replace(/^\/api/, ''),
      },
    },
  },
});