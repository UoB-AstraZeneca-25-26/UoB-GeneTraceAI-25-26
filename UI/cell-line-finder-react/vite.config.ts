import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

// Neither API sends CORS headers yet, so in dev we route calls through these
// proxies: the browser talks to the same-origin dev server, and Vite forwards
// to the API server-side (where CORS doesn't apply). Remove each entry once
// its API sends CORS headers (or configure CORS on the Lambda Function URL).
export default defineConfig({
  plugins: [react()],
  server: {
    proxy: {
      '/api': {
        // Local override for local-to-local testing (uncommitted -- see api.ts's
        // USE_MOCK for the matching override): points at `uvicorn
        // final_pipeline.api.app:app` on :8000 instead of the deployed Lambda.
        // Local routes have no /prod stage prefix, so strip that too.
        // Revert to the deployed target below before committing.
        target: 'http://localhost:8000',
        changeOrigin: true,
        secure: false,
        rewrite: (path) => path.replace(/^\/api\/prod/, ''),
        // target: 'https://9368clqa34.execute-api.eu-north-1.amazonaws.com',
        // secure: true,
        // rewrite: (path) => path.replace(/^\/api/, ''),
      },
      '/agent-api': {
        target: 'https://mhyvpwmjma3gqy44dk4wskkxje0gqfud.lambda-url.eu-west-2.on.aws',
        changeOrigin: true,
        secure: true,
        rewrite: (path) => path.replace(/^\/agent-api/, ''),
      },
    },
  },
});