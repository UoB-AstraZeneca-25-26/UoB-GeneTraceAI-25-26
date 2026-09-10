// Central place for backend base URLs. Import from here — no hardcoded URLs in components.

// The agent Lambda Function URL's CORS is app-controlled (see AIAgent/api/middleware.py),
// so in dev we route through Vite's same-origin proxy (see vite.config.ts) to avoid needing
// localhost added to ALLOWED_ORIGINS on every redeploy. In a production build, call it
// directly — this requires the deployed origin to be in the app's ALLOWED_ORIGINS.
export const AGENT_API_URL = import.meta.env.DEV
  ? '/agent-api'
  : 'https://mhyvpwmjma3gqy44dk4wskkxje0gqfud.lambda-url.eu-west-2.on.aws';

// Scoring API is consumed through lib/api.ts, which has its own dev-proxy vs.
// prod switch. Kept here for reference/consistency; not wired into that file.
export const SCORING_API_URL = 'https://9368clqa34.execute-api.eu-north-1.amazonaws.com';
