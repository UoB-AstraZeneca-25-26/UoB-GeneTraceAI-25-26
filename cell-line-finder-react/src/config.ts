// Central place for backend base URLs. Import from here — no hardcoded URLs in components.

// The agent Lambda Function URL has no CORS headers yet, so in dev we route
// through Vite's same-origin proxy (see vite.config.ts) to sidestep it. In a
// production build, call it directly — this requires the Lambda Function
// URL's CORS config to allow the deployed origin.
export const AGENT_API_URL = import.meta.env.DEV
  ? '/agent-api'
  : 'https://kpqr75ugvaeq5oijot3b27d5sq0xagsw.lambda-url.eu-west-2.on.aws';

// Scoring API is consumed through lib/api.ts, which has its own dev-proxy vs.
// prod switch. Kept here for reference/consistency; not wired into that file.
export const SCORING_API_URL = 'https://9368clqa34.execute-api.eu-north-1.amazonaws.com';
