/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      fontFamily: {
        // `sans` is the app default; `display` for headings; `mono` for IDs/numbers.
        sans: ['"IBM Plex Sans"', 'system-ui', 'sans-serif'],
        display: ['"Space Grotesk"', 'system-ui', 'sans-serif'],
        mono: ['"IBM Plex Mono"', 'ui-monospace', 'monospace'],
      },
      colors: {
        // Optional teal accent (GeneTraceAI identity). Not applied by default —
        // swap indigo-* for accent-* to adopt it. Kept here so it's one edit away.
        accent: {
          light: '#E6F1F1',
          DEFAULT: '#0B7C86',
          dark: '#075A62',
        },
      },
    },
  },
  plugins: [],
}