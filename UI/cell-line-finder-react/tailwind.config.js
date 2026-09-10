/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      fontFamily: {
        // AstraZeneca house type: IBM Plex Sans throughout, IBM Plex Mono for IDs/numbers.
        sans: ['"IBM Plex Sans"', 'system-ui', 'sans-serif'],
        display: ['"IBM Plex Sans"', 'system-ui', 'sans-serif'],
        mono: ['"IBM Plex Mono"', 'ui-monospace', 'monospace'],
      },
      colors: {
        // AstraZeneca mulberry (brand ≈ #830051 at 600). Primary accent across the app.
        mulberry: {
          50: '#fbeef4',
          100: '#f7dbe9',
          200: '#eeb3d0',
          300: '#e186b3',
          400: '#cc5391',
          500: '#ad2870',
          600: '#830051',
          700: '#6d0044',
          800: '#580037',
          900: '#43012a',
        },
        // Dark plum used for the sidebar / rails.
        plum: {
          DEFAULT: '#2a1b28',
          accent: '#3a2838',
          border: '#41303e',
          muted: '#9c8494',
        },
        // Confidence tiers (AstraZeneca-aligned hues).
        tier: {
          high: '#2f8f5b',
          medium: '#9a7b23',
          low: '#78767f',
          context: '#3a72a8',
        },
        // Optional teal accent (GeneTraceAI identity). One edit away if needed.
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
