/// <reference types="vitest/config" />
import { defineConfig } from 'vite';

// base './' makes every asset URL relative, so the built site works from a GitHub Pages
// subpath (https://<user>.github.io/<repo>/) as well as from the domain root.
export default defineConfig({
  base: './',
  build: {
    target: 'es2022',
    chunkSizeWarningLimit: 4000,
    rollupOptions: {
      output: {
        manualChunks: {
          plotly: ['plotly.js-cartesian-dist-min'],
          three: ['three'],
        },
      },
    },
  },
  test: {
    include: ['tests/**/*.test.ts'],
    environment: 'node',
    maxWorkers: 2,
  },
});
