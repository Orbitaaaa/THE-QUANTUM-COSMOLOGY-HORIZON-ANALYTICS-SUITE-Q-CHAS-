/// <reference types="vitest/config" />
import { defineConfig, type Plugin } from 'vite';

// Some hosts reject a literal U+FFFD in served scripts; Plotly's bundled UTF-8 decoder has a few
// inside string literals. Emit the equivalent escape so every chunk is plain, portable text.
const escapeReplacementChar = (): Plugin => ({
  name: 'escape-replacement-char',
  enforce: 'post',
  // after minification, which would otherwise turn the escape back into the literal character
  generateBundle(_options, bundle) {
    for (const file of Object.values(bundle)) {
      if (file.type === 'chunk' && file.code.includes('\uFFFD')) file.code = file.code.replace(/\uFFFD/g, '\\ufffd');
    }
  },
});

// base './' makes every asset URL relative, so the built site works from a GitHub Pages
// subpath (https://<user>.github.io/<repo>/) as well as from the domain root.
export default defineConfig({
  base: './',
  plugins: [escapeReplacementChar()],
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
