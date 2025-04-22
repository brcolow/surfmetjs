// vitest.config.js
import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    environment: 'node', // If you want to test DOM stuff, first install `npm install --save-dev jsdom` and use 'jsdom' instead of 'node'
    globals: true, // Enables global APIs like `describe`, `it`, etc.
    // Should you want to generate coverage reports, first install `npm install --save-dev @vitest/coverage-v8`
    // and uncomment the following:
    // coverage: {
    //   provider: 'v8', // Use V8 coverage
    //   reporter: ['text', 'json', 'html'], // Coverage report formats
    // },
  },
});
