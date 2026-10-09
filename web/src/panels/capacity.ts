/** Crossed-Product Capacity: placeholder until the tool is implemented (see web/README.md). */

import { mountStub } from './stub';

export const title = 'Crossed-Product Capacity';

export const mount = mountStub({
  sources: ['[BH] Theorem G', '[QG] Theorem C', '[QG] Theorem U'],
  plan: [
    'The algebra ladder: type III₁ (no trace) → II∞ (crossed product) → II₁ (energy window)',
    'Window capacity ln tr Π against A/4G for any mass and window width, computed live',
    'De Sitter entropy 3π/GΛ of the measured cosmological constant and the area gap',
  ],
  physics: [
    { module: 'physics/capacity.ts', provides: 'capacityWindow, capacityCurve, deSitterEntropy, minimalArea, bitsPerArea' },
  ],
});
