/** Page Curve: placeholder until the tool is implemented (see web/README.md). */

import { mountStub } from './stub';

export const title = 'Page Curve';

export const mount = mountStub({
  sources: ['[BH] Theorems C-F', '[BH] Fig. 1'],
  plan: [
    "Hawking's curve β(1 − m²) against the capacity m² and the Page law min{S_H, A/4G}",
    'The island minimal entropy, the Page time and the 2-design lower bound for finite S₀',
    'Interior volume growing while the capacity shrinks (Theorem F)',
  ],
  physics: [
    { module: 'physics/page.ts', provides: 'hawkingEntropy, pageCurve, pageCurveLower, islandMinEntropy, interiorVolume, curves' },
  ],
});
