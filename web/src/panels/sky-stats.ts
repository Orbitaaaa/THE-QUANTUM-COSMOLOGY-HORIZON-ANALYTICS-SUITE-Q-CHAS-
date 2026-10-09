/** Sky Statistics: placeholder until the tool is implemented (see web/README.md). */

import { mountStub } from './stub';

export const title = 'Sky Statistics';

export const mount = mountStub({
  sources: ['[TOP] Sec. 8', '[TOP] Sec. 11.2'],
  plan: [
    'Observed D_ℓ (WMAP, Planck when available) against ΛCDM and torus models, 2 ≤ ℓ ≤ 30',
    'Quadrupole and S1/2 statistics with their ΛCDM p-values',
    'KL information versus d/D with the Pinsker bound, and the causal-ball ceiling',
    'Survey footprints (SPHEREx frames) on the sky',
  ],
});
