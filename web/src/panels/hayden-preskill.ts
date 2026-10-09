/** Hayden-Preskill Recovery: placeholder until the tool is implemented (see web/README.md). */

import { mountStub } from './stub';

export const title = 'Hayden-Preskill Recovery';

export const mount = mountStub({
  sources: ['[BH] Sec. 11', '[BH] Theorem E(iv)', '[QG] Theorem R'],
  plan: [
    'Exactly unitary evaporation of an n-qubit black hole: S(R) against the Page envelope',
    'A diary thrown in early or late and the step at which it is recovered from the radiation',
    'Decoding infidelity against the 2^(1 + k − s) bound for Haar and ETH dynamics',
  ],
  physics: [
    { module: 'physics/page.ts', provides: 'decodingInfidelityBound, haydenPreskillEmissions' },
  ],
});
