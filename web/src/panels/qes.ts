/** QES & Atmosphere: placeholder until the tool is implemented (see web/README.md). */

import { mountStub } from './stub';

export const title = 'QES & Atmosphere';

export const mount = mountStub({
  sources: ['[BH] Sec. 9', '[BH] Secs. 10, 15'],
  plan: [
    'Quantum extremal surface location for any mass, central charge and collection radius, solved live',
    'The [BH] Sec. 9 table reproduced row by row, with the depth and lag laws',
    'The entangled quantum atmosphere: lattice shares and the continuum profile',
  ],
  physics: [
    { module: 'physics/qes.ts', provides: 'solveQes, depthFit, lagFit, PAPER_TABLE' },
  ],
});
