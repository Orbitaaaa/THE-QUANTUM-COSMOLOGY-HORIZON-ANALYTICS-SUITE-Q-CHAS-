/** Manifold Selector: placeholder until the tool is implemented (see web/README.md). */

import { mountStub } from './stub';

export const title = 'Manifold Selector';

export const mount = mountStub({
  sources: ['[TOP] Sec. 1', '[TOP] Sec. 4', '[TOP] Sec. 9.3'],
  plan: [
    'The 18 complete flat 3-manifolds E1-E18: closed/open, orientability, holonomy',
    'Principles P1-P4 toggled on and off, with the theorem that excludes each manifold (Theorems A, C′, B-loc, D)',
    'The surviving topology (the 3-torus) and why each principle is needed',
    'Curved candidates (Poincaré dodecahedral, octahedral, Weeks) with the quoted Bayes factors',
  ],
});
