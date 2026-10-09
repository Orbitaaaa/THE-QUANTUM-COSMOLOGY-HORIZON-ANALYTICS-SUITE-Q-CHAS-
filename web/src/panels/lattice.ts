/** Lattice Shape Optimizer: placeholder until the tool is implemented (see web/README.md). */

import { mountStub } from './stub';

export const title = 'Lattice Shape Optimizer';

export const mount = mountStub({
  sources: ['[TOP] Sec. 6', '[TOP] Prop. E.2', '[TOP] Theorem H.2'],
  plan: [
    'Z(2) at unit covolume along the rhombohedral family (Fig. 2), computed live',
    'SC, BCC, FCC critical points with Hessian spectra versus the certified values',
    'Shape-flow trajectories and basin counts (FCC vs BCC) from the pipeline',
    '3D view of the fundamental domain for any basis',
  ],
  physics: [
    { module: 'physics/lattice.ts', provides: 'epsteinZ2, rhombohedralBasis, cubicBasis, lllReduce, latticeVectors' },
  ],
});
