/** Ghost-Twin Tracker: placeholder until the tool is implemented (see web/README.md). */

import { mountStub } from './stub';

export const title = 'Ghost-Twin Tracker';

export const mount = mountStub({
  sources: ['[TOP] Sec. 7', '[TOP] Corollary F.1'],
  plan: [
    'Matched-circle pairs for any torus basis, size d/D and orientation, on the real CMB map',
    'Click a point on a circle to see its twin n′ = n + g/χ* and the chord through the observer',
    'Ghost images of a source at any distance, inside and beyond the horizon',
    'The [TOP] Sec. 9 best-fit FCC orientation (not significant) as a preset',
  ],
  physics: [
    { module: 'physics/circles.ts', provides: 'circlePairs, matchedPoint, circlePoints, ghostImages' },
    { module: 'lib/sky.ts', provides: 'globe and Mollweide views with circle, point, chord and footprint overlays' },
  ],
});
