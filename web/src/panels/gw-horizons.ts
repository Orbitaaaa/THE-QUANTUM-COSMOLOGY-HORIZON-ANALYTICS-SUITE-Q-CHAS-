/** GW Horizons (LIGO): placeholder until the tool is implemented (see web/README.md). */

import { mountStub } from './stub';

export const title = 'GW Horizons (LIGO)';

export const mount = mountStub({
  sources: ['[QG] Sec. 17'],
  plan: [
    'Whitened, band-passed H1 and L1 strain of GW150914, LVT151012, GW151226, GW170104',
    'Ringdown fits against the remnant mass and spin; echo delays',
    'Remnant entropy, Hawking temperature, evaporation and Page times',
  ],
});
