"""Registry of the public datasets Q-CHAS ingests, with mirrors in priority order.

Nothing here is simulated: every entry is an official data product of a NASA, ESA or
LIGO-Virgo-KAGRA release. Mirrors are tried in order, so the same registry works on a
laptop, in GitHub Actions, and in sandboxes where only some archives are reachable.
"""

from __future__ import annotations

from dataclasses import dataclass, field

LAMBDA_S3 = "https://nasa-lambda.s3.amazonaws.com/wmap/dr5"
LAMBDA_WEB = "https://lambda.gsfc.nasa.gov/data/map/dr5"
IRSA_PLANCK = "https://irsa.ipac.caltech.edu/data/Planck/release_3"
PLA_AIO = "https://pla.esac.esa.int/pla/aio/product-action"
LOSC_TUTORIAL = "https://raw.githubusercontent.com/losc-tutorial/LOSC_Event_tutorial/master"
GWOSC = "https://gwosc.org"
SPHEREX_S3 = "https://nasa-irsa-spherex.s3.amazonaws.com"


@dataclass(frozen=True)
class Dataset:
    key: str
    filename: str
    mirrors: tuple[str, ...]
    mission: str
    description: str
    citation: str
    min_bytes: int = 1
    tags: tuple[str, ...] = field(default_factory=tuple)


DATASETS: dict[str, Dataset] = {}


def _register(ds: Dataset) -> None:
    DATASETS[ds.key] = ds


# --- WMAP 9-year (NASA LAMBDA) -----------------------------------------------------
_register(Dataset(
    key="wmap_ilc",
    filename="wmap_ilc_9yr_v5.fits",
    mirrors=(f"{LAMBDA_S3}/dfp/ilc/wmap_ilc_9yr_v5.fits", f"{LAMBDA_WEB}/dfp/ilc/wmap_ilc_9yr_v5.fits"),
    mission="WMAP",
    description="WMAP 9-year Internal Linear Combination CMB temperature map (Nside 512, 1 deg FWHM, mK)",
    citation="Bennett et al. 2013, ApJS 208, 20",
    min_bytes=20_000_000,
    tags=("cmb", "map"),
))
_register(Dataset(
    key="wmap_tt",
    filename="wmap_tt_spectrum_9yr_v5.txt",
    mirrors=(f"{LAMBDA_S3}/dcp/spectra/wmap_tt_spectrum_9yr_v5.txt", f"{LAMBDA_WEB}/dcp/spectra/wmap_tt_spectrum_9yr_v5.txt"),
    mission="WMAP",
    description="WMAP 9-year combined TT power spectrum D_l (uK^2)",
    citation="Hinshaw et al. 2013, ApJS 208, 19",
    min_bytes=50_000,
    tags=("cmb", "spectrum"),
))
_register(Dataset(
    key="wmap_kq75",
    filename="wmap_temperature_kq75_analysis_mask_r9_9yr_v5.fits",
    mirrors=(
        f"{LAMBDA_S3}/ancillary/masks/wmap_temperature_kq75_analysis_mask_r9_9yr_v5.fits",
        f"{LAMBDA_WEB}/ancillary/masks/wmap_temperature_kq75_analysis_mask_r9_9yr_v5.fits",
    ),
    mission="WMAP",
    description="WMAP 9-year KQ75 temperature analysis mask (Galactic plane + point sources)",
    citation="Bennett et al. 2013, ApJS 208, 20",
    min_bytes=20_000_000,
    tags=("cmb", "mask"),
))

# --- Planck PR3 (ESA, mirrored at NASA/IPAC IRSA) ----------------------------------------
_register(Dataset(
    key="planck_commander",
    filename="COM_CMB_IQU-commander_2048_R3.00_full.fits",
    mirrors=(
        f"{IRSA_PLANCK}/all-sky-maps/maps/component-maps/cmb/COM_CMB_IQU-commander_2048_R3.00_full.fits",
        f"{PLA_AIO}?MAP.MAP_ID=COM_CMB_IQU-commander_2048_R3.00_full.fits",
    ),
    mission="Planck",
    description="Planck 2018 (PR3) Commander CMB map, inpainted temperature (Nside 2048, K_CMB), 1.6 GB",
    citation="Planck Collaboration 2020, A&A 641, A4",
    min_bytes=1_000_000_000,
    tags=("cmb", "map", "large"),
))
_register(Dataset(
    key="planck_tt",
    filename="COM_PowerSpect_CMB-TT-full_R3.01.txt",
    mirrors=(
        f"{IRSA_PLANCK}/ancillary-data/cosmoparams/COM_PowerSpect_CMB-TT-full_R3.01.txt",
        f"{PLA_AIO}?COSMOLOGY.FILE_ID=COM_PowerSpect_CMB-TT-full_R3.01.txt",
    ),
    mission="Planck",
    description="Planck 2018 TT power spectrum D_l with asymmetric errors (Commander for l < 30)",
    citation="Planck Collaboration 2020, A&A 641, A5",
    min_bytes=50_000,
    tags=("cmb", "spectrum"),
))

# --- LIGO-Virgo strain (official LIGO open-data tutorial release of GWOSC files) --------
_register(Dataset(
    key="losc_events",
    filename="BBH_events_v3.json",
    mirrors=(f"{LOSC_TUTORIAL}/BBH_events_v3.json",),
    mission="LIGO",
    description="LIGO open-data tutorial event list (GW150914, LVT151012, GW151226, GW170104)",
    citation="LIGO Scientific Collaboration, GWOSC tutorial release",
    min_bytes=1_000,
    tags=("gw", "catalog"),
))
_LOSC_FILES = {
    "GW150914": ("H-H1_LOSC_4_V2-1126259446-32.hdf5", "L-L1_LOSC_4_V2-1126259446-32.hdf5"),
    "LVT151012": ("H-H1_LOSC_4_V2-1128678884-32.hdf5", "L-L1_LOSC_4_V2-1128678884-32.hdf5"),
    "GW151226": ("H-H1_LOSC_4_V2-1135136334-32.hdf5", "L-L1_LOSC_4_V2-1135136334-32.hdf5"),
    "GW170104": ("H-H1_LOSC_4_V1-1167559920-32.hdf5", "L-L1_LOSC_4_V1-1167559920-32.hdf5"),
}
for _event, _files in _LOSC_FILES.items():
    for _det, _fname in zip(("H1", "L1"), _files):
        _register(Dataset(
            key=f"strain_{_event.lower()}_{_det.lower()}",
            filename=_fname,
            mirrors=(f"{LOSC_TUTORIAL}/{_fname}",),
            mission="LIGO",
            description=f"{_event} {_det} strain, 32 s at 4096 Hz (GWOSC release via the LIGO open-data tutorial)",
            citation="LIGO Scientific Collaboration and Virgo Collaboration; GWOSC",
            min_bytes=500_000,
            tags=("gw", "strain"),
        ))
_register(Dataset(
    key="gwtc",
    filename="gwtc_allevents.json",
    mirrors=(f"{GWOSC}/eventapi/json/allevents/", f"{GWOSC}/eventapi/json/GWTC/"),
    mission="LIGO-Virgo-KAGRA",
    description="GWOSC event API: every released compact-binary event with source parameters",
    citation="LIGO-Virgo-KAGRA, GWOSC (Abbott et al. 2023, ApJS 267, 29)",
    min_bytes=10_000,
    tags=("gw", "catalog"),
))


def datasets_with_tag(tag: str) -> list[Dataset]:
    return [d for d in DATASETS.values() if tag in d.tags]
