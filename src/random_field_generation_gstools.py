# generate random numbers
import gstools as gs
import numpy as np

########################################################################################################################
# shared utilities
def update_randomfield_use_autocc(auto_cc, old_rf, new_rf):
    # update using temporal auto (i.e., lag-1) correlation
    new_rf = old_rf * auto_cc + new_rf * np.sqrt(1 - auto_cc ** 2)
    return new_rf


def update_randomfield_use_crosscc(cross_cc, source_rf, target_rf):
    # source_rf: this random field will be used to adjust the target_rf
    target_rf = source_rf * cross_cc(1) + np.sqrt(1 - cross_cc**2) * target_rf
    return target_rf

########################################################################################################################
# random field generation using gstools
# not so good for large clen

def generate_Gaussian_latlon_random_field(lat, lon, clen, randseed):
    # generate a random field on geographical coordinates
    # Reference: https://geostat-framework.readthedocs.io/projects/gstools/en/stable/examples/08_geo_coordinates/00_field_generation.html
    # len_scale parameter: https://gmd.copernicus.org/articles/15/3161/2022/

    # Example from Reference:
    # model = gs.Gaussian(latlon=True, var=1, len_scale=777, rescale=gs.EARTH_RADIUS)
    # lat = lon = range(-80, 81)
    # srf = gs.SRF(model, seed=1234)
    # field = srf.structured((lat, lon))
    # srf.plot()

    model = gs.Gaussian(latlon=True, var=1, len_scale=clen, rescale=gs.EARTH_RADIUS) # var: sigma
    srf = gs.SRF(model, seed=randseed)
    field = srf.structured((lat, lon))

    return field


def gs_better_seeds(num, masterseed=0):
    seed0 = gs.random.MasterRNG(masterseed)
    seeds = np.zeros(num, dtype=int)
    for i in range(num):
        seeds[i] = seed0()
    return seeds

# random field generation using gstools
########################################################################################################################





