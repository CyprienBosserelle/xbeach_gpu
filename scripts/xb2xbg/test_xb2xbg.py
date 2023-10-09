"""
Unit tests for xb2xbg.py

The core functions responsible for XB & XBG parameter/variable extraction
and transformation have test cases, using an example XB parameter file
as input. The emphasis here is for anyone changing the script to be
alerted to the effects of a change to transformation code, and to be
able to refactor code and make other non-functionality changing modifications
with confidence.

Note that the generate() function is not tested here, nor are parameter and
output variables document strings.
"""

import os

from xb2xbg import manual_url
from xb2xbg import extract_xbg_params, ParamObjectBuilder
from xb2xbg import read_xbeach_params
from xb2xbg import transform_outvars
from xb2xbg import ParamTransformer

# Example XB parameters file name
EXAMPLE_XB_PARAMS = "example-XBeach-params.txt"

# Expected results
expected_xb_params = {
    'xori': 381150.0, 'yori': 6414200.0, 'alfa': 33.0, 'dx': 5.0, 'dy': 5.0,
    'nx': 192, 'ny': 256,
    'projection=+proj=utm +zone=50 +south +ellps=GRS80 +units': 'm', 
    'posdwn': -1, 'depfile': 'DoT_2009Lidar_Crop_MGAz50_mAHD.z',
    'front': 1, 'back': 2, 'left': 0, 'right': 0, 'rugdepth': 0.011, 
    'tunits': 'seconds since 2000-01-01 00:00:00', 'instat': 41, 'break': 1, 
    'scheme': 1, 'order': 1, 'leftwave': 'wavecrest', 'rightwave': 'wavecrest', 
    'random': 1, 'zs0file': 'zs0input.txt', 'tidelen': 121, 'tideloc': 1, 
    'zs0': 0.0, 'hmin': 0.01, 'wci': 0, 'alpha': 1, 'delta': 0.0, 'n': 10, 
    'roh': 1025, 'g': 9.81, 'thetamin': -80, 'thetamax': 80, 'dtheta': 10.0, 
    'beta': 0.1, 'roller': 1, 'gamma': 0.55, 'gammax': 1.0, 'bcfile': 'jonswap1.inp', 
    'sedtrans': 0, 'morfac': 0.0, 'morphology': 0, 'cf': 0.01, 'paulrevere': 0,
    'eps': 0.01, 'epsi': 0.001, 'tstart': 0, 'tint': 3600.0, 'tstop': 820800,
    'tintm': 3600.0, 'CFL': 0.8, 'nuh': 0.5, 'nuhfac': 0, 'umin': 0.1, 'smag': 0,
    'oldhu': 1, 'outputformat': 'netcdf', 'ncfilename': 'outputs/xboutput.nc',
    'nglobalvar': 15, 'nmeanvar': 10,
    'global_vars': ['H', 'thetamean', 'hh', 'u', 'v', 'D', 'R', 'k', 'E', 'ue',
                    've', 'urms', 'zb', 'zs', 'Qb'],
    'mean_vars': ['H', 'u', 'v', 'ue', 've', 'R', 'D', 'urms', 'k', 'Qb']
}

expected_initial_xbg_param_values = { \
    "General Parameters":
        [('GPUDEVICE', 0), ('dx', 0.0), ('flow', 1), ('grdalpha', 0.0),
         ('morphology', 0), ('name', 'GeneralParameters'), ('nx', 0), ('ny', 0),
         ('posdown', 1), ('sedtrans', 0), ('swave', 1)],
    "Flow Parameters":
        [('Cd', 0.002), ('cf', 0.01), ('cfreef', 0.01), ('cfsand', 0.01),
         ('eps', 0.01), ('fc', 0.0), ('g', 9.81), ('hwci', 0.1), ('lat', 0.0),
         ('name', 'FlowParameters'), ('nuh', 1.0), ('nuhfac', 1.0), ('rho', 1025.0),
         ('smag', 0.3), ('usesmago', 1), ('wci', 0.0)],
    "Wave Parameters":
        [('alpha', 1.0), ('beta', 0.15), ('breakmodel', 1), ('dtheta', None),
         ('fw', 0.001), ('fwreef', 0.001), ('fwsand', 0.001), ('gamma', 0.6),
         ('gammax', 2.0), ('n', 8.0), ('name', 'WaveParameters'), ('ntheta', 1),
         ('roller', 1), ('thetamax', 90), ('thetamin', -90)],
    "Sediment Parameters":
        [('D50', 0.00038), ('D90', 0.00053), ('bed', 1.0), ('drydzmax', 1.0),
         ('facas', 0.2), ('facsk', 0.2), ('maxslpchg', 0.01), ('morfac', 1.0),
         ('name', 'SedimentParameters'), ('por', 0.4), ('rhosed', 2650.0),
         ('sus', 1.0), ('wetdzmax', 2.0), ('wws', 0.0509)],
    "Files":
        [('SedThkfile', ''), ('TSOfile', ''), ('TSnode', ''), ('bathy', ''),
         ('depfile', ''), ('name', 'Files'),
         ('outfile', ''), ('outvars', ''),
         ('slbndfile', ''), ('wavebndfile', ''),
         ('windbndfile', '')],
    "Time Keeping":
        [('CFL', 0.7), ('endtime', 0.0), ('name', 'TimeKeeping'),
         ('outputtimestep', 0.0), ('sedstart', 3600.0)],
    "Wave boundary variables":
        [('dtbc', 1.0), ('fcutoff', 0.0), ('name', 'Waveboundaryvariables'),
         ('nmax', 0.8), ('random', 0), ('rtlength', 3600.0), ('sprdthr', 0.8),
         ('wavebndtype', 2)]
}

expected_final_parameter_values = { \
    "General Parameters":
        [('GPUDEVICE', 0), ('dx', 5.0), ('flow', 1), ('grdalpha', 33.0),
         ('morphology', 0), ('name', 'GeneralParameters'), ('nx', 193), ('ny', 257),
         ('posdown', 0), ('sedtrans', 0), ('swave', 1)],
    "Flow Parameters":
        [('Cd', 0.002), ('cf', 0.01), ('cfreef', 0.01), ('cfsand', 0.01),
         ('eps', 0.01), ('fc', 0.0), ('g', 9.81), ('hwci', 0.1), ('lat', 0.0),
         ('name', 'FlowParameters'), ('nuh', 0.5), ('nuhfac', 0), ('rho', 1025.0),
         ('smag', 0.3), ('usesmago', 0), ('wci', 0)],
    "Wave Parameters":
        [('alpha', 1), ('beta', 0.1), ('breakmodel', 1), ('dtheta', 10.0),
         ('fw', 0.001), ('fwreef', 0.001), ('fwsand', 0.001), ('gamma', 0.55),
         ('gammax', 1.0), ('n', 10), ('name', 'WaveParameters'), ('ntheta', 1),
         ('roller', 1), ('thetamax', 80), ('thetamin', -80)],
    "Sediment Parameters":
        [('D50', 0.00038), ('D90', 0.00053), ('bed', 1.0), ('drydzmax', 1.0),
         ('facas', 0.2), ('facsk', 0.2), ('maxslpchg', 0.01), ('morfac', 0.0),
         ('name', 'SedimentParameters'), ('por', 0.4), ('rhosed', 2650.0),
         ('sus', 1.0), ('wetdzmax', 2.0), ('wws', 0.0509)],
    "Files":
        [('SedThkfile', ''), ('TSOfile', ''), ('TSnode', ''), ('bathy', ''),
         ('depfile', 'DoT_2009Lidar_Crop_MGAz50_mAHD.dep'), ('name', 'Files'),
         ('outfile', 'outputs/xboutput.nc'), ('outvars', ''),
         ('slbndfile', 'zs0input.txt'), ('wavebndfile', 'jonswap1-xbg.inp'),
         ('windbndfile', '')],
    "Time Keeping":
        [('CFL', 0.8), ('endtime', 820800), ('name', 'TimeKeeping'),
         ('outputtimestep', 3600.0), ('sedstart', 3600.0)],
    "Wave boundary variables":
        [('dtbc', 1.0), ('fcutoff', 0.0), ('name', 'Waveboundaryvariables'),
         ('nmax', 0.8), ('random', 1), ('rtlength', 3600.0), ('sprdthr', 0.8),
         ('wavebndtype', 2)]
}

expected_transformed_parameter_name_mappings = { \
    "General Parameters":
        {'alfa': 'grdalpha', 'dx': 'dx', 'nx': 'nx', 'ny': 'ny',
         'posdwn': 'posdown', 'sedtrans': 'sedtrans', 'morphology': 'morphology'},
    "Flow Parameters":
        {'wci': 'wci', 'g': 'g', 'cf': 'cf', 'eps': 'eps', 'nuh': 'nuh',
         'nuhfac': 'nuhfac', 'smag': 'usesmago'},
    "Wave Parameters":
        {'alpha': 'alpha', 'n': 'n', 'thetamin': 'thetamin',
         'thetamax': 'thetamax', 'dtheta': 'dtheta', 'beta': 'beta', 'roller': 'roller', 'gamma': 'gamma', 'gammax': 'gammax'},
    "Sediment Parameters":
        {'morfac': 'morfac'},
    "Files":
        {'depfile': 'depfile', 'zs0file': 'slbndfile', 'bcfile': 'wavebndfile', 'ncfilename': 'outfile'},
    "Time Keeping":
        {'tint': 'outputtimestep', 'tstop': 'endtime', 'CFL': 'CFL'},
    "Wave boundary variables":
        {'random': 'random'}
}

expected_initial_xbg_outvar_names = \
    ['C', 'Cmean', 'D', 'DR', 'E', 'Fx', 'Fy', 'H', 'Hmean', 'R', 'c', 'ceqsg',
     'cfm', 'cg', 'cgx', 'cgy', 'ctheta', 'cx', 'cy', 'dhdx', 'dhdy', 'dudx',
     'dudy', 'dvdx', 'dvdy', 'dzb', 'dzsdt', 'dzsdx', 'dzsdy', 'ee', 'fwm', 'hh',
     'hhmean', 'hu', 'hum', 'hv', 'hvm', 'k', 'kh', 'kturb', 'rolthick', 'rr',
     'sigm', 'sinh2kh', 'stdep', 'thetamean', 'ududx', 'udvdx', 'ueu', 'urms',
     'ust', 'uu', 'uumean', 'uv', 'vdudy', 'vdvdy', 'vev', 'vmageu', 'vmagev',
     'vu', 'vv', 'vvmean', 'wci', 'zb', 'zs', 'zsmax', 'zsmean']

expected_transformed_xbg_outvar_names = \
    ['D', 'E', 'H', 'Hmean', 'hh', 'k', 'thetamean', 'urms', 'zb', 'zs']

expected_transformed_xbg_outvar_names_all_known = \
    ['D', 'E', 'H', 'Hmean', 'R', 'hh', 'k', 'thetamean', 'urms', 'zb', 'zs']


def test_read_xbeach_params():
    xb_params = read_xbeach_params(full_path_from_filename(EXAMPLE_XB_PARAMS))

    assert len(xb_params) == 65
    assert xb_params == expected_xb_params
    

def test_extract_xbg_params():
    builder = ParamObjectBuilder()
    extract_xbg_params(manual_url(), builder)

    for xbg_params in builder.param_objs:
        assert xbg_params.values() == \
            expected_initial_xbg_param_values[xbg_params.section_name]

    assert sorted(builder.outvars.keys()) == expected_initial_xbg_outvar_names


def test_transform_outvars():
    xb_params = read_xbeach_params(full_path_from_filename(EXAMPLE_XB_PARAMS))
    builder = ParamObjectBuilder()
    extract_xbg_params(manual_url(), builder)
    xbg_outvars = transform_outvars(xb_params, builder, verbose=False,
                                    use_all_known_outvars=False)

    assert sorted(xbg_outvars) == expected_transformed_xbg_outvar_names


def test_transform_all_known_outvars():
    xb_params = read_xbeach_params(full_path_from_filename(EXAMPLE_XB_PARAMS))
    builder = ParamObjectBuilder()
    extract_xbg_params(manual_url(), builder)
    xbg_outvars = transform_outvars(xb_params, builder, verbose=False,
                                    use_all_known_outvars=True)

    assert sorted(xbg_outvars) == expected_transformed_xbg_outvar_names_all_known


def test_transform_params():
    builder = ParamObjectBuilder()
    extract_xbg_params(manual_url(), builder)
    xb_params = read_xbeach_params(full_path_from_filename(EXAMPLE_XB_PARAMS))

    # transform each initial XBG parameter group in the presence
    # of XB input parameters 
    for xbg_params in builder.param_objs:
        transformer = ParamTransformer(xb_file_root="data", xbg_output_dir=None,
                                       use_defaults=False,
                                       xbg_params=xbg_params, xb_params=xb_params,
                                       directional_spread_coefficient=400,
                                       peak_enhancement_factor=3.3,
                                       verbose=False)
        
        xb_params_transformed = transformer.transform_params()

        assert xb_params_transformed == \
            expected_transformed_parameter_name_mappings[xbg_params.section_name]

        assert xbg_params.values() == \
            expected_final_parameter_values[xbg_params.section_name]

    # clean up    
    os.unlink("jonswap1-xbg.inp")
    os.unlink("zs0input.txt")
    os.unlink("DoT_2009Lidar_Crop_MGAz50_mAHD.dep")


# Helpers
     
def full_path_from_filename(filename: str) -> str:
    path = None
    for root, dir, files in os.walk("."):
        if filename in files:
            path = os.path.join(root, filename)
    return path
