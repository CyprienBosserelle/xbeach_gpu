"""
Convert XBeach params.txt file to XBG_param.txt

Assumes Python 3.7 or higher and param library
"""

import argparse
import numpy as np
import os
import pandas as pd
import re
import shutil
import sys
import textwrap

from collections import OrderedDict

from dataclasses import dataclass

from io import TextIOWrapper
from urllib.request import urlopen

from typing import Any, Dict, List, Tuple, Union

###############################
# XBeach-GPU parameters class #
###############################

@dataclass
class XBGParam:
    """XBG parameter data class"""
    name: str
    value: Union[int, float, str] = None
    default: Union[int, float, str] = None
    doc: str = ""

    def name_value(self) -> Tuple[str, Union[int, float, str]]:
        """Return name-value pair"""
        return (self.name, self.value)
        

class XBGParams:
    """
    Collection of XBG parameters corresponding to
    a particular XBG user manual section.
    """
    def __init__(self, section_name: str):
        self.section_name = section_name
        self.params = OrderedDict()
        # add section name as a parameter
        self.params["name"] = XBGParam("name", self.section_name.replace(" ", ""))

    @staticmethod
    def title():
        """Header comment lines for the XBG output file"""
        return "{}\n{}".format("\n### ".\
                               join(["#"*73,
                                     "XBGPU parameter settings input file" + \
                                     " "*30 + " ###",
                                     "Note: for binary parameters such as GPUDEVICE & "
                                     "flow, 1=YES, 0=NO ###"]),
                               "#"*73)

    def decorated_section(self):
        return "{0} {1} {2}".format("#"*3, self.section_name, "#"*(68-len(self.section_name)))

    def put(self, name: str, value: XBGParam):
        """Associate a XBG parameter object with a parameter name"""
        self.params[name] = value

    def set(self, name: str, value: Union[int, float, str]):
        """Set the value of an existing parameter"""
        self.params[name].value = value

    def get(self, name: str) -> XBGParam:
        """Get the XBG parameter correspoding to a parameter name"""
        return self.params[name]

    def values(self) -> List[Tuple[str, XBGParam]]:
        return [self.params[key].name_value() for key in sorted(self.params)]


#############################
# Parameter builder classes #
#############################

class ParamObjectBuilder:
    """
    This class allows an object to be built that captures
    XBG parameters and output variables.
    """
    def __init__(self):
        self.types = {int: "Integer", float: "Number", str: "String"}
        self.param_name = None
        self.param_type = None
        self.param_default = None
        self.param_doc = None
        self.param_objs: List[XBGParams] = []
        self.outvars = {}

    def add_param(self, section: str):
        self.param_objs.append(XBGParams(section_name=section))

    def set_param_name(self, name: str):
        self.param_name = name

    def set_param_type_and_default(self, ptype: str, val: Union[int, float, str]):
        self.param_type = ptype
        self.param_default = val

    def set_param_doc(self, doc: str):
        self.param_doc = doc

    def end_param(self):
        self.param_objs[-1].put(self.param_name,
                                XBGParam(self.param_name, value=self.param_default,
                                default=self.param_default, doc=self.param_doc))

    def add_outvar(self, var: str, doc: str):
        self.outvars[var] = doc


def extract_xbg_params(source: str, builder: ParamObjectBuilder):
    """
    Extract XBeach-GPU parameter information from XBeach GPU User
    Manual HTML file or URL. See also manual_url().

    :param source: XBeach GPU User Manual HTML file or URL
    :param builder: subclass of ParamObjectBuilder that creates a XBG parameter
                    product or generates output
    """
    def xbg_param_lines(source: str):
        if source[0:4] == "http":
            return TextIOWrapper(urlopen(source), encoding="utf-8")
        else:
            return open(source, encoding="utf-8").readlines()

    section_pattern = re.compile(r"^\s*<h3>([^<]+)</h3>\s*$")
    param_col_pattern = re.compile(r"^\s*<td>([^<]+)</td>\s*$")
    outvar_pattern = re.compile(r"^\s*<li><p><strong>\s*(\w+)</strong>:(.*)</p></li>\s*$")

    state = "SECTION"
    for line in xbg_param_lines(source):
        if state == "SECTION":
            if line.find("<h1>Output</h1>") != -1:
                state = "OUTVARS"
            match = section_pattern.match(line)
            if match is not None:
                section = match.group(1)
                builder.add_param(section)
                state = "PARAM"
                td_count = 1
        elif state == "PARAM":
            if line.find("</tbody>") != -1:
                state = "SECTION"
            else:
                match = param_col_pattern.match(line)
                if match is not None:
                    text = match.group(1)
                    if td_count == 1:
                        builder.set_param_name(text.strip())
                        td_count += 1
                    elif td_count == 2:
                        try:
                            builder.set_param_type_and_default(int, int(text))
                        except ValueError:
                            try:
                                builder.set_param_type_and_default(float, float(text))
                            except ValueError:
                                if text in ["[]", '""']:
                                    text = ""
                                    builder.set_param_type_and_default(str, text)
                                elif text == "N/A":
                                    text = None
                                    # assume number; see dtype in user manual
                                    builder.set_param_type_and_default(float, text)
                        td_count += 1
                    else:
                        builder.set_param_doc(text.strip())
                        td_count += 1

                    if td_count > 3:
                        td_count = 1
                        builder.end_param()
        elif state == "OUTVARS":
            match = outvar_pattern.match(line)
            if match is not None:
                builder.add_outvar(match.group(1), match.group(2))
        else:
            raise ValueError("Unknown state encountered while processing XBG parameter file")


##############################################
# XBeach to XBeach-GPU parameter transformer #
##############################################

class ParamTransformer:
    """
    Transforms XB to XBG parameters for a particular XBG parameter
    section, e.g. Wave Parameters, captured by xbg_params
    """
    def __init__(self, xb_file_root: str, xbg_output_dir: str, use_defaults: bool,
                 xbg_params: XBGParams, xb_params: Dict[str, Any],
                 directional_spread_coefficient: int,
                 peak_enhancement_factor: float, verbose: bool):
        """
        :param xb_file_root: root directory for XB input files; may be None
        :param xbg_output_dir: directory for XBG output files; may be None
        :param use_defaults: use XBG parameter defaults instead of XB-derived values
        :param xbg_params: object containing a section of XBG parameters
        :param xb_params: dictionary of XB parameters
        :param directional_spread_coefficient: directional spread coeffient for jonswap
        :param peak_enhancement_factor: Peak enhancement factor for jonswap
        :param verbose: verbosity flag
        """
        self.xb_file_root = xb_file_root
        self.xbg_output_dir = xbg_output_dir
        self.use_defaults = use_defaults
        self.xbg_params = xbg_params
        self.xb_params = xb_params
        self.directional_spread_coefficient = directional_spread_coefficient
        self.peak_enhancement_factor = peak_enhancement_factor
        self.verbose = verbose

        # simple XB -> XBG parameter name mapping...
        self.xb2xbg_names = {
            "alfa": "grdalpha", "break": "breakmodel", "bcfile": "wavebndfile",
            "dzmax": "maxslpchg", "dryslp": "drydzmax", "morstart": "sedstart",
            "ncfilename": "outfile", "ne_layer": "SedThkfile", "posdwn": "posdown",
            "rt": "rtlength", "smag": "usesmago", "tint": "outputtimestep",
            "tintg": "outputtimestep", "tstop": "endtime", "wetslp": "wetdzmax",
            "windfile": "windbndfile", "ws": "wws", "zs0file": "slbndfile"
        }
        # where something more is needed to transform from XB to XBG variable...
        self.xb2xbg_transformers = {
            "depfile": self.xb2xbgfile,
            "breakmodel": self.breakermodel,
            "nx": self.inc_by_1,
            "ny": self.inc_by_1,
            "posdown": self.posdown,
            "slbndfile": self.xb2xbgfile,
            "wavebndfile": self.xb2xbgfile,
            "windbndfile": self.xb2xbgfile
        }

        self.params_transformed = {}

    def transform_params(self) -> Dict[str, str]:
        """Transform XB parameters for use with XBG.
           :return: a mapping from XB to XBG parameters transformed
        """
        if not self.use_defaults:
            if self.verbose:
                print("\nProcessing '{}'...".format(self.xbg_params.section_name))
            xbg_names = [param[0] for param in self.xbg_params.values()]
            for xb_name in self.xb_params:
                transformed = self.transform_param(xb_name, xbg_names)
                if transformed:
                    if xb_name in self.xb2xbg_names:
                        target_name = self.xb2xbg_names[xb_name]
                    else:
                        target_name = xb_name
                    self.params_transformed[xb_name] = target_name
        return self.params_transformed

    def transform_param(self, xb_name: str, xbg_names: List[str]) -> bool:
        """Transform a single XB parameter for use with XBG.
           This may involve a name change and/or a value change.
           :return: whether a transformation happened
        """
        try:
            transformed = False
            xbg_name = None

            if xb_name in self.xb2xbg_names:
                # XB and XBG names are different.
                xbg_name = self.xb2xbg_names[xb_name]
                # The XBG name may not be in this particular parameters section.
                # If not, ignore.
                if xbg_name not in xbg_names:
                    xbg_name = None
            elif xb_name in xbg_names:
                # XB and XBG names are the same.
                xbg_name = xb_name

            if xbg_name is not None:
                # XB parameter corresponds to an XBG parameter.
                # If a transformer function exists, call that.
                # Otherwise set the XBG parameter to the XB parameter's
                # value.
                if xbg_name in self.xb2xbg_transformers:
                    transformed = self.xb2xbg_transformers[xbg_name](xb_name, xbg_name)
                else:
                    self.xbg_params.set(xbg_name, self.xb_params[xb_name])
                    transformed = True
        except Exception as ex:
            print(str(ex), file=sys.stderr)

        return transformed

    def breakermodel(self, xb_name: str, xbg_name: str) -> bool:
        transformed = False
        breakmodel_val_map = {1: "roelvink (same as XB roelvink2)", 2: "breakmodel"}
        xb_val = self.xb_params[xb_name]
        if xb_val in [2, "baldock"]:
            self.xbg_params.set(xbg_name, 2)
            transformed = True
        elif xb_val in [3, "roelvink2"]:
            self.xbg_params.set(xbg_name, 1)
            transformed = True
        elif self.verbose:
            xbg_default = self.xbg_params.get(xbg_name).default
            print("** XBG '{}' can only use XB '{}' values 2 (baldock) and 3 "
                  "(roelvink2)\n** => XB '{}' is {}, so using XBG '{}' "
                  "default of {}: {}".
                  format(xbg_name, xb_name, xb_name, xb_val, xbg_name,
                         xbg_default, breakmodel_val_map[xbg_default]),
                  file=sys.stderr)

        return transformed

    def inc_by_1(self, xb_name: str, xbg_name: str) -> bool:
        xb_val = self.xb_params[xb_name]
        self.xbg_params.set(xbg_name, xb_val+1)
        return True

    def posdown(self, xb_name: str, xbg_name: str) -> bool:
        transformed = False
        xb_val = self.xb_params[xb_name]
        if xb_val == -1:
            self.xbg_params.set(xbg_name, 0)
            transformed = True
        elif xb_val == 1:
            self.xbg_params.set(xbg_name, 1)
            transformed = True
        elif self.verbose:
            print("** Expected -1 or 1 for XB parameter '{}'."\
                  .format(xb_name), file=sys.stderr)
        return transformed

    def xb2xbgfile(self, xb_name: str, xbg_name: str) -> bool:
        def full_path(file_dir, file_path):
            if file_dir is not None:
                file_path = os.path.join(file_dir, file_path)
            return file_path

        xb_file = self.xb_params[xb_name]

        try:
            xbg_file_path = None
            result = True

            if xbg_name == "wavebndfile":
                suffix_start = xb_file.find(".")
                if suffix_start != -1:
                    xbg_file = "{}{}{}".format(xb_file[0:suffix_start], "-xbg",
                                            xb_file[suffix_start:])
                else:
                    xbg_file = "{}{}".format(xb_file[0:suffix_start], "-xbg")

                xbg_file_path = full_path(self.xbg_output_dir, xbg_file)

                xb_file_path = full_path(self.xb_file_root, xb_file)

                is_jonwap = False
                # XB manual is inconsistent re: permitted values of wbctype
                if "wbctype" in self.xb_params and "jons" in self.xb_params["wbctype"]:
                    is_jonwap = True
                else:
                    is_jonwap = "jons" in xb_file_path.lower() or "jswap" in xb_file_path.lower()

                if is_jonwap:
                    # convert the XB to XBG jonswap file and set the wavebndtype
                    convert_jonswap(xb_file_path, xbg_file_path,
                                    self.directional_spread_coefficient, 
                                    self.peak_enhancement_factor,
                                    self.verbose)
                    if self.verbose:
                        print("** ensure that wavebndtype is set to 4 to use this JONSWAP file",
                              file=sys.stderr)

                result = is_jonwap

            elif xbg_name == "depfile":
                xb_file_path = full_path(self.xb_file_root, xb_file)
                xbg_file_path = full_path(self.xbg_output_dir, xb_file)
                xbg_file_path = xbg_file_path.replace(".z", ".dep")
                xbg_file_path = xbg_file_path.replace(".txt", ".dep")
                if xb_file_path != xbg_file_path:
                    shutil.copy(xb_file_path, xbg_file_path)
                if self.verbose:
                    print("** check posdown setting (1 = bathymetry is positive up, "
                        "0 = bathymetry is positive down)", file=sys.stderr)

            elif xbg_name in ["slbndfile", "windbndfile"]:
                xb_file_path = full_path(self.xb_file_root, xb_file)
                xbg_file_path = full_path(self.xbg_output_dir, xb_file)
                if xb_file_path != xbg_file_path:
                    shutil.copy(xb_file_path, xbg_file_path)

            # set the XBG parameter value to the final path
            if xbg_file_path is not None:
                self.xbg_params.set(xbg_name, xbg_file_path)
                if self.verbose and xb_file_path != xbg_file_path:
                    print(">> copied {} to {} ({})".\
                          format(xb_file_path, xbg_file_path, xbg_name))

        except IOError as err:
            result = False
            print(str(err), file=sys.stderr)

        return result


def transform_outvars(xb_params: Dict[str, Any],
                      builder: ParamObjectBuilder,
                      verbose: bool, use_all_known_outvars: bool) -> List[str]:
    """Transform XB output variables for use with XBG.

       :param xb_params: XBeach parameter dictionary
       :param builder: XBeach-GPU parameter/variable object builder
       :param use_all_outvars: allow all output variables to be used?
       :param verbose: verbosity flag
       :return: a list of included output variables"""
    xbg_outvars = []

    # Mean output variables for XBG

    mean_vars_xb2xbg = {
        "H": {"name": "Hmean", "msg": None},
        "hh": {"name": "hhmean", "msg": None},
        "thetamean": {"name": "thetamean",
                      "msg": "units differ (radians in XB, degrees in XBG)"},
        "u": {"name": None, "msg": "XBG has no umean variable; consider uumean instead"},
        "v": {"name": None, "msg": "XBG has no vmean variable; consider vvmean instead"},
        "uu": {"name": "uumean", "msg": None},
        "vv": {"name": "vvmean", "msg": None},
        "zs": {"name": "zsmean", "msg": None}
    }

    # Just the mean XBG variable names for use in excluding from globals
    mean_vars = [mean_vars_xb2xbg[key]["name"] for key in mean_vars_xb2xbg
                 if mean_vars_xb2xbg[key]["name"] is not None]

    # Global output variables for XBG

    # Global output variables with no documentation
    global_var_no_doc = ["C", "DR", "R", "dhdx", "dhdy", "dzsdt", "dzsdx", "dzsdy",
                         "dudx", "dudy", "dvdx", "dvdy", "fwm", "kturb",
                         "rolthick", "sinh2kh", "ududx", "udvdx", "ust",
                         "vdudy", "vdvdy", "vmageu", "vmagev", "wci"]

    # Global output variables about which we are currently uncertain.
    # When we have some certainty, we can remove variables from this list.
    # Note: Cmean & zsmean are mean variables, but are included here because there
    #       does not appear to be any corresponding XB variable and we want to
    #       exclude it for now
    global_var_uncertain = ["Cmean", "cfm", "cgx", "cgy", "ctheta", "dzb", "kh",
                            "sigm", "stdep", "zsmax", "zsmean"]

    # Global output variables that we are including
    global_vars_xb2xbg = OrderedDict()
    for outvar in sorted(builder.outvars):
        if outvar not in global_var_uncertain and \
            (outvar not in global_var_no_doc or use_all_known_outvars) and \
                outvar not in mean_vars:
            global_vars_xb2xbg[outvar] = {"name": outvar, "msg": None}

    global_vars_xb2xbg["thetamean"] = mean_vars_xb2xbg["thetamean"]  # needed here?
    global_vars_xb2xbg["u"] = {"name": None,
                               "msg": "XBG has no u variable; consider uu or uv instead"}
    global_vars_xb2xbg["v"] = {"name": None,
                               "msg": "XBG has no v variable; consider vu or vv instead"}

    def maybe_add_outvar(xb_var_kind: str, xb2xbg_vars: Dict[str, Dict[str, str]]):
        kind = xb_var_kind[0:xb_var_kind.find("_")]
        if verbose:
            print("\nProcessing {} output variables...".format(kind))
        if xb_var_kind in xb_params:
            for xb_var in xb_params[xb_var_kind]:
                if xb_var in xb2xbg_vars:
                    xbg_var = xb2xbg_vars[xb_var]["name"]
                    if xbg_var is not None:
                        xbg_outvars.append(xbg_var)
                    msg = xb2xbg_vars[xb_var]["msg"]
                    if msg is not None:
                        if verbose:
                            if xbg_var is not None:
                                print("** {} (XB) -> {} (XBG): {}".
                                      format(xb_var, xbg_var, msg), file=sys.stderr)
                            else:
                                print("** {} (XB): {}".
                                      format(xb_var, msg), file=sys.stderr)
                else:
                    if verbose:
                        print("** No equivalent of XB's {} {} output variable in XBG".
                              format(kind, xb_var), file=sys.stderr)

    maybe_add_outvar("global_vars", global_vars_xb2xbg)
    maybe_add_outvar("mean_vars", mean_vars_xb2xbg)

    return xbg_outvars


def convert_jonswap(jonswap_in_path: str, jonswap_out_path: str,
                    directional_spread_coefficient: int,
                    peak_enhancement_factor: float,
                    verbose: bool):
    """
    Convert XB jonswap file for use with XBG.

    See also:
    - https://xbeach.readthedocs.io/en/latest/xbeach_manual.html?highlight=jonswap#partable-jonswap
    - Need for XBG: time (s), significant wave height (HS in m), peak period (Tp in s),
      peak direction (Dp in degrees), gamma (peak enhancement factor, usually 3.3),
      directional spread coefficient (cos^2s) in this order per row:
      
        time, wave height, peak period, peak direction, directional spread, gamma

    :param jonswap_in_path: XB jonswap input file path
    :param jonswap_out_path: XBG jonswap output file path
    :param directional_spread_coefficient: Directional spread coefficient
    :param peak_enhancement_factor: Peak enhancement factor
    :param verbose: verbosity flag

    Raises an IOError if the file is not found
    """
    # Open XB jonswap file
    # if os.path.exists(jonswap_in_path):
    jswap = pd.read_csv(jonswap_in_path, delim_whitespace=True, header=None)

    # Create array of timestep values
    timestep = jswap.iloc[0, 5]
    times = [t for t in np.arange(0, len(jswap)*timestep, timestep)]

    # Add time column to a copy of the jonswap dataframe
    new_jswap = jswap.copy()
    new_jswap.insert(0, "time", times)

    # Drop colums 3, 4, 5, and 6. xbeach_gpu does not use 5 and 6
    # Columns 3 and 4 will be replaced below.
    new_jswap = new_jswap.drop([3, 4, 5, 6], axis=1)

    # Set column 3 to directional spread coefficient
    new_jswap[3] = [directional_spread_coefficient for n in range(len(jswap))]

    # Add gamma (peak enhancement factor) at end
    new_jswap.insert(5, "4", [peak_enhancement_factor for n in range(len(jswap))])

    # Write new jonswap file as CSV
    new_jswap.to_csv(jonswap_out_path, index=None, header=None)


def handle_keyval_pair(s: str, params: Dict[str, Any], state: str) -> str:
    """
    Extract the key-value pair from a string, store in the XB
    parameters dictionary, and return the next state.

    :param s: string (hopefully) containing key=val
    :param params: XBeach parameters dictionary
    :param state: current state
    :returns: next state
    """
    matcher = re.match(r"\s*(.+)=(.+)\s*$", s)
    if matcher is not None:
        key = matcher.group(1).strip()
        val = matcher.group(2).strip()
        try:
            val = int(val)
        except ValueError:
            try:
                val = float(val)
            except ValueError:
                pass

        if key is not None:
            params[key] = val
            if key in ["nglobalvar", "nmeanvar"]:
                state = key

    return state


def read_xbeach_params(path: str) -> Dict[str, Any]:
    """
    Read XBeach parameters file given a path to it, and return a dictionary
    of parameter values.

    :param path: XBeach parameter file path
    :return: dictionary of XBeach parameters to values
    """
    params = {}
    global_vars = []
    mean_vars = []
    state = "params"
    with open(path) as xbeach:
        for line in xbeach.readlines():
            if state == "params":
                if len(line.strip()) == 0 or line.lstrip()[0] in ["#", "%"]:
                    continue
                elif re.match(r"^\s*projection=+$", line) is not None:
                    continue
                elif line.find("=") != -1:
                    state = handle_keyval_pair(line, params, state)
            elif "var" in state:
                if "=" in line:
                    # nglobalvar=N or nmeanvar=M
                    state = handle_keyval_pair(line, params, state)
                elif state == "nglobalvar":
                    global_vars.append(line.strip())
                elif state == "nmeanvar":
                    mean_vars.append(line.strip())
            else:
                raise ValueError("Unknown state encountered while processing XB parameter file")

        if len(global_vars) != 0:
            params["global_vars"] = global_vars

        if len(mean_vars) != 0:
            params["mean_vars"] = mean_vars

        return params


def generate(args: argparse.Namespace):
    """Convert XBeach parameters file for use with XBG, generating a
       parameter file.

    :param args: command-line arguments
    """
    def include_param(name: str, value: Any):
        return value not in ["\"\"", "x,y", None]

    xb_params = read_xbeach_params(args.xb_params_path)

    xbg_params_path = "XBG_param.txt"
    if args.xbg_output_dir is not None:
        xbg_params_path = os.path.join(args.xbg_output_dir, xbg_params_path)

    with open(xbg_params_path, "w") as xbg_out:
        print("{}\n".format(XBGParams.title()), file=xbg_out)
        builder = ParamObjectBuilder()
        extract_xbg_params(args.xbg_user_manual, builder)
        for xbg_params in builder.param_objs:
            print("{}".format(xbg_params.decorated_section()), file=xbg_out)

            transformer = ParamTransformer(args.xb_file_root, args.xbg_output_dir,
                                           args.gen_defaults, xbg_params, xb_params,
                                           args.directional_spread_coefficient,
                                           args.peak_enhancement_factor,
                                           args.verbose)
            xb_params_transformed = transformer.transform_params()

            if args.verbose:
                transformed_params = []
                for xb_name in xb_params_transformed:
                    if xb_name != xb_params_transformed[xb_name]:
                        transformed_params.append("{} -> {}".
                                                  format(xb_name,
                                                         xb_params_transformed[xb_name]))
                    else:
                        transformed_params.append(xb_name)
                print(">> parameters: {}".format(", ".join(transformed_params)))

            for xbg_param in xbg_params.values():
                name = xbg_param[0]
                if name in ["name", "outvars"]:
                    continue

                value = xbg_param[1]
                par = xbg_params.get(name)

                if include_param(name, value):
                    if args.gen_doc:
                        if par.default is None or str(par.default).strip() == "":
                            default = "?"
                        else:
                            default = par.default

                        wrapped_doc = "\n# ".join(textwrap.wrap(par.doc, 60))
                        print("\n# {} (default: {})".format(wrapped_doc, default),
                              file=xbg_out)
                    print("{} = {}".format(name.ljust(20, " "), value), file=xbg_out)
            print(file=xbg_out)

        section_name = "Output variables"
        print("{0} {1} {2}".format("#"*3, section_name, "#"*(68-len(section_name))), file=xbg_out)

        transformed_outvars = transform_outvars(xb_params, builder, args.verbose,
                                                args.use_all_known_output_variables)
        if len(transformed_outvars) == 0:
            # no specified/transformed variables, so generate all
            transformed_outvars = builder.outvars.keys()

        if args.gen_doc:
            for outvar in sorted(builder.outvars):
                if outvar in transformed_outvars:
                    wrapped_doc = "\n# ".join(textwrap.wrap(builder.outvars[outvar], 65))
                    print("# {}: {}".format(outvar, wrapped_doc),
                          file=xbg_out)

        print("{} = {}".format("outvars".ljust(20, " "), ",".join(sorted(transformed_outvars))), 
              file=xbg_out)
        print("\n# Generated from {} by {} with arguments:".
                format(args.xb_params_path, os.path.basename(sys.argv[0])), file=xbg_out)
        args_dict = args.__dict__
        print("{}".format("\n".join(["# - {}: {}".
                            format(key.replace("_", "-"),
                                   args_dict[key]) for key in sorted(args_dict)])),
                        file=xbg_out)


def manual_url():
    """
    This XBG URL will need to change if the location of the HTML manual moves.
    Note that it can be overridden via the --xbg-user-manual command-line argument.
    """
    return "https://raw.githubusercontent.com/CyprienBosserelle/xbeach_gpu/gh-pages/Manual.html"


def create_arg_parser() -> argparse.ArgumentParser:
    """
    Create and return the command-line parser.
    """
    parser = argparse.ArgumentParser(
            description="Convert XB params file for use with XBG.")

    parser.add_argument("--xb-params-path", "-b",
                        dest="xb_params_path",
                        default="params.txt",
                        help="XB params input file path")

    parser.add_argument("--xb-file-root", "-r",
                        dest="xb_file_root",
                        default=None,
                        help="XB input file (e.g. bcfile) root. "
                             "We allow XB parameter file location and input "
                             "file locations to vary; if specified, will "
                             "apply to all input files within XB parameters file")

    parser.add_argument("--xbg-output-dir", "-p",
                        dest="xbg_output_dir",
                        default=None,
                        help="XBG params output file directory")

    parser.add_argument("--xbg-user-manual", "-m",
                        dest="xbg_user_manual",
                        default=manual_url(),
                        help="XBG user manual HTML input file path or URL")

    parser.add_argument("--directional-spread-coefficient", "-s",
                        dest="directional_spread_coefficient",
                        type=int,
                        default=400,
                        help="Directional spread coefficient (default: 400)")

    parser.add_argument("--peak-enhancement-factor", "-e",
                        dest="peak_enhancement_factor",
                        type=float,
                        default=3.3,
                        help="Peak enhancement factor (default: 3.3)")

    parser.add_argument("--gen-doc", "-g",
                        dest="gen_doc",
                        default=False,
                        action="store_true",
                        help="Generate documentation strings along with parameters")

    parser.add_argument("--gen-defaults", "-d",
                        dest="gen_defaults",
                        default=False,
                        action="store_true",
                        help="Generate all XBG default parameters "
                             "(this argument overrides --xb-params)")

    parser.add_argument("--verbose", "-v",
                        dest="verbose",
                        default=False,
                        action="store_true",
                        help="Verbose output mode")

    parser.add_argument("--use-all-known-output-variables", "-u",
                        dest="use_all_known_output_variables",
                        action="store_true",
                        help="Allow all known XBG output variables to be used, "
                             "irrespective of uncertainty of relationship to XB")

    return parser


def check_args(parser: argparse.ArgumentParser, args: argparse.Namespace):
    """
    Check args and exit if any problem is found.

    :param parser: Command-line parser.
    :param args: Command-line options.
    """
    msgs = []

    def exists(path, can_be_empty=False):
        return can_be_empty and path == "" or \
            (path is not None and os.path.exists(path))

    if not exists(args.xb_params_path):
        msgs.append("No path to XB parameter file")

    if args.xb_file_root is not None and not exists(args.xb_file_root):
        msgs.append("XB input files root not found")

    if args.xbg_output_dir is not None and \
       not exists(args.xbg_output_dir, can_be_empty=True):
        os.makedirs(args.xbg_output_dir)
        print("** creating XBG output directory: {}".format(args.xbg_output_dir), file=sys.stdout)

    if len(msgs) != 0:
        usage(parser, msgs)


def usage(parser: argparse.ArgumentParser, msgs: List[str] = None):
    """
    Print usage message and exit.

    :param parser: The command-line parser whose help we seek.
    :param msgs: Optional messages to print along with usage.
    """
    if msgs is not None:
        for msg in msgs:
            print("{0}{1}".format(os.linesep, msg), file=sys.stderr)
    print(os.linesep)

    parser.print_help()

    sys.exit(1)


def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    check_args(parser, args)
    generate(args)


if __name__ == "__main__":
    main()
