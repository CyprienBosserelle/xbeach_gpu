### Overview
Python script to convert XBeach parameter files for use with XBeach-GPU. It can also be used to generate default XBeach-GPU parameters.

An ordinary (non-"reuse") JONSWAP file listed in an XBeach parameters file will also be converted. A directional spread coefficient can be specified to be used in the converted JONSWAP file. Any bathymetry file with `.z` suffix will be replaced with `.dep`. In addition, sea-level and wind boundary files will also be copied to the specified output location.

Requires Python 3.7 or higher and the param library.

#### Usage

    python3 xb2xbg.py -h
    usage: xb2xbg.py [-h] [--xb-params-path XB_PARAMS_PATH] [--xb-file-root XB_FILE_ROOT] [--xbg-output-dir XBG_OUTPUT_DIR]
                    [--xbg-user-manual XBG_USER_MANUAL] [--directional-spread-coefficient DIRECTIONAL_SPREAD_COEFFICIENT]
                    [--peak-enhancement-factor PEAK_ENHANCEMENT_FACTOR] [--gen-doc] [--gen-defaults] [--verbose]
                    [--use-all-known-output-variables]

    Convert XB params file for use with XBG.

    optional arguments:
    -h, --help            show this help message and exit
    --xb-params-path XB_PARAMS_PATH, -b XB_PARAMS_PATH
                            XB params input file path
    --xb-file-root XB_FILE_ROOT, -r XB_FILE_ROOT
                            XB input file (e.g. bcfile) root. We allow XB parameter file location and input file locations to
                            vary; if specified, will apply to all input files within XB parameters file
    --xbg-output-dir XBG_OUTPUT_DIR, -p XBG_OUTPUT_DIR
                            XBG params output file directory
    --xbg-user-manual XBG_USER_MANUAL, -m XBG_USER_MANUAL
                            XBG user manual HTML input file path or URL
    --directional-spread-coefficient DIRECTIONAL_SPREAD_COEFFICIENT, -s DIRECTIONAL_SPREAD_COEFFICIENT
                            Directional spread coefficient (default: 400)
    --peak-enhancement-factor PEAK_ENHANCEMENT_FACTOR, -e PEAK_ENHANCEMENT_FACTOR
                            Peak enhancement factor (default: 3.3)
    --gen-doc, -g         Generate documentation strings along with parameters
    --gen-defaults, -d    Generate all XBG default parameters (this argument overrides --xb-params)
    --verbose, -v         Verbose output mode
    --use-all-known-output-variables, -u
                            Allow all known XBG output variables to be used, irrespective of uncertainty of relationship to XB

#### Example Invocations

    # convert XB parameters, generating XBG_params.txt
    python xb2xbg.py --xb-params-path=data/example-XBeach-params.txt

    # as above, but with verbose output
    python xb2xbg.py --xb-params-path=data/example-XBeach-params.txt --verbose

    # also add docs per parameter and output variable where they exist
    python xb2xbg.py --xb-params-path=data/example-XBeach-params.txt --verbose --gen-doc

    # XB input files will have "data" prepended to path found in XB parameters file
    # and XBG files (e.g. XBG parameters, jonswap) will go into the `xbg` directory
    python xb2xbg.py --xb-params-path=data/example-XBeach-params.txt --xb-file-root=data --xbg-output-dir=xbg

    # default XBG parameters will be generated instead of those values that would normally
    # be converted from XB
    python xb2xbg.py --xb-params-path=data/example-XBeach-params.txt --gen-defaults

    # directional spread coefficient for use in jonswap file conversion; defaults to 400
    python xb2xbg.py --xb-params-path=data/example-XBeach-params.txt --directional-spread-coefficient=300

    # peak enhancement factor for use in jonswap file conversion; defaults to 4.2
    python3 xb2xbg.py --xb-params-path=data/example-XBeach-params.txt --peak-enhancement-factor=4.2

    # any output variables in XBG where there is some uncertainty re: relationship to XB
    # will be generated anyway; use with caution
    python xb2xbg.py --xb-params-path=data/example-XBeach-params.txt --verbose --use-all-known-output-variables

    # explicit reference to XBG user manual HTML
    python xb2xbg.py --xb-params-path=data/example-XBeach-params.txt --xbg-user-manual=https://raw.githubusercontent.com/CyprienBosserelle/xbeach_gpu/gh-pages/Manual.html

#### Run Unit Tests (optional)
Requires `pytest` library.

    pytest -v

or, if output from library warnings is too noisy (Python/library version dependent)

    pytest -v --disable-warnings

#### Possible Future Work
* handle TSnode and TSOfile parameters and XBeach equivalents (see npoints, npointvar, tintp)
* refinement/addition of output variable and parameter documentation
* --no-empty-values: don't include parameters with no value in output
* --doc-line-wrap (default to None); generated doc strings currently wrapped at or before column 60