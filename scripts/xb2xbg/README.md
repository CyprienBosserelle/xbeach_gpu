### convert_xbeach_params.py
Python script to convert XBeach parameter files for use with XBeach-GPU. It can also be used to generate default XBeach-GPU parameters.

#### Usage
    usage: convert_xbeach_params.py [-h] [--xb-params-path XB_PARAMS_PATH]
                                    [--xb-file-root XB_FILE_ROOT]
                                    [--xbg-params-root XBG_PARAMS_ROOT]
                                    [--xbg-user-manual XBG_USER_MANUAL]
                                    [--directional-spread-coefficient DIRECTIONAL_SPREAD_COEFFICIENT]
                                    [--gen-doc] [--gen-defaults] [--verbose]
                                    [--use-all-known-output-variables]

    Convert XB params file for use with XBG.

    optional arguments:
    -h, --help            show this help message and exit
    --xb-params-path XB_PARAMS_PATH, -b XB_PARAMS_PATH
                            XB params input file path
    --xb-file-root XB_FILE_ROOT, -r XB_FILE_ROOT
                            XB input file (e.g. bcfile) root. We allow XB
                            parameter file location and input file locations to
                            vary; if specified, will apply to all input files
                            within XB parameters file.
    --xbg-params-root XBG_PARAMS_ROOT, -p XBG_PARAMS_ROOT
                            XBG params output file root
    --xbg-user-manual XBG_USER_MANUAL, -m XBG_USER_MANUAL
                            XBG user manual HTML input file path or URL
    --directional-spread-coefficient DIRECTIONAL_SPREAD_COEFFICIENT, -s DIRECTIONAL_SPREAD_COEFFICIENT
                            Directional spread coefficient (default: 400)
    --gen-doc, -g         Generate documentation strings along with parameters
    --gen-defaults, -d    Generate all XBG default parameters (this argument
                            overrides --xb-params)
    --verbose, -v         Verbose output mode
    --use-all-known-output-variables, -u
                            Allow all known XBG output variables to be used,
                            irrespective of uncertainty of relationship to XB

#### Example Invocations

    # convert XB parameters, generating XBG_params.txt
    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt

    # as above, but with verbose output
    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt --verbose

    # also add docs per parameter and output variable where they exist
    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt --verbose
    --gen-doc

    # XB input files will have "data" prepended to what is found in XB parameters file
    # and XBG output files (XBG parameters, jonswap) will go into the data directory
    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt
    --xbg-params-root=data --xbg-params-root=data

    # default XBG parameters will be generated instead of those values that would normally
    # be converted from XB
    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt
    --gen-defaults

    # directional spread coefficient for use in jonswap file conversion; defaults to 400
    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt
    --directional-spread-coefficient=300 --verbose

    # any output variables in XBG where there is some uncertainty re: relationship to XB
    # will be generated anyway; use with caution
    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt --verbose
    --use-all-known-output-variables

#### Run Unit Tests

    pytest -v --disable-warnings

or

    nosetests -v