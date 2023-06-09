### convert_xbeach_params.py
Python script to convert XBeach parameter files for use with XBeach-GPU. It can also be used to generate default XBeach-GPU parameters.

#### Usage
    python convert_xbeach_params.py --help
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
                            XB input file (e.g. bcfile) root. We want to allow XB
                            params and input file paths to vary; if specified,
                            will apply to all input files within a XB parameters
                            file.
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
                            Allow all known XBG output variables to be used

#### Example Invocations

    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt

    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt \
    --xbg-params-root=xbg-params-out-dir  # . is the default directory to which XBG_param.txt is written

    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt \
    --gen-defaults  # default XBG parameters will be generated instead of those converted from XB

    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt --verbose

    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt --verbose --gen-doc

    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt --verbose --use-all-known-output-variables

    python convert_xbeach_params.py --xb-params-path=data/example-XBeach-params.txt 
    \--directional-spread-coefficient=300--verbose  # directional spread coefficient for use in jonswap file conversion defaults to 400

#### Run Unit Tests

    pytest -v --disable-warnings

    nosetests -v