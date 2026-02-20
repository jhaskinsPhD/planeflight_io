"""
Functions to make readable Planelight.dat.YYYYMMDD input files for GEOS-Chem
and read in outputted planeflight.log files in as pandas dataframes or xarray
datasets with attached meta data and option to concatenate output files.

Created on Sun Mar 21 14:21:14 2021

@author: Dr. Jessica D. Haskins
GitHub: @jhaskinsPhD
Email: jessica.haskins@utah.edu
"""
import os
import sys
import warnings
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime
from planeflight_utils import (_check_simtype, _read_planeflight_diags_yml, _display_diags,
                               _parse_gc_config, _check_str_arr_inputs, _read_planelog_to_df,
                               _build_output_meta_dict, _get_tracer_name_num_mapping,
                               _convert_moleccm3_to_mr, _save_outputs)


# See example folder for usage examples...

###############################################################################
#  Functions for creating input plane.dat files for GEOS-Chem simulations
###############################################################################
def remove_nan_rows(*arrays):
    """Function to remove rows containing NaN values from multiple data collections,
    while maintaining their original types. Provided to help clean campaign data
    inputs for lat/long/time/pressure. Will not change datatype from input type.

    INPUT:
    ------
    arrays : Variable number of collections (such as lists, numpy arrays, or pandas Series)
             or xarray.DataArray that need to be processed for NaN removal.
             The function accepts:
               - list
               - numpy.ndarray
               - pandas.Series
               - pandas.Timestamp
               - xarray.DataArray

    OUTPUT:
    -------
    cleaned_arrays : tuple
        A tuple containing the cleaned collections. Each element has the same
        data type as its corresponding input, with all rows removed where any
        input collection had a NaN.

    BEHAVIOR:
    ---------

    1. Converts each data collection to a pandas DataFrame to handle NaNs effectively.
    2. Removes any rows from this DataFrame where any column contains NaN, ensuring that
       all correlated data across the collections is consistent.
    3. Converts the DataFrame back into the original data type and format for each input
       before returning, thus preserving the type integrity of the inputs.
    """
    # Determine input types and convert inputs to a DataFrame
    input_types = []
    data_dict = {}

    for i, array in enumerate(arrays):
        input_types.append(type(array))

        if isinstance(array, (list, np.ndarray, pd.Series)):
            data_dict[f'col{i}'] = pd.Series(array)
        elif isinstance(array, xr.DataArray):
            data_dict[f'col{i}'] = array.to_series()
        elif isinstance(array, pd.Timestamp):
            data_dict[f'col{i}'] = pd.Series([pd.Timestamp(array) if not pd.isnull(array)
                                             else pd.NaT for _ in range(len(arrays[0]))])
        else:
            raise ValueError(f"Unsupported data type: {type(array)}")

    df = pd.DataFrame(data_dict)

    # Drop rows with any NaN values
    df_cleaned = df.dropna()

    # Convert back to the original input types
    cleaned_data = []
    for i, typ in enumerate(input_types):
        if typ == list:
            cleaned_data.append(df_cleaned[f'col{i}'].tolist())
        elif typ == np.ndarray:
            cleaned_data.append(df_cleaned[f'col{i}'].to_numpy())
        elif typ == pd.Series:
            cleaned_data.append(pd.Series(df_cleaned[f'col{i}'].to_list()))
        elif typ == xr.DataArray:
            cleaned_series = df_cleaned[f'col{i}'].to_series()
            cleaned_data.append(xr.DataArray(cleaned_series.values, dims=array.dims))

    return tuple(cleaned_data)


def get_compatible_input_diags(simtype: str = '', these_collections: list = [],
                               display: bool = False):
    """Function to load all optional planeflight diagnotics in an
    output dicitonary as keys with info on diagnostic's FULLNAME/UNITS stored
    in values that are compatiable with a user's simulation type & print results
    to screen if requested.

    INPUT:
    -------

    (1) simtype           - STRING containing user's simulation type.
                            Used to determine what available planflight diagnostics are
                            compatiable with your simulation type.

    (2) these_collections - (OPTIONAL): LIST of STRINGS containing diagnostic
                            collection names that you would like to retrieve.
                            Useful if you only want to get specific collections
                            of diags, not just ALL diags compatiable with simtype:
                            Valid options to pass include:
                              ['aer_uptake', 'aodb', 'aodc', 'aq_aer', 'broken',
                               'chem_fams', 'defaults', 'gmao_ice', 'gmao_met',
                               'hg', 'htep', 'tomas']

    (2) display           - (OPTIONAL): BOOL indicating if compatible requested
                            diagnotics should be also be printed to screen.
    OUTPUT:
    -------

    (1) diag_dict - DICTIONARY with keys corresponding to diagnostic names
                    used in planeflight.dat input files / plane.log output
                    files with VALUES being a dictionary with "FULLNAME" and
                    "UNITS" keys describing what the diagnostic is / its units.
    """
    # Check that user passed an allowed/recognized simtype or throw warning.
    _check_simtype(simtype, blank_allowed=False)

    # Read in all valid input planeflight diags (not including advected species) as dictionary.
    diags = _read_planeflight_diags_yml(inputs_only=True)

    # Initialize dict holding strings of all available diagnostic collections:
    ok_diags = list(diags.keys())

    ###########################################################################
    # Remove invalid collections from "ok_diags" based on user's simulation type
    ###########################################################################
    if simtype.lower() != 'hg':
        ok_diags.remove('hg')
    if simtype.lower() != 'tomas':
        ok_diags.remove('tomas')
    if simtype.lower() == 'carbon':
        for item in ['aer_uptake', 'htep', 'chem_fams']:
            ok_diags.remove(item)

    # The "broken" collection includes diags that do not work with any simtype b/c
    # they need source code updates in GEOS-Chem... Only return if "all" is requested.
    ok_diags.remove('broken')

    ###########################################################################
    # If user asked for specific compatible diags by listing collections:
    ###########################################################################
    if len(these_collections) > 0:

        # Check that user only passed valid args for "collections" or throw error:
        invalid = [item for item in these_collections if item.lower() not in list(diags.keys())]
        if len(invalid) > 0:
            raise ValueError('Invalid planeflight diagnostic collection name(s):\n\t' +
                             ','.join(invalid) + '\nValid collection names include: \n\t' +
                             ','.join(list(diags.keys()))).with_traceback(sys.exc_info()[2])

        # Check for any requested collections that aren't compatiable with simtype:
        not_compat = [item for item in these_collections if item not in ok_diags]
        if len(not_compat) > 0:
            # If all requested collections aren't copmatiable, throw error:
            if len(not_compat) == len(these_collections):
                info = _display_diags(diags, ok_diags)
                raise ValueError('None of the requested collections are compatiable with your ' +
                                 'simulation type. The following options are those compatiable:' + info)

            else:
                # Otherwise just warn we're dropping some...
                warnings.warn('The following requested collections are not compatiable ' +
                              'with your simulation type and will NOT be included in the ' +
                              'outputted diagnostic collection:\n\t' + ','.join(not_compat) +
                              '\nUse pln.get_compatible_input_diags(simtype,display=True) to display ' +
                              'all diagnostics compatiable with your simulation type.',
                              UserWarning, stacklevel=2)

        # Define collections so it only includes compatiable collections:
        collections = [item for item in these_collections if item in ok_diags]

    else:  # If specific collections weren't requested, return all compatiable options:
        collections = ok_diags

    ###########################################################################
    # Parse the now curated collections list and build output diags dict:
    ###########################################################################
    diag_dict = dict({})
    for c in collections:
        for key in diags[c]['Diagnostics']:
            diag_dict[key] = diags[c]['Diagnostics'][key]

    ###########################################################################
    # Print compatiable diags to screen if asked
    ###########################################################################
    if display is True:
        # Format info about diags in parse_this to print to screen:
        info = _display_diags(diags, collections)
        print(info)

    return diag_dict


def make_planeflight_inputs(savedir: str,
                            gc_config,
                            datetimes,
                            lat_arr,
                            lon_arr,
                            vert_arr,
                            vert_is_pres: bool,
                            tracers,
                            tracers_minus: list = [],
                            diags=[],
                            diags_minus: list = [],
                            simtype: str = '',
                            username: str = 'user',
                            overwrite: bool = False,
                            use_tracer_names: bool = False):
    """Function to create planeflight.dat files in correct format for GEOS-Chem input.

    INPUT:
    -------
    Note: For all inputs where "ARRAY" is accepted, vars can be 1-D lists, np.ndarray,
          pd.Series, or xr.dataarray.

       (1) savedir      - STRING containing path to directory in which to save
                          the output planeflight.dat files at.

       (2) gc_config   - STRING containing path to your geoschem_config.yml file,
                          OR None if you do not have the file available.
                          When a path is given, it is used to read the advected
                          species list, detect the simulation type, and map tracer
                          names to numbers (so output is in mol/mol).
                          When None, you must supply simtype= explicitly, tracers
                          must be an explicit list (not '?ALL?'), and tracer names
                          are always written (output will be in molec/cm3 unless
                          you pass convert2_molmol=True when reading the output).

       (3) datetimes    - PANDAS SERIES of datetimes in UTC stored as Timestamps at which
                          to sample the model.Type(datetimes[0]) should return:
                          <class 'pandas._libs.tslibs.timestamps.Timestamp'>

                          To create input in the correct format do:
                               date_range = pd.date_range(start='2017-01-01', end='2017-01-03', freq='60s')
                               datetimes=pd.Series(date_range)

       (4) lat_arr      - ARRAY of latitudes at which to sample the model. (range: -90 to 90 deg)

       (5) lon_arr      - ARRAY of longitudes at which to sample the model (range: -180 to 180 deg)

       (6) vert_arr     - EITHER an ARRAY of pressures (hPa) OR altitudes above
                          the ground (meters) at which to sample the model. See
                          https://github.com/geoschem/geos-chem/issues/320
                          for discussion on whether input altitudes should be
                          "above ground" or "above sea level".

       (7) vert_is_pres - BOOLEAN indicating if "vert_arr" containined pressures or not.
                          When TRUE,  values are assumed to be pressures (hPa).
                          When FALSE, values are assumed to be altitudes (meters).

       (8) tracers     -  Either (1) an ARRAY of specific advected tracers you want to sample
                              OR (2) a STRING equal to '?ALL?' to sample all advected species
                          listed in your geoschem_config.yml file.

       (9) tracers_minus - (OPTIONAL) ARRAY containing STRINGS with all advected species
                          tracers you don't want to include (only relevant if you passed
                          tracers='?ALL?').

       (10) diags       - (OPTIONAL) Either (1) an ARRAY containing STRINGS of any additional
                          diagnostics to sample from model (beyond tracers) OR (2) a
                          STRING equal to '?ALL?' to sample all available additional
                          diagnostics compatiable with your simulation type. Default is
                          to include the grid-box indexes planflight pulled from &
                          Pres/Temp/RH at center of grid box (e.g. ['GMAO_IIEV',
                          'GMAO_JJEV', 'GMAO_LLEV', 'GMAO_PRES','GMAO_RELH', 'GMAO_TEMP']).

       (11) diags_minus - (OPTIONAL) ARRAY containing STRINGS with all additional
                          diagnostics you don't want to include (only relevant
                          if you passed diags='?ALL?').

       (12) simtype    - (OPTIONAL) STRING identifying the simulation type. Only
                          required when gc_config=None; ignored otherwise (the
                          simulation type is read from the config file instead).
                          Valid options: 'fullchem', 'aerosol', 'carbon', 'Hg',
                          'POPs', 'tagO3', 'TransportTracers', 'metals', 'CH4',
                          'CO2', 'tagCO'.

       (13) username   - (OPTIONAL) STRING contaiing name of user who created files.
                         This gets listed in header of resulting planedat input files.

       (14) overwrite  - (OPTIONAL) BOOLEAN of whether to overwrite any existing files
                          at 'savedir' with this name or not. If FALSE, & any files
                          under 'savedir' would be overwritten, a new sub-directory
                          under 'savedir' named "NEW_YYYYMMDD_HHMMSS" is created to
                          hold all the new output files. If TRUE, only files
                          with conflicting names under 'savedir' are overwritten.
                          Default is set to FALSE (not to overwrite files).

       (15) use_tracer_names - (OPTIONAL) BOOL indicating if you want to write the file
                          with tracer names rather than with tracer numbers. Default is
                          set to FALSE (to write file using tracer numbers) since outputs
                          like this have units of mol/mol dry. If set to true, input
                          files will be written with tracer names instead, which will
                          result in the output files have units of molec/cm3. While
                          using tracer names is more readable, it adds an extra step of
                          potential error to compare directly to observations. Thus,
                          writing files with tracer numbers is reccommended. This code
                          can accomodate either though when reading in output files.

    OUTPUT:
    ------
       (1) One file for each day listed in 'datetimes' named "planeflight.dat.YYYYMMDD"
           written to "savedir"that can be passed to GEOS-Chem as input files for
           the planeflight diagnostic.

    """
    # -------------------------------------------------------------------------
    #                   Perform User Input Checks:
    # --------------------------------------------------------------------------

    # Check that the output "savedir" directory exists:
    if not os.path.isdir(savedir):
        raise OSError("The following directory passed for savedir='" + savedir +
                      "' does not exist.").with_traceback(sys.exc_info()[2])

    # Check that datetimes, lat_arr, lon_arr, and vert_arr don't contain NaNs/NaT/Inf:
    if (pd.isna(datetimes)).any():
        raise ValueError('Input "datetimes" contains np.NaT or np.nan.').with_traceback(sys.exc_info()[2])
    if ((np.isnan(np.array(lat_arr)).any()) or (np.isinf(np.array(lat_arr)).any())):
        raise ValueError('Input "lat_arr" contains np.nan or np.inf.').with_traceback(sys.exc_info()[2])
    if ((np.isnan(np.array(lon_arr)).any()) or (np.isinf(np.array(lon_arr)).any())):
        raise ValueError('Input "lon_arr" contains np.nan or np.inf.').with_traceback(sys.exc_info()[2])
    if ((np.isnan(np.array(lat_arr)).any()) or (np.isinf(np.array(lat_arr)).any())):
        raise ValueError('Input "vert_arr" contains np.nan or np.inf.').with_traceback(sys.exc_info()[2])

    # Check that all lats are in range -90 to 90 or raise error:
    if np.any(np.abs(lat_arr) > 90):
        raise ValueError('Some |Latitudes| are > 90 (e.g. out of range)' +
                         'in input "lat_arr".').with_traceback(sys.exc_info()[2])

    # Check that all lons are in range -180 to 180 or raise error:
    if np.any(np.abs(lon_arr) > 180):
        raise ValueError('Some |Longitudes| are > 180 (e.g. out of range)' +
                         'in input "lon_arr".').with_traceback(sys.exc_info()[2])

    # Check that pressures/alts in vert_arr are not < 0 or raise error:
    if np.any(np.abs(vert_arr) < 0):
        raise ValueError('Some values in "vert_arr" are < 0 (e.g. out of range).'
                         ).with_traceback(sys.exc_info()[2])

    # Check that no pressures are > 1100 hPa or ask if they want to proceed:
    if ((vert_is_pres is True) and (np.any(np.abs(vert_arr) > 1100))):
        print('WARNING: Some |Pressures| in "vert_arr" are > 1100 hPa.' +
              'Did you input them in a unit other than hPa by mistake?')
        while True:
            user_input = input('Enter "1" to proceed if this was intended.\n' +
                               'Enter "2" to exit now.')
            if user_input == 2:
                exit(0)
            elif user_input != 1:
                print("Invalid input. Please enter '1' to proceed or '2' to exit.")

    # Check that NOT all alts <15 m, ask them if they want to proceed or not.
    elif ((vert_is_pres is True) and (np.any(np.abs(vert_arr) < 15))):
        print('WARNING: All altitude inputs in "vert_arr" are < 15 m.' +
              'Did you input them in kilometers by mistake?')
        user_input = input('Enter "1" to proceed if this was intended.\n' +
                           'Enter "2" to exit now.')
        if user_input == 2:
            exit(0)
        elif user_input != 1:
            print("Invalid input. Please enter '1' to proceed or '2' to exit.")

    # Check that 'lat_arr', 'lon_arr', 'vert_arr', and 'datetimes' are all the same length:
    input_arr_lens = [len(var) for var in [datetimes, lat_arr, lon_arr, vert_arr]]
    if any([len_i != input_arr_lens[0] for len_i in input_arr_lens]):
        raise ValueError('Inputs for "datetimes","lat_arr","long_arr" and "vert_arr" ' +
                         'must all be the same length.' +
                         '\n\tlen(datetimes)=' + str(input_arr_lens[0]) +
                         '\n\tlen(lat_arr)=' + str(input_arr_lens[1]) +
                         '\n\tlen(long_arr)=' + str(input_arr_lens[2]) +
                         '\n\tlen(vert_arr)=' + str(input_arr_lens[3])).with_traceback(sys.exc_info()[2])

    if gc_config is not None:
        # Parse config, validate tracers against adv_species, map names→nums when needed.
        adv_species, simtype = _parse_gc_config(gc_config)

        # Check that input for "tracers" is either our wild card string ('?ALL?') or
        # is a list-like type of tracers that doesn't contain any tracers not tracked in model:
        tracer_err, tracer_err_msg, tracers = _check_str_arr_inputs(
            tracers, 'tracers', adv_species, 'the list of advected species in geoschem_config.yml', wildcard_str=True)

        # If any checks on "tracers" input were failed, throw descriptive error
        if tracer_err is True:
            raise ValueError('Invalid input for "tracers".' + tracer_err_msg +
                             ' Input for "tracers" must either be\n:' +
                             '\t(1) A STRING equal to "?ALL?" to sample all' +
                             'advected species listed in geoschem_config.yml' +
                             '(e.g. tracers="?ALL/").\n''\t(2) A LIST/ARRAY of ' +
                             'STRINGS containing tracers you wish to sample ' +
                             '(e.g. tracers=["NO","NO2"]').with_traceback(sys.exc_info()[2])

        if len(tracers_minus) > 0:  # If the has user asked to exclude specific advected species tracers...

            # Check that input for "tracers_minus" is a list/array of diagnostics
            # and that it only includes diagnositics we recognize from list of diags:
            tracers_minus_err, tracers_minus_err_msg, tracers_minus = _check_str_arr_inputs(
                tracers_minus, 'tracers_minus', adv_species,
                'the list of advected species in geoschem_config.yml',
                wildcard_str=False)

            # If any checks on "tracers_minus" were failed, throw descriptive error
            if tracers_minus_err is True:
                raise ValueError('Invalid input for "tracers_minus".' + tracers_minus_err_msg +
                                 ' Input for "tracers_minus" must be a LIST/ARRAY ' +
                                 'of STRINGS containing advected species you do NOT wish ' +
                                 'to sample (e.g. tracers_minus=["NO2","NO"].'
                                 ).with_traceback(sys.exc_info()[2])

            # Remove any diagnostics in list "tracers" if they apprear in "tracers_minus":
            tracers = [item for item in tracers if item not in tracers_minus]

        # Convert tracer names to tracer numbers so output is given in mol/mol not molec/cm3 if asked.
        if use_tracer_names is False:
            names2nums = _get_tracer_name_num_mapping(gc_config, keys_are_nums=False)
            tracers = [names2nums[name] for name in tracers]

    else:
        # gc_config=None path: user must supply simtype and an explicit tracer list.
        if not simtype:
            raise ValueError(
                "When gc_config=None, you must supply simtype= explicitly.\n"
                "Valid options: 'fullchem', 'aerosol', 'carbon', 'Hg', 'POPs', "
                "'tagO3', 'TransportTracers', 'metals', 'CH4', 'CO2', 'tagCO'.")
        _check_simtype(simtype, blank_allowed=False)

        # '?ALL?' requires the advected species list from gc_config — not possible here
        if isinstance(tracers, str):
            raise ValueError(
                "tracers='?ALL?' requires gc_config to retrieve the full advected species "
                "list. Provide an explicit list of tracer names, or supply gc_config.")

        # Validate tracers is a list-like of strings (can't check against species list)
        if not isinstance(tracers, (list, np.ndarray, pd.Series)):
            raise TypeError("'tracers' must be a list of tracer name strings when gc_config=None.")
        if not all(isinstance(t, (str, np.str_)) for t in tracers):
            raise TypeError("All items in 'tracers' must be strings (tracer names).")
        tracers = list(tracers)

        # tracers_minus: just subtract, no need to validate against species list
        if len(tracers_minus) > 0:
            tracers = [t for t in tracers if t not in list(tracers_minus)]

        # Name→number mapping requires gc_config — warn and continue with names
        if use_tracer_names is False:
            warnings.warn(
                "use_tracer_names=False has no effect when gc_config=None — the name→number "
                "mapping requires geoschem_config.yml. Tracer names will be used, so advected "
                "species concentrations in the output plane.log will be in molec/cm³. "
                "Pass convert2_molmol=True when reading the output to convert to mol/mol.",
                UserWarning, stacklevel=2)

    if len(diags) > 0:  # If the has user asked to include specific optional diagnostics...

        # Get dictionary of ALL optional diagnostics compatiable with user's simtype.
        diag_dict = get_compatible_input_diags(simtype, display=False)

        # Check that input for requested "diags" is either our wild card string ('?ALL?') or
        # is a list/array of diagnostics including ONLY those compatiable with user's simulation type.
        diag_err, diag_err_msg, diags = _check_str_arr_inputs(diags, 'diags', list(diag_dict.keys()),
                                                              'the list of optional diagnostics compatible' +
                                                              'with simtype="' + simtype + '" outputted from ' +
                                                              'get_diag_info(simtype,return_all=False)',
                                                              wildcard_str=True)

        # If any checks on "diag" were failed, throw descriptive error
        if diag_err is True:
            raise ValueError('Invalid input for "diags".' + diag_err_msg +
                             ' Input for "diags" must either be\n:\t(1) A STRING ' +
                             'equal to "?ALL?" to include all optional diagnostics ' +
                             'compatible with your simtype="' + simtype + '"(e.g. dias="?ALL/").\n' +
                             '\t(2) A LIST/ARRAY of STRINGS containing diags you wish ' +
                             'to sample (e.g. diags=["GMAO_RELH","RO2"].\n To see a list of ' +
                             'all optional diagnostics compatible with your simulation ' +
                             'try: pln.get_diag_info("' + simtype + '",print_options=True)'
                             ).with_traceback(sys.exc_info()[2])

        if isinstance(diags, str):
            # If input is string (which will match wildcard if made it this far),
            # then grab all valid options from diag_dict to define list of diags to include.
            diags = [*diag_dict]

        if len(diags_minus) > 0:  # If the has user asked to exclude specific optional diagnostics...

            # Check that input for "diags_minus" is a list/array of diagnostics
            # and that it only includes diagnositics we recognize from list of diags:
            diag_minus_err, diag_minus_err_msg, diags_minus = _check_str_arr_inputs(
                diags_minus, 'diags_minus',
                list(diag_dict.keys()),
                'the list of optional diagnostics compatible'
                'with simtype="' + simtype + '" outputted from '
                'get_diag_info(simtype,baddies=False)',
                wildcard_str=False)

            # If any checks on "diags_minus" were failed, throw descriptive error
            if diag_minus_err is True:
                raise ValueError('Invalid input for "diags_minus".' + diag_minus_err_msg +
                                 ' Input for "diags_minus" must be a LIST/ARRAY ' +
                                 'of STRINGS containing diags you do NOT wish ' +
                                 'to sample (e.g. diags_minus=["GMAO_RELH","RO2"].\n To see a list of ' +
                                 'all optional diagnostics compatible with your simulation ' +
                                 'try: pln.get_diag_info("' + simtype + '",print_options=True)'
                                 ).with_traceback(sys.exc_info()[2])

            # Remove any diagnostics in list "diags" if they apprear in "diags_minus":
            diags = [item for item in diags if item not in diags_minus]

    # -------------------------------------------------------------------------
    # Assign "TYPE" column in output planedat files appropriately for vert inputs:
    # --------------------------------------------------------------------------
    # GEOSCore/planeflight_mod.F90 decides input is pres/alt based on the string
    # assigned as "TYPE". But, it only converts alts to pressures to pull sample from
    # model if one of the following strings are passed as the "TYPE":
    valid_alt_types = [  # noqa: F841
        'Aacg', 'Aabne', 'Acar', 'Acma', 'Acrv', 'Adnd', 'Aesp', 'Aetl',
        'Ahil', 'Ahip', 'Alef', 'Anha', 'Apfa', 'Arta', 'Asca', 'Asgp', 'Atgc', 'Athd', 'Awbi',
        'Tbao', 'Tamt', 'Tcrb', 'Tlef', 'Tlew', 'Tmbo', 'Tmvy',
        'Tmwo', 'Tnwr', 'Tsct', 'Tsgp', 'Tstr', 'Twbi', 'Twgc', 'Twkt',
    ]

    if vert_is_pres:
        # If user indicates values in vert_arr are PRESSURE, assign TYPE as "Pinp"
        # which is NOT one of the strings that will trigger GEOS-Chem to assume input
        # is altiutude and convert it to pressure:
        type2assign = 'Pinp'
    else:
        # If user indicates values in vert_arr are ALTITUDES, assign TYPE to be
        # one of the strings that triggers GEOS-Chem to convert alts to pressures:
        type2assign = 'Tbao'  # valid_alt_types[0]

    # Make list of all tracers/diagnostics that planedat files should include:
    tracer_list = diags + tracers
    ntracers = str(int(len(tracer_list)))  # count the number of quantities.

    # Designate a few variables that we'll use to make the header lines of the files
    today = str(datetime.now().strftime('%m-%d-%Y %H:%M:%S'))
    spacer = '-------------------------------------------------------------------------------'
    title = '  Now give the times and locations of the flight'

    # Parse the passed dates, get list of individual days we need planeflight files for:
    all_dates = datetimes.dt.date
    unq_dates = np.unique(datetimes.dt.date)

    # Make list of filenames we'll create & check to see if any existing files
    # under savedir match that pattern/ would be overwritten:
    new_files = ['Planeflight.dat.' + unq_dates[i].strftime('%Y%m%d') for i in range(0, len(unq_dates))]
    old_files = [f for f in os.listdir(savedir) if os.path.isfile(os.path.join(savedir, f)) and 'Planeflight.dat' in f]
    dupe_files = [file for file in new_files if file in old_files]

    # Create a new subdirectory under savedir to hold files if overwrite=False.
    if len(dupe_files) > 0 and overwrite is False:
        print('Planeflight.dat files found that would be overwritten under: ' + savedir)
        new_savedir = os.path.join(savedir, 'NEW_' + today.replace('-', '').replace(':', ''))
        print('Creating new directory to hold new planeflight.dat files at:' + new_savedir)
        os.mkdir(new_savedir)
        savedir = new_savedir

    # Loop over each unique day in date range we need to make a planeflight.dat file for:
    for i in range(0, len(unq_dates)):

        # Create a filename based on the date:
        filename = 'Planeflight.dat.' + unq_dates[i].strftime('%Y%m%d')

        if overwrite:  # Delete any existing files under savedir with this name.
            if os.path.isfile(os.path.join(savedir, filename)):
                os.remove(os.path.join(savedir, filename))

        # Get indexes in larger array of flight obs that took place on this date.
        inds = np.where(all_dates == unq_dates[i])[0]

        if len(inds) > 0:
            # ==========================================================================
            # Build arrays of all the data we need in Planeflight.dat for this date.
            # Rount to # of decmial places allowed in GEOS-Chem input for these vars.
            # ==========================================================================

            # Define array with that will fill POINTS column (indexing var)
            points = np.arange(1, len(inds) + 1).astype(str)
            # Define array holding data to fill OBS column full of 9999.000
            obs = np.full(len(inds), 9999.000)
            # Define/Format array that will hold TYPE col (allowex max of 7 chars)
            typez = np.full(len(inds), '{typ: >6}'.format(typ=type2assign))
            # Define array holding DATE (GMT) for each obs point in DD-MM-YYYY format
            day = datetimes.dt.strftime('%d-%m-%Y')[inds]
            # Define array holding TIME (GMT) for each obs point in HH:MM format:
            tms = datetimes.dt.strftime('%H:%M')[inds]
            # Define array holding LAT, LON, and ALT/PRE info formatted to max 2 decimals:
            lats = np.around(lat_arr[inds], decimals=2)
            lons = np.around(lon_arr[inds], decimals=2)
            verts = np.around(vert_arr[inds], decimals=2)  # GC doesn't allow more than 2 decimals

            # Pack all this into a dataframe, because pandas writes tab delimited files nicely!
            df = pd.DataFrame({'POINT': points, 'TYPE': typez, 'DD-MM-YYYY': day,
                              'HH:MM': tms, 'LAT': lats, 'LON': lons, 'ALT/PRE': verts,
                               'OBS': obs})

            # Append a line at the bottom of the DF that says its the end!
            new_row = pd.DataFrame([{'POINT': 99999, 'TYPE': '{typ: >6}'.format(typ='END'),
                                     'DD-MM-YYYY': '00-00-0000', 'HH:MM': '00:00',
                                     'LAT': 0.00, 'LON': 0.00, 'ALT/PRE': 0.00,
                                     'OBS': 0.000}])

            # Use pd.concat to append the new row to the existing DataFrame
            df = pd.concat([df, new_row], ignore_index=True)

            # =========================================================================
            # Format everything so they're in GEOS-Chem's expected format!
            # ==========================================================================

            # Note: I got these vals used for each var from planeflight_mod.f90 i/o checks.
            df.POINT = df.POINT.map('{: >5}'.format)  # POINTS has len 5 char strings
            df.LAT = df.LAT.map('{:7.2f}'.format)    # LAT has max 6 total digits, w/ 2 decimals
            df.LON = df.LON.map('{:7.2f}'.format)    # LON has max 6 total digits, w/ 2 decimals
            df['ALT/PRE'] = df['ALT/PRE'].map('{:7.2f}'.format)  # ALT/PRE has max 6 total digits, w/ 2 decimals
            df.OBS = df.OBS.map('{:10.3f}'.format)   # OBS has max 10 total digits, w/ 3 decimals

            # Write the header lines that GEOS Chem expects to the file, considering format.
            header = '{strr: >6}'.format(strr='POINT ') + \
                '{strr: >7}'.format(strr='TYPE ') + \
                '{strr: >11}'.format(strr='DD-MM-YYYY ') + \
                '{strr: >6}'.format(strr='HH:MM ') +\
                '{strr: >8}'.format(strr='LAT ') +\
                '{strr: >8}'.format(strr='LON ') +\
                '{strr: >7}'.format(strr='ALT/PRE') +\
                '{strr: >11}'.format(strr='OBS')

            # ======================================================================
            #            Begin writing the planeflight.dat text file
            # =====================================================================
            # This is just a list containing our header lines in the right order...
            textList = [filename, username, today, spacer, ntracers, spacer] + \
                tracer_list + [spacer, title, spacer, header]

            # Write the data to a temporary file using ASCII encoding. Pandas default
            # in Python 3 uses UTF-8 encoding, which GEOS-Chem can't read.
            df.to_csv(savedir + filename + '_0', header=False, index=None, sep=' ', mode='a',
                      encoding='ascii')

            # Open the output file and write headers, then append the temp file contents.
            # Because we use a space delimiter, pandas wraps strings in quotes; read the
            # temp file line by line and strip those quotes before writing to the output file.
            with open(savedir + filename, "w") as outF:
                for line in textList:
                    outF.write(line)
                    outF.write("\n")
                with open(savedir + filename + '_0', 'r') as fin:
                    for line in fin:
                        outF.write(line.replace('"', ''))

            os.remove(savedir + filename + '_0')  # Delete the temp file.

            print('Output saved at: ' + savedir + filename)  # Tell where output is saved.

    return

###############################################################################
#  Functions for reading in output plane.log files from GEOS-Chem simulations
###############################################################################


def read_planelog(planelog_file: str, spdb_yaml: str, config_yaml: str,
                  as_xarr: bool = False, convert2_molmol: bool = False,
                  output_dir: str = '', output_file: str = '', overwrite: bool = False):
    """Function to read a single plane.log ouput file into either a pandas dataframe
    or xarray dataset. Handles data split across multiple lines in a plane.log file
    if a large number of output diagnostics are requested and auto converts
    the "YYYYMMDD" and "HHMM" date/time diagnostics into a pandas datetime object
    used as either the pandas dataframe index or the xarray coordinate. Automatically
    corrects output units to mol/mol if tracer names rather than tracer numbers used
    if asked. Will always output metadata telling you about the vars in the outputs.

    INPUT:
    ------
        (1) planelog_file  -  STRING containing absolute path to a plane.log file.

        (2) spdb_yaml      -  STRING containing absolute path to the species_database.yml
                              file associated with your run. Used to add metadata.

        (3) config_yaml    -  STRING containing path to geoschem_config.yml file
                              associated with your run. Used to add metadata.

        (4) as_xarr         - (OPTIONAL) BOOLEAN indicating if you want to output the data as
                              an xarray dataset instead of as a pandas dataframe.
                              If this is set to TRUE, the data arrays within the
                              output dataset will have attributes containing metadata
                              assigned to them. Otherwise, output in the xarray
                              dataset is the same as that in the df. Default is set
                              to FALSE (to return a pandas df instead).

        (5) convert2_molmol - (OPTIONAL) BOOLEAN indicating if you wish to convert all
                              outputted concentrations into units of mol/mol.
                              Only variables with concentration units are
                              converted & others are left in their native unit.
                              Default is set to FALSE (not to do conversion).
                              If set to "TRUE" and "TRA_###" appears in plane.log headers,
                              no conversion is actually performed to prevent against
                              improper usage.

        (6) output_dir       - STRING containing absolute path to a directory where
                              output file will be saved. If none is passed, then
                              output saved in same dir as input planelog file.

        (7) output_file     - STRING containing desired name of ouptput file. If none
                              is passed, the output file will be named:
                                  "planelog_[startYYYYMMDD]_[endYYYYMMDD]"

        (8) overwrite       - BOOL indicating that if there is already a file named
                              "output_file" in the output directory, "out_dir", whether
                              that file should be overwritten or if a new unique name
                              should be used.

        NOTE: Outputted advected species concentrations are in molec/cm3 if Tracer
              Names are used, but are in mol/mol if Tracer Numbers are used in input files.
              For more details on why this happens, see GitHub issue #796:
                  https://github.com/geoschem/geos-chem/issues/796

    OUTPUT:
    -------
    If output_xarr==FALSE (Default):

        (1) df        - PANDAS DATAFRAME indexed along a pd.datetime object index named 'time_UTC'.
                         Saved at: output_dir+output_file+'.csv'

        (2) info_dict - DICTIONARY containing info about output vars & units.
                        Saved as yaml at: output_dir+output_file+'.yaml'

    If output_xarr == TRUE:

        (1) ds - XARRAY DATASET indexed along a np.datetime64 type coordinate
                 & dimension named 'time_UTC'. All data vars within the ds have
                 individual attributes for "long_name" and "units" assigned to them
                 and will reflect unit conversion changes to mol/mol from molec/cm3
                 if convert2_molmol is set to TRUE. If this is done an additional
                 attribute 'note' is included indicating this occurred.
                 Saved at: output_dir+output_file+'.nc'

    Regardless of the output type selected, the resulting df/ds contains columns/data arrays
    corresponding to all output diagnostics except "YYYYMMDD" and "HHMM" which are
    dropped after constructing the time index/coordinate. Column/data_var names are
    nearly identical to the diagnostics in plane.log, with the exception that any
    '-' characters are replaced with '_' (e.g. 'P-I' in plane.log becomes 'P_I').
    """

    # Read raw file in as dataframe (no unit conversions / renaming)
    df = _read_planelog_to_df(planelog_file)

    # Build meta data dictionary with info about all vars in output dataframe.
    #   This Assumes advected species units are mol/mol if tracer numbers used, but
    #    molec/cm3 if tracer numbers used.
    #
    # **NOTE**: This line MUST remain before renaming df columns with tracer numbers
    #           to tracer names for output units to be correct in info_dict since
    #           assigning them properly requires us to know if the column names
    #           were originally tracer numbers or tracer names.
    info_dict = _build_output_meta_dict(df, config_yaml, spdb_yaml)

    # Get dictionary to map tracer numbers to tracer names:
    nums2names = _get_tracer_name_num_mapping(config_yaml, keys_are_nums=True)

    # Rename all columns in dataframe so columns with tracer numbers are converted
    # to tracer names. Won't affect anything if cols are already named with tracer names.
    df.rename(columns=nums2names, inplace=True)

    # Convert concentrations to mol/mol from molec/cm3 & update units in info_dict
    # accordingly. Also add note if conversion performed or not ...
    if convert2_molmol is True:
        # NOTE: Function will not actually do anything if conc_units!='molec/cm3'
        # which is the case if tracer numbers were used in analogous input files.
        df, info_dict = _convert_moleccm3_to_mr(df, info_dict)

    # =========================================================================
    # Configure output filename/directory:
    # =========================================================================
    # Verify that output dirpath exists, otherwise assign output_dir to dir of input file.
    if not os.path.isdir(output_dir):
        warnings.warn('The output directory passed to read_planelog() could not '
                      'be found:\n\t' + output_dir + '\n Output will be saved at '
                      'input planelog file directory instead.',
                      UserWarning, stacklevel=2)
        output_dir = os.path.dirname(planelog_file)

    if len(output_file) == 0:
        # Set to name of input file witout extension...
        output_file = os.path.splitext(os.path.basename(planelog_file))[0]

    # Combine output directory and output filename for abs path to output files
    savefile = os.path.join(output_dir, output_file)

    # Save output / convert to xarray if asked & print to screen where outpput saved:
    out = _save_outputs(df, info_dict, savefile, as_xarr=as_xarr, overwrite=overwrite)

    return out


def read_and_concat_planelogs(planelog_dir: str, spdb_yaml: str, config_yaml: str,
                              as_xarr: bool = False, convert2_molmol: bool = False,
                              output_dir: str = '', output_file: str = '', overwrite: bool = False):
    """
    Function to read in all plane.log files within a directory and concatenate them
    into either a single pandas dataframe or into an xarray dataset indexed in UTC
    time. If you only want to read in a single file, use read_planelog() instead.

    INPUT:
    -------
    (1) planelog_dir   -  STRING containing absolute path to a directory containing
                          GEOSChem output plane.log files.

    (2) spdb_yaml      -  STRING containing absolute path to the species_database.yml
                          file associated with your run. Used to add metadata.

    (3) config_yaml    -  STRING containing path to geoschem_config.yml file
                          associated with your run. Used to add metadata.

    (4) as_xarr        -  (OPTIONAL) BOOLEAN indicating if you want to output data as
                          an xarray dataset instead of as a pandas dataframe.
                          If this is set to TRUE, the data arrays within the
                          output dataset will have attributes containing metadata
                          assigned to them. Otherwise, output in the xarray
                          dataset is the same as that in the df. Default is set
                          to FALSE (to return a pandas df instead).

    (5) convert2_molmol - (OPTIONAL) BOOLEAN indicating if you wish to convert all
                          outputted concentrations into units of mol/mol.
                          Only variables with concentration units are converted
                          & others are left in their native unit. Default is set
                          to FALSE (not to do conversion).If set to "TRUE" and
                          "TRA_###" appears in plane.log headers, no conversion
                          is actually performed to prevent against improper usage.

    (6) output_dir       - STRING containing absolute path to a directory where
                          concatenated output file will be saved. If none is
                          passed, then output saved in planelog_dir.

    (7) output_file     - STRING containing desired name of ouptput file. If none
                          is passed, the output file will be named:
                              "planelog_concat_[startYYYYMMDD]_[endYYYYMMDD]"

    (8) overwrite       - BOOL indicating that if there is already a file named
                          "output_file" in the output directory, "output_dir", whether
                          that file should be overwritten or if a new unique name
                          should be used.

    NOTE: Outputted advected species concentrations are in molec/cm3 if Tracer
          Names are used, but are in mol/mol if Tracer Numbers are used in input files.
          For more details on why this happens, see GitHub issue #796:
              https://github.com/geoschem/geos-chem/issues/796


    OUTPUT:
    -------
    If output_xarr==FALSE (Default):

        (1) df        - PANDAS DATAFRAME indexed along a pd.datetime object index named 'time_UTC'.
                         Saved at: output_dir+output_file+'.csv'

        (2) info_dict - DICTIONARY containing info about output vars & units.
                        Saved as yaml at: output_dir+output_file+'.yaml'

    If output_xarr == TRUE:

        (1) ds - XARRAY DATASET indexed along a np.datetime64 type coordinate
                 & dimension named 'time_UTC'. All data vars within the ds have
                 individual attributes for "long_name" and "units" assigned to them
                 and will reflect unit conversion changes to mol/mol from molec/cm3
                 if convert2_molmol is set to TRUE. If this is done an additional
                 attribute 'note' is included indicating this occurred.
                 Saved at: output_dir+output_file+'.nc'

    Regardless of the output type selected, the resulting df/ds contains columns/data arrays
    corresponding to all output diagnostics except "YYYYMMDD" and "HHMM" which are
    dropped after constructing the time index/coordinate. Column/data_var names are
    nearly identical to the diagnostics in plane.log, with the exception that any
    '-' characters are replaced with '_' (e.g. 'P-I' in plane.log becomes 'P_I').
    """
    # =========================================================================
    # Check User arguments & raise errors if files/dirs not found
    # =========================================================================
    # Normalize the directory path (remove trailing slash if present)
    planelog_dir = planelog_dir.rstrip('/')

    # Verify that input dirpath exists, otherwise thrown an error:
    if not os.path.isdir(planelog_dir):
        raise FileNotFoundError('The following input directory passed to' +
                                'read_and_concat_planelogs() could not ' +
                                'be found:\n\t' + planelog_dir).with_traceback(sys.exc_info()[2])
    else:
        # Verify that at least 1 file matching plane.log.* exists
        # within the dirpath to read in, otherwise throw an error:
        files = os.listdir(planelog_dir)
        planelog_files = [f for f in files if f.startswith('plane.log.')]
        if len(planelog_files) == 0:
            raise FileNotFoundError('No files matching "plane.log.*" were found' +
                                    'under the following directory to read in:\n\t' +
                                    planelog_dir).with_traceback(sys.exc_info()[2])

    # =========================================================================
    # Loop over all plane.log files in directory & concat data into single df
    # =========================================================================
    for i, file in enumerate(planelog_files):

        # Read raw file in as dataframe (no unit conversions / renaming)
        df_i = _read_planelog_to_df(os.path.join(planelog_dir, file))

        # Add a new column to df indicating the file this data came from:
        df_i['PlaneLog_File'] = np.full(len(df_i), os.path.join(planelog_dir, file))

        # Add a new column to df indicating the index of the file this data came from:
        df_i['File_Index'] = np.full(len(df_i), i + 1)

        if i == 0:
            # W/ first file, define concatenated "df" w/ deep copy of df_i.
            df = df_i.copy()
        else:
            # For all subsequent loops, concatenate a deep copy of df_i with
            # existing data in df_all from other plane.log files. Deep copy is
            # precaution to ensure we're not inadvertently referencing a view.
            df = pd.concat([df, df_i.copy()])

    # =========================================================================
    # Build meta data dictionary with info about all vars in output dataframe.
    # =========================================================================
    # This Assumes advected species units are mol/mol if tracer numbers used, but
    # molec/cm3 if tracer numbers used.
    #
    # **NOTE**: This line MUST remain before renaming df columns with tracer numbers
    #           to tracer names for output units to be correct in info_dict since
    #           assigning them properly requires us to know if the column names
    #           were originally tracer numbers or tracer names.

    info_dict = _build_output_meta_dict(df, config_yaml, spdb_yaml)

    # Get dictionary to map tracer numbers to tracer names:
    nums2names = _get_tracer_name_num_mapping(config_yaml, keys_are_nums=True)

    # Rename all columns in dataframe so columns with tracer numbers are converted
    # to tracer names. Won't affect anything if cols are already named with tracer names.
    df.rename(columns=nums2names, inplace=True)

    # =========================================================================
    # (OPTIONALLY) Convert units to mol/mol from molec/cm3 if tracer names were used.
    # =========================================================================
    # Convert concentrations to mol/mol from molec/cm3 & update units in info_dict
    # accordingly. Also add note if conversion performed or not ...
    if convert2_molmol is True:
        # NOTE: Function will not actually do anything if conc_units!='molec/cm3'
        # which is the case if tracer numbers were used in analogous input files.
        df, info_dict = _convert_moleccm3_to_mr(df, info_dict)

    # =========================================================================
    # Configure output filename/directory:
    # =========================================================================
    # Verify that output dirpath exists, otherwise assign out_dir to planelog_dir.
    if len(output_dir) == 0 or not os.path.isdir(output_dir):
        output_dir = planelog_dir
        if not os.path.isdir(output_dir):
            warnings.warn('The output directory passed to read_and_concat_planelogs() could not '
                          'be found:\n\t' + planelog_dir + '\n Output will be saved at '
                          'input planelog_dir instead.',
                          UserWarning, stacklevel=2)

    if len(output_file) == 0:
        # Create filename from start/end date of planelogs if none passed by extracting
        # the dates from parsed filenames & find the min/max dates
        dates = []
        good_dates = True  # Make empty list to save dates from file
        for file in planelog_files:
            date_str = file.split('.')[-1]  # Get the date part (expecting 'plane.log.YYYYMMDD')
            if len(date_str) == 8:  # Ensure it's in YYYYMMDD format
                dates.append(date_str)
            else:
                good_dates = False  # Update bool to ID that not all dates could be parsed.
                warnings.warn('Dates extracted from plane.log.YYYYMMDD filenames '
                              'were not 8 chars long as expected, which can result '
                              'if you changed the format of the planeflight output '
                              'filenames from their default format in your '
                              'geoschem_config.yml file. So, the output filename of '
                              'the concatenated data will not contain the date range '
                              'of all data concatenated within it.',
                              UserWarning, stacklevel=2)
        if good_dates is True:
            # If we were able to get len(8) dates from all file names, take min/max.
            earliest = min(dates)
            latest = max(dates)
            # And use them to define the name of the output file:
            output_file = 'planelog_concat_' + earliest + '_' + latest
        else:
            # If we couldn't get dates from all file names, make dateless outputfile:
            output_file = 'planelog_concat'

    # Combine output directory and output filename for abs path to output files
    savefile = os.path.join(output_dir, output_file)

    # Save output / convert to xarray if asked & print to screen where outpput saved:
    out = _save_outputs(df, info_dict, savefile, as_xarr=as_xarr, overwrite=overwrite)

    return out
