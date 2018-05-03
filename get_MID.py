#!/usr/bin/env python
# basics
import numpy as np
import pandas as pd
import datetime
import argparse

# for reading mzml files
from pyteomics import mzml

# for parsing ion structures
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops

# global vars
C13_NEUTRON = 13.0033548378 - 12
H2_NEUTRON = 2.01410177811 - 1.00782503224

def parse_arguments():
    """
    Parses arguments from the command line.
    """
    try:    
        parser = argparse.ArgumentParser(description='Gets an MID out of an mzML file given an ion structure.')
        parser.add_argument(
            '-f',
            "--mzml_file", 
            action='store',
            help="path to mzML file from which to extract MID",
            required=True,
            type=str
                        )
        parser.add_argument(
            '-s',
            "--structure", 
            action='store',
            help="SMILES string or path to .sdf file of a single ion structure",
            required=True,
            type=str
                        )
                        
        parser.add_argument(
            '-o',
            "--output_file", 
            action='store',
            help="filename of .csv to write MID data",
            required=False,
            default='out.csv',
            type=str
                        )
        
        parser.add_argument(
            '-r',
            "--rt_start", 
            action='store',
            help="Retention time in minutes at which to start MID extraction",
            required=False,
            default=0,
            type=float
                        )
        
        parser.add_argument(
            '-t',
            "--rt_stop", 
            action='store',
            help="Retention time in minutes at which to end MID extraction",
            required=False,
            default=6,
            type=float
                        )
        
        parser.add_argument(
            '-m',
            "--mz_tol", 
            action='store',
            help="Mass tolerance in parts per million used for MID extraction",
            required=False,
            default=35,
            type=float
                        )
        
        parser.add_argument(
            '-n',
            "--max_neutrons", 
            action='store',
            help="maximum heavy neutrons to consider when extracting MIDs",
            required=False,
            default=5,
            type=int
                        )
        
        parsed = parser.parse_args()
        
        return(parsed)
    except ValueError:
        print('Unable to parse arguments.')


def calc_mz_windows(monoiso_mz, span, ppm, charge=1):
    """
    Finds a set of m/z ranges that contain 13C and 2H isotopologues of a given monoisotopic m/z.
    :param monoiso_mz:   float   the monoisotopic m/z value in Da
    :param span:         int     the maximum number of heavy neutrons to consider
    :param ppm:          float   mass accuracy in parts per million
    :param charge:       int     the assumed (absolute value of) the charge on the monoisotopic m/z
    :return mz_windows:  array   (span+1-by-2) matrix with columns mz_min and mz_max for rows M0, M1 ... 
    """    
    # initialize output matrix
    mz_windows = np.zeros(shape=(span+1, 2),
                          dtype='float')
    
    row_idxs = range(span+1)
    
    # calc ppm factor
    low_ppm = 1-ppm/1e6
    high_ppm = 1+ppm/1e6
    
    # loop through rows
    for m in row_idxs:
        mz_m = np.array((monoiso_mz + m*C13_NEUTRON, monoiso_mz + m*H2_NEUTRON))
        mz_windows[m, :] = np.min(mz_m * low_ppm), np.max(mz_m * high_ppm)
    
    return mz_windows


def get_monoisotopic_mz_and_z(structure):
    """
    Determines the monoisotopic m/z value and charge of an ion provided as a SMILES string or .sdf file.
    :param structure:    str     a valid SMILES string OR a path to an .sdf file containg a single ion structure.
    :return out_dict:    dict    w/ entries "charge" (int) and "monoiso_mz" (float in Daltons) and rdkit mol obj.
    """
    # parse input
    try:
        mol = Chem.MolFromSmiles(structure)
        if mol is None:
            raise TypeError('The provided structure was not a valid SMILES, assuming it is a path to an .sdf file...')
    except TypeError:
        try:
            lst = [mol for mol in Chem.SDMolSupplier(structure)]
            mol = lst[0]
        except OSError:
            raise TypeError('The provide structure was neither a valid SMILES string nor a path to an .sdf file.')
    
    # ensure mol exists
    if not mol:
        raise NotImplementedError('For unknown reasons, the provided structure could not be analyzed.')
    
    # determine properties of mol
    monoiso_mz = rdMolDescriptors.CalcExactMolWt(mol)
    charge = rdmolops.GetFormalCharge(mol)
    
    # ensure provided structure is of an ion
    if not charge:
        raise ValueError('Provided structures must be of ions, not neutral molecules.')
    
    charge = int(charge)
    out_dict = {'charge': charge, 'monoiso_mz': monoiso_mz, 'mol':mol}
    return out_dict


def extract_mid_from_file(mzml_path, mz_windows, rt_window, ppm):
    """
    Pulls out MID info from mzML file, looking only in given mz_windows and rt_window.
    :param mzml_path:    str    path to mzml_file for extraction
    :param mz_windows:   array  (span+1-by-2) matrix with columns mz_min and mz_max for rows M0, M1 ... 
    :param rt_window:    array  (1-by-2) array containing rt_min and rt_max in minutes for peak of interest
    :param ppm:          float  parts-per-million mass accuracy 
    :return out:         dict w/ keys:
                                           m, np.array of ints, the number of heavy neutrons
                                           mean_mz, the intensity-weighted m/z of the m peak
                                           total_i, the total intensity of isotopologue m
    """
    with mzml.read(mzml_path) as reader:
        # identify rt window
        rt_min, rt_max = tuple(rt_window)

        # initialize output matrix of mean mz and intensities
        n_rows = mz_windows.shape[0]
        m_span = range(n_rows)
        
        # initialize lists to store relevant raw data points
        mzs = [[] for m in m_span]
        intensities = [[] for m in m_span]

        # loop through scans, keeping data points in mz_windows and rt_window only
        for spec in reader:
            try:
                rt = spec['scanList']['scan'][0]['scan start time']
            except (KeyError, IndexError):
                continue
            if rt >= rt_min and rt <= rt_max:
                # get raw scan data
                these_mzs = spec['m/z array']
                these_intensities = spec['intensity array']

                # index into mz_windows to find relevant data in scan
                index_mat = np.searchsorted(these_mzs, mz_windows)
                start = index_mat[:, 0]
                stop = index_mat[:, 1]

                for m in m_span:
                    # if scan has no mz values of interest, skip it
                    if start[m] != stop[m]:
                        mzs[m].extend(list(these_mzs[start[m]:stop[m]]))
                        intensities[m].extend(list(these_intensities[start[m]:stop[m]]))

        mean_mz = np.asarray([np.average(mzs[m], weights=intensities[m]) for m in m_span])
        total_i = np.asarray([np.sum(intensities[m]) for m in m_span])

        return({'m': np.asarray(m_span), 
               'mean_mz': mean_mz, 
               'total_i': total_i})

def main():
    """
    From an ion structure, mzML file, and RT window, determine an MID.
    Outputs a .csv file and prints a pandas dataframe.  
    The dataframe has columns:
        date,         str,   date this script was run
        ppm,          float, ppm used for mass accuracy
        ion_smiles,   str,   analyzed ion structure
        m,            int,  integer mass isotopologue of interest
        c13_theo_mz,  float, m/z of ion with only m 13C atoms
        h2_theo_mz,   float, m/z of ion with only m 2H 
        mean_mz,      float, mean m/z observed in data
        raw_intensity,float, total intensity in mz_windows and rt_window
        mid,          float, normalized raw intensities
    """
    #unpack arguments
    args = parse_arguments()
    
    # determine monoisotopic m/z and charge from structure
    structure = args.structure
    ion_info = get_monoisotopic_mz_and_z(structure)
    monoiso_mz = ion_info['monoiso_mz']
    charge = ion_info['charge']
    
    # determine mz_windows from monoiso_mz
    max_neutrons = args.max_neutrons
    ppm = args.mz_tol
    
    mz_windows = calc_mz_windows(monoiso_mz, 
                                 span=max_neutrons, 
                                 ppm=ppm, 
                                 charge=charge)
        
    # calculate MID raw data
    rt_window = (args.rt_start, args.rt_stop)
    mzml_file = args.mzml_file
    raw_mid = extract_mid_from_file(mzml_file, 
                                    mz_windows=mz_windows,
                                    rt_window=rt_window,
                                    ppm=ppm)
    
    # preprare output
    today = datetime.datetime.today()
    m = np.arange(max_neutrons + 1)
    C13_theo_mz = monoiso_mz + m*C13_NEUTRON
    H2_theo_mz = monoiso_mz + m*H2_NEUTRON
    raw_ints = raw_mid['total_i']
    n_rows = len(m)
    
    
    out_df = pd.DataFrame({'m': m,
                           'mzml_file': mzml_file,
                           'structure': structure,
                           'analysis_date': today,
                           'C13_theo_mz': C13_theo_mz,
                           'H2_theo_mz': H2_theo_mz,
                           'mean_mz': raw_mid['mean_mz'],
                           'raw_intensity': raw_ints,
                           'mid': raw_ints / raw_ints.sum()
                           })
    print(out_df)
    out_df.to_csv(args.output_file)

if __name__ == '__main__':
    main()