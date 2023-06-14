#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Stephen Szwiec, 2023 June 14
#
# invoke on a folder with .mol files to get rdkit descriptors in a CSV
# fragments optional with -f or --fragments flag
# example: python ./free_descriptors.py -i /path/to/molfiles -o /path/to/output.csv
#
# several options exist to compute additional descriptors:
# - MACCS keys: -M or --MACCS flag
# - ECFP6 fingerprints: -E or --ECFP6 flag
# - Mordred descriptors: -m or --mordred flag
# - Macrocycle descriptors: -c or --macrocycle flag
#
# Optional descriptors have a tendency to take a long time to compute and to fail on some molecules
#
# This script is a modification of the rdkit_descriptors.py script by Petr Skoda
#
# This script requires the RDKit and Mordred packages to be installed
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
# The full text of the GNU General Public License is available at <https://www.gnu.org/licenses/gpl-3.0.html>.

"""Compute RDKit descriptors for molecules/fragments.
Usage:
    python free_descriptors.py
        -i {path to molecules for reading}
        -o {path to output CSV file}
        --fragments Use fragments else use molecules.
                    Default is to use molecules.
This file can be also imported as a python script. In such case please
use the extract_fragments method.
"""
# region Imports
import os
import argparse
import logging
import glob
import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem, MACCSkeys
from mordred import Calculator, descriptors
from mordred.RingCount import RingCount
# endregion Imports

# region Metadata
__author__ = 'Stephen Szwiec'
__license__ = 'GPLv3'
__email__ = 'stephen.szwiec@ndus.edu'
# endregion Metadata

# region Classes
class MACCS:
    """
    Compute MACCS keys for molecules. The output is a CSV file.
    """
    def __init__(self, smiles):
        """
        Initialize the class.

        :param smiles: list of SMILES strings
        :type smiles: list

        :return: None
        """
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles

    def compute_MACCS(self, name):
        """
        Compute MACCS keys for molecules. The output is a CSV file.

        :param name: name of the output file
        :type name: str

        :return: None
        """
        MACCS_list = []
        header = ['bit' + str(i) for i in range(167)]
        for i in range(len(self.mols)):
            ds = list(MACCSkeys.GenMACCSKeys(self.mols[i]).ToBitString())
            MACCS_list.append(ds)
        df = pd.DataFrame(MACCS_list,columns=header)
        df.insert(loc=0, column='smiles', value=self.smiles)
        df.to_csv(name[:-4]+'_MACCS.csv', index=False)

class ECFP6:
    """
    Compute ECFP6 fingerprints for molecules. The output is a CSV file.
    """
    def __init__(self, smiles):
        """
        Initialize the class.

        :param smiles: list of SMILES strings
        :type smiles: list

        :return: None
        """
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles

    def mol2fp(self, mol, radius = 3):
        """
        Output array of ECFP6 fingerprints for a molecule for a given radius.

        :param mol: molecule
        :type mol: rdkit.Chem.rdchem.Mol

        :param radius: radius of the fingerprint
        :type radius: int

        :return: array of ECFP6 fingerprints
        :rtype: numpy.ndarray
        """
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius = radius)
        array = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, array)
        return array

    def compute_ECFP6(self, name):
        """
        Compute ECFP6 fingerprints for molecules. The output is a CSV file.

        :param name: name of the output file
        :type name: str

        :return: None
        """
        bit_headers = ['bit' + str(i) for i in range(2048)]
        arr = np.empty((0,2048), int).astype(int)
        for i in self.mols:
            fp = self.mol2fp(i)
            arr = np.vstack((arr, fp))
        df_ecfp6 = pd.DataFrame(np.asarray(arr).astype(int),columns=bit_headers)
        df_ecfp6.insert(loc=0, column='smiles', value=self.smiles)
        df_ecfp6.to_csv(name[:-4]+'_ECFP6.csv', index=False)

class Macrocycle_Descriptors:
    """
    Compute Macrocycle descriptors for molecules. The output is a CSV file.
    """
    def __init__(self, smiles):
        """
        Initialize the class.

        :param smiles: list of SMILES strings
        :type smiles: list

        :return: None
        """
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles
        self.mordred = None


    def compute_ringsize(self, mol):
        """
        check for macrolides of RS 3 to 99, return a  list of ring counts.
        [RS3,RS4,.....,RS99]
        [0,0,0,...,1,...,0]

        :param mol: molecule
        :type mol: rdkit.Chem.rdchem.Mol

        :return: list of ring counts
        :rtype: list
        """
        RS_3_99 = [i+3 for i in range(97)]
        RS_count = []
        for j in RS_3_99:
            RS = RingCount(order=j)(mol)
            RS_count.append(RS)
        return RS_count

    def macrolide_ring_info(self):
        """
        Compute ring sizes for macrolides; outputs a pandas dataframe.

        :return: pandas dataframe
        """
        headers = ['n'+str(i+13)+'Ring' for i in range(87)]+['SmallestRS','LargestRS']
        # up to nR12 is already with mordred, start with nR13 to nR99
        ring_sizes = []
        for i in range(len(self.mols)):
            RS = self.compute_ringsize(self.mols[i])  # nR3 to nR99
            RS_12_99 = RS[9:]    # start with nR12 up to nR99
            ring_indices = [i for i,x in enumerate(RS_12_99) if x!=0]  # get index if item isn't equal to 0
            # if there is a particular ring present, the frequency won't be zero. Find those indexes.
            if ring_indices:
                # Add 12 (starting ring count) to get up to the actual ring size
                smallest_RS = ring_indices[0]+12     # Retrieve the first index (for the smallest core RS - note the list is in ascending order)
                largest_RS = ring_indices[-1]+12	 # Retrieve the last index (for the largest core RS)
                RS_12_99.append(smallest_RS)  # Smallest RS
                RS_12_99.append(largest_RS)  # Largest RS
            else:
                RS_12_99.extend(['',''])
            ring_sizes.append(RS_12_99[1:]) # up to nR12 is already with mordred, start with nR13 to nR99
        df = pd.DataFrame(ring_sizes, columns=headers)
        return df

    def sugar_count(self):
        """
        Compute sugar counts for macrolides; outputs a pandas dataframe.

        :return: pandas dataframe
        """
        sugar_patterns = [
        '[OX2;$([r5]1@C@C@C(O)@C1),$([r6]1@C@C@C(O)@C(O)@C1)]',
        '[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]',
        '[OX2;$([r5]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C1),$([r6]1@C(!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C(O)@C(O)@C1)]',
        '[OX2;$([r5]1@C(!@[OX2H1])@C@C@C1),$([r6]1@C(!@[OX2H1])@C@C@C@C1)]',
        '[OX2;$([r5]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]',
        '[OX2;$([r5]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C1),$([r6]1@[C@](!@[OX2,NX3,SX2,FX1,ClX1,BrX1,IX1])@C@C@C@C1)]',
        ]
        sugar_mols = [Chem.MolFromSmarts(i) for i in sugar_patterns]
        sugar_counts = []
        for i in self.mols:
            matches_total = []
            for s_mol in sugar_mols:
                raw_matches = i.GetSubstructMatches(s_mol)
                matches = list(sum(raw_matches, ()))
                if matches not in matches_total and len(matches) !=0:
                    matches_total.append(matches)
            sugar_indices = set((list(itertools.chain(*matches_total))))
            count = len(sugar_indices)
            sugar_counts.append(count)
        df = pd.DataFrame(sugar_counts, columns=['nSugars'])
        return df

    def core_ester_count(self):
        """
        Returns pandas frame containing the count of esters in core rings of >=12 membered macrocycles.

        :return: pandas dataframe
        """
        ester_smarts = '[CX3](=[OX1])O@[r;!r3;!r4;!r5;!r6;!r7;!r8;!r9;!r10;!r11]'
        core_ester = []
        ester_mol = Chem.MolFromSmarts(ester_smarts)
        for i in self.mols:
            ester_count = len(i.GetSubstructMatches(ester_mol))
            core_ester.append(ester_count)
        df = pd.DataFrame(core_ester, columns=['core_ester'])
        return df

    def mordred_compute(self, name):
        """
        Compute mordred descriptors for a list of molecules. Outputs a CSV file.

        :param name: name of the file to be saved
        :type name: str


        :return: none
        """
        calc = Calculator(descriptors, ignore_3D=True)
        df = calc.pandas(self.mols)
        self.mordred = df
        df.insert(loc=0, column='smiles', value=self.smiles)
        df.to_csv(name[:-4]+'_mordred.csv', index=False)

    def compute_mordred_macrocycle(self, name):
        """
        Compute mordred descriptors for a list of molecules. Outputs a CSV file.

        :param name: name of the file to be saved
        :type name: str

        :return: none
        """
        if not isinstance(self.mordred, pd.DataFrame):
            self.mordred = self.mordred_compute(name)
        ring_df = self.macrolide_ring_info()
        sugar_df = self.sugar_count()
        ester_df = self.core_ester_count()
#        self.mrc = pd.concat([ring_df,sugar_df, ester_df], axis=1)
        mordred_mrc = pd.concat([self.mordred, ring_df,sugar_df, ester_df], axis=1)
        mordred_mrc.to_csv(name[:-4]+'_mordred_mrc.csv', index=False)
# endregion Classes

# region Descriptors definition
_names = [
    'MolWt',
    'HeavyAtomMolWt',
    'ExactMolWt',
    'NumValenceElectrons',
    'NumRadicalElectrons',
    'MaxEStateIndex',
    'MinEStateIndex',
    'MaxAbsEStateIndex',
    'MinAbsEStateIndex',
    'BalabanJ',
    'BertzCT',
    'Chi0',
    'Chi0n',
    'Chi0v',
    'Chi1',
    'Chi1n',
    'Chi1v',
    'Chi2n',
    'Chi2v',
    'Chi3n',
    'Chi3v',
    'Chi4n',
    'Chi4v',
    'EState_VSA1',
    'EState_VSA10',
    'EState_VSA11',
    'EState_VSA2',
    'EState_VSA3',
    'EState_VSA4',
    'EState_VSA5',
    'EState_VSA6',
    'EState_VSA7',
    'EState_VSA8',
    'EState_VSA9',
    'FractionCSP3',
    'HallKierAlpha',
    'HeavyAtomCount',
    'Ipc',
    'Kappa1',
    'Kappa2',
    'Kappa3',
    'LabuteASA',
    'MolLogP',
    'MolMR',
    'NHOHCount',
    'NOCount',
    'NumAliphaticCarbocycles',
    'NumAliphaticHeterocycles',
    'NumAliphaticRings',
    'NumAromaticCarbocycles',
    'NumAromaticHeterocycles',
    'NumAromaticRings',
    'NumHAcceptors',
    'NumHDonors',
    'NumHeteroatoms',
    'NumRotatableBonds',
    'NumSaturatedCarbocycles',
    'NumSaturatedHeterocycles',
    'NumSaturatedRings',
    'PEOE_VSA1',
    'PEOE_VSA10',
    'PEOE_VSA11',
    'PEOE_VSA12',
    'PEOE_VSA13',
    'PEOE_VSA14',
    'PEOE_VSA2',
    'PEOE_VSA3',
    'PEOE_VSA4',
    'PEOE_VSA5',
    'PEOE_VSA6',
    'PEOE_VSA7',
    'PEOE_VSA8',
    'PEOE_VSA9',
    'RingCount',
    'SMR_VSA1',
    'SMR_VSA10',
    'SMR_VSA2',
    'SMR_VSA3',
    'SMR_VSA4',
    'SMR_VSA5',
    'SMR_VSA6',
    'SMR_VSA7',
    'SMR_VSA8',
    'SMR_VSA9',
    'SlogP_VSA1',
    'SlogP_VSA10',
    'SlogP_VSA11',
    'SlogP_VSA12',
    'SlogP_VSA2',
    'SlogP_VSA3',
    'SlogP_VSA4',
    'SlogP_VSA5',
    'SlogP_VSA6',
    'SlogP_VSA7',
    'SlogP_VSA8',
    'SlogP_VSA9',
    'TPSA',
    'VSA_EState1',
    'VSA_EState10',
    'VSA_EState2',
    'VSA_EState3',
    'VSA_EState4',
    'VSA_EState5',
    'VSA_EState6',
    'VSA_EState7',
    'VSA_EState8',
    'VSA_EState9',
    'fr_Al_COO',
    'fr_Al_OH',
    'fr_Al_OH_noTert',
    'fr_ArN',
    'fr_Ar_COO',
    'fr_Ar_N',
    'fr_Ar_NH',
    'fr_Ar_OH',
    'fr_COO',
    'fr_COO2',
    'fr_C_O',
    'fr_C_O_noCOO',
    'fr_C_S',
    'fr_HOCCN',
    'fr_Imine',
    'fr_NH0',
    'fr_NH1',
    'fr_NH2',
    'fr_N_O',
    'fr_Ndealkylation1',
    'fr_Ndealkylation2',
    'fr_Nhpyrrole',
    'fr_SH',
    'fr_aldehyde',
    'fr_alkyl_carbamate',
    'fr_alkyl_halide',
    'fr_allylic_oxid',
    'fr_amide',
    'fr_amidine',
    'fr_aniline',
    'fr_aryl_methyl',
    'fr_azide',
    'fr_azo',
    'fr_barbitur',
    'fr_benzene',
    'fr_benzodiazepine',
    'fr_bicyclic',
    'fr_diazo',
    'fr_dihydropyridine',
    'fr_epoxide',
    'fr_ester',
    'fr_ether',
    'fr_furan',
    'fr_guanido',
    'fr_halogen',
    'fr_hdrzine',
    'fr_hdrzone',
    'fr_imidazole',
    'fr_imide',
    'fr_isocyan',
    'fr_isothiocyan',
    'fr_ketone',
    'fr_ketone_Topliss',
    'fr_lactam',
    'fr_lactone',
    'fr_methoxy',
    'fr_morpholine',
    'fr_nitrile',
    'fr_nitro',
    'fr_nitro_arom',
    'fr_nitro_arom_nonortho',
    'fr_nitroso',
    'fr_oxazole',
    'fr_oxime',
    'fr_para_hydroxylation',
    'fr_phenol',
    'fr_phenol_noOrthoHbond',
    'fr_phos_acid',
    'fr_phos_ester',
    'fr_piperdine',
    'fr_piperzine',
    'fr_priamide',
    'fr_prisulfonamd',
    'fr_pyridine',
    'fr_quatN',
    'fr_sulfide',
    'fr_sulfonamd',
    'fr_sulfone',
    'fr_term_acetylene',
    'fr_tetrazole',
    'fr_thiazole',
    'fr_thiocyan',
    'fr_thiophene',
    'fr_unbrch_alkane',
    'fr_urea'
]

# http://www.rdkit.org/docs/api/rdkit.Chem.Descriptors-module.html
_functions = [
    Descriptors.MolWt,
    Descriptors.HeavyAtomMolWt,
    Descriptors.ExactMolWt,
    Descriptors.NumValenceElectrons,
    Descriptors.NumRadicalElectrons,
    Descriptors.MaxEStateIndex,
    Descriptors.MinEStateIndex,
    Descriptors.MaxAbsEStateIndex,
    Descriptors.MinAbsEStateIndex,
    Descriptors.BalabanJ,
    Descriptors.BertzCT,
    Descriptors.Chi0,
    Descriptors.Chi0n,
    Descriptors.Chi0v,
    Descriptors.Chi1,
    Descriptors.Chi1n,
    Descriptors.Chi1v,
    Descriptors.Chi2n,
    Descriptors.Chi2v,
    Descriptors.Chi3n,
    Descriptors.Chi3v,
    Descriptors.Chi4n,
    Descriptors.Chi4v,
    Descriptors.EState_VSA1,
    Descriptors.EState_VSA10,
    Descriptors.EState_VSA11,
    Descriptors.EState_VSA2,
    Descriptors.EState_VSA3,
    Descriptors.EState_VSA4,
    Descriptors.EState_VSA5,
    Descriptors.EState_VSA6,
    Descriptors.EState_VSA7,
    Descriptors.EState_VSA8,
    Descriptors.EState_VSA9,
    Descriptors.FractionCSP3,
    Descriptors.HallKierAlpha,
    Descriptors.HeavyAtomCount,
    Descriptors.Ipc,
    Descriptors.Kappa1,
    Descriptors.Kappa2,
    Descriptors.Kappa3,
    Descriptors.LabuteASA,
    Descriptors.MolLogP,
    Descriptors.MolMR,
    Descriptors.NHOHCount,
    Descriptors.NOCount,
    Descriptors.NumAliphaticCarbocycles,
    Descriptors.NumAliphaticHeterocycles,
    Descriptors.NumAliphaticRings,
    Descriptors.NumAromaticCarbocycles,
    Descriptors.NumAromaticHeterocycles,
    Descriptors.NumAromaticRings,
    Descriptors.NumHAcceptors,
    Descriptors.NumHDonors,
    Descriptors.NumHeteroatoms,
    Descriptors.NumRotatableBonds,
    Descriptors.NumSaturatedCarbocycles,
    Descriptors.NumSaturatedHeterocycles,
    Descriptors.NumSaturatedRings,
    Descriptors.PEOE_VSA1,
    Descriptors.PEOE_VSA10,
    Descriptors.PEOE_VSA11,
    Descriptors.PEOE_VSA12,
    Descriptors.PEOE_VSA13,
    Descriptors.PEOE_VSA14,
    Descriptors.PEOE_VSA2,
    Descriptors.PEOE_VSA3,
    Descriptors.PEOE_VSA4,
    Descriptors.PEOE_VSA5,
    Descriptors.PEOE_VSA6,
    Descriptors.PEOE_VSA7,
    Descriptors.PEOE_VSA8,
    Descriptors.PEOE_VSA9,
    Descriptors.RingCount,
    Descriptors.SMR_VSA1,
    Descriptors.SMR_VSA10,
    Descriptors.SMR_VSA2,
    Descriptors.SMR_VSA3,
    Descriptors.SMR_VSA4,
    Descriptors.SMR_VSA5,
    Descriptors.SMR_VSA6,
    Descriptors.SMR_VSA7,
    Descriptors.SMR_VSA8,
    Descriptors.SMR_VSA9,
    Descriptors.SlogP_VSA1,
    Descriptors.SlogP_VSA10,
    Descriptors.SlogP_VSA11,
    Descriptors.SlogP_VSA12,
    Descriptors.SlogP_VSA2,
    Descriptors.SlogP_VSA3,
    Descriptors.SlogP_VSA4,
    Descriptors.SlogP_VSA5,
    Descriptors.SlogP_VSA6,
    Descriptors.SlogP_VSA7,
    Descriptors.SlogP_VSA8,
    Descriptors.SlogP_VSA9,
    Descriptors.TPSA,
    Descriptors.VSA_EState1,
    Descriptors.VSA_EState10,
    Descriptors.VSA_EState2,
    Descriptors.VSA_EState3,
    Descriptors.VSA_EState4,
    Descriptors.VSA_EState5,
    Descriptors.VSA_EState6,
    Descriptors.VSA_EState7,
    Descriptors.VSA_EState8,
    Descriptors.VSA_EState9,
    Descriptors.fr_Al_COO,
    Descriptors.fr_Al_OH,
    Descriptors.fr_Al_OH_noTert,
    Descriptors.fr_ArN,
    Descriptors.fr_Ar_COO,
    Descriptors.fr_Ar_N,
    Descriptors.fr_Ar_NH,
    Descriptors.fr_Ar_OH,
    Descriptors.fr_COO,
    Descriptors.fr_COO2,
    Descriptors.fr_C_O,
    Descriptors.fr_C_O_noCOO,
    Descriptors.fr_C_S,
    Descriptors.fr_HOCCN,
    Descriptors.fr_Imine,
    Descriptors.fr_NH0,
    Descriptors.fr_NH1,
    Descriptors.fr_NH2,
    Descriptors.fr_N_O,
    Descriptors.fr_Ndealkylation1,
    Descriptors.fr_Ndealkylation2,
    Descriptors.fr_Nhpyrrole,
    Descriptors.fr_SH,
    Descriptors.fr_aldehyde,
    Descriptors.fr_alkyl_carbamate,
    Descriptors.fr_alkyl_halide,
    Descriptors.fr_allylic_oxid,
    Descriptors.fr_amide,
    Descriptors.fr_amidine,
    Descriptors.fr_aniline,
    Descriptors.fr_aryl_methyl,
    Descriptors.fr_azide,
    Descriptors.fr_azo,
    Descriptors.fr_barbitur,
    Descriptors.fr_benzene,
    Descriptors.fr_benzodiazepine,
    Descriptors.fr_bicyclic,
    Descriptors.fr_diazo,
    Descriptors.fr_dihydropyridine,
    Descriptors.fr_epoxide,
    Descriptors.fr_ester,
    Descriptors.fr_ether,
    Descriptors.fr_furan,
    Descriptors.fr_guanido,
    Descriptors.fr_halogen,
    Descriptors.fr_hdrzine,
    Descriptors.fr_hdrzone,
    Descriptors.fr_imidazole,
    Descriptors.fr_imide,
    Descriptors.fr_isocyan,
    Descriptors.fr_isothiocyan,
    Descriptors.fr_ketone,
    Descriptors.fr_ketone_Topliss,
    Descriptors.fr_lactam,
    Descriptors.fr_lactone,
    Descriptors.fr_methoxy,
    Descriptors.fr_morpholine,
    Descriptors.fr_nitrile,
    Descriptors.fr_nitro,
    Descriptors.fr_nitro_arom,
    Descriptors.fr_nitro_arom_nonortho,
    Descriptors.fr_nitroso,
    Descriptors.fr_oxazole,
    Descriptors.fr_oxime,
    Descriptors.fr_para_hydroxylation,
    Descriptors.fr_phenol,
    Descriptors.fr_phenol_noOrthoHbond,
    Descriptors.fr_phos_acid,
    Descriptors.fr_phos_ester,
    Descriptors.fr_piperdine,
    Descriptors.fr_piperzine,
    Descriptors.fr_priamide,
    Descriptors.fr_prisulfonamd,
    Descriptors.fr_pyridine,
    Descriptors.fr_quatN,
    Descriptors.fr_sulfide,
    Descriptors.fr_sulfonamd,
    Descriptors.fr_sulfone,
    Descriptors.fr_term_acetylene,
    Descriptors.fr_tetrazole,
    Descriptors.fr_thiazole,
    Descriptors.fr_thiocyan,
    Descriptors.fr_thiophene,
    Descriptors.fr_unbrch_alkane,
    Descriptors.fr_urea
]
# endregion Descriptors definition

# region Functions definition
def _read_configuration():
    """Get and return application settings.
    :return:
    """
    parser = argparse.ArgumentParser(
        description='Compute RDKit descriptors for given'
                    'molecules/fragments.')
    parser.add_argument('-i', type=str, dest='input',
                        help='input directory', required=True),
    parser.add_argument('-o', type=str, dest='output',
                        help='output CSV file', required=True)
    parser.add_argument('-f', dest='fragments',
                        help='use fragments instead of molecules',
                        action='store_true', required=False)
    parser.add_argument('--fragments', dest='fragments',
                        help='use fragments instead of molecules',
                        action='store_true', required=False)
    parser.add_argument('-M', dest='MACCS',
                        help="compute MACCS keys",
                        action='store_true', required=False)
    parser.add_argument('--MACCS', dest='MACCS',
                        help="compute MACCS keys",
                        action='store_true', required=False)
    parser.add_argument('-E', dest='ECFP6',
                        help="compute ECFP6 fingerprints",
                        action='store_true', required=False)
    parser.add_argument('--ECFP6', dest='ECFP6',
                        help="compute ECFP6 fingerprints",
                        action='store_true', required=False)
    parser.add_argument('-m', dest='mordred',
                        help="compute Mordred descriptors",
                        action='store_true', required=False)
    parser.add_argument('--mordred', dest='mordred',
                        help="compute Mordred descriptors",
                        action='store_true', required=False)
    parser.add_argument('-c', dest='macrocycle',
                        help="compute macrocycle descriptors",
                        action='store_true', required=False)
    parser.add_argument('--macrocycle', dest='macrocycle',
                        help="compute macrocycle descriptors",
                        action='store_true', required=False)

    return vars(parser.parse_args())


def compute_descriptors(input_path, output_file, use_fragments, use_MACCS, use_ECFP6, use_mordred, use_macrocycle):
    """Compute descriptors for molecules/fragments in given input file.

    :param input_path: Input path: directory to *.mol files.
    :param output_file: Output file.
    :param use_fragments: If true use fragments instead of molecules.

    :return: Summary object.
    """
    data = []
    os.chdir(input_path)
    for file in glob.glob("*.mol"):
        data.append(Chem.MolFromMolFile(file))
    # Gather data.
    smiles_set = set()
    if use_fragments:
        for molecule in data:
            for fragment in molecule['fragments']:
                if not fragment['smiles'] in smiles_set:
                    smiles_set.add(fragment['smiles'])
    else:
        for molecule in data:
                smiles_set.add(Chem.MolToSmiles(molecule))
    # Compute and write descriptors.
    sanitize_operation = rdkit.Chem.SanitizeFlags.SANITIZE_ALL ^ \
                         rdkit.Chem.SanitizeFlags.SANITIZE_KEKULIZE
    number_of_invalid = 0
    output_file = output_file + '.csv'
    with open(output_file, 'w') as stream:
        stream.write('smiles,')
        stream.write(','.join(_names))
        stream.write('\n')
        counter = 0
        counter_step = int(len(smiles_set) / 10)
        if counter_step == 0:
            counter_step = 1
        for smiles in smiles_set:
            if counter % counter_step == 0:
                logging.info('%d/%d', counter, len(smiles_set))
            counter += 1
            # SMILES.
            stream.write('"')
            stream.write(smiles)
            stream.write('",')
            # Construct molecule, compute and write properties.
            molecule = rdkit.Chem.MolFromSmiles(str(smiles), sanitize=False)
            # Do not kekulize molecule.
            rdkit.Chem.SanitizeMol(molecule, sanitizeOps=sanitize_operation)
            # Check if molecule is valid.
            if molecule is None:
                logging.error('Invalid molecule detected: %s', smiles)
                number_of_invalid += 1
                continue
            print([str(fnc(molecule)) for fnc in _functions])
            stream.write(','.join([str(fnc(molecule)) for fnc in _functions]))
            stream.write('\n')
    # Optional descriptors computation and writing.
    if use_MACCS or use_ECFP6 or use_mordred or use_macrocycle:
        smiles_list = list(smiles_set)
    if use_MACCS:
        M = MACCS(smiles_list)
        try:
            M.compute_MACCS(output_file)
        except:
            logging.info('MACCS failed on one or more molecules')
    if use_ECFP6:
        E = ECFP6(smiles_list)
        try:
            E.compute_ECFP6(output_file)
        except:
            logging.info('ECFP6 failed on one or more molecules')
    if use_mordred:
        m = Macrocycle_Descriptors(smiles_list)
        try:
            m.mordred_compute(output_file)
        except:
            logging.info('Mordred failed on one or more molecules')
    if use_macrocycle:
        c = Macrocycle_Descriptors(smiles_list)
        try:
            c.compute_mordred_macrocycle(output_file)
        except:
            logging.info('Macrocycle failed on one or more molecules')
    # Log and return summary.
    logging.info('Invalid molecules: %d/%d', number_of_invalid, len(smiles_set))
    return {
        'number_of_invalid': number_of_invalid,
        'total': len(smiles_set)
    }
# endregion Functions definition


# region Main
def _main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(module)s - %(message)s',
        datefmt='%H:%M:%S')
    configuration = _read_configuration()
    use_fragments = 'fragments' in configuration and configuration['fragments']
    use_MACCS = 'MACCS' in configuration and configuration['MACCS']
    use_ECFC6 = 'ECFP6' in configuration and configuration['ECFP6']
    use_mordred = 'mordred' in configuration and configuration['mordred']
    use_macrocycle = 'macrocycle' in configuration and configuration['macrocycle']
    compute_descriptors(configuration['input'], configuration['output'], use_fragments, use_MACCS, use_ECFC6, use_mordred, use_macrocycle)

if __name__ == '__main__':
    _main()
# endregion Main
