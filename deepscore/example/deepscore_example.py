################################################################################
#   Copyright (C) 2016 Neelesh Soni <neeleshsoni03@gmail.com>, 
#   <neelesh.soni@alumni.iiserpune.ac.in>
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################


import sys
from pathlib import Path
import os.path
from random import random
from src.DeepScoreClass import DeepScore

deepscore_root = str(Path(__file__).parent.parent)


def example_run( data_directory: str = deepscore_root+'/DATASET',
                     rand_id: str = "tmp_"+str(int(random() * 100000000)),
                     example_pdb: str = deepscore_root+'/example/2grn.pdb',
                     example_chain: str = 'A',
                     skip_depth_run = False,
                     overwrite_directory = False):
    """
    Example function to demonstrate the DeepScore functionality.

    Args:
        data_directory (str): Data directory of new files
        rand_id (str): Unique string for identifying the output
        example_pdb (str): PDB format input file for scoring
        example_chain (str): Chain id of the PDB file for consideration
        skip_depth_run (bool): Compute depth or Not
        overwrite_directory (bool): Overwrite output directory or not

    Returns:
        The return value. True for success, False otherwise.

    """

    # GARLIC class instantiation
    ds_obj = DeepScore(rand_id, data_directory, example_pdb,
                        example_chain, overwrite_directory)
    

    if skip_depth_run==False:
        # Get solvated model
        ds_obj.Get_Solvated_Models()


    ds_obj.Read_Solvated_Models_Depths()

    #DEPTH values are stored in ds_obj.ResidueDepths python dictionary.
    #Key-value pairs: key: (chain,residuenumber), value: [DEPTH value , Standard deviation]
    
    #for k,v in ds_obj.ResidueDepths.items():
    #    print(k,v)

    #for k,v in ds_obj.AtomicDepths.items():
    #    print(k,v)

    for Residue_ID, Depth_vals in ds_obj.ResidueDepths.items():
        chain,residuenumber = Residue_ID
        DEPTH_value , Depth_SD =  Depth_vals
    
    for Residue_ID, Depth_vals in ds_obj.AtomicDepths.items():
        chain,residuenumber,atomtype = Residue_ID
        DEPTH_value , Depth_SD =  Depth_vals
        

    ds_obj.Plot_Residue_Depths()

    ds_obj.Plot_Atomic_Depths()


    #TODO:
    #Add modules for analysis, Extract from GARLIC code


    return True






