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



"""Main DeepScore class to analyze depth values.

"""


import subprocess
from pathlib import Path
from sys import exit
from os.path import basename
import os


deepscore_root = str(Path(__file__).parent.parent)

class DeepScore:

    def __init__(self,
                random_id,
                dataset_root=deepscore_root+'/DATASET',
                pdb_user=deepscore_root+'/example/3tw2.pdb',
                overwrite_directory=False):

        from src.version import __version__
        print("\nDeepScore Version:", __version__)


        # Following variables provided to the file depending upon the input
        self.__pdb_user_org = pdb_user

        self.__pdb_user = basename(pdb_user)
        #self.__pdb_chain = pdb_chain
        #self.__pdb_name = self.__pdb_user[:-4]+self.__pdb_chain+'.pdb'
        self.__pdb_name = self.__pdb_user

        # Get a 8 digit random number as parameter
        self.__rnd_num = str(random_id)

        # Sub-dircetories
        self.__root = deepscore_root

        self.__dataset_root = dataset_root

        self.__input_dir = os.path.join(self.__dataset_root, self.__rnd_num )

        self.__depth_exe = os.path.join(self.__root ,'Externals/DEPTH/bin/DEPTH')
        self.__depth_iter = '10'

        self.__models_path = os.path.join(self.__root ,'/src/Models/')

        self.__pdb_file = os.path.join(self.__input_dir , self.__pdb_name)

        self.__sol_file_path = os.path.join(self.__input_dir , 'Solvated_Files/')

        self.__strc_info_file_path = os.path.join(self.__input_dir, 'StrcInfo_Files/')

        self.__strc_info_file = os.path.join(self.__strc_info_file_path, self.__pdb_name, '_Model_Stats.dat')

        if overwrite_directory==False:

            if Path(self.__input_dir).is_dir():
                print("Temporary Directory Exists! Exiting...")
                exit()
            self.Generate_SubDirectories()


    def Generate_SubDirectories(self):
        """
        This function creates the subdirectories in the project DATASET folder.
        The sub-directories will contain a copy of all input and output files
        generated during the program run. Final scores, plots will also be
        stored in this subdirectory. This function is called by default by the
        constructor of the DeepScore class
        
        
        Args:
            DeepScore class instance

        Returns:
            None

        """
        # try creating DATASET directory
        try:
            # Makeing Root Dataset Directory
            subprocess.run("mkdir " + self.__dataset_root, shell=True, check=True)
        except:
            pass

        # Make tmp subdirectories
        print("\nGenerating temporary Sub-directories in: ", self.__input_dir)
        subprocess.run("mkdir " + self.__input_dir, shell=True, check=True)

        print("\nGenerating sub-directories for Depth files")
        subprocess.run("mkdir " + self.__input_dir + "/Depth_Files",
                       shell=True,
                       check=True)

        print(
            "\nGenerating sub-directories for Solvation files and Structural Info files"
        )
        subprocess.run("mkdir " + self.__input_dir + "/Solvated_Files",
                       shell=True,
                       check=True)
        subprocess.run("mkdir " + self.__input_dir + "/StrcInfo_Files",
                       shell=True,
                       check=True)

        print("\nCopying Input PDB file to temporary directory")
        # Copy org pdb file in the tmp location
        subprocess.run("cp " + self.__pdb_user_org + " " +
                       self.__input_dir,
                       shell=True,
                       check=True)

        return


    def Get_Solvated_Models(self):
        """
        This function executes the DEPTH software in the background and measures
        the DEPTH vaues of the protein atoms. The module also generates the
        solvated protein coordinates with one layer of solvent atoms. These
        coordinates are used in the stats modules.

        Args:
            DeepScore class instance

        Returns:
            None

        """

        print("\nSolvating the Protein")
        # Execute DEPTH for the given PDB file
        self.DEPTH_CMD  = self.__depth_exe
        self.DEPTH_CMD += " -i " + os.path.join(self.__input_dir,self.__pdb_name) 
        self.DEPTH_CMD += " -o " + os.path.join(self.__input_dir , self.__pdb_name + "_out")
        self.DEPTH_CMD += " -keep " + os.path.join(self.__input_dir , self.__pdb_name + "_sol")
        self.DEPTH_CMD += " -n " + self.__depth_iter

        subprocess.run(self.DEPTH_CMD, shell=True, check=True)

        # mv log files generated in the current working directory by Depth
        subprocess.run("mv " + str(os.getcwd()) + "/*.log " + self.__input_dir,
                       shell=True)

        return



    def Read_Solvated_Models_Depths(self):
        """
        This function extract the depth values form the depth files.
        
        Args:
            DeepScore class instance

        Returns:
            None

        """
        from src.read_solvated_models import Read_Models_Residue_Depths
        from src.read_solvated_models import Read_Models_Atomic_Depths

        depthfile =  os.path.join(self.__input_dir , "Depth_Files/", self.__pdb_user)

        depscore_score_file_residue = os.path.splitext(depthfile)[0] + '-residue.depth'
        depscore_score_file_atomic = os.path.splitext(depthfile)[0] + '-atom.depth'

        if os.path.isfile(depscore_score_file_residue):
            self.ResidueDepths = Read_Models_Residue_Depths(depscore_score_file_residue)
            
        else:
            print("File doesn't Exists!. Exiting",depscore_score_file_residue)
            exit()


        if os.path.isfile(depscore_score_file_atomic):
            self.AtomicDepths = Read_Models_Atomic_Depths(depscore_score_file_atomic)
            
        else:
            print("File doesn't Exists!. Exiting",depscore_score_file_atomic)
            exit()

        return


    def Plot_Residue_Depths(self):

        import matplotlib.pyplot as plt
        import matplotlib as mpl

        residue_depths = []
        residue_depths_sd = []
        residues = []
        for k,v in self.ResidueDepths.items():
            
            residue_depths.append(v[0])
            residue_depths_sd.append(v[1])
            residues.append(int(k[1]))
            
        plt.scatter(residues,residue_depths,c=residue_depths, s=5 ,cmap=mpl.colormaps['cool'],label='residue depth')

        plt.hlines(5.0, xmin=min(residues), xmax=max(residues), linestyles='dashed')

        depthfile =  os.path.join(self.__input_dir , self.__pdb_user+'_residuedepth.png')
        
        plt.colorbar()
        plt.title("Residue Depths PLot")
        plt.xlabel("Residue Numbers")
        plt.ylabel("Residue Depths")
        plt.savefig(depthfile)
        #plt.show()
        print("Plots are saved in ",depthfile)

        return

    def Plot_Atomic_Depths(self):

        import matplotlib.pyplot as plt
        import matplotlib as mpl

        atomic_depths = []
        atomic_depths_sd = []
        atoms = []
        for k,v in self.AtomicDepths.items():
            
            atomic_depths.append(v[0])
            atomic_depths_sd.append(v[1])
            atoms.append(int(k[1]))
            
        plt.scatter(atoms, atomic_depths,c=atomic_depths,s=1, cmap=mpl.colormaps['cool'],label='atomic depth')

        plt.hlines(5.0, xmin=min(atoms), xmax=max(atoms), linestyles='dashed')

        depthfile =  os.path.join(self.__input_dir , self.__pdb_user+'_atomicdepth.png')
        plt.colorbar()
        plt.title("Atomic Depths Plot")
        plt.xlabel("Residue Number")
        plt.ylabel("Atomic Depths")
        plt.savefig(depthfile)
        #plt.show()
        print("Plots are saved in ",depthfile)

        return






