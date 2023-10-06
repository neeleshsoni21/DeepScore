
**********************
Software Information
**********************

DeepScore is a software to analyze DEPTH values generated from DEPTH software.

**********************
Installation
**********************
Please refer to INSTALL file in the distribution.

**********************
Usage and Requirements
**********************
DEPTH code: if depth-generated files are not provided as input


**********************
Compiling DEPTH: See Also; http://cospi.iiserpune.ac.in/depth/download/
**********************

cd deepscore/Externals/DEPTH/bin
make


To run the DeepScore with an example file, use the following:

python simple_API.py

This will plot the residue-wise depth values and atomwise depth values. The actual values can be accessed using DeepScore object as follows:

#DEPTH values are stored in ds_obj.ResidueDepths python dictionary.
#Key-value pairs: key: (chain,residuenumber), value: [DEPTH value , Standard deviation]

    for Residue_ID, Depth_vals in ds_obj.ResidueDepths.items():
        chain,residuenumber = Residue_ID
        DEPTH_value , Depth_SD =  Depth_vals
    
    for Residue_ID, Depth_vals in ds_obj.AtomicDepths.items():
        chain,residuenumber,atomtype = Residue_ID
        DEPTH_value , Depth_SD =  Depth_vals