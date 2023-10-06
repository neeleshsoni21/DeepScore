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


import deepscore as ds


#Run DeepScore by first generating input files using DEPTH and then processing
ds.example_run(data_directory = './output/', 
    rand_id='pdb_2grn', 
    skip_depth_run = False, 
    overwrite_directory=True)


#Run DeepScore by skipping DEPTH runs (This works only if you provide depth files seperately)
ds.example_run(data_directory = './output/', 
    rand_id='pdb_2grn', 
    skip_depth_run = True, 
    overwrite_directory=True)


