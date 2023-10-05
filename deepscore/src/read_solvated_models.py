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


'''
This function parse the PDB file and return the coordinates of all atoms
'''
from collections import OrderedDict

def Read_Models_Residue_Depths(fname):

    ResidueDepths = OrderedDict()

    with open(fname,'r') as inf:
        lines = inf.readlines()
        for l in lines:

            if l[0]=='#':
                continue

            if len(l)<=1:
                continue

            toks = l.strip().split()

            chain , residue = toks[0].split(':')

            ResidueDepths[(chain,residue)]=[float(toks[2]), float(toks[3])]

    return ResidueDepths


def Read_Models_Atomic_Depths(fname):

    ResidueDepths = OrderedDict()

    with open(fname,'r') as inf:
        lines = inf.readlines()
        for l in lines:

            if l[0]=='#':
                continue

            if len(l)<=1:
                continue

            toks = l.strip().split()

            try:
                chain , residue = toks[1].split(':')
            except:
                continue
            atomtype = toks[2]

            ResidueDepths[(chain,residue,atomtype)]=[float(toks[3]), float(toks[4])]

    return ResidueDepths