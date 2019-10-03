import numpy as numpy
import os

class FEFF4Molecular(coords_file):

	def __init__(self):
		self.path = os.getcwd()

		print('current running path: ',self.path)
		params = self.read_parameters(self.path+'input.dat')
		self.office = params['WORKING_STATION']
		self.sec_atom = params['SEC_ATOM']
		self.atom_edge = params['ABSORPTION_ATOM']
		self.KL23_edge = params['ABSORPTION_EDGE']
		self.lattice_group = params['LATTICE_GROUP']
		self.cd_num = int(params['Cd_NUM'])
		self.side = params['SIDE']
		if self.sec_atom == 'Se':
		    self.sec_atomic_num = 34
		elif self.sec_atom == 'S':
		    self.sec_atomic_num = 16

		if self.atom_edge == 'Cd':
		    self.atom_atomic_num = 48
		elif self.atom_edge == 'Se':
		    self.atom_atomic_num = 34
		elif self.atom_edge == 'S':
		    self.atom_atomic_num = 16

		return

	def read_parameters(self,filename):
	    config = configparser.ConfigParser()
	    config.read(filename)
	    parameters = config['Parameters']
	    return parameters

	def coords_reform(self):

		with open(self.path+coords_file,'r') as f:
			original_data = []
		    for line in f.readlines():
		        original_data.append(line.split())
		for item in original_data:
		    element = []
		    if item == []:
		        continue
		    if len(item) == 1:
		        continue
		    if item[0] != 'Cd' and item[0] != 'O' and item[0] != 'C' and item[0] != '{:s}'.format(self.sec_atom):
		        continue

		    if len(item) > 7:
		        element += item[5:8]
		    else:
		        element += item[2:5]

		    if item[0] == 'Cd':
		        element += ['1']
		    elif item[0] == '{:s}'.format(sec_atom):
		        element += ['2']
		    elif item[0] == 'O':
		        element += ['3']
		    elif item[0] == 'C':
		        element += ['4']
		    element.append(item[0])
		    coordinate_list.append(element)

		g = open(self.path+coords_file, 'w')

		for i in coordinate_list:
		    g.write("%s  " % i[0])
		    g.write("%s  " % i[1])
		    g.write("%s  " % i[2])
		    g.write("%s  " % i[3])
		    g.write("%s" % i[4])
		    g.write("\n")
		g.close()

		print('Coordinates withdrawed')

	def check_atom(self):
	    datafile = open(self.path+coords_file).read()
	    if 'O' in datafile:
	        O_potential = True
	    else: O_potential = False
	    if 'C ' in datafile:
	        C_potential = True
	    else: C_potential = False
	    return O_potential,C_potential

	def genfeffinp(self):

		O_potential,C_potential = self.check_atom()
		print(' Existence of O and C : ', O_potential,C_potential)

		absorb_atom_num = 0
		for filenumber, i in enumerate(coordinate_data):
	    if i[4] == '{:s}'.format(atom_edge):
	        absorb_atom = list(i)
	        absorb_atom[3] = 0 #change absorbing atom index into 0
	        #       print(i[3])
	        #    print(absorb_atom)
	        #    filenumber =+1
	        absorb_atom_num = absorb_atom_num + 1
	        print(absorb_atom_num)
	        g = open('{:d}feff.inp'.format(absorb_atom_num), 'w')
	        g.write('''TITLE Cd{:s}_nano
                
CONTROL	{:s}
PRINT	{:s}
EXCHANGE 0 -1 0.0 -1
SCF 4.5 0 100 0.2 1  
EDGE {:s}	{:s}
COREHOLE {:s}
LDOS	{:s}
RPATH {:s}
EXAFS	{:s}
DEBYE {:s}



POTENTIALS
*	ipot	z	label	lmax1	lmax2
0 {:d} {:s}
1 48 Cd
2 {:d} {:s} \n'''.format(self.sec_atom,self.CONTROL, self.PRINT, self.KL23_edge,self.SO2,
	self.COREHOLE, self.LDOS, self.RPATH, self.EXAFS,self.DEBYE, self.atom_atomic_num, 
	self.atom_edge, self.sec_atomic_num, self.sec_atom))
	        if O_potential:
	            g.write('''3 8 O \n''')
	        if C_potential:
	            g.write('''4 6 C \n''')

	        g.write('ATOMS \n')

	        for ii in absorb_atom:
	            g.write('%s' % ii)
	            g.write('  ')
	        g.write('\n')

	        for j in coordinate_data:
	            if j[0:3] == i[0:3]:
	                continue
	            else:
	                for jj in j:
	                    g.write('%s' % jj)
	                    g.write('  ')
	                g.write('\n')

	        g.write('END')
	        g.close()

		print('FEFF file generated')




