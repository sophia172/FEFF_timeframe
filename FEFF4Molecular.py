import configparser

import numpy as np
import os
from scipy.interpolate import interp1d
import glob
import configparser


class FEFF4Molecular():
	def __init__(self, coords_file):
		self.coords_file = coords_file
		self.path = os.getcwd()
		self.out_dir = self.path
		print('current running path: ',self.path)
		params = self.read_parameters(self.path+'/input.dat')
		self.office = params['WORKING_STATION']
		self.sec_atom = params['SEC_ATOM']
		self.atom_edge = params['ABSORPTION_ATOM']
		self.KL23_edge = params['ABSORPTION_EDGE']
		self.lattice_group = params['LATTICE_GROUP']
		self.cd_num = int(params['Cd_NUM'])
		self.side = params['SIDE']
		##feff.inp file parameter
		self.CONTROL = params['CONTROL']
		self.PRINT = params['PRINT']
		self.COREHOLE = params['COREHOLE']
		self.SO2 = params['SO2']
		self.SCF = params['SCF']
		self.RPATH = params['RPATH']
		self.EXAFS = params['EXAFS']
		self.AFOLP = params['AFOLP']
		self.EXCHANGE = params['EXCHANGE']
		self.LDOS = params['LDOS']
		self.DEBYE = params['DEBYE']

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
		coordinate_list = []
		with open(self.path+'/'+self.coords_file,'r') as f:
			original_data = []
			for line in f.readlines():
				original_data.append(line.split())
		for item in original_data:
			element = []
			if item == []:
				continue
			elif len(item) == 1:
				continue
			elif item[0] != 'Cd' and item[0] != 'O' and item[0] != 'C' and item[0] != '{:s}'.format(self.sec_atom):
				continue

			elif len(item) > 7:
				element += item[5:8]
			else:
				element += item[2:5]

			if item[0] == 'Cd':
				element += ['1']
			elif item[0] == '{:s}'.format(self.sec_atom):
				element += ['2']
			elif item[0] == 'O':
				element += ['3']
			elif item[0] == 'C':
				element += ['4']
			element.append(item[0])
			coordinate_list.append(element)

		g = open(self.path+'/'+self.coords_file, 'w')

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
		datafile = open(self.path+'/'+self.coords_file).read()
		if 'O' in datafile:
			O_potential = True
		else: O_potential = False
		if 'C ' in datafile:
			C_potential = True
		else: C_potential = False
		return O_potential,C_potential

	def genfeffinp(self):

		with open(self.path+'/'+self.coords_file, 'r') as f:
			coordinate_data = []
			for line in f.readlines():
				coordinate_data.append(line.split())

		O_potential,C_potential = self.check_atom()
		print(' Existence of O and C : ', O_potential,C_potential)

		absorb_atom_num = 0
		for filenumber, i in enumerate(coordinate_data):
			if i[4] == '{:s}'.format(self.atom_edge):
				absorb_atom = list(i)
				absorb_atom[3] = 0 #change absorbing atom index into 0
				#       print(i[3])
				#    print(absorb_atom)
				#    filenumber =+1
				absorb_atom_num = absorb_atom_num + 1
				print(absorb_atom_num)
				g = open('{:d}feff.inp'.format(absorb_atom_num), 'w')
				g.write('''TITLE Cd{:s}_nano\n
EDGE {:s}	{:s}	
CONTROL	{:s}
PRINT	{:s}
EXCHANGE 0 -1 0.0 -1
SCF 4.5 0 100 0.2 1 
COREHOLE {:s}
LDOS	{:s}
RPATH {:s}
EXAFS	{:s}
DEBYE {:s}
	
POTENTIALS
*	ipot	z	label	lmax1	lmax2
0 {:d} {:s}
1 48 Cd
2 {:d} {:s} \n'''.format(self.sec_atom,self.KL23_edge,self.SO2,self.CONTROL, self.PRINT,
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
		return

	def newfolder_runFEFF(self):
		os.system('mkdir {:s}'.format(self.coords_file[:-4]))
		os.system('mv ' + self.path + '/' + self.coords_file + ' '+self.path + '/' + self.coords_file[:-4] + '/')
		os.system('cp ' + self.path + '/autorun.sh' + ' ' + self.path + '/' + self.coords_file[:-4] + '/')
		os.system('mv ' + self.path + '/*feff.inp' + ' ' + self.path + '/' + self.coords_file[:-4] + '/')
		os.chdir(self.path+'/'+self.coords_file[:-4]+'/')
		print('check current directory in subfolder: ',self.path+'/'+ self.coords_file[:-4]+'/')
		os.system('bash ./autorun.sh')
		return

	def average(self,file_list, data_cols_to_average):
		'''
	    #reads files and calcs averages
	    #note - big question is: do we want to average on each column, or read energy, ff, interpolate and then average ion the interpolated result??
	    #for now, just copying the 1st file energy column, then averaging experimental points on each data 'row'
	    #input parameters: a list containing the file names to be processed and a list of integers with the columns to average out
	    #output : data matrix (array) with energy on 1st column, e (energy - Fermi) on 2nd column, k on 3rd column then averaged data columns
	    '''
		print('Averaging:')
		print(file_list)
		print('data_cols_to_average:', data_cols_to_average)
		print('#########data in each file will be interpolate with same amount of sample points...#######')
		run_once = 0
		nrows_min = np.inf
		nrows_max = -np.inf
		for ifile in file_list:
			data_in = np.genfromtxt(ifile)
			irows = np.shape(data_in)[0]  # shape gives row number and column number
			print('row number recorded yeehha  ',ifile)
			nrows_min = min(nrows_min, irows)
			nrows_max = max(nrows_max, irows)
		if (nrows_min == nrows_max):
			print('ok, all data have same size,', nrows_min)
		else:
			print('Warning: data have different sizes! min:', nrows_min, '    max:', nrows_max)

		print('record the energy data points from last file ......')
		raw_e_value = data_in[:,:data_cols_to_average[0]]
		energy_avg = np.zeros((nrows_min,data_cols_to_average[0]))
		#    for icol in data_cols_to_average[0]:      # This column is for energy spaced out in log. K is not in part
	    #        energy_avg[:,icol] = np.logspace(raw_e_value[0,icol], raw_e_value[-1,icol], num=nrows_min)
		#        #use third column k to calculate energy and e
		energy_avg[:,0] = np.linspace(raw_e_value[0,0], raw_e_value[-1,0], num=nrows_min)
		energy_avg[:,1] = np.linspace(raw_e_value[0,1], raw_e_value[-1,1], num=nrows_min)
		energy_avg[:,2] = np.linspace(raw_e_value[0,2], raw_e_value[-1,2], num=nrows_min)
		print('calculating averages on min no. of rows with interpolate data sets')
		# now do the averages
		ncols = np.size(data_cols_to_average) + data_cols_to_average[0]  # keep all columns before first averaged column
		data_out = np.zeros((nrows_min, ncols))
		data_out[:, 0:data_cols_to_average[0]] = energy_avg[:,:] # The array shape should be same here
		run_once = 0
		# first sum data
		for ifile in file_list:
			data_in = np.genfromtxt(ifile)
			#data_in[:, 5] = 1 / np.exp(data_in[:, 5])
			for icol in range(data_cols_to_average[0]):
				raw_spec = interp1d(data_in[:,icol],data_in[:,data_cols_to_average[icol]],fill_value="extrapolate") # interp1d(energy,mu)/(e,mu0)/(k,chi). However, This would be a problem if there are more y columns than x columns.
				data_in[:nrows_min,data_cols_to_average[icol]] = raw_spec(energy_avg[:,icol])
				print('interpolation successful in column  ',icol)
			for icol in data_cols_to_average:
				# data_cols_to_average need to be in sequence
				##print(icol, data_out[:, icol])
				data_out[:, icol] = data_out[:, icol] + data_in[:nrows_min, icol]

		# now divide by n
		nfiles = np.size(file_list)
		for icol in data_cols_to_average:
			data_out[:, icol] = data_out[:, icol] / nfiles
			print('another one column done')

		return data_out

	def average_allcols(self,group):
		print('inside average_allcols function ')
		print('output xmu filename: ',group)
		data_in = np.genfromtxt(group[0])
		print('numpy function successful')
		ncols = np.shape(data_in)[1]
		print('Average_allcols - ncols:', ncols)
		if (ncols == 6):
			avg = self.average(group, np.arange(3,6))
			print('average calculation successful')
		else:
			print('column number not right,change code or change file')
		if (self.out_dir == None):
			file_out = os.path.basename('{:s}simulation_'+self.coords_file[:4]+'.dat'.format(self.side))  # save in present directory
		else:
			file_out = self.out_dir + '/' + os.path.basename('{:s}simulation_'+self.coords_file[:4]+'.dat'.format(self.side))
		print('file_out:', file_out)
		np.savetxt(file_out, avg,fmt='%g')
		return

	def excute(self):
		self.coords_reform()
		self.genfeffinp()
		self.newfolder_runFEFF()
		file_list = glob.glob(self.path+'/'+self.coords_file[:-4]+'/*xmu.dat')
		self.average_allcols(file_list)
		return

if __name__ == '__main__':
	
	file = FEFF4Molecular('CdS37_Coordinates.txt')
	file.excute()

