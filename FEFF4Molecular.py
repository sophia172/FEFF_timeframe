import configparser

import numpy as np
import os
from scipy.interpolate import interp1d
import glob
import configparser
import multiprocess as mp
import threading
import time
import fileinput


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
		self.feff_path = ''

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
		file_end_line = 0
		with open(self.path+'/'+self.coords_file,'r') as f:
			original_data = []
			for line in f.readlines():
				original_data.append(line.split())

		self.file_num = 0
		for n,item in enumerate(original_data):

			element = []
			print(item)
			if item == []:
				continue
			elif item[0].isspace():
				continue
			elif len(item) == 1:
				self.atom_num = int(item[0])
				file_end_line += self.atom_num + 2
				self.file_num += 1
				continue
			elif item[0].startswith('i'):
				try:
					self.file_num = int(item[2].replace(',',''))
				except:
					print('no i in this line')
				continue


			elif len(item) > 7:
				element += item[5:8]
			else:
				element += item[1:4]

			if item[0].startswith('Cd'):
				element += ['1']
				element.append('Cd')
				print('element attached to coordinate_list', element)
			elif item[0].startswith('{:s}'.format(self.sec_atom)):
				element += ['2']
				element.append(self.sec_atom)
			elif item[0].startswith('O'):
				element += ['3']
				element.append('O')
			elif item[0].startswith('C'):
				element += ['4']
				element.append('C')
			coordinate_list.append(element)
			print('element added to coordination list', element)
			if n == file_end_line -1:
				print('n == line_x:', file_end_line ,'in file: ',self.file_num)
				g = open(self.path + '/' + 'Cd{:s}'.format(self.sec_atom) + str(self.file_num) + '.dat', 'w')

				for i in coordinate_list:
					g.write("%s  " % i[0])
					g.write("%s  " % i[1])
					g.write("%s  " % i[2])
					g.write("%s  " % i[3])
					g.write("%s" % i[4])
					g.write("\n")
				g.close()
				coordinate_list = []


		print('Coordinates withdrawed')

	def check_atom(self):
		datafile = open(self.path+'/' + self.coords_file[:-4] + '/' +self.coords_file).read()
		if 'O' in datafile:
			O_potential = True
		else: O_potential = False
		if 'C ' in datafile:
			C_potential = True
		else: C_potential = False
		return O_potential,C_potential

	def genfeffinp(self):
		os.system('mkdir {:s}'.format(self.coords_file[:-4]))
		os.system('mv ' + self.path + '/' + self.coords_file + ' ' + self.path + '/' + self.coords_file[:-4] + '/')
		with open(self.path + '/' + self.coords_file[:-4] + '/'+ self.coords_file, 'r') as f:
			print('open coordination file ',self.coords_file)
			self.coordinate_data = []
			for line in f.readlines():
				self.coordinate_data.append(line.split())
		O_potential,C_potential = self.check_atom()
		print(' Existence of O and C : ', O_potential,C_potential)

		absorb_atom_num = -1
		for filenumber, i in enumerate(self.coordinate_data):
			if i[4] == '{:s}'.format(self.atom_edge):
				absorb_atom = list(i)
				absorb_atom[3] = 0 #change absorbing atom index into 0
				#       print(i[3])
				#    print(absorb_atom)
				#    filenumber =+1
				absorb_atom_num += 1
				g = open(self.path + '/' + self.coords_file[:-4] + '/' + '{:d}feff.inp'.format(absorb_atom_num), 'w')
				g.write('''TITLE Cd{:s}_nano\n
EDGE {:s}	{:s}	
CONTROL	{:s}
PRINT	{:s}
EXCHANGE  {:s}
COREHOLE {:s}
SCF {:s}
RPATH {:s}
EXAFS	{:s}
DEBYE {:s}
	
POTENTIALS
*	ipot	z	label	lmax1	lmax2
0 {:d} {:s}
1 48 Cd 
2 {:d} {:s} \n'''.format(self.sec_atom,self.KL23_edge,self.SO2,self.CONTROL, self.PRINT, self.EXCHANGE,
					 self.COREHOLE, self.SCF, self.RPATH, self.EXAFS,self.DEBYE, self.atom_atomic_num,
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

				for j in self.coordinate_data:
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

	def feff_folder(self):
		self.cd_num = len(glob.glob(self.path + '/' + self.coords_file[:-4] + '/*feff.inp'))
		for i in range(self.cd_num):
			os.system('mkdir ' + self.path + '/' + self.coords_file[:-4] + '/' + '{:d}feff'.format(i))
			os.system('mv ' + self.path + '/' + self.coords_file[:-4] + '/' + '{:d}feff.inp'.format(i) + ' '
					  + self.path + '/' + self.coords_file[:-4] + '/' + '{:d}feff/feff.inp'.format(i))
		return

	def runFEFF(self):
		print('check current directory in subfolder: ',self.feff_path)
		os.chdir(self.feff_path)
		os.system('~/Packages/JFEFF/feff feff.inp')
		return

	def xmu_exist(self):
		_xmu_exist = os.path.exists(self.feff_path + '/xmu.dat')
		return _xmu_exist
	
	def change_listdat(self):
		print('change list.dat file in : ',self.feff_path)
		os.chdir(self.feff_path)
		with open(self.feff_path + '/list.dat') as file:
			lines = []
			for line in file.readlines():
				line = line.split()
				try:
					r = float(line[-1])
				except:
					lines.append(line)
					continue
				if r < 2.45 :
					line[1] = '   0.01'  # bulk CdO not known
				elif r < 3:
					line[1] = '  -0.0057' # -0.0057 bulk CdS
				elif r < 5 :
					line[1] = '   0.0045'  # 0.0045 bulk CdS
				elif r < 6 :
					line[1] = '  0.006' # a lot so can be ignored
				elif r < 7 :
					line[1] = '   0.008' # a lot so can be ignored
				else:
					line[1] = '0.01'
				lines.append(line)
		with open(self.feff_path + '/list.dat','w') as file:
			for line in lines:
				for item in line:
					file.write('%s  '%item)
				file.write('\n')
		return

	def change_feffinp(self):
		print('change feff input file: ',self.feff_path)
		os.chdir(self.feff_path)
		for line in fileinput.FileInput(self.feff_path + '/feff.inp',inplace=1,backup='.bak'):
			line = line.replace('1 1 1 1 1 1', '1 1 1 1 0 1')
			print(line)

		return

	def del_xmu(self):
		print('deleting ', self.feff_path,'/xmu.dat')
		os.system('rm ' + self.feff_path + '/xmu.dat')
		return

	def delete_listdat(self):
		print('change feff input file: ',self.feff_path)
		os.chdir(self.feff_path)
		os.system('rm list.dat')
		return


	def group_xmus(self):
		for i in range(self.cd_num):
			if not os.path.exists(self.coords_file[:-4]+'/{:d}feff/list.dat'.format(i))\
					and os.path.exists(self.coords_file[:-4]+'/{:d}feff/xmu.dat'.format(i)):
				os.system('mv ' + self.path + '/' + self.coords_file[:-4] + '/{:d}feff/xmu.dat'.format(i) + ' '
						  + self.path + '/' + self.coords_file[:-4] + '/{:d}xmu.dat'.format(i))
				os.system(('rm -r ' + self.path + '/' + self.coords_file[:-4] + '/{:d}feff'.format(i)))
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
			file_out = os.path.basename('EXAFSsimulation_'+self.coords_file)  # save in present directory
		else:
			file_out = self.out_dir + '/' + os.path.basename('EXAFSsimulation_'+self.coords_file)
		print('file_out:', file_out)
		np.savetxt(file_out, avg,fmt='%g')
		return

	def feff_running_workflow(self):
		self.runFEFF()
		while not self.xmu_exist():
			time.sleep(60)
		print('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<time to change list.dat')
		self.change_listdat()
		print('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<finished changing list.dat, time to delete xmu.dat')
		self.del_xmu()
		print('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<finished delete xmu.dat, time to change feff.inp')
		self.change_feffinp()
		self.runFEFF()
		while not self.xmu_exist():
			time.sleep(60)
		self.delete_listdat()

def assign_task(file,feff_path):
	
	print('start calculating coordinates frame  >>>>>>>>>>>>>>>>>>>',i )
	file.feff_path = feff_path
	file.feff_running_workflow()

	return



def loop_del_folders():
	while True:
		for i in range(1,int(file_num)+1):

			
			_folder_exist = os.path.exists(files[str(i)].path + '/' + files[str(i)].coords_file[:-4])
			if _folder_exist:
				print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ', files[str(i)].coords_file[:-4], 'folder exists')
				# print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>CdS',i,' cleaning')

				files[str(i)].group_xmus()
				_feff_folder_exist = any([os.path.exists(files[str(i)].path + '/' + files[str(i)].coords_file[:-4]
												+ '/'+str(j) + 'feff') for j in range(files[str(i)].cd_num)])
				if not _feff_folder_exist:
					print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ', files[str(i)].coords_file[:-4], '/',str(j), 'feff folder not exist')
					files[str(i)].group_xmus()
					file_list = glob.glob(files[str(i)].path+'/'+files[str(i)].coords_file[:-4]+'/*xmu.dat')
					files[str(i)].average_allcols(file_list)
					os.system('rm -r '+ files[str(i)].path + '/' + files[str(i)].coords_file[:-4])

		_last_feff_folder_exist = os.path.exists(files[str(int(file_num))].path + '/' + files[str(int(file_num))].coords_file[:-4]
											+ '/'+str(files[str(int(file_num))].cd_num-1) + 'feff')
		print(files[str(int(file_num))].path + '/' + files[str(int(file_num))].coords_file[:-4]
											+ '/'+str(files[str(int(file_num))].cd_num-1) + 'feff', _last_feff_folder_exist)
		if _last_feff_folder_exist:

			time.sleep(100)
			continue
		else:
			print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> last feff folder not exist')
			for i in range(1,int(file_num)+1):
				files[str(i)].group_xmus()
				_folder_exist = os.path.exists(files[str(i)].path + '/' + files[str(i)].coords_file[:-4])
				if _folder_exist:
					file_list = glob.glob(files[str(i)].path+'/'+files[str(i)].coords_file[:-4]+'/*xmu.dat')
					print(file_list)
					files[str(i)].average_allcols(file_list)
			os.system('rm -r CdS*')
			break


if __name__ == '__main__':
	xyz_file = FEFF4Molecular('CdS_betaSnV1_add_oxygen.xyz')
	xyz_file.coords_reform()
	file_num = xyz_file.file_num
	files={}
	for i in range(1,int(file_num)+1):
		files[str(i)]=FEFF4Molecular('CdS'+str(i)+'.dat')
		files[str(i)].genfeffinp()
		files[str(i)].feff_folder()

	t = threading.Thread(target=loop_del_folders)
	t.start()
	#
	pool = mp.Pool(mp.cpu_count())
	print('CPU number : ',mp.cpu_count())
	for i in range(1,int(file_num)+1):
		for j in range(files[str(i)].cd_num):
			feff_path = files[str(i)].path + '/' + files[str(i)].coords_file[:-4] + '/{:d}feff'.format(j)
			# assign_task(files[str(i)],feff_path)
			pool.apply_async(assign_task,(files[str(i)],feff_path,))
	pool.close()
	pool.join()

