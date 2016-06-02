from os.path import join

import xlrd

from src.utils import PROJECT_DIR


class InsulinMassSpecAndModelDataReader:

	_default_datafile_for_nmeth_2703 = join(PROJECT_DIR,'data','insulin','mass_spec','nmeth.2703-S5.xlsx')
	_default_datafile_for_model_time_series = join(PROJECT_DIR,'data','insulin','time_series_BIOMD0000000223.csv')

	_relevant_proteins_nmeth_2703 = {'IR': 'P08069', 'RasGAP': 'Q96PV0',  'RAF': 'P04049',
	                     'MEK2': 'P36507', \
	                     'ERK': 'Q99759', 'IRS 2': 'Q9Y4H2',  'mTOR': 'P42345'} #'AKT1-substrate': 'Q96B36', 'SHP2 [ptpn13]': 'Q9UQC2',

	def __init__(self):
		pass

	def loadnmeth_2703(self,file_=_default_datafile_for_nmeth_2703):
		'''
		Loads the mass spec data about the insulin pathway from the paper
			Quantifying protein interaction dynamics by SWATH mass spectrometry: application to the 14-3-3 system, Ben et al. Nature Methods 2013
		:param file_:where to find the xlsx data file, the default location should work fine unless you did something weird
		:return: a tuple of
			1) dictionary with key=protein names in uniprot and value=array of log2 fold change compared to 0 minutes at t=[-60,1,10,30,100]mins respectively
			2) time indexes of measurements in minutes ["-60", "1", "10", "30", "100"]
			3) a map of relevant proteins with key=human friendly name, value=unpirot id
		'''
		sheet = xlrd.open_workbook(file_).sheet_by_index(0)
		res = {}
		for row_idx in range(6,500):
			name = sheet.cell_value(row_idx,0)
			res[name] = []
			for i,col_idx in enumerate([1,3,5,7,9]):
				try:
					k = sheet.cell_value(row_idx, col_idx)
					res[name].append(self._transform(k))
				except ValueError:
					raise ValueError('Values row {0}, cell {1} is not a number.'.format(row_idx,col_idx))

		return res,["-60", "1", "10", "30", "100"],self._relevant_proteins_nmeth_2703

	def _transform(self,x):
		return float(x)

	def load_model_data(self,file_=_default_datafile_for_model_time_series):
		line_num = 2
		t,y = [],[]
		with open(file_) as pin:
			pin.readline()
			for line in pin:
				s_line = line.split(',')
				try:
					t.append(self._transform(s_line[0]))
					y.append([self._transform(e) for e in s_line[1:]])
				except ValueError:
					raise ValueError('Values in row {0} are not a number.'.format(line_num))
				line_num += 1
		return t,y
