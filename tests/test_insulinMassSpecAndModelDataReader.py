from unittest import TestCase

from src.readers.mass_spec_data_reader import InsulinMassSpecAndModelDataReader

class TestInsulinMassSpecAndModelDataReader(TestCase):


	def test_loadnmeth_2703(self):
		self.reader = InsulinMassSpecAndModelDataReader()
		res,times,relevant_proteins_nmeth_2703 = self.reader.loadnmeth_2703()
		self.assertTrue(res)
		self.assertTrue(times)
		self.assertTrue(relevant_proteins_nmeth_2703)

	def test_load_model_data(self):
		self.reader = InsulinMassSpecAndModelDataReader()
		[t,y] = self.reader.load_model_data()
		self.assertTrue(t)
		self.assertTrue(y)
