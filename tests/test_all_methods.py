from unittest import TestCase

from utils import getMockInstanceSingleTimeSeriesView3Nodes
from src.methods.methods import methods

class TestAllMethods(TestCase):

	def setUp(self):
		self.instance = getMockInstanceSingleTimeSeriesView3Nodes()

	def test_OutputFormat(self):
		for method in methods:
			print method
			self.singleMethodOutputFormatTest(methods[method])

	def singleMethodOutputFormatTest(self,method):
		res = method(self.instance)
		self.assertEqual(res.shape,(self.instance.n_nodes,self.instance.n_nodes))
		self.assertTrue('float' in res.dtype.name)


