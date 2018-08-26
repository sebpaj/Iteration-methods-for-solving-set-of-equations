import unittest
from main import *

class Tests(unittest.TestCase):
    def test_jacobi_x1(self):
        self.assertAlmostEqual(jacobi_method([[-2, 1, 0, 0, 0],
               [1, -2, 1, 0, 0],
               [0, 1, -2, 1, 0],
               [0, 0, 1, -2, -1]])[0],
               0.200000, delta=0.000001)
    def test_jacobi_x2(self):
        self.assertAlmostEqual(jacobi_method([[-2, 1, 0, 0, 0],
               [1, -2, 1, 0, 0],
               [0, 1, -2, 1, 0],
               [0, 0, 1, -2, -1]])[1],
               0.400000, delta=0.000001)
    def test_jacobi_x3(self):
        self.assertAlmostEqual(jacobi_method([[-2, 1, 0, 0, 0],
               [1, -2, 1, 0, 0],
               [0, 1, -2, 1, 0],
               [0, 0, 1, -2, -1]])[2],
               0.600000, delta=0.000001)
    def test_jacobi_x4(self):
        self.assertAlmostEqual(jacobi_method([[-2, 1, 0, 0, 0],
               [1, -2, 1, 0, 0],
               [0, 1, -2, 1, 0],
               [0, 0, 1, -2, -1]])[3],
               0.800000, delta=0.000001)
    def test_gauss_siedl_x1(self):
        self.assertAlmostEqual(gauss_siedl([[-2, 1, 0, 0, 0],
               [1, -2, 1, 0, 0],
               [0, 1, -2, 1, 0],
               [0, 0, 1, -2, -1]])[0],
               0.200000, delta=0.000001)
    def test_gauss_siedl_x2(self):
        self.assertAlmostEqual(gauss_siedl([[-2, 1, 0, 0, 0],
               [1, -2, 1, 0, 0],
               [0, 1, -2, 1, 0],
               [0, 0, 1, -2, -1]])[1],
               0.400000, delta=0.000001)
    def test_gauss_siedl_x3(self):
        self.assertAlmostEqual(gauss_siedl([[-2, 1, 0, 0, 0],
               [1, -2, 1, 0, 0],
               [0, 1, -2, 1, 0],
               [0, 0, 1, -2, -1]])[2],
               0.600000, delta=0.000001)
    def test_gauss_siedl_x4(self):
        self.assertAlmostEqual(gauss_siedl([[-2, 1, 0, 0, 0],
               [1, -2, 1, 0, 0],
               [0, 1, -2, 1, 0],
               [0, 0, 1, -2, -1]])[3],
               0.800000, delta=0.000001)
    def test_sor_x1(self):
        self.assertAlmostEqual(sor([[-2, 1, 0, 0, 0],
               [1, -2, 1, 0, 0],
               [0, 1, -2, 1, 0],
               [0, 0, 1, -2, -1]])[0],
               0.200000, delta=0.000001)
    def test_sor_x2(self):
        self.assertAlmostEqual(sor([[-2, 1, 0, 0, 0],
               [1, -2, 1, 0, 0],
               [0, 1, -2, 1, 0],
               [0, 0, 1, -2, -1]])[1],
               0.400000, delta=0.000001)
    def test_sor_x3(self):
        self.assertAlmostEqual(sor([[-2, 1, 0, 0, 0],
               [1, -2, 1, 0, 0],
               [0, 1, -2, 1, 0],
               [0, 0, 1, -2, -1]])[2],
               0.600000, delta=0.000001)
    def test_sor_x4(self):
        self.assertAlmostEqual(sor([[-2, 1, 0, 0, 0],
               [1, -2, 1, 0, 0],
               [0, 1, -2, 1, 0],
               [0, 0, 1, -2, -1]])[3],
               0.800000, delta=0.000001)