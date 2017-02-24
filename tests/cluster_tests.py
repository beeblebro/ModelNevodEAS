from numpy import array
import unittest

from cluster import Cluster

class TestCluster(unittest.TestCase):

    def setUp(self):
        """Создание кластера"""
        self.cluster = Cluster()

        self.center = array(0, 0, 0)
