from numpy import array
import unittest

from core.utils import get_distance


class TestUtils(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_ger_distance(self):
        res = get_distance(st_coord=[0, 0, 0], n=[1, 0, 1], x0=1, y0=0)
        self.assertEqual(res, 1.0)


if __name__ == '__main__':
    unittest.main()
