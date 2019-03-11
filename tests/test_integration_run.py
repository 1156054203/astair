from __future__ import print_function

import pdb
import csv
import sys
import unittest
import subprocess
from os import path

import scripts.run as run

current = path.abspath(path.dirname(__file__))


class AsTaiRunTest(unittest.TestCase):
    """"""

    def test_run(self):
        """"""
        self.assertEqual(run.cli.help, 'asTair (tools for processing cytosine modification sequencing data)')
        self.assertEqual(run.cli.epilog, '\n__________________________________About__________________________________\nasTair was written by Gergana V. Velikova and Benjamin Schuster-Boeckler.\nThis code is made available under the GNU General Public License, see \nLICENSE.txt for more details.\n\n                                                         Version: 3.0.0\n')
        
        
if __name__ == '__main__':
    unittest.main()
