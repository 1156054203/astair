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
        
if __name__ == '__main__':
    unittest.main()
