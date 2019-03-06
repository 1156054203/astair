#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import print_function

def astair():
    print(
    """Program: asTair (tools for processing cytosine modification sequencing data)
    Version: 3.0

    Usage: Command [command specific options/arguments]

    Commands:
  
    astair_align         Align reads
    
    astair_call          Call methylation
    
    astair_simulate      Simulate TAPS or WGBS data
    
    astair_mbias         Visualise modification bias
    
    astair_phred         Visualise Phred scores
    

    __________________________________About__________________________________
    
    asTair was written by Gergana V. Velikova and Benjamin Schuster-Boeckler.
    This code is made available under the GNU General Public License, see 
    LICENSE.txt for more details.""")
    
if __name__ == "__main__":
   astair()
