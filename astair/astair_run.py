#!/usr/bin/env python
#-*- coding: utf-8 -*-

import click

@click.group()
def cli():
    """
    asTair (tools for processing cytosine modification sequencing data)
   
    Version: 3.0
    __________________________________About__________________________________
    
    asTair was written by Gergana V. Velikova and Benjamin Schuster-Boeckler.
    This code is made available under the GNU General Public License, see 
    LICENSE.txt for more details."""
pass

import aligner
import caller 
import phred
import simulator
import mbias


cli.add_command(aligner.align)
cli.add_command(caller.call)
cli.add_command(simulator.simulate)
cli.add_command(phred.phred)
cli.add_command(mbias.mbias)

if __name__ == '__main__':
	cli()
