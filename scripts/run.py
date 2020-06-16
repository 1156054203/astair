#!/usr/bin/env python
#-*- coding: utf-8 -*-

import click
import logging

import astair
import astair.mbias as mbias
import astair.phred as phred
import astair.caller as caller
import astair.finder as finder
import astair.filter as filter 
import astair.aligner as aligner
import astair.summary as summary
import astair.idbiaser as idbiaser
import astair.simulator as simulator


# TODO make this config properly configurable using command line options
# For example, we could use a global option -v to change the log level to
# DEBUG and produce verbose output for all commands
logging.basicConfig(level=logging.WARNING)
logs = logging.getLogger(__name__)

@click.group()
def cli():
    """
    asTair (tools for processing cytosine modification sequencing data)
    """
    pass

cli.epilog = """
__________________________________About__________________________________
asTair was written by Gergana V. Velikova and Benjamin Schuster-Boeckler.
This code is made available under the GNU General Public License, see 
LICENSE.txt for more details.

                                                         Version: __version__
"""
cli.epilog = cli.epilog.replace('__version__', astair.__version__)


cli.add_command(caller.call)
cli.add_command(phred.phred)
cli.add_command(mbias.mbias)
cli.add_command(finder.find)
cli.add_command(aligner.align)
cli.add_command(simulator.simulate)
cli.add_command(filter.filter)
cli.add_command(summary.summarise)
cli.add_command(idbiaser.idbias)


if __name__ == '__main__':
    cli()
