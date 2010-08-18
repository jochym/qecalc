#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Alex Dementsov
#                      California Institute of Technology
#                        (C) 2009  All Rights Reserved
#
# {LicenseText}
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

PROJECT = qecalc
PACKAGE = qetask/qeparser

BUILD_DIRS = \
	inputs \
        outputs \
	qestructureparser \


OTHER_DIRS = \

RECURSE_DIRS = $(BUILD_DIRS) $(OTHER_DIRS)

#--------------------------------------------------------------------------
#

all: export-package-python-modules  #export
	BLD_ACTION="all" $(MM) recurse

tidy::
	BLD_ACTION="tidy" $(MM) recurse

clean::
	BLD_ACTION="clean" $(MM) recurse

distclean::
	BLD_ACTION="distclean" $(MM) recurse

#--------------------------------------------------------------------------
# export

EXPORT_PYTHON_MODULES = \
    block.py \
    card.py \
    d3input.py \
    d3qpoints.py \
    __init__.py \
    matdyninput.py \
    matdynqpoints.py \
    namelist.py \
    orderedDict.py \
    phinput.py \
    phqpoints.py \
    pwinput.py \
    cpinput.py \
    pwkpoints.py \
    qeinput.py \
    qesinput.py \
    qelattice.py \
    qeoutput.py \
    qeparser.py \
    qestructure.py \
    qeatom.py \



__date__ = "$Jan 17, 2010 4:13:02 PM$"
