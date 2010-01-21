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
PACKAGE = qetask/qeparser/inputs

BUILD_DIRS = \


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
    __init__.py \
    inputbands.py \
    inputcppp.py \
    inputd3.py \
    inputdos.py \
    inputdynmat.py \
    inputgipaw.py \
    inputinitial_state.py \
    inputld1.py \
    inputmatdyn.py \
    inputph.py \
    inputpp.py \
    inputprojwfc.py \
    inputpwcond.py \
    inputpw.py \
    inputq2r.py \



__date__ = "$Jan 17, 2010 4:13:02 PM$"
