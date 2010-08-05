# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2005  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROJECT = vnfb
PACKAGE = qeutils/qeparser/inputs


#--------------------------------------------------------------------------
#

EXPORT_PYTHON_MODULES = \
        __init__.py \
        inputbands.py \
        inputcp.py \
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

BUILD_DIRS = \

OTHER_DIRS = \

RECURSE_DIRS = $(BUILD_DIRS) $(OTHER_DIRS)

all: export

tidy::
	BLD_ACTION="tidy" $(MM) recurse

clean::
	BLD_ACTION="clean" $(MM) recurse

distclean::
	BLD_ACTION="distclean" $(MM) recurse

export:: export-package-python-modules
	BLD_ACTION="export" $(MM) recurse


# version
# $Id$

# End of file