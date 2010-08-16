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
PACKAGE = qetask

BUILD_DIRS = \
    qeparser \


OTHER_DIRS = \

RECURSE_DIRS = $(BUILD_DIRS) $(OTHER_DIRS)

#--------------------------------------------------------------------------
#

all: export-package-python-modules #export-python-modules
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
    cptask.py \
    d3task.py \
    dostask.py \
    projwfctask.py \
    dynmattask.py \
    __init__.py \
    matdyntask.py \
    phtask.py \
    pwphmerger.py \
    pwtask.py \
    q2rtask.py \
    qetask.py \
    qetorque.py \
    setting.py \
    taskmerger.py \



__date__ = "$Jan 17, 2010 4:13:02 PM$"
