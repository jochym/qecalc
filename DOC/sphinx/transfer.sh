#!/usr/bin/env sh

ROOT_UID=0   # Root has $UID 0.

if [ "$UID" -eq 0 ]  # Will the real "root" please stand up?
then
  su jbk #become someone with permission to move the docs
fi


svn up
make html
scp -r _build/html jbrkeith@login.cacr.caltech.edu:qecalc
