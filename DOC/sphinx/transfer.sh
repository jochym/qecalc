#!/usr/bin/env sh

ROOT_UID=0   # Root has $UID 0.

if [[ $USER -eq root ]]
then 
	su jbk #become someone with permission to move the docs 
fi

if [[ $USER -eq jbk ]]
then 
	svn up
	make html
	scp -r _build/html/* jbrkeith@login.cacr.caltech.edu:qecalc
fi

# hi nikolay...just do a script like the one above...