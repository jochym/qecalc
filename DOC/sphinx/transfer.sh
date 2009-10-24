#!/usr/bin/env bash

if [[ $USER -eq jbk ]] then 
	svn up
	make html
	scp -r _build/html/* jbrkeith@login.cacr.caltech.edu:projects/danse/docs.danse.us/docroot/AbInitio/qecalc
	# hi nikolay...just do a script like the one above...
elif [[ $USER -eq root ]] then 
	su jbk #become someone with permission to move the docs 
fi