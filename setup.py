#!/usr/bin/env python


"""QECalc2- wrapper for QE.

"""

from setuptools import setup, find_packages

# define distribution
setup(
        name = "qecalc",
        version = "1.0",
        #namespace_packages = ['matter'],
        packages = find_packages(exclude=['DOC','examples']),
        test_suite = 'tests',
#        install_requires = [
#            'PyCifRW',
#        ],
#        dependency_links = [
#            'http://www.diffpy.org/packages/',
#        ],

        author = 'Nikolay',
        author_email = 'jbrkeith@gmail.edu',
        description = "wrapper for QE",
        license = 'BSD',
        keywords = "quantum espresso",
        classifiers = [
            # List of possible values at
            # http://pypi.python.org/pypi?:action=list_classifiers
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'Operating System :: MacOS',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python :: 2.5',
            'Topic :: Scientific/Engineering :: Physics',
        ],
)

# End of file
