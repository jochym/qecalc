#!/usr/bin/env python


"""QECalc - wrapper for Quantum Espresso.

"""

from setuptools import setup, find_packages

# define distribution
setup(
        name = "qecalc",
        version = "0.2-rc1",
        #namespace_packages = ['matter'],
        packages = find_packages(exclude=['DOC','examples']),
        test_suite = 'examples',
        install_requires = [
            'diffpy.Structure',
            'numpy',
        ],
        dependency_links = [
            'http://www.diffpy.org/packages/',
        ],

        author = 'Nikolay Markovskiy',
        author_email = 'markovskiy@gmail.com',
        description = "wrapper for Quantum Espresso",
        license = 'BSD',
        keywords = "quantum espresso",
        url = "http://docs.danse.us/AbInitio/espresso/qecalc/index.html",
        download_url = 'http://danse.cacr.caltech.edu/packages/dev_danse_us/',
        classifiers = [
            # List of possible values at
            # http://pypi.python.org/pypi?:action=list_classifiers
            'Development Status :: 2 - Pre-Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            #'Operating System :: MacOS',
            #'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python :: 2.5',
            'Topic :: Scientific/Engineering :: Physics',
        ],
)

# End of file
