import numpy as np
import string
#from idf import Polarizations, Omega2
import qe_io_dict as io_dict
def matdyn(fname):
    matdynDict = io_dict.read_file(fname)
    qKeys = io_dict.find_all_keys_from_string(matdynDict, 'q =')
    qP = []
    for i in qKeys:
        qP.append( [ float(qi) for qi  in string.split( matdynDict[ i ] )[2:5] ] )
    qPoints = np.array(qP)
    
    
# find number of atoms  per unit cell and dimensionality       
# get frequencies keys for the last q-point:    
    fKeys = io_dict.find_all_keys_from_string_afterkey( matdynDict, qKeys[-1], 'omega')
#  get sample displacement vector for the first frequency of the last q point and find its dimensionality
#  here omegaShift = 1 means no space between displacements and frequencies (2 means 1 extra line ...)
    omegaShift = 1
    nDim =  len( matdynDict[ fKeys[0] + omegaShift  ].split()[1:-1] )/2
#  get number of atoms in unit cell    
    nAtom =  fKeys[1] - fKeys[0] - omegaShift

# qShift = 2 means 1 exra line between q line and frequencies line
    qShift = 2
    
    
#  create numpy array in accordance with idf format specification    

    
#    Pol = [ [ [ [ ] ] ] ]
    Pol = []
    Omegas = []
    for i in qKeys:
#        Mode = [ [ [ ] ] ]
        Mode = []
        Omega = []
        for iOmega in range( i + qShift, nDim*nAtom*(nAtom + omegaShift) +  i + qShift, nAtom+omegaShift):
	    # get omegas in THz:

#            print float( matdynDict[ iOmega].split('=')[1].split()[0] )
            Omega.append( float( matdynDict[ iOmega].split('=')[1].split()[0] ) )
#	    Atom = [ [ ] ]
            Atom = []
            for iAtom in range( iOmega + omegaShift, iOmega + omegaShift + nAtom ):
                vecStr = matdynDict[ iAtom ].split()[1:-1]
#		        vec = [  ]
                vec = [ float(v1) + 1.0j*float(v2) for v1,v2 in zip( vecStr[:-1:2], vecStr[1::2] ) ]
                Atom.append(vec)
            Mode.append( Atom )
        Pol.append( Mode )
        Omegas.append( Omega )
    npOmega = np.array(Omegas)
    npPol = np.array(Pol)
    THz2meV = 4.1357 # meV * s
    return npPol, npOmega*THz2meV, qPoints
    
    
if __name__ == '__main__':
    Pol, Omega, qPoints = matdyn( 'matdyn.modes' )
    Omega2.write( ((Omega/4.1357)**2)*1.e+24,'Omega2.idf','')
#    Polarizations.write(Pol, 'Polarizations.idf','Polarizations')
#    idfdata = Polarizations.read('Polarizations.idf')
#    idfdataOmega = Omega2.read('Omega2.idf')
#    print idfdata
#    print idfdataOmega
#    print qPoints
#    print len(qPoints[:,1])
