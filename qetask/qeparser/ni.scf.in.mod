&CONTROL
    calculation = 'nscf',
    wf_collect = .true.,
    tstress = .true.,
    tprnfor = .true.,
    verbosity = 'high',
    prefix = 'mgb2',
    pseudo_dir = 'mgb2',
    lkpoint_dir = .false.,
    outdir = 'temp',
    title = Ni,
/
&SYSTEM
    ibrav = 4,
    celldm(1) = 5.78739785,
    celldm(2) = 5.78739785,
    celldm(3) = 1.135794331,
    nat = 3,
    ntyp = 2,
    nspin = 1,
    nbnd = 12,
    occupations = 'smearing',
    degauss = 0.025,
    smearing = 'methfessel-paxton',
    ecutwfc = 32.0,
    ecutrho = 256.0,
    la2f = .false.,
/
&ELECTRONS
    conv_thr = 1.0d-12,
    diago_full_acc = .TRUE.,
/
ATOMIC_POSITIONS (alat)
 Say Hi! :)
K_POINTS (automatic)
 24 24 24 0 0 0
