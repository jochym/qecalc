from qe_io_dict import read_file, find_key_from_string, save_dic

def varnameValue(fname, var, val):
    'Dummy! Assumes no comments for a variable in question and one variable per string'
    import string
    pwscfDic = read_file(fname)
    key = find_key_from_string(pwscfDic, var)
    pwscfDic[key] = '    ' + var + ' = ' + str(val) + ',\n'
    save_dic(pwscfDic, fname)

def k_points(fname, k_points):
    'Dummy! Assumes no comments and implies k-point values are in the next line after K_POINTS AUTOMATIC'
    import string
    pwscfDic = read_file(fname)
    key = find_key_from_string(pwscfDic, 'K_POINTS')
    k_string = string.join([str(v) + ' ' for v in k_points])
    pwscfDic[key+1] = k_string + '\n'
    save_dic(pwscfDic, fname)
    
def atomic_positions(fname, geometry):
    'Dummy! Assumes geometry description is right after "ATOMIC_POSITIONS" (no empty lines)'
    import string
    pwscfDic = read_file(fname)
    key = find_key_from_string(pwscfDic, 'ATOMIC_POSITIONS')
    for i in range(len(geometry)):
        pwscfDic[key+1+i] = geometry[i] + '\n'
    save_dic(pwscfDic, fname)