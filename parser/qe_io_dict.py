# Dictionary is assumed to be created from a text file. Keys are
# integers, corresponding to line numbers 

def read_file(fname):
    fileIn = open(fname, "r")
    line = fileIn.readline()
    dic = {}
    keycounter = 1
    while line:
        key = keycounter
        dic[key] = line
        keycounter = keycounter + 1
        line = fileIn.readline()
    return dic
    
def save_dic(dic, fname):
    fileOut = open(fname, 'w')
    for k,v in dic.iteritems(): fileOut.write(dic[k])    
    fileOut.close()
    
def find_key_from_marker_string(dic, marker, val):
    """return the key of dictionary dic given the string marker and value"""
    return [k for k, v in dic.iteritems() if v[0] == marker and val in v][0]
    
def find_all_keys_from_marker_string(dic, marker, val):
    """return the key of dictionary dic given the string marker and value"""
    return [k for k, v in dic.iteritems() if v[0] == marker and val in v]

def find_all_keys_from_string(dic, val):
    """return the key of dictionary dic given the string marker and value"""
    return [k for k, v in dic.iteritems() if val in v]

def find_key_from_string(dic, val):
    """return the key of dictionary dic given the value"""
    return [k for k, v in dic.iteritems() if val in v][0]
    
def find_key_from_string_afterkey(dic, key0, val):
    """return the key of dictionary dic given the value and starting line number"""
    nmaxKey = len(dic)
    for ikey in range(key0, nmaxKey+1):
        if val in dic[ikey]: return ikey
    assert True, "find_key_from_string_afterkey: The value is not in the dictionaty"
    
def find_all_keys_from_string_afterkey(dic, key0, val):
    """return the key of dictionary dic given the value and starting line number"""
    nmaxKey = len(dic)
    return [ ikey for ikey in range(key0, nmaxKey+1) if val in dic[ikey] ]
    assert True, "find_key_from_string_afterkey: The value is not in the dictionaty"
  
