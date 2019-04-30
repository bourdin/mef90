def Dictwritetxt(D,filename,overwrite=True):
    if overwrite:
        f = open(filename,'w')
    else:
        f = open(filename,'a')
    K = D.keys()
    K.sort()
    for key in K:
        f.write('%s \t\t %s\n'%(key, str(D[key])))
    f.close()
  
def DictwriteJSON(D,filename,overwrite=True):
    try:
        import json
        if overwrite:
            jsonfile = open(filename,'w')
        else:
            jsonfile = open(filename,'a')

        json.dump(D,jsonfile, sort_keys = True, indent = 4)
        jsonfile.close()
    except ImportError:
        print ('JSON module not available, skipping DictJSONwrite')


def Dictreadtxt(filename):
    D = {}
    for l in open(filename,'r').readlines():
        l = l.strip()
        k, v = l.split(' ', 1) if l.count(' ') > 0 else (l, '')
        v = v.strip()
        try:
            v = float(v)            
            if int(v) == v:
                v = int(v)
        except ValueError:
            pass
        D[k] = v
    # special case: make sure that prefix is read as a string or all hell will break loose...
    try:
        D['prefix']+'oulala'
    except TypeError:
        D['prefix'] = repr(D['prefix'])
    return D
    

