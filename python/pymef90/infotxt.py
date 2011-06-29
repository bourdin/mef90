def mygetenv(Dict,key,defaultvalue=None):
  ### I could do better 
  import os
  tmp = os.getenv(key)
  if tmp == None:
    Dict[key] = defaultvalue
  else:
    try:
      Dict[key] = float(tmp)
    except ValueError:
      Dict[key] = tmp

def Dictwritetxt(D,filename):
  f = open(filename,'a')
  K = D.keys()
  K.sort()
  for key in K:
    f.write('%s \t\t %s\n'%(key, str(D[key])))
  f.close()
  
def DictwriteJSON(D,filename):
	try:
		import json
	
		jsonfile = open(filename,'w')
		jsonfile.write(json.encoder.JSONEncoder().encode(D))
		jsonfile.flush()
		jsonfile.close()
	except ImportError:
		print 'JSON module not available, skipping DictJSONwrite'


def Dictreadtxt(filename):
    D = {}
    for l in open(filename,'r').readlines():
        l = l.strip()
        k, v = l.split(' ', 1) if l.count(' ') > 0 else (l, '')
        v = v.strip()
        try:
            v = float(v)
        except ValueError:
            pass
        D[k] = v
    return D
    

