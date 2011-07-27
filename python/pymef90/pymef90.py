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


def getnumproc(prefix,pattern='%s-%04i.gen'):
    import os
    import os.path
    
    np = 0
    while os.path.isfile(pattern%(prefix,np)):
        np += 1
    return np

def argsWrite(infilename,outfilename,Dict):
    infile=open(infilename,'r')
    outfile=open(outfilename,'w')
    outfile.write(infile.read()%Dict)
    outfile.close()
    infile.close()

