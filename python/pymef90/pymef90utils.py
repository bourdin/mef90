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
    outfile.write(infile.read().format(**Dict))
    outfile.close()
    infile.close()

def argsWriteOldStyle(infilename,outfilename,Dict):
    infile=open(infilename,'r')
    outfile=open(outfilename,'w')
    outfile.write(infile.read()%Dict)
    outfile.close()
    infile.close()

def confirm(prompt=None, resp=False):
    """prompts for yes or no response from the user. Returns True for yes and
    False for no.

    'resp' should be set to the default value assumed by the caller when
    user simply types ENTER.

    >>> confirm(prompt='Create Directory?', resp=True)
    Create Directory? [y]|n: 
    True
    >>> confirm(prompt='Create Directory?', resp=False)
    Create Directory? [n]|y: 
    False
    >>> confirm(prompt='Create Directory?', resp=False)
    Create Directory? [n]|y: y
    True

    """
    from sys import version_info

    if prompt is None:
        prompt = 'Confirm'

    if resp:
        prompt = '%s [%s]|%s: ' % (prompt, 'y', 'n')
    else:
        prompt = '%s [%s]|%s: ' % (prompt, 'n', 'y')
        
    while True:
        if version_info.major == 3:
            ans = input(prompt)
        else:
            ans = raw_input(prompt)
        if not ans:
            return resp
        if ans not in ['y', 'Y', 'n', 'N']:
            print ('please enter y or n.')
            continue
        if ans == 'y' or ans == 'Y':
            return True
        if ans == 'n' or ans == 'N':
            return False
