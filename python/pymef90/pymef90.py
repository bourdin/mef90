def getnumproc(prefix,pattern='%s-%04i.gen'):
    import os
    import os.path
    
    np = 0
    while os.path.isfile(pattern%(prefix,np)):
        np += 1
    return np
