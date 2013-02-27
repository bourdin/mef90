import fnmatch
import json
import os

from getpass import getuser
from time import asctime, localtime
from decimal import getcontext, localcontext, Decimal, Context, InvalidOperation

# Monkeypatch Decimal's __str__ to handle engineering notation better
def str_eng(self, eng=False, context=None):
    """Return string representation of the number in scientific notation.

    Captures all of the information in the underlying representation.
    """

    sign = ['', '-'][self._sign]
    if self._is_special:
        if self._exp == 'F':
            return sign + 'Infinity'
        elif self._exp == 'n':
            return sign + 'NaN' + self._int
        else: # self._exp == 'N'
            return sign + 'sNaN' + self._int

    # number of digits of self._int to left of decimal point
    leftdigits = self._exp + len(self._int)

    # dotplace is number of digits of self._int to the left of the
    # decimal point in the mantissa of the output string (that is,
    # after adjusting the exponent)
    if self._exp <= 0 and leftdigits > -2:
        # no exponent required
        dotplace = leftdigits
    elif not eng:
        # usual scientific notation: 1 digit on left of the point
        dotplace = 1
    elif self._int == '0':
        # engineering notation, zero
        dotplace = (leftdigits + 1) % 3 - 1
    else:
        # engineering notation, nonzero
        dotplace = (leftdigits - 1) % 3 + 1

    if dotplace <= 0:
        intpart = '0'
        fracpart = '.' + '0'*(-dotplace) + self._int
    elif dotplace >= len(self._int):
        intpart = self._int+'0'*(dotplace-len(self._int))
        fracpart = ''
    else:
        intpart = self._int[:dotplace]
        fracpart = '.' + self._int[dotplace:]
    if leftdigits == dotplace:
        exp = ''
    else:
        if context is None:
            context = getcontext()
        exp = ['e', 'E'][context.capitals] + "%+d" % (leftdigits-dotplace)

    return sign + intpart + fracpart + exp

Decimal.__str__ = str_eng

def eng(count, prec=6):
  """ Return the count or an approximately equal value in engr. notation. """
  with localcontext(Context(prec=prec)):
    return (0 + Decimal(count)).to_eng_string().lower()

def scan(outdir, paths, excludedirs, excludefiles, project, save_json=False):
  """ Scan paths for jobs, output json file containing parsed data. """
  data = scan_dirs(paths, excludedirs, excludefiles)
  data['outdir'] = outdir
  data['project'] = project
  data['date'] = asctime(localtime())
  data['user'] = getuser()

  if save_json:
    with open(os.path.join(outdir, 'data.json'), 'w') as outdata:
      json.dump(data, outdata, sort_keys=True, indent=2)

  return data

def scan_dirs(paths, exdirs, exfiles):
  """ Walk the directory specified by path and find all computations.
      Exclude files specified by the wildcard excludes.  """

  db = {}
  jobs = db['jobs'] = {}
  headers = db['headers'] = []

  for path in paths:
    for root, dirs, files in os.walk(path):
      # Don't walk directories we don't care about
      for d in [x for ex in exdirs for x in fnmatch.filter(dirs, ex)]:
        dirs.remove(d)

      if not dirs and '00_INFO.txt' in files:
        job = os.path.split(root)[1]
        jobdict = jobs[job] = {}

        with open(os.path.join(root, '00_INFO.txt')) as jobinfo:
          jobdata = {}
          for l in jobinfo:
            l = l.strip()
            k, v = l.split(' ', 1) if l.count(' ') > 0 else (l, '')
            if k not in headers:
              headers.append(k)
            v = v.strip()
            # We check to see if the value is a number and try to cast it
            if '.' in v:
              try:
                v = eng(v)
              except InvalidOperation:
                pass
            jobdata[k] = v
          jobdict['info'] = jobdata

        # Remove files from output that we excluded in command line
        for f in [x for ex in exfiles for x in fnmatch.filter(files, ex)]:
          files.remove(f)
        files.remove('00_INFO.txt')

        jobdict['thumbable'] = [x for ext in ('*.png', '*.pdf', '*.mov', '*.avi')
                    for x in fnmatch.filter(files, ext)]

        # Set the file listing for the job
        jobdict['files'] = files
        jobdict['path'] = root

  return db

def check_mtime(infile, outfile):
  if not os.path.exists(infile):
    return True;
  if os.path.exists(outfile):
    istat, ostat = os.stat(infile), os.stat(outfile)
    if ostat.st_mtime > istat.st_mtime:
      return True
  return False

