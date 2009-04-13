#!/usr/bin/env python

import sys

def run():
  import re

  newline = re.compile(r'(\'[\\a-zA-Z\_,\(\) ]+\')c')
  for line in sys.stdin.readlines():
    sys.stdout.write(re.sub(newline, r'\1', line))
  return

if __name__ == '__main__':
  run()
