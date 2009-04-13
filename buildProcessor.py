#!/usr/bin/env python

import sys

def run():
  for line in sys.stdin.readlines():
    sys.stdout.write(line)
  return

if __name__ == '__main__':
  run()
