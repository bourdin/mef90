import argparse

from sanson import Action
from scanner import scan, check_mtime

class MEF90Action(Action):
  """ Default action to help define common args"""
  def __init__(self):
    self.parser = argparse.ArgumentParser(add_help=False)
    self.parser.add_argument('paths', metavar='path', nargs='+',
                             help='paths to include in the database')
    self.parser.add_argument('--outdir', help='path to output directory, default is current dir',
                             default='.')
    self.parser.add_argument('--project', help='name of project',
                             default='')
    self.parser.add_argument('--excludedir', action='append',
                             dest='excludedirs', metavar='DIRS', default=[],
                             help='May be used multiple times to exclude directories')
    self.parser.add_argument('--excludefile', action='append',
                             dest='excludefiles', metavar='FILES', default=[],
                             help='May be used multiple times to exclude files')

    self.parser.add_argument('--plotener', default=False,action='store_true',
                             help='Generates an energy plot')
    self.parser.add_argument('--plotpng', default=False,action='store_true',
                             help='Plot last step [req. VisIt]')
    self.parser.add_argument('--movtransient', default=False,action='store_true',
                             help='generates a transient movie [req. VisIt]')
    self.parser.add_argument('--mov3d', default=False,action='store_true',
                             help='generates a 3D movie [req. VisIt]')
    self.parser.add_argument('--all', default=False,action='store_true',
                             help='Do all of the above except movies')

__all__ = [scan, check_mtime, MEF90Action]
