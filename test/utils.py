import os
import sys

def set_search_paths(topdir):
    """Set search paths so that we can import Python modules"""
    pydir = os.path.join(topdir, 'pyext', 'src')
    os.environ['PYTHONPATH'] = pydir + os.pathsep \
                               + os.environ.get('PYTHONPATH', '')
    sys.path.insert(0, pydir)
