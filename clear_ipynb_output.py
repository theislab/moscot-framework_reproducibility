
import sys
import io
import os
from nbformat import read, write, NO_CONVERT


def remove_outputs(nb):
    """remove the outputs from a notebook"""
    cell_rm = []
    for cell in nb.cells:
        cell.outputs = []
    nb.metadata.history = []


if __name__ == '__main__':
    fname = sys.argv[1]
    with io.open(fname, 'r') as f:
        nb = read(f, NO_CONVERT)
    remove_outputs(nb)
    base, ext = os.path.splitext(fname)
    new_ipynb = "%s_removed%s" % (base, ext)
    with io.open(new_ipynb, 'w', encoding='utf8') as f:
        write(nb, f, NO_CONVERT)
    print("wrote %s"% new_ipynb)
