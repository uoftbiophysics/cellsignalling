import csv
import numpy as np
import os

from settings import DIR_INPUT


def load_input(filename, verbose=True):
    with open(os.path.join(os.getcwd(), DIR_INPUT, filename + '.csv'), 'r') as f:
        reader = csv.reader(f)
        read_list = list(reader)
        if verbose:
            print(filename)
        numeric_read_list = [(read_list[0][0],read_list[0][1])] + [(float(el[0]), float(el[1])) for el in read_list[1:]]
    return numeric_read_list


def get_csvnames():
    files = os.listdir(DIR_INPUT)
    valid_names = [f[:-4] for f in files if f[-4:] == '.csv']
    return valid_names


def gen_datadict(verbose=True):
    """
    Stores all input csv files in dict of the form
        {csvfilename -> {xlab, ylab, xpts, ypts} }
    """
    print("Loading expressions exported from Mathematica")
    datadict = {}
    csvnames = get_csvnames()
    for csvname in csvnames:
        list_of_tuple = load_input(csvname, verbose=verbose)
        header = list_of_tuple[0]
        arr = np.array(list_of_tuple[1:])
        datadict[csvname] = {'xlab': header[0], 
                             'ylab': header[1],
                             'xpts': arr[:, 0],
                             'ypts': arr[:, 1]}
    return datadict

if __name__ == '__main__':
    DATADICT = gen_datadict()
