import os

def write_function(textfile, equationfile, path):

    filename = os.path.splitext(textfile)[0]
    filename = filename.replace(path + os.path.sep,'')

    with open(textfile, 'rU') as f:
        eqn = f.read().strip()

    with open(equationfile,'a+') as g:
        g.write("\n\ndef %s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    return %s" % (filename, eqn) )

    return 0


if __name__ == '__main__':
    python_equation_file = 'equations.py'
    path = os.getcwd() + os.path.sep + 'input' + os.path.sep + 'mathematica_equations'

    with open(python_equation_file,'w+') as f:
        f.write('import numpy as np\nimport os\nfrom heatmaps import KON, KP, T, KF')

    for r, d, f in os.walk(path):
        for file in f:
            if '.txt' in file:
                write_function(path + os.path.sep + file, python_equation_file, path)
