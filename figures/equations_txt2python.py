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

    # trace equation and eigenvalue equations
    with open(python_equation_file,'a+') as g:
        g.write("\n\ndef traceSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    return (SigmacrlbC%(name)s(c, koff, kon, T, KF, KP)+SigmacrlbK%(name)s(c, koff, kon, T, KF, KP))" %{'name' : '2'})
        g.write("\n\ndef traceSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    return (SigmacrlbC%(name)s(c, koff, kon, T, KF, KP)+SigmacrlbK%(name)s(c, koff, kon, T, KF, KP))" %{'name' : '3'})
        g.write("\n\ndef traceSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    return (SigmacrlbC%(name)s(c, koff, kon, T, KF, KP)+SigmacrlbK%(name)s(c, koff, kon, T, KF, KP))" %{'name' : '4'})
        g.write("\n\ndef traceSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    return (SigmacrlbC%(name)s(c, koff, kon, T, KF, KP)+SigmacrlbK%(name)s(c, koff, kon, T, KF, KP))" %{'name' : '2NoTrace'})
        g.write("\n\ndef traceSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    return (SigmacrlbC%(name)s(c, koff, kon, T, KF, KP)+SigmacrlbK%(name)s(c, koff, kon, T, KF, KP))" %{'name' : '3NoTrace'})
        g.write("\n\ndef traceSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    return (SigmacrlbC%(name)s(c, koff, kon, T, KF, KP)+SigmacrlbK%(name)s(c, koff, kon, T, KF, KP))" %{'name' : '4NoTrace'})
        g.write("\n\ndef evalplusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : '2'})
        g.write("\n\ndef evalminusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr - 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : '2'})
        g.write("\n\ndef evalplusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : '3'})
        g.write("\n\ndef evalminusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr - 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : '3'})
        g.write("\n\ndef evalplusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : '4'})
        g.write("\n\ndef evalminusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr - 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : '4'})
        g.write("\n\ndef evalplusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : '2NoTrace'})
        g.write("\n\ndef evalminusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr - 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : '2NoTrace'})
        g.write("\n\ndef evalplusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : '3NoTrace'})
        g.write("\n\ndef evalminusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr - 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : '3NoTrace'})
        g.write("\n\ndef evalplusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : '4NoTrace'})
        g.write("\n\ndef evalminusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr - 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : '4NoTrace'})
