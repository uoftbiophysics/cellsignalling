import os

def write_function(textfile, equationfile, path):

    filename = os.path.splitext(textfile)[0]
    filename = filename.replace(path + os.path.sep,'')

    with open(textfile, 'rU') as f:
        eqn = f.read().strip()

    with open(equationfile,'a+') as g:
        g.write("\n\ndef %(filename)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    return %(eqn)s" %{'filename' : filename, 'eqn': eqn} )

    return 0

def write_function_2ligands(textfile, equationfile, path):

    filename = os.path.splitext(textfile)[0]
    filename = filename.replace(path + os.path.sep,'')

    with open(textfile, 'rU') as f:
        eqn = f.read().strip()

    with open(equationfile,'a+') as g:
        g.write("\n\ndef %(filename)s(c1, koff, c2, koff2, kon=KON, T=T, KP=KP, ALPHA=ALPHA):\n    return %(eqn)s" %{'filename' : filename, 'eqn': eqn} )

    return 0



if __name__ == '__main__':
    python_equation_file = 'equations2ligands.py'

    if python_equation_file == 'equations.py':
        path = os.getcwd() + os.path.sep + 'input' + os.path.sep + 'mathematica_equations'
        func_write = write_function
    elif python_equation_file == 'equations2ligands.py':
        path = os.getcwd() + os.path.sep + 'input' + os.path.sep + 'mathematica_equations_2ligands'
        func_write = write_function_2ligands

    with open(python_equation_file,'w+') as f:
        f.write('import numpy as np\nimport os\nfrom heatmaps import KON, KP, T, KF, ALPHA\n\ndef Sqrt(x):\n    return np.sqrt(x)')

    for r, d, f in os.walk(path):
        for file in f:
            if '.txt' in file:
                    func_write(path + os.path.sep + file, python_equation_file, path)

    # extra functions to write for the 1 ligand case
    if python_equation_file == 'equations.py':
        # writing function Det/(c^2 koff^2)
        for r, d, f in os.walk(path):
            for file in f:
                if 'DetSigmacrlb' in file:
                    with open(python_equation_file,'a+') as g:
                        filename = os.path.splitext(file)[0]
                        g.write("\n\ndef Rel%(filename)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    return ( %(filename)s(c, koff, kon, T, KF, KP) )/( c**2 * koff**2 )" %{'filename' : filename} )


        # trace equation and eigenvalue equations
        with open(python_equation_file,'a+') as g:

            # dedim relative errors
            for model in ['2', '3', '4', '2NoTrace','3NoTrace','4NoTrace']:
                # trace
                g.write("\n\ndef traceSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    return (SigmacrlbC%(name)s(c, koff, kon, T, KF, KP)+SigmacrlbK%(name)s(c, koff, kon, T, KF, KP))" %{'name' : model})
                # high and low eigenvalue
                g.write("\n\ndef evalplusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr + 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : model})
                g.write("\n\ndef evalminusSigmacrlb%(name)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    tr = traceSigmacrlb%(name)s(c, koff, kon, T, KF, KP); det = DetSigmacrlb%(name)s(c, koff, kon, T, KF, KP);\n    return ( 0.5*tr - 0.5*np.sqrt( tr**2-4*det ) )" %{'name' : model})

                for estimate in ['C', 'K']:
                    # dedimensionalized error (scaled by kp t)
                    g.write("\n\ndef dedimRelErr%(estimate)s%(model)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    return KP*T*RelErr%(estimate)s%(model)s(c, koff, kon, T, KF, KP)" %{'estimate' : estimate, 'model' : model})

            for model in ['1NoTrace']:
                for estimate in ['X']:
                    # edimensionalized error for X1 (scaled by kp t)
                    g.write("\n\ndef dedimRelError%(estimate)s%(model)s(c, koff, kon=KON, T=T, KF=KF, KP=KP):\n    return KP*T*RelError%(estimate)s%(model)s(c, koff, kon, T, KF, KP)" %{'estimate' : estimate, 'model' : model})

    # extra functions to write for the 1 ligand case
    if python_equation_file == 'equations2ligands.py':

        # matrix equations
        with open(python_equation_file,'a+') as g:
            for matrix in ["dmudthetaInv","sigmadata"]:
                g.write("\n\ndef matrix_%(matrix)s(c1, koff, c2, koff2, kon=KON, T=T, KP=KP,  ALPHA=ALPHA):\n    return np.array([[%(matrix)s11(c1, koff, c2, koff2, kon, T, KP, ALPHA), %(matrix)s12(c1, koff, c2, koff2, kon, T, KP, ALPHA), %(matrix)s13(c1, koff, c2, koff2, kon, T, KP, ALPHA), %(matrix)s14(c1, koff, c2, koff2, kon, T, KP, ALPHA)],[%(matrix)s21(c1, koff, c2, koff2, kon, T, KP, ALPHA), %(matrix)s22(c1, koff, c2, koff2, kon, T, KP, ALPHA), %(matrix)s23(c1, koff, c2, koff2, kon, T, KP, ALPHA), %(matrix)s24(c1, koff, c2, koff2, kon, T, KP, ALPHA)],[%(matrix)s31(c1, koff, c2, koff2, kon, T, KP, ALPHA), %(matrix)s32(c1, koff, c2, koff2, kon, T, KP, ALPHA), %(matrix)s33(c1, koff, c2, koff2, kon, T, KP, ALPHA), %(matrix)s34(c1, koff, c2, koff2, kon, T, KP, ALPHA)],[%(matrix)s41(c1, koff, c2, koff2, kon, T, KP, ALPHA), %(matrix)s42(c1, koff, c2, koff2, kon, T, KP, ALPHA), %(matrix)s43(c1, koff, c2, koff2, kon, T, KP, ALPHA), %(matrix)s44(c1, koff, c2, koff2, kon, T, KP, ALPHA)]])" %{'matrix' : matrix})
