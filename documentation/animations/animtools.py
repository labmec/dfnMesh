# import os,sys,inspect
# currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# parentdir = os.path.dirname(currentdir)
# sys.path.insert(0,parentdir) 

# import animtools



def MatrixToString(fmatrix):
    rows = len(fmatrix)
    cols = len(fmatrix[0])
    fstring = ""
    for i in range(rows):
        for j in range(cols):
            if(p0[i][j] >= 0.0): fstring += ' '
            fstring += str(p0[i][j])+' '
        fstring += "\n"
    # print(fstring)
    return fstring

