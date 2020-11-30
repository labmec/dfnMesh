def MatrixToString2(fmatrix):
    rows = len(fmatrix)
    cols = len(fmatrix[0])
    fstring = ""
    for i in range(rows):
        for j in range(cols):
            if(p0[i][j] >= 0.0): fstring += ' '
            fstring += str(p0[i][j])+' '
        fstring += "\n"
    print(fstring)
    return fstring