
def addVec(df, raCol='ditheredRA', decCol='ditheredDec'):

    thetas  = - df[decCol] + np.pi /2.
    phis = df[raCol]
    df['vec'] = list(hp.ang2vec(thetas, phis))
