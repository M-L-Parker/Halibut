
from numpy import where, zeros, floor, mean, std, newaxis
import numpy as np
from matplotlib import pyplot as p

def dcfirr (ta,ra,tb,rb,lag=None,minpt=100,minlag=-30.,maxlag=30.,numf=0,numpt=0, indices=[-1], indfinal=[-1]):

    # Lag convention is lag = ta-tb, so positive lag values show the A lightcurve lagging behind the B.

    #;; avoid overflows
    minpt=long(minpt)

    #;;; --- Lag terms --- ;;;

    #;; Only calculate lag terms if array of lag index pairs isn't given (for speed)
    if len(indices) != (ta.size*tb.size) :

        #;; calculate lags between all time pairs
        lagarr = ta - tb[newaxis].T 

        #;; sort lag values
        lagarr = lagarr.ravel()
        indices = lagarr.argsort()

        laglist=lagarr[indices]#.ravel()

        #pl=p.plot(laglist)
        #;; bin lags by number of points
        numf=int(floor( laglist.size/minpt))
        lag=zeros(numf)
        for i in range(numf) :
            lag[i]=mean( laglist[ (i*minpt) : ((i+1)*minpt - 1) ] )

        #;; extract lags between minlag and maxlag
        (indfinal,) = where( ((minlag < lag)&(lag < maxlag)) )

        if indfinal.size > 0 :
            lag=lag[indfinal]
        else :
            print('no lags in given range')
            return float('NaN'),float('NaN')

        numf=lag.size

        numpt=zeros(numf)+minpt


    #;;; --- End lag terms --- ;;



    #;;; --- rate terms --- ;;;

    #;; Normalise lightcurve distribution - mean subtract here, correct variance at end
    ranorm=ra-mean(ra)
    rbnorm=rb-mean(rb)

    #;; calculate pairwise correlations
    corarr=ranorm * rbnorm[newaxis].T

    #;; sort correlations by lag
    corarr=corarr.ravel()
    corlist=corarr[indices]


    #;; bin correlation pairs
    cor=zeros(numf)
    for i in range(numf) :
        cor[i]=mean( corlist [ (indfinal[i]*minpt):((indfinal[i]+1)*minpt-1) ] )


    #;; correct for variance of original distribution
    cor=cor/std(ra)/std(rb)

    #;;; --- End rate terms --- ;;;

    return lag,cor,numf,indices,indfinal

    # END
