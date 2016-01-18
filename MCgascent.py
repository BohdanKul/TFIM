import pylab    as py
import numpy    as np
import numpy.ma as ma
import kevent, argparse, tfim, MCstat, Hbuilder, Hfile

#------------------------------------------------------------------------------
# Inversly-proportional annealing schedule.
# LL - left limit
# RL - right limit
# nsteps - number of steps
#------------------------------------------------------------------------------
def IPschedule(LL, RL, nsteps):
    b    = (1.0*nsteps*RL)/(LL - RL)
    a    = LL*b

    return a, b

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def checkAutoCorr(data, dim=0):
    return False
    
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def getStats(data,dim=0):
    ''' Get the average and error of all columns in the data matrix. '''

    dataAve = py.average(data,dim) 
    bins    = MCstat.bin(data) 
    dataErr = py.amax(bins,axis=0)

    return dataAve, dataErr


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def detTotalBins(errs, weights, r, bins):
    ''' Estimate the optimal number of extra bins to take for each vector
        in the dataset in order to converge the errors below the tolerance 
        treshold. The later is controlled by r. It is a fraction < 1 and 
        is computed for each average considered in the simulation. The 
        estimate is obtained with the help of Lagrange multipliers 
        method. The function to minize is sum(bins/r2s) that is the estimate
        of the total number of bins. '''
    print 'bins:    ', bins
    print 'weights: ', weights
    print 'errors:  ', errs
    print 'r:       ', r
    bins = np.transpose(np.atleast_2d(np.array(bins)))           # convert vector to column-shaped
    iLambda = np.sum((errs*errs*r*r*weights*weights),   axis=0)  # compute Lagrange multiplier
    iLambda/= np.sum(errs*np.abs(weights)*(np.sqrt(bins)), axis=0) 
    r2s     = np.sqrt(bins)/(np.abs(weights)*errs)*iLambda               
    
    print 'r2s:    ', r2s
    maxr2s = np.amin(r2s, axis=1)    # find the smallest ones
    print 'maxr2s: ', maxr2s
    # compute
    tbins  = np.ceil(np.transpose(np.atleast_2d(bins))/maxr2s)
    print 'tbins: ', tbins
    return np.squeeze(tbins).tolist()


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def simpleMovingAverage(period,data):
    assert period == int(period) and period > 0, "Period must be an integer >0"
  
    #avoid averaging over a larger number of elements than there is 
    period = min(size(data,0),period)

    import numpy as np
    weightings = np.repeat(1.0, period) / period 
    return np.convolve(data, weightings, 'valid')

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--data',       help='Data file',                  type=str)
    parser.add_argument('--seed',       help='Seed used to generate data', type=int, default=0)
    parser.add_argument('--inter','-I', help='Interactions file',          type=str)
    parser.add_argument('--beta', '-b', help='Inverse temperature ',       type=float)
    args = vars(parser.parse_args())

    
    # Load data generating Hamiltonian ----------------------------------------
    beta = args['beta']
    Ns   = 0
    (gHs, gDs, gJs, bonds) = (np.array([]), np.array([]), np.array([]), np.array([]))  
    if (args['inter'] is not None):
        gHs, gDs, gJs, bonds = Hfile.LoadInters(args['inter'])
        Ns = len(gHs)
    else:
        print "Error: enter interaction file"
        return 1

    # Convert bonds datastructure to the one compatible with TFIM class
    t = []
    for l in bonds.tolist():
        t.append((l[0],l[1]))
    bonds = t
    Nb = len(bonds) 
    
    # Generate initial guess for the couplings---------------------------------
    Hs = np.random.randn(Ns)*max([max(gHs),0.5])+gHs
    Ds = np.abs(np.random.randn(Ns)*max([max(gDs),0.5]))+gDs
    Js = np.random.randn(Nb)*max([max(gJs),0.5])+gJs
    
    # Obtain or generate a training dataset -----------------------------------
    data, weights = Hfile.GetData(args['data'], Ns, args['seed'], 10000, bonds, gHs, gDs, gJs, beta) 
    data    = data[:2]
    weights = weights[:2]
    data += [[]]          # add a vector with no bits clamped
    Nd = len(data)
    weights = np.append(weights, -beta) # and its weight
    weights = weights.reshape((Nd,1))  # take a transpose
    #weights = np.transpose(np.atleast_2d(weights))# take a transpose
    
    # Set up MC constants------------------------------------------------------ 
    mSlice = 1025  # minimum number of bins to take in attempt to converge the autocorrelation
    tol = 0.1      # fractional error tolerance in estimation of the gradient
    
    # Set up gradient descent constants----------------------------------------
    a, b = IPschedule(0.9/beta, 0.1/beta, 200)
    Norm = 100.0 # gradient norm 
    step = 0     # step in gradient descent
    nsteps = 1000

    # the gradient
    gval = np.zeros(len(gHs)+len(gDs)+len(gJs)) 
    subset = range(Nd)

    # Gradient descent --------------------------------------------------------
    while ((step < nsteps) and (Norm > 0.00001)):
    
        # Initiate MC solvers -----------------------------------------------------
        TFIMs = []                      # list of MC objects, one for each data vector
        (Zs,   Xs,  ZZs) = ([], [], []) # lists of raw measurements for each data vector
        for vector in data:
            TFIMs.append(tfim.TFIM(vector, bonds, Hs.tolist(), Ds.tolist(), Js.tolist(), Ns,beta, 1))
            Zs.append([]); Xs.append([]); ZZs.append([])
         
        # Reset necessary datastructures ------------------------------------------
        (bins, mtargets) = (np.zeros(Nd), mSlice*np.ones(Nd)) # bins taken and total number of bins to take 
        (aZs, aXs, aZZs) = (np.zeros((Nd, Ns)), np.zeros((Nd, Ns)), np.zeros((Nd, Nb)))  # estimators' averages
        (eZs, eXs, eZZs) = (np.zeros((Nd, Ns)), np.zeros((Nd, Ns)), np.zeros((Nd, Nb)))  # and errors   
        (isAutoCorr, isUncertain) = (True, True) # flags controlling simulation runtime

        # run MC until the error is converged and is below the error tolerance level
        while isAutoCorr or isUncertain:
            # for each vector in dataset
            for i in range(Nd):
                # run a pre-determined number of measurements
                while (bins[i] < mtargets[i]):
                    for sweep in range(100):
                        # diagonal update
                        while  TFIMs[i].DMove() == 1:  # if the operator list is too small
                               TFIMs[i].Adjust()       # adjust its length 
                               bins[i] = 0             # set bin count to zero
                               (Zs[i], Xs[i], ZZs[i]) = ([], [], []) # reset timeseries
                               isAutoCorr = True       # signal possible autocorrelation 
                        # off-diagonal update
                        TFIMs[i].ODMove()
                    
                    # take a measurement with the help of auxilary variables
                    bins[i] += 1
                    (tZs, tXs, tZZs) = ([],[],[]) # reset auxilary variables 
                    TFIMs[i].Measure(tZs, tXs, tZZs)
                    # accumulate them in lists
                    Zs[i].append(tZs); Xs[i].append(tXs); ZZs[i].append(tZZs)
                    
            # check whether the measurements are autocorrelated
            if  isAutoCorr: 
                subset  = [i for i in subset if (checkAutoCorr(Zs[i]) or checkAutoCorr(Xs[i]) or checkAutoCorr(ZZs[i]))]
                isAutoCorr = len(subset) > 0
                # assign more measurements to take for samples with high autocorrelation
                if isAutoCorr: mtargets[np.array(subset)] += mSlice
                # if all samples are converged, reset the subset of samples to the full set 
                
            # start taking averages only when all chains have converged errors
            if  not isAutoCorr:
                for i in range(Nd):
                    ((aZs[i], eZs[i]), (aXs[i], eXs[i]), (aZZs[i], eZZs[i])) = (getStats(np.array(Zs[i])), getStats(np.array(Xs[i])), getStats(np.array(ZZs[i])))
                # stack measurements for an easy manipulation
                aves = np.hstack([aZs, aXs, aZZs])
                errs = np.hstack([eZs, eXs, eZZs])
                
                # compute the gradient components averages and errors
                gval = np.sum((aves*weights),  axis=0)
                gerr = np.sqrt(np.sum((errs*errs*weights*weights), axis=0))
                
                # determine the indices of averages falling below the tolerance treshold 
                #indices = ma.getmask(ma.masked_where(gval*tol-gerr < 0, ))
                print 'bins:    ', bins
                print 'targets: ', mtargets 
                print 'errors:  ', errs[errs<0]
                indices = np.squeeze(np.where(gval*tol-gerr < 0))
                # if there are no such averages, signal that the errors are within tolerance bounds
                if len(indices) == 0: isUncertain = False
                # otherwise, estimate how many more bins to take
                else: mtargets = detTotalBins(errs[:, indices], weights, (np.abs(gval)*tol/gerr)[indices], bins)
                
                return 0
        # Destroy MC solvers -----------------------------------------------------
        for i in range(len(TFIMs)):
            del TFIMs[i]
        
        # Annealing constant 
        eta   = a/(b+1.0*step)

        # Follow the negative gradient 
        (dZ, dX, dZZ) = np.split(gval, [Ns, Ns+Ns])    
        Zs   += -eta*dZ
        Xs   += -eta*dX
        ZZs  += -eta*dZZ

    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
    #fig = py.figure(1,figsize=(13,6))
    #ax  = py.subplot(111)
    #py.connect('key_press_event',kevent.press)
    #ax.plot(Zs[:,0],color=colors[0],linewidth=1,marker='None',linestyle='-')
    #ax.plot(Zs[:,1],color=colors[1],linewidth=1,marker='None',linestyle='-')
    #ax.plot(Zs[:,2],color=colors[2],linewidth=1,marker='None',linestyle='-')
    #py.show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
