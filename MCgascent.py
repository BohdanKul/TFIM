import pylab    as py
import numpy    as np
import numpy.ma as ma
import numpy.random as rm
import kevent, argparse, tfim, MCstat, Hbuilder, Hfile, bmachine

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
        method. The function to minize is sum(bins/r2s) - an estimate
        of the total number of bins required to achieve convergence. '''
    print 'bins:    ', bins
    print 'weights: ', weights.T
    print 'errors:  ', errs
    print 'r:       ', r
    bins = np.transpose(np.atleast_2d(np.array(bins)))           # convert vector to column-shaped
    iLambda = np.sum((errs*errs*r*r*weights*weights),   axis=0)  # compute Lagrange multiplier
    iLambda/= np.sum(errs*np.abs(weights)*np.sqrt(bins), axis=0) 
    r2s     = np.sqrt(bins)/(np.abs(weights)*errs)*iLambda               
    
    #print 'r2s:    ', r2s
    maxr2s = np.amin(r2s, axis=1)    # find the smallest ones
    #print 'maxr2s: ', maxr2s
    # compute
    bins   = np.squeeze(np.transpose(np.atleast_2d(bins)))
    tbins  = np.squeeze(np.ceil(bins/maxr2s))
    print 'projected bins: ', tbins
    maxbins = 10000
    tbins = np.asarray(np.minimum(tbins, bins+np.ones_like(tbins)*maxbins))
    tbins = tbins.tolist()
    if isinstance(tbins, float): tbins = [tbins]
    return tbins


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

    modes = ['aquant', 'equant'] 
    # setup the command line parser options 
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--data',       help='Data file',                  type=str)
    parser.add_argument('--inter','-I', help='Interactions file',          type=str)
    parser.add_argument('--seed',       help='Seed used to generate data', type=int,   default=0)
    parser.add_argument('--beta', '-b', help='Inverse temperature ',       type=float, default=1.0)
    parser.add_argument('--mode', '-m', help='Training mode',           choices=modes, default='aquant')
    args = vars(parser.parse_args())

    
    # Load data generating Hamiltonian ----------------------------------------
    beta = args['beta']
    Ns   = 0
    (Hs, Ds, Js, bonds) = (np.array([]), np.array([]), np.array([]), np.array([]))  
    if (args['inter'] is not None):
        Hs, Ds, Js, bonds = Hfile.LoadInters(args['inter'])
        Ns = len(Hs)
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
    np.random.seed(args['seed'])
    #Hs = np.random.randn(Ns)*max([max(Hs),0.5])+Hs
    #Ds = np.abs(np.random.randn(Ns)*max([max(Ds),0.1]))+Ds
    #Js = np.random.randn(Nb)*max([max(Js),0.5])+Js
    
    Hs[:,-1] = (rm.rand(Ns)-0.5)*0.1
    Js[:,-1] = (rm.rand(Nb)-0.5)*0.1
    Ds[:,-1] = np.ones(Ns) * 2.0
    
    # Obtain or generate a training dataset -----------------------------------
    # Vectors in a dataset are assumed to be a series of 0,1 (Python convention)
    data, weights = Hfile.GetData(args['data'], Ns, args['seed'], 10000, bonds, Hs, Ds, Js, beta) 
    data    = data
    weights = weights
    data += [[]]          # add a vector with no bits clamped
    Nd = len(data)
    weights = np.append(weights, -beta) # and its weight
    weights = weights.reshape((Nd,1))  # take a transpose
    #weights = np.transpose(np.atleast_2d(weights))# take a transpose
    
    # Keep a copy for the ED solver
    pdata    = data[:]
    pweights = np.copy(weights)
   
    gval   = np.zeros(Ns+Ns+Nb) # the gradient 
    gpos   = np.zeros(Ns+Ns+Nb) # positive contribution to the gradient
    
    # It can be computed exactly in the classical and inequality quantum trainings
    if args['mode'] == 'aquant':
        dstats = np.zeros((Nd-1, Ns+Ns+Nb))
        for i, cbits in enumerate(data[:-1]):
            dstats[i, :Ns] = -1.0*(np.array(cbits)*2-1)*beta
            for j,bond in enumerate(bonds): 
                dstats[i, 2*Ns+j] = (cbits[bond[0]]*2-1)*(cbits[bond[1]]*2-1)*beta
        gpos = np.sum(dstats*weights[:-1, :], axis=0)

        # Get rid of data we don't need to work with anymore,
        # leave only the identity projector 
        data = [data[-1]]
        weights = np.array(weights[-1, 0]) 
        Nd = 1
        
    # Set up MC constants------------------------------------------------------ 
    mSlice = 1025   # minimum number of bins to take in attempt to converge the autocorrelation
    gtol  = 0.10    # fractional error tolerance in estimation of the gradient
    etol  = 0.08    # tolerance  
    prec  = 0.01    # gradient resolution below which it is considered to be converged 
    

    # Set up gradient descent constants----------------------------------------
    a, b = IPschedule(0.9/beta, 0.1/beta, 200)
    Norm = 100.0 # gradient norm 
    step = 0     # step in gradient descent
    nsteps = 1000

    subset  = range(Nd)
    
    # Identify the indices of trained parameters to nulify in the gradient
    if   args['mode']=='equant': tparams = np.ones(Ns+Ns+Nb)
    elif args['mode']=='aquant': tparams = np.hstack([np.ones(Ns), np.zeros(Ns), np.ones(Nb)])

    print 'gpos: ', gpos
    print 'grad: ', gval
    print 'mask: ', tparams


    # Gradient descent --------------------------------------------------------
    while ((step < nsteps) and (Norm > 0.00001)):
        eta   = a/(b+1.0*step)
        
        print 'Hs:   ', Hs[:,-1]
        print 'Ds:   ', Ds[:,-1]
        print 'Js:   ', Js[:,-1]
    
        # Initiate ED solver --------------------------------------------------
        kwargs = {'X': Ds, 'Z1': Hs, 'Z2': Js}
        BM = bmachine.BoltzmannMachine(Ns, beta, **kwargs)
    
        # Evaluate log-likelihood
        vLL = 0
        for k,cbits in enumerate(pdata[:-1]):
            BM.setProjector(cbits)
            vLL -= np.log(np.real(BM.evaluateProjector())) * pweights[k]
        
        print 'ED LL:   ', vLL 
        BM.setProjector([])
        edAves = BM.computeLocalAverages()*pweights[-1]
        print 'ED grad: ', gpos + edAves
        del BM

        # Initiate MC solvers -------------------------------------------------
        TFIMs = []                      # list of MC objects, one for each data vector
        (Zs,   Xs,  ZZs) = ([], [], []) # lists of raw measurements for each data vector
        for vector in data:
            TFIMs.append(tfim.TFIM(vector, bonds, Hs[:,-1].tolist(), Ds[:,-1].tolist(), Js[:,-1].tolist(), Ns,beta, 1))
            Zs.append([]); Xs.append([]); ZZs.append([])
         
        # Reset necessary datastructures --------------------------------------
        (bins, mtargets) = (np.zeros(Nd), mSlice*np.ones(Nd)) # bins taken and total number of bins to take 
        (aZs, aXs, aZZs) = (np.zeros((Nd, Ns)), np.zeros((Nd, Ns)), np.zeros((Nd, Nb)))  # estimators' averages
        (eZs, eXs, eZZs) = (np.zeros((Nd, Ns)), np.zeros((Nd, Ns)), np.zeros((Nd, Nb)))  # and errors   
        (isAutoCorr, isUncertain) = (True, True) # flags controlling simulation runtime

        print 'bins:    ', bins
        print 'targets: ', mtargets
        # run MC until the error is converged and is below the error tolerance level
        while isAutoCorr or isUncertain:
            
            
            # for each vector in dataset
            print "\n\n------------------------------------------------------------" 
            for i in range(Nd):
                # run a pre-determined number of measurements
                while (bins[i] < mtargets[i]):
                    for sweep in range(100):
                        # diagonal update
                        while  TFIMs[i].DMove() == 1:  # if the operator list is too small
                               TFIMs[i].Adjust()       # adjust its length 
                               bins[i] = 0             # set bin count to zero
                               (Zs[i], Xs[i], ZZs[i]) = ([], [], []) # reset timeseries
                               (bins, mtargets) = (np.zeros(Nd), mSlice*np.ones(Nd))
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
                mcAves = np.sum((aves*weights),  axis=0)*tparams
                gval = gpos + mcAves
                gerr = np.sqrt(np.sum((errs*errs*weights*weights), axis=0))*tparams
                print "eta:      ", eta
                print "LL:       ", vLL
                print "gradient: ", gval
                print "exact:    ", (gpos+edAves)*tparams
                #print "MC aves:  ", mcAves
                #print 'ED aves:  ', edAves 
                print "gerror:   ", gerr

                # determine the indices of averages falling below the tolerance treshold 
                #indices = ma.getmask(ma.masked_where(gval*gtol-gerr < 0, ))
                #indices = np.squeeze(np.where(np.abs(gval)*gtol-gerr < 0))
                indices  = np.squeeze(np.where(np.all([
                                                        np.abs(gval)*gtol-gerr < 0,
                                                        np.abs(gval)     -prec > 0
                                                      ], axis=0)
                                               )
                                     )  
                print "indices:   ", indices
                if  indices.shape == (): indices = np.array([indices])

                # if there are no such averages, signal that the errors are within gtolerance bounds
                if (len(indices) == 0) or (np.amax(gerr)<etol): isUncertain = False
                # otherwise, estimate how many more bins to take
                else: mtargets = detTotalBins(errs[:, indices], weights, np.abs(gval)[indices]*gtol/gerr[indices], bins)
                
        # Destroy MC solvers -----------------------------------------------------
        for j in reversed(range(len(TFIMs))):
            del TFIMs[j]
        
        # Annealing constant 

        # Follow the negative gradient 
        (dH, dD, dJ) = np.split(gval, [Ns, Ns+Ns])    
        Hs[:,-1] += -eta*dH
        Ds[:,-1] += -eta*dD
        Js[:,-1] += -eta*dJ
        step += 1

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
