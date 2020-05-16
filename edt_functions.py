# EDT functions
# Written by José Ignacio Hernández
# May 2020

# Changelog:
# v0.1: Rewrite all the code from zero

# Load packages
import pandas as pd
import numpy as np
import time
import datetime
import itertools as itt
import re
from scipy.stats import chi2_contingency

# EffDesign main routine
def effdesign(ATTLIST,NALT,NCS,OPTS):

    # Set defaults in case some options were not defined
    OPTS = effdesign_options(OPTS)

    # Check inputs integrity
    chkmess,chkcode = checkatts(ATTLIST,NALT,NCS,OPTS)

    # If some error is detected, stop
    if chkcode == 1:
        # print(chkmess)
        raise ValueError(chkmess)

    # Define scalars
    NATT = len(ATTLIST)
    NRUNS = NALT*NCS

    # Create arrays that will be used from now on
    OPTOUT = OPTS['optout']
    ASC = OPTS['asc']
    ASCPAR = OPTS['ascpar']
    SEED = OPTS['seed']
    ITERLIM = OPTS['iterlim']
    NOIMPROVLIM = OPTS['noimprovlim']
    TIMELIM = OPTS['timelim']
    NBLOCKS = OPTS['nblocks']
    COND = OPTS['cond']

    NAMES = []
    LEVS = []
    CODS = []
    PAR = []

    for k in range(0,NATT):
        NAMES = NAMES + [ATTLIST[k]['name']]
        LEVS = LEVS + [len(ATTLIST[k]['levels'])]
        CODS = CODS + [ATTLIST[k]['coding']]
        PAR = PAR + ATTLIST[k]['par']

    # Add the ASC parameter if present
    if len(ASC) > 0:
        PAR = PAR + ASCPAR


    # Set random seed if defined
    if SEED >= 0:
        np.random.seed(SEED)

    ############################################################
    ########## Step 1: Generate initial design matrix ##########
    ############################################################

    # Generate 10 random designs
    # Keep the one with lowest D-error as initial design
    print('Generating the initial design matrix and D-error')
    desmat = []
    initd = np.inf

    for _ in range(0,10):
        desmat0 = initdesign(ATTLIST,LEVS,NATT,NRUNS,NAMES,COND)
        derr0 = imat_ubalance(desmat0,CODS,PAR,NATT,NALT,NCS,NRUNS,OPTOUT,ASC)
        if derr0 < initd:
            desmat = desmat0.copy()
            initd = derr0

    ############################################################
    ############## Step 2: Initialize algorighm ################
    ############################################################

    # Execute Swapping algorithm
    bestdes, bestderr, best_t, elapsed_time = swapalg(desmat,initd,NATT,NALT,NCS,NRUNS,CODS,PAR,OPTOUT,ASC,NAMES,COND,ITERLIM,NOIMPROVLIM,TIMELIM)

    # Compute utility balance ratio
    ub = imat_ubalance(desmat0,CODS,PAR,NATT,NALT,NCS,NRUNS,OPTOUT,ASC,ubalance=True)

    ############################################################
    ############## Step 3: Arange final design #################
    ############################################################

    # If opt-out option is present, add to design
    if OPTOUT == 1:
        bestdes.shape = (NCS,NALT,NATT)
        bestdes = np.concatenate((bestdes,np.zeros((NCS,1,NATT))),axis=1)
        bestdes.shape = (NRUNS+NCS,NATT)
        
        grprow = np.repeat(np.arange(1,NCS+1),NALT+1)
        grprow.shape = (NRUNS+NCS,1)
        altrow = np.tile(np.arange(1,NALT+1+1),NCS)
        altrow.shape = (NRUNS+NCS,1)

    else:
        grprow = np.repeat(np.arange(1,NCS),NALT)
        grprow.shape = (NRUNS,1)
        altrow = np.tile(np.arange(1,NALT+1),NCS)
        altrow.shape = (NRUNS,1)

    bestdes = np.concatenate((grprow,altrow,bestdes),1).copy()

    # Generate blocks
    if NBLOCKS > 0:
        print('Generating ' + str(NBLOCKS) + ' blocks...')
        blocksrow = blockgen(bestdes,NBLOCKS,NCS,1000)
        bestdes = np.concatenate((bestdes,blocksrow),1).copy()

    # Create Pandas DataFrame
    if NBLOCKS > 0:
        exportdes = pd.DataFrame(bestdes,columns=['CS','Alt'] + NAMES + ['Block'],dtype='int64')
    else:
        exportdes = pd.DataFrame(bestdes,columns=['CS','Alt'] + NAMES,dtype='int64')

    # Return a dictionary
    di = {'final.it': best_t, 'init.derr': initd, 'final.derr': bestderr, 'balance.ratio': ub, 'optimal.des': exportdes}

    print('Optimization complete')
    print('Elapsed time: ' + str(datetime.timedelta(seconds=elapsed_time))[:7])
    print('D-error of initial design: ',round(initd,6))
    print('D-error of last stored design: ',round(bestderr,6))
    print('Utility Balance ratio: ',round(ub,2),'%')
    print('Algorithm iterations: ',best_t)
    print('')
    return(di)

# Swapping algorithm function
def swapalg(DES,INITD,NATT,NALT,NCS,NRUNS,CODS,PAR,OPTOUT,ASC,NAMES,COND,ITERLIM,NOIMPROVLIM,TIMELIM):
    
    # Lock design matrix
    desmat = DES.copy()

    # If conditionals are set, create list of feasibility
    if len(COND) > 0:
        condlist = []
        # Generate condition list
        for i in range(0,len(COND)):
            condlist = condlist + condgen(NAMES,COND[i],init=False)

    # Generate matrix of all possible permutations
    combmat = np.array(list(itt.combinations(range(1,NRUNS+1), 2))) - 1

    # Start stopwatch (falta)
    a = time.time()
    b = time.time()

    # Initialize algorithm parameters
    i = 0
    t = 0
    ni = 0
    iterd = INITD

    # Start algorithm
    while True:
        
        # Iteration No.
        t = t+1
        
        # If one stopping criterion is satisfied, break!
        if ni >= NOIMPROVLIM or t >= ITERLIM or (b-a)/60 >= TIMELIM:
            break
        
        # Take a random swap
        pairswap = combmat[np.random.choice(combmat.shape[0],1)][0]
        
        # If attribute levels differ, do the swap.
        if desmat[pairswap[0],i] != desmat[pairswap[1],i]:
            swapdes = desmat.copy()
            swapdes[pairswap[0],i] = desmat[pairswap[1],i]
            swapdes[pairswap[1],i] = desmat[pairswap[0],i]
            
            # If conditions are defined, check if the swap satisfies them
            if len(COND) > 0:
            
                # Generate check vector
                check = []
                for co in range(0,len(condlist)):
                    check = check + [np.any(np.logical_and(eval(condlist[co][0]),np.logical_not(eval(condlist[co][1]))))]
                
                # If all conditions are satisfied, compute D-error
                if not(np.any(check)):
                    newd = imat_ubalance(swapdes,CODS,PAR,NATT,NALT,NCS,NRUNS,OPTOUT,ASC)
                    
                #...else, keep original D-error
                else:
                    newd = iterd.copy()
            
            #...else if conditions are not defined, pass directly to D-error computation
            else:
                newd = imat_ubalance(swapdes,CODS,PAR,NATT,NALT,NCS,NRUNS,OPTOUT,ASC)

        # ...else if they do not differ, keep the D-error
        else:
            newd = iterd.copy()
            
        # If the swap made an improvement, keep the design.
        # Else, pass to the next attribute
        if newd < iterd:
            desmat = swapdes.copy()
            iterd = newd.copy()
            ni = 0
            
            # Add this lines for GUI version (update text bars)
            # self.line_DERR.setText(str(round(iterd,6)))
            # self.line_ITER.setText(str(t))
            # self.line_NOIMPROV.setText(str(0))
            # self.line_ELAPSED.setText(str(datetime.timedelta(seconds=b-a))[:7])
            
            # Update progress bar
            print('Optimizing. Press ESC to stop. / ' + 'Elapsed: ' + str(datetime.timedelta(seconds=b-a))[:7] + ' / D-Error: ' + str(round(iterd,6)),end='\r')
            
        else:
            i = i+1
            ni = ni+1
        
        # If the algorithm reach the attribute limit, reset to attribute 1
        if i > NATT-1:
            i = 0
        
        # Update progress bar each second
        if (b-a)%1 == 0:
            print('Optimizing. Press ESC to stop. / ' + 'Elapsed: ' + str(datetime.timedelta(seconds=b-a))[:7] + ' / D-Error: ' + str(round(iterd,6)),end='\r')
            # Add this lines for GUI version (update text bars)
            # self.line_DERR.setText(str(round(iterd,6)))
            # self.line_ITER.setText(str(t))
            # self.line_NOIMPROV.setText(str(0))
            # self.line_ELAPSED.setText(str(datetime.timedelta(seconds=b-a))[:7])
                
        # Cancel if ESC is pressed...
        # if msvcrt.kbhit():
        #     if ord(msvcrt.getch()) == 27:
        #         break
        
        b = time.time()
        difftime = b - a
        # In GUI version, this part prevents window freezing
        # QtWidgets.qApp.processEvents()
    
    # Return optimal design plus D-error
    print('\n')
    return(desmat,iterd,t,difftime)

# Block generation function
def blockgen(DES,NBLOCKS,NCS,REPS):
    
    blocks = np.repeat(np.arange(1,NBLOCKS+1),NCS/NBLOCKS).copy()
    blocks.shape = (NCS,1)
    bestcorr = np.inf
    bestblock = blocks.copy()
    
    for _ in range(0,REPS):
        np.random.shuffle(blocks)
        blockmat = np.repeat(blocks,int(np.max(DES[:,1]))).copy()
        blockmat.shape = (blockmat.shape[0],1)
        sumcorr = 0
        
        for a in range(2,DES.shape[1]):
            d = DES[:,a].copy()
            d.shape = (d.shape[0],1)
            c = cross(blockmat,d)
            corr = chi2_contingency(c)[1]
            sumcorr = sumcorr + corr
        
        if sumcorr < bestcorr:
            bestblock = blockmat.copy()
            bestcorr = sumcorr
        
    bestblock.shape = (bestblock.shape[0],1)

    return(bestblock)

# Function for dummy generation
def dummygen(var):

    NLEVS = len(np.unique(var))

    dm = []
    if NLEVS > 1:
        for l in range(1,NLEVS):
            dm = np.hstack((dm,(var == l)*1))

    return(dm)

# Information matrix and utility balance function
def imat_ubalance(DES,CODS,PAR,NATT,NALT,NCS,NRUNS,OPTOUT,ASC,ubalance=False):

    # Initialize scalars
    NPAR = len(PAR)
    DES = DES.copy()

    # If opt-out option is present, add to design
    if OPTOUT == 1:
        DES.shape = (NCS,NALT,NATT)
        NALT = NALT+1
        NRUNS = NRUNS+NCS
        DES = np.concatenate((DES,np.zeros((NCS,1,NATT))),axis=1)
        DES.shape = (NRUNS,NATT)

    # Populate the estimable design matrix
    estdes = []
    
    for k in range(0,NATT):
        if CODS[k] > 1:
            dum = dummygen(DES[:,k])
            estdes = np.hstack((estdes,dum))
        else:
            estdes = np.hstack((estdes,DES[:,k]))

    estdes.shape = ((NPAR-len(ASC)),NRUNS)

    # Add alternative-specific constants
    if len(ASC)> 0:
        altrow = np.tile(np.arange(1,NALT+1),NCS)
        
        for i in range(0,len(ASC)):
            aa = (altrow == ASC[i])*1
            estdes = np.vstack((estdes,aa))
        
    estdes = estdes.T

    # Calculate Probability
    v = estdes.dot(PAR)
    ev = np.exp(v)
    ev = np.reshape(ev,(NALT,NCS),order = 'F').copy()
    sev = np.sum(ev,0)
    p = ev/sev

    # Calculate utility balance if required
    if ubalance:
        B = (p/(1/NALT)).copy()
        B = np.prod(B,axis=0)*100
        B = np.sum(B,axis=0)/NCS

        return(B)

    else:
        p = np.reshape(p,(NRUNS,1),order = 'F').copy()
        
        # Calculate Information Matrix
        ia = np.diag(p.flat).dot(estdes)
        ia = ia.T.dot(estdes)
        
        ib = np.diag(p.flat).dot(estdes)
        ib = np.reshape(ib,(NALT,NCS,NPAR),order = 'F').copy()
        ib = np.sum(ib,0)
        ib = ib.T.dot(ib)
        
        im = ia - ib
        
        # Calculate D-error
        if np.linalg.det(im) != 0:
            vce = np.linalg.solve(im,np.eye(im.shape[0]))

            if len(ASC) > 0:
                vce = vce[:vce.shape[0]-1,:vce.shape[0]-1].copy()

            detvce = np.linalg.det(vce)
            dr = detvce**(1/vce.shape[0])

        else:
            dr = np.inf

        return(dr)

# Generate initial design matrix
def initdesign(ATTLIST,LEVS,NATT,NRUNS,NAMES,COND):

    desmat = np.zeros((NRUNS,NATT))

    # Populate the initial design matrix
    for k in range(0,NATT):
        dd = np.repeat(ATTLIST[k]['levels'],NRUNS/LEVS[k])
        np.random.shuffle(dd)
        desmat[:,k] = dd

    # If there exists some condition, apply to initial design
    if len(COND) > 0:
        condlist = []
        
        # Generate condition list
        for i in range(0,len(COND)):
            condlist = condlist + [condgen(NAMES,COND[i])]
        
        # Apply conditions
        for k in range(0,len(COND)):
            for p in range(0,len(condlist[k])):
                samplevec = np.array(ATTLIST[condlist[k][p][1]]['levels'])
                condition = 'samplevec' + str(condlist[k][p][2]) + str(condlist[k][p][3])
                samplevec = samplevec[eval(condition)].copy()
                
                for i in range(0,NRUNS):
                    if eval(condlist[k][p][0][0]) and not(eval(condlist[k][p][0][1])):
                        desmat[i,int(condlist[k][p][1])] = np.random.choice(samplevec,1)

    return(desmat)

# Function to set defaults in case some options are not defined
def effdesign_options(OPTS):
    
    # Take the original options dictionary
    complete_opts = OPTS.copy()
    
    # Add default values if missing
    complete_opts.setdefault('optout',0)
    complete_opts.setdefault('asc',[])

    if len(complete_opts['asc']) == 0:
        complete_opts['ascpar'] = []

    complete_opts.setdefault('seed',-1)
    complete_opts.setdefault('iterlim',np.inf)
    complete_opts.setdefault('noimprovlim',np.inf)
    complete_opts.setdefault('timelim',np.inf)
    complete_opts.setdefault('nblocks',0)
    complete_opts.setdefault('cond',[])

    # Return new dictionary
    return(complete_opts)

# Integrity check function
def checkatts(ATTLIST,NALT,NCS,OPTS):
    NATT = len(ATTLIST)
    NRUNS = NALT*NCS

    code = 0
    mess = "All OK"

    # Start check loop among attributes
    for k in range(0,NATT):

        # Check if the element k is a dictionary
        if type(ATTLIST[k]) is not(dict):
            mess = "Error: attribute " + str(k+1) + " is not a dictionary."
            code = 1
            return(mess,code)

        # Check if all elements are present
        for e in ['name','levels','coding','par']:
            if e not in ATTLIST[k]:
                mess = "Error: " + e + " is not defined in attribute " + str(k+1)
                code = 1
                return(mess,code)

        # Check if 'names' is a string
        if type(ATTLIST[k]['name']) is not(str):
            mess = "Error: name of attribute " + str(k+1) + " is not a string."
            code = 1
            return(mess,code)

        # Check if 'levels' and 'par' are lists
        for e in ['levels','par']:
            if type(ATTLIST[k][e]) is not(list):
                mess = "Error: " + e + " of attribute " + str(k+1) + " is not a list."
                code = 1
                return(mess,code)

        # Check if 'coding' is an integer
        if type(ATTLIST[k]['coding']) is not(int):
            mess = "Error: coding of attribute " + str(k+1) + " is not an integer."
            code = 1
            return(mess,code)

        # Check if NRUNS is divisible by number of attribute levels
        if NRUNS%len(ATTLIST[k]['levels']) != 0:
            mess = "Error: No. of Choice sets times No. of alternatives is not divisible by number of levels of attribute " + str(k+1) + "."
            code = 1
            return(mess,code)


        # Check if number of prior parameters are well-defined (one less than number of levels)
        if ATTLIST[k]['coding'] != 1:
            if len(ATTLIST[k]['par']) != len(ATTLIST[k]['levels']) - 1:
                mess = "Error: number of prior parameters of attribute " + str(k) + " must be one less than its corresponding number of attribute levels."
                code = 1
                return(mess,code)
        
        else:
            if len(ATTLIST[k]['par']) > 1:
                mess = "Error: number of prior parameters of attribute " + str(k) + " must be of length one."
                code = 1
                return(mess,code)

    # Check if there are enough alternatives defined.
    if (NALT < 2 and OPTS['optout']==0):
        mess = "Error: at least two alternatives (including opt-out) must be specified."
        code = 1
        return(mess,code)

    # Check if enough choice sets are defined.
    if NCS <2:
        mess = "Error: at least two choice sets must be specified."
        code = 1
        return(mess,code)

    # Check if ASC and its priors are well-defined
    if len(OPTS['asc']) > 0:
        if OPTS['optout'] == 1:
            if len(OPTS['asc']) - 1 >= NALT:
                mess = "Number of ASC must be less than number of alternatives."
                code = 1
                return(mess,code)
            
            if max(OPTS['asc']) > NALT + 1 or min(OPTS['asc']) <= 0:
                mess = "An ASC is defined outside the number of alternatives."
                code = 1
                return(mess,code)

        else:
            if len(OPTS['asc']) >= NALT:
                mess = "Number of ASC must be less than number of alternatives."
                code = 1
                return(mess,code)
            
            if max(OPTS['asc']) > NALT or min(OPTS['asc']) <= 0:
                mess = "An ASC is defined outside the number of alternatives."
                code = 1
                return(mess,code)

        if len(OPTS['asc']) != len(OPTS['ascpar']):
            mess = "Number of ASC prior parameters must be equal to number of ASC."
            code = 1
            return(mess,code)

    return(mess,code)

# Crosstab function
def cross(x,y):
    tab = []
    
    for i in np.unique(x):
        cols = []
        
        for j in np.unique(y):
            c = np.count_nonzero((x==i) & (y==j))
            cols = cols + [c]
        
        tab = tab + [cols]

    return(np.array(tab))

# Condition generation function
def condgen(NAMES,cc,init=True):
    
    # Split if > then conditions from the comma expression
    a0 = cc.split(',')
    b0 = a0[0]
    c0 = a0[1]
    
    # Match names with attribute columns
    a1 = np.where(np.isin(NAMES,re.sub('[^A-Za-z]+','',b0)))[0]
    b1 = re.sub('[^<>=]+','',b0)
    c1 = re.sub('[^0-9]+','',b0)
    
    c0 = c0.split('and')
    
    for i in range(0,len(c0)):
        locals()['a' + str(i+2)] = np.where(np.isin(NAMES,re.sub('[^A-Za-z]+','',c0[i])))[0]
        locals()['b' + str(i+2)] = re.sub('[^<>=]+','',c0[i])
        locals()['c' + str(i+2)] = re.sub('[^0-9]+','',c0[i])

    if init:
        clist = []
        
        for i in range(0,len(c0)):
            clist = clist + [	[	[	'desmat[i' + ',' + str(a1[0]) + ']' + str(b1) + str(c1) ,
                                        'desmat[i' + ',' + str(locals()['a' + str(i+2)][0]) + ']' + str(locals()['b' + str(i+2)]) + str(locals()['c' + str(i+2)])
                                    ],
                                    locals()['a' + str(i+2)][0],
                                    locals()['b' + str(i+2)],
                                    int(locals()['c' + str(i+2)])
                                ]
                            ]
        
    else:
        clist = []
        
        for i in range(0,len(c0)):
            clist = clist + [	[	'swapdes[:' + ',' + str(a1) + ']' + str(b1) + str(c1),
                                    'swapdes[:' + ',' + str(locals()['a' + str(i+2)]) + ']' + str(locals()['b' + str(i+2)]) + str(locals()['c' + str(i+2)])
                                ]
                            ]
    
    return(clist)