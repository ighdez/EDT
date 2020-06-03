# EDT functions
# Written by José Ignacio Hernández
# June 2020

# Changelog:
# -----------------------------------------------
# v0.2: - Several changes in condition generation
#       - Bug fixes and code improvements
# -----------------------------------------------
# v0.1: - Rewrite all the code from zero
# -----------------------------------------------

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

    target_atts = np.arange(0,len(NAMES))

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
    bestdes, bestderr, best_t, elapsed_time = swapalg(desmat,initd,target_atts,NATT,NALT,NCS,NRUNS,CODS,PAR,OPTOUT,ASC,NAMES,COND,ITERLIM,NOIMPROVLIM,TIMELIM)

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
def swapalg(DES,INITD,TARGET_ATTS,NATT,NALT,NCS,NRUNS,CODS,PAR,OPTOUT,ASC,NAMES,COND,ITERLIM,NOIMPROVLIM,TIMELIM):
    
    # Lock design matrix
    desmat = DES.copy()

    # If conditionals are set, create list of feasibility
    if len(COND) > 0:
        condlist = []
        # Generate condition list
        for i in range(0,len(COND)):
            condlist = condlist + [condgen(NAMES,"swapdes",COND[i],init=False)]

    # Generate matrix of all possible permutations
    combmat = np.array(list(itt.combinations(range(1,NRUNS+1), 2))) - 1

    # Start stopwatch
    t0 = time.time()
    t1 = time.time()

    difftime = 0

    # Initialize algorithm parameters
    i = np.random.choice(TARGET_ATTS,1)[0]
    t = 0
    ni = 0
    improv = 0
    iterd = INITD
    newd = INITD

    # Start algorithm
    while True:
        
        # Iteration No.
        t = t+1
        
        # If one stopping criterion is satisfied, break!
        if ni >= NOIMPROVLIM or t >= ITERLIM or (difftime)/60 >= TIMELIM:
            break
        
        # Take a random swap
        pairswap = combmat[np.random.choice(len(combmat),1)][0]
        
        # Check if attribute levels differ
        check_difflevels = desmat[pairswap[0],i] != desmat[pairswap[1],i]

        # If attribute levels differ, do the swap and check for conditions (if defined)
        if check_difflevels:
            swapdes = desmat.copy()
            swapdes[pairswap[0],i] = desmat[pairswap[1],i]
            swapdes[pairswap[1],i] = desmat[pairswap[0],i]
            
            # Check if conditions are satisfied after a swap
            check_satisfied_conds = True  
            
            # If conditions are defined, this section will check that are satisfied, and rewrite 'check_satisfied_conds' if neccesary
            if len(COND) > 0:
            
                # Generate check vector
                condcheck = []

                # Loop among conditions
                for j in range(0,len(condlist)):
                    for c in range(0,len(condlist[j][1])):
                        test = np.logical_or(np.logical_not(eval(condlist[j][0])),eval(condlist[j][1][c]))
                        condcheck = condcheck + [np.all(test)]
                
                # Rewrite 'check_satisfied_conds'
                check_satisfied_conds = np.all(condcheck)

            # If all conditions are satisfied, compute D-error
            if check_satisfied_conds:
                newd = imat_ubalance(swapdes,CODS,PAR,NATT,NALT,NCS,NRUNS,OPTOUT,ASC)

        # ...else if they do not differ, keep the D-error
        else:
            newd = iterd.copy()
            
        # If the swap made an improvement, keep the design and update progress bar
        if newd < iterd:
            desmat = swapdes.copy()
            iterd = newd.copy()
            ni = 0
            improv = improv + 1
            
            # Add this lines for GUI version (update text bars)
            # self.line_DERR.setText(str(round(iterd,6)))
            # self.line_ITER.setText(str(t))
            # self.line_NOIMPROV.setText(str(0))
            # self.line_ELAPSED.setText(str(datetime.timedelta(seconds=b-a))[:7])
            
            # Update progress bar
            print('Optimizing. Press ESC to stop. / ' + 'Elapsed: ' + str(datetime.timedelta(seconds=difftime))[:7] + '/ Improv.: ' + str(improv) + ' / D-Error: ' + str(round(iterd,6)),end='\r')

        # ...else, pass to a random attribute and increment the 'no improvement' counter by 1.
        else:
            i = np.random.choice(TARGET_ATTS,1)[0]
            ni = ni+1
        
        # Update progress bar each second
        if (difftime)%1 < 0.1:
            print('Optimizing. Press ESC to stop. / ' + 'Elapsed: ' + str(datetime.timedelta(seconds=difftime))[:7] + '/ Improv.: ' + str(improv) + ' / D-Error: ' + str(round(iterd,6)),end='\r')

            # Add this lines for GUI version (update text bars)
            # self.line_DERR.setText(str(round(iterd,6)))
            # self.line_ITER.setText(str(t))
            # self.line_NOIMPROV.setText(str(0))
            # self.line_ELAPSED.setText(str(datetime.timedelta(seconds=b-a))[:7])
                
        # Cancel if ESC is pressed...
        # if msvcrt.kbhit():
        #     if ord(msvcrt.getch()) == 27:
        #         break
        
        t1 = time.time()
        difftime = t1-t0
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
            condlist = condlist + [condgen(NAMES,"desmat",COND[i])]

        # Apply conditions
        for j in range(0,len(condlist)):
            for i in range(0,len(desmat)):

                # Create test object
                test = not(eval(condlist[j][0])) or eval(condlist[j][1])
                
                # If the statement is true, then replace the value that doesn't satisfy the statement
                if not(test):
                    
                    # For multiple "then" statements
                    if len(condlist[j][2])>1:

                        # Create sample vector
                        samplevec = []

                        for k in range(0,len(condlist[j][2])):
                            samp = np.array(ATTLIST[condlist[j][2][k]]['levels'])
                            condition = 'samp' + str(condlist[j][3][k]) + str(condlist[j][4][k])
                            samplevec = samplevec + [samp[eval(condition)]]

                        # Replace values until the statement is true
                        while not(test):
                            for k in range(0,len(samplevec)):
                                desmat[i,int(condlist[j][2][k])] = np.random.choice(samplevec[k],1)
                            
                            test = not(eval(condlist[j][0])) or eval(condlist[j][1])
                    
                    # For single statement
                    else: 
                        samp = np.array(ATTLIST[condlist[j][2][0]]['levels'])
                        condition = 'samp' + str(condlist[j][3][0]) + str(condlist[j][4][0])
                        samplevec = samp[eval(condition)]

                        while not(test):
                            desmat[i,int(condlist[j][2][0])] = np.random.choice(samplevec,1)                                
                            test = not(eval(condlist[j][0])) or eval(condlist[j][1])

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
def condgen(NAMES,d,cc,init=True):
    
    delimiters = '<=|>=|==|<|>'

    # Split conditions
    if_part = cc.split(',')[0].replace('if ','')
    if_pos = re.split(delimiters,if_part)[0].replace(' ','')
    if_pos = np.where(np.isin(NAMES,if_pos))[0][0]
    if_sign = ''.join(re.findall(delimiters,if_part))
    if_value = int(re.split(delimiters,if_part)[1].replace(' ',''))

    then_part = cc.split(',')[1].split('and')

    if init:
        # Create part 1 of condition expression
        part1 = str(d) + '[i,' + str(if_pos) + ']' + if_sign + str(if_value)

        # Create part 2 of condition expression
        then_att = re.split(delimiters,then_part[0])[0].replace(' ','')
        then_pos = [np.where(np.isin(NAMES,then_att))[0][0]]
        then_sign = [''.join(re.findall(delimiters,then_part[0]))]
        then_value = [int(re.split(delimiters,then_part[0])[1].replace(' ',''))]

        part2 = str(d) + '[i,' + str(then_pos[0]) + ']' + then_sign[0] + str(then_value[0])

        # If the second part of the "then" expression contains more than 1 condition, then join
        if len(then_part) > 1:
            for i in range(1,len(then_part)):
                then_att = re.split(delimiters,then_part[i])[0].replace(' ','')
                then_pos = then_pos + [np.where(np.isin(NAMES,then_att))[0][0]]
                then_sign = then_sign + [''.join(re.findall(delimiters,then_part[i]))]
                then_value = then_value + [int(re.split(delimiters,then_part[i])[1].replace(' ',''))]

                part2 = part2 + ' and ' + str(d) + '[i,' + str(then_pos[i]) + ']' + then_sign[i] + str(then_value[i])

        clist = [part1,part2,then_pos,then_sign,then_value]

    else:
        # Create part 1 of condition expression
        part1 = str(d) + '[:,' + str(if_pos) + ']' + if_sign + str(if_value)

        # Create part 2 of condition expression
        part2 = []
        for i in range(0,len(then_part)):
            then_att = re.split(delimiters,then_part[i])[0].replace(' ','')
            then_pos = np.where(np.isin(NAMES,then_att))[0][0]
            then_sign = ''.join(re.findall(delimiters,then_part[i]))
            then_value = int(re.split(delimiters,then_part[i])[1].replace(' ',''))

            part2 = part2 + [str(d) + '[:,' + str(then_pos) + ']' + then_sign + str(then_value)]
        clist = [part1,part2]

    return(clist)