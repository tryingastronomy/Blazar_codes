#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 22:31:35 2022

@author: joffresd
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 16:21:08 2022

@author: joffresd

Do same streamlined process but now instead of 4FGL Gal. sources do it for the FAVA detections
"""

from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astroquery.ipac.ned import Ned
from astroquery.ipac.irsa import Irsa


import astropy.coordinates as coord
import astropy.units as u

from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import copy

#Table creation
from astropy.table import Table, Column

#Polygon
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

import pandas as pd

import matplotlib.path as mpltPath

from astropy.table import vstack

import bisect

#Note: df is the initial flaring sources, v_result is the result of pulling AW srcs from Vizier


#DEFINE FUNCTION FOR FUTURE USE
def pointing (x,y,x1,y1,x2,y2):
    n_ = (x-x1)*(y2-y1) + (y-y1)*(-x2+x1)
    return n_
### Lot of files used for testing ###
df = pd.read_csv() # 4FGL blazars
                    # 4FGL Galactic sources
                    #Novae
                    #unclassified 4FGL

########################################

#insert file from FAVA here
df = pd.read_csv('')

#Note: if want to check 4FGL sources uncomment code for ra, dec, and search_radius

all_3_index = []
all_3_ul = []
blz_assoc_list = []

#NOTE: can adjust ra, dec, and search_radius definitions for whatever data file you are looking at with
for k in range(len(df['flareID'])): #len(df['flareID'])
    #k = 21
    #k=0
    #First do just as one and then convert to loop through all of the sources
    ra = df['best_ra'][k] #from 4FGL #df['RAJ2000'][k]
    dec =  df['best_dec'][k] #from 4FGL #df['DEJ2000][k]
    search_radius = df['best_r95'][k] #from 4FGL #df['Conf_95_SemiMajor'][k]
    Vizier.ROW_LIMIT= -1  #feel free to change, especially if pulling a lot fewer
    #when go to 10k have these troublesome ones in final [0, 1, 2, 4, 8, 10, 11, 16, 19, 21]
    # List of catalogs to search. If empty, all catalogs will be checked.
    cat_list = ['AllWISE']
    
    vizier_table_list = Vizier.query_region(coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg),
                                                frame='icrs' ), radius=search_radius*u.deg, catalog=cat_list)
    
    
    #instead of removing sources, try putting in exception with continue command
    try:
        v_result1 = copy.deepcopy(vizier_table_list[0])
    except IndexError:
        continue
        
    
    
    
    #Need to be able to detect upper limits
    #This can be done by checking the e_mag field, if nan then it will be an UL - see later
    
    
    
    
    # RA and Dec column names will need to be looked up in the table
    ra_col_name = 'RAJ2000'
    dec_col_name = 'DEJ2000'
    
    c1 = SkyCoord(v_result1[ra_col_name], v_result1[dec_col_name], unit=(u.hourangle, u.deg), frame='icrs')
    c2 = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    sep = c1.separation(c2)
    
    # Add coordinate separation to table
    v_result1['sep'] = sep.arcminute
    
    proper_sources = v_result1['sep'] < 60*search_radius
    
    v_result = v_result1[proper_sources]
    
    # get nearest source
    #oh shoot, when I test max, I am getting wayyyyy too many sources...
    v_result[v_result['sep'] == v_result['sep'].min()]#['texto']
    
    #delete all that have sep>search radius
    
    
    
    #See how many pulled
    print(Vizier.ROW_LIMIT)
    #See how many are actually here
    actual_row_number = len(v_result['e_W3mag'])
    print("actual number of AW sources pulled: " + str(actual_row_number))
    
    v_result['e_W3mag'].value.mask
    
    '''
    Need to adapt below so that we don't just get rid of the rows where there are upper limits
        Results with deletion was -> 19 non blazar srcs (need to investigate how blazars slipped in, 3FGL vs 4FGL uncert. regions )
    
    We need to keep upper limits depending on the color combination and where it lies
    
    eg. W1-W2 vs W2-W3: if W2 is an upper limit, W1W2 -> a *lower* limit. So if the point is to the
    left of the blazar strip, it will still be counted in. CAN'T  JUST GET RID OF ALL
    
    Clarification, when I say UPPER or LOWER I am saying the total value/magnitude - if W1>W2 and W1 is an upper limit (aka smallest value possible so can be bigger) that means
    that W1-W2 can also get bigger, therefore, this is a Lower limit of magnitude, since it is bigger. 
    
    KEY: LOWER LIMIT = SMALLEST VALUE of magnitude
         UPPER LIMIT = LARGEST VALUE of magnitude
    
    If W1 is UL, then LOWER limit on W1-W2  (smallest value, can get bigger)
    
    If W2 is UL for W1-W2 then UPPER limit (largest value calculatd, can get smaller) - 
    
    If W2 is UL for W2-W3 then LOWER limit (smallest mag value calculated, can get bigger)
    
    If W3 is UL, for W3-W4 then LOWER limit (smallest mag value calculated, can get bigger)
    
    If W3 is UL then W2W3 is UPPER limit (largest magnitude of mag value calculated, can get smaller)
    
    If W4 is UL, then UPPER limit (largest magnitude of mag value calculated, can get smaller)
    
    '''
    
    #see full table
    #v_result
    #still need to get rid of rows that have values that are bad. 
    #oh we don't even need the hole turn into -99 thing. Can remove with indices
    #But will update index values each time, so need to define these after each interation - I hate coding
    
    #Here we go, this is what matters
    #Find and remove, need to do each time since index will change
    '''
    Instead of removing, create flag system, W1, W2, W3, W4 
    0-> no limit
    1 -> upper limit
    4 columns, one for each color
    
    lol will also have to demask it because pyfits is stupid 
    
    If it is an upper limit, then it does not have an error associated with it... (just UL)
    '''
    

    
    #v_result.add_column(v_result['e_W1mag'] + v_result['e_W2mag'], name='w1w2_err')
    v_result.add_column(0,name='W1_lim')
    v_result.add_column(0,name='W2_lim')
    v_result.add_column(0,name='W3_lim')
    v_result.add_column(0,name='W4_lim')
    
    w1_maskarr = np.where(v_result['W1mag'].value.mask)[0] #-> you get the mask in  in the uncert. 
    #now replace the column values at the above index with 1 for UL
    v_result['W1_lim'][w1_maskarr] = 1
    

    ew1_maskarr = np.where(v_result['e_W1mag'].value.mask)[0]
    v_result['W1_lim'][ew1_maskarr] = 1
    
    #v_result.remove_rows(ew1_maskarr)
    #print(len(v_result['AllWISE']))
    
    w2_maskarr = np.where(v_result['W2mag'].value.mask)[0]
    v_result['W2_lim'][w2_maskarr] = 1
    
    ew2_maskarr = np.where(v_result['e_W2mag'].value.mask)[0]
    #v_result.remove_rows(ew2_maskarr)
    #print(len(v_result['AllWISE']))
    v_result['W2_lim'][ew2_maskarr] = 1
    
    w3_maskarr = np.where(v_result['W3mag'].value.mask)[0]
    v_result['W3_lim'][w3_maskarr] = 1
    
    ew3_maskarr = np.where(v_result['e_W3mag'].value.mask)[0]
    #v_result.remove_rows(ew3_maskarr)
    #print(len(v_result['AllWISE']))
    v_result['W3_lim'][ew3_maskarr] = 1
   
    w4_maskarr = np.where(v_result['W4mag'].value.mask)[0]
    v_result['W4_lim'][w4_maskarr] = 1
    
    ew4_maskarr = np.where(v_result['e_W4mag'].value.mask)[0]
    #v_result.remove_rows(ew4_maskarr)
    #print(len(v_result['AllWISE']))
    v_result['W4_lim'][ew4_maskarr] = 1
    
    v_result = Table(v_result, masked=False)
    
    ul_list = list(ew1_maskarr) + list(ew2_maskarr) + list(ew3_maskarr) + list(ew4_maskarr)
    UL_index = list(set(ul_list))
    
    
    #Create new table of srces with UL in their colors
    t = Table(v_result, masked=False) #so t_ul is of the AW sources with UL. Therefore, if the non-upper limits already meet this requirement of being in all 3, we don't need to check again!
    t_ul = t[UL_index]
    
    #now, delete these from v_result, and 
    v_result.remove_rows(UL_index)
    
    
    '''
    w1_maskarr = np.where(v_result['W1mag'].value.mask)[0]
    v_result.remove_rows(w1_maskarr)
    print(len(v_result['AllWISE']))
    
    
    w2_maskarr = np.where(v_result['W2mag'].value.mask)[0]
    v_result.remove_rows(w2_maskarr)
    print(len(v_result['AllWISE']))
    
    w3_maskarr = np.where(v_result['W3mag'].value.mask)[0]
    v_result.remove_rows(w3_maskarr)
    print(len(v_result['AllWISE']))
    
    w4_maskarr = np.where(v_result['W4mag'].value.mask)[0]
    v_result.remove_rows(w4_maskarr)
    print(len(v_result['AllWISE']))
    '''
    
    
    #Can see how many remain after each step
    
    #There should be none left
    #print(w3_maskarr)
    
    post_uL_row_num = len(v_result['W1mag'])
    print("No upper limits: " + str(post_uL_row_num))
    
    
    #Now actually plot Wise Blazar Strips
    #We are currently only doing W1-W2 and W2-W3 but not W3-W4
    v_result.add_column(v_result['W1mag'] - v_result['W2mag'], name='w1w2')
    v_result.add_column(v_result['W2mag'] - v_result['W3mag'], name='w2w3')
    v_result.add_column(v_result['W3mag'] - v_result['W4mag'], name='w3w4')
    
    
    w1w2 = v_result['w1w2']
    w2w3 = v_result['w2w3']
    w3w4 = v_result['w3w4']
    
    #Create Columns for errors - errors are positive
    #so this is going to create nan's for where one (or both?) are an upper limit
    v_result.add_column(v_result['e_W1mag'] + v_result['e_W2mag'], name='w1w2_err')
    v_result.add_column(v_result['e_W2mag'] + v_result['e_W3mag'], name='w2w3_err')
    v_result.add_column(v_result['e_W3mag']+v_result['e_W4mag'], name='w3w4_err')
    
    
    w1w2_err = v_result['w1w2_err']
    w2w3_err = v_result['w2w3_err']
    w3w4_err = v_result['w3w4_err']
    
    
    #Also do for UL table
    t_ul.add_column(t_ul['W1mag'] - t_ul['W2mag'], name='w1w2')
    t_ul.add_column(t_ul['W2mag'] - t_ul['W3mag'], name='w2w3')
    t_ul.add_column(t_ul['W3mag'] - t_ul['W4mag'], name='w3w4')
    
    
    #Want to check for error bars as well. 
    #w1w2 and w2w3 are like y0 and x0. Will need a double for checking x0,y0 with x0y0  - +/- for both x and y directions
    
    
    #so now check (y_up,w2w3) (y_low,w2w3), (w1w2,x_up), (w1w2,x_low),
    #write as loop to go through each dimension of the blazar strip
    
    #Would probably work better if had two separate tables, one with errors and one without, so can keep
    #the old and make something for the new. I am hoping that we can find a technique that does not 
    #include looping, and instead just operate on all UL, but I find that less likely...
    
    for h in range(3):
        if h ==0: #W1-W2 vs. W2-W3
            #if v_result['W1_lim'][]:
                
            y_up = w1w2 + w1w2_err
            y_low = w1w2 - w1w2_err
            x_up = w2w3 + w2w3_err
            x_low = w2w3 - w2w3_err
            
            data_points = np.array([[0.0]*2]*len(w2w3))
            dp_yup = np.array([[0.0]*2]*len(w2w3))
            dp_xup = np.array([[0.0]*2]*len(w2w3))
            dp_ylow = np.array([[0.0]*2]*len(w2w3))
            dp_xlow = np.array([[0.0]*2]*len(w2w3))
    
            #Now do the same for the rest of the error range combos
            for i in range(len(y_up)):
                dp_yup[i,0] = w2w3[i] #x value
                dp_yup[i,1] = y_up[i] #y value
                
            for i in range(len(y_low)):
                dp_ylow[i,0] = w2w3[i] #x value
                dp_ylow[i,1] = y_low[i] #y value
                
            for i in range(len(x_up)):
                dp_xup[i,0] = x_up[i] #x value
                dp_xup[i,1] = w1w2[i] #y value 
    
            for i in range(len(x_low)):
                dp_xlow[i,0] = x_low[i] #x value
                dp_xlow[i,1] = w1w2[i] #y value 
    
    
    
            #W1W2 vs. W2W3 (first one)
            #Create polygons of blazar WISE strip
            
            poly_bll = Polygon([(2.01,0.37),(3.30,1.17),(2.59,1.20),(1.52,0.51)]) #had wrong y for 3rd point
            poly_fsrq = Polygon([(2.90,0.85),(3.81,1.17),(3.29,1.67),(2.29,1.08)])
            #path = mpltPath.Path(polygon)
            #poly_bll
            #poly_fsrq
            
        elif h ==1: #W2-W3 vs. W3-W4
            data_points = np.array([[0.0]*2]*len(w2w3))
            dp_yup = np.array([[0.0]*2]*len(w2w3))
            dp_xup = np.array([[0.0]*2]*len(w2w3))
            dp_ylow = np.array([[0.0]*2]*len(w2w3))
            dp_xlow = np.array([[0.0]*2]*len(w2w3))
            
            y_up = w2w3 + w2w3_err
            y_low = w2w3 - w2w3_err
            x_up = w3w4 + w3w4_err
            x_low = w3w4 - w3w4_err
            #Add in the upper/lower limits 
    
            dp_yup = np.array([[0.0]*2]*len(w2w3))
            dp_xup = np.array([[0.0]*2]*len(w2w3))
            dp_ylow = np.array([[0.0]*2]*len(w2w3))
            dp_xlow = np.array([[0.0]*2]*len(w2w3))
    
            for i in range(len(y_up)):
                dp_yup[i,0] = w3w4[i] #x value
                dp_yup[i,1] = y_up[i] #y value
                
            for i in range(len(y_low)):
                dp_ylow[i,0] = w3w4[i] #x value
                dp_ylow[i,1] = y_low[i] #y value
                
            for i in range(len(x_up)):
                dp_xup[i,0] = x_up[i] #x value
                dp_xup[i,1] = w2w3[i] #y value 
    
            for i in range(len(x_low)):
                dp_xlow[i,0] = x_low[i] #x value
                dp_xlow[i,1] = w2w3[i] #y value 
            
            poly_bll = Polygon([(2.20,1.65),(2.72,2.57),(2.29,3.30),(1.20,1.96)])
            poly_fsrq = Polygon([(2.25,2.22),(3.04,3.05),(2.67,3.70),(1.68,2.85)])
                
                
        elif h ==2: #W1-W2 vs. W3-W4
            y_up = w1w2 + w1w2_err
            y_low = w1w2 - w1w2_err
            x_up = w3w4 + w3w4_err
            x_low = w3w4 - w3w4_err
            #Add in the upper/lower limits 
    
            dp_yup = np.array([[0.0]*2]*len(w2w3))
            dp_xup = np.array([[0.0]*2]*len(w2w3))
            dp_ylow = np.array([[0.0]*2]*len(w2w3))
            dp_xlow = np.array([[0.0]*2]*len(w2w3))
    
            for i in range(len(y_up)):
                dp_yup[i,0] = w3w4[i] #x value
                dp_yup[i,1] = y_up[i] #y value
                
            for i in range(len(y_low)):
                dp_ylow[i,0] = w3w4[i] #x value
                dp_ylow[i,1] = y_low[i] #y value
                
            for i in range(len(x_up)):
                dp_xup[i,0] = x_up[i] #x value
                dp_xup[i,1] = w1w2[i] #y value 
    
            for i in range(len(x_low)):
                dp_xlow[i,0] = x_low[i] #x value
                dp_xlow[i,1] = w1w2[i] #y value 
                
            #W1-W2 vs. W3-W4
            poly_bll = Polygon([(2.05,0.33),(2.83,1.07),(2.28,1.21),(1.20,0.73)])
            poly_fsrq = Polygon([(2.48,0.78),(3.05,1.17),(2.55,1.50),(1.72,1.12)])
        
        
        #Done with if statements but still have h=0,1,2 loop values, now add to the running list of truthers
        #This portion of the code remains the same for each h value
        #convert data_points to points
        bll_area_list = []
        fsrq_area_list = []
        #Identify if the WISE blazar strip contains our data points
        for i in range(post_uL_row_num):
            point = Point(data_points[i])
            bll_area_list.append(poly_bll.contains(point))
            fsrq_area_list.append(poly_fsrq.contains(point))
    
        #points = Point(data_points)
    
        #print(poly_bll.contains(points))
        len(bll_area_list)
    
        bal_yup = []
        fal_yup = []
    
        for i in range(post_uL_row_num):
            point = Point(dp_yup[i])
            bal_yup.append(poly_bll.contains(point))
            fal_yup.append(poly_fsrq.contains(point))
    
        bal_ylow = []
        fal_ylow = []
    
        for i in range(post_uL_row_num):
            point = Point(dp_ylow[i])
            bal_ylow.append(poly_bll.contains(point))
            fal_ylow.append(poly_fsrq.contains(point))
            
            
        bal_xup = []
        fal_xup = []
    
        for i in range(post_uL_row_num):
            point = Point(dp_xup[i])
            bal_xup.append(poly_bll.contains(point))
            fal_xup.append(poly_fsrq.contains(point))
            
        bal_xlow = []
        fal_xlow = []
    
        for i in range(post_uL_row_num):
            point = Point(dp_xlow[i])
            bal_xlow.append(poly_bll.contains(point))
            fal_xlow.append(poly_fsrq.contains(point))
            
    
        #Identify which sources fall within this region
        truthers_bll = list(filter(lambda i: bll_area_list[i], range(len(bll_area_list))))
        truthers_fsrq = list(filter(lambda i: fsrq_area_list[i], range(len(fsrq_area_list))))
        #truthers_bll
        #tbpf = truthers_bll +truthers_fsrq
        #print(tbpf)
    
    
        
        
        #Do for up/lower limits
        truthers_bll_yup = list(filter(lambda i: bal_yup[i], range(len(bal_yup))))
        truthers_fsrq_yup = list(filter(lambda i: fal_yup[i], range(len(fal_yup))))
    
        truthers_bll_ylow = list(filter(lambda i: bal_ylow[i], range(len(bal_ylow))))
        truthers_fsrq_ylow = list(filter(lambda i: fal_ylow[i], range(len(fal_ylow))))
    
        truthers_bll_xlow = list(filter(lambda i: bal_xlow[i], range(len(bal_xlow))))
        truthers_fsrq_xlow = list(filter(lambda i: fal_xlow[i], range(len(fal_xlow))))
    
        truthers_bll_xup = list(filter(lambda i: bal_xup[i], range(len(bal_xup))))
        truthers_fsrq_xup = list(filter(lambda i: fal_xup[i], range(len(fal_xup))))
    
    
        truthers_bll = truthers_bll + truthers_bll_yup + truthers_bll_ylow + truthers_bll_xup + truthers_bll_xlow
        #bll_name = v_result['AllWISE'][truthers_bll]
        #bll_name
    
    
        truthers_fsrq = truthers_fsrq +truthers_fsrq_yup + truthers_fsrq_ylow + truthers_fsrq_xup +truthers_fsrq_xlow
        #fsrq_name = v_result['AllWISE'][truthers_fsrq]
        #fsrq_name
    
        tbpf = truthers_bll + truthers_fsrq
        no_dup = list(set(tbpf)) #have to be careful using the same indice for for loops (i-> h)
        print(tbpf)#so we can see if duplicates
        print('\n')
        if h ==0: #W12 vs. W23
            A = no_dup
        elif h ==1:#W23 vs. W34
            B = no_dup
        elif h == 2:#W12 vs W34
            C = no_dup
    #Pull this out of h loop
    #now, the h loop is done, but still have the k value

    '''
    So, we can 
    '''    
    A_B = set(A) & set(B)
    ABC = set(A_B) & set(C)
    in3 = list(ABC)
    print('Number of srcs in all 3 strips: ' + str(len(in3)) + '\n')
    #Running list with indices of sources that are in all 3 blazar strips
    if len(in3) > 0:#aka is there a source in all 3D
        all_3_index.append(k) #this is the value that will be reported
        blz_assoc_list.append(df['assoc'][k])
        #################
        #Now address the sources that have upper limits for their Wmags (good news, don't have to bother with color error)
        #Note: this is inside the loop to do each blazar strip
        '''
        If W1 is UL, then UPPER limit
        
        If W2 is UL for W1-W2 then LOWER limit
        
        If W2 is UL for W2-W3 then UPPER limit
        
        If W3 is UL, for W3-W4 then UPPER limit
        
        If W3 is UL then W2W3 is LOWER limit
        
        If W4 is UL, then UPPER limit
        '''
        
        
        #yes need to check both dimensions 
        
        #skip this step if k is in all_3_index
        '''
        ADD THIS^^^
        pull out of h loop, and then disregard if already present in the all_3_index list
        then do similar, each p value has its own h0, h1, h2 list, and compare these, if in all 3 get p value, and reflect back to
        t_ul, which is completely different since 2 different Tables!
        
        What I have rn is ok, but you need blazars_srcs_AW_3d to compare to df
        
        READ BELOW!!!
        I need to go, reogranize your thoughts
        df - 536
        
        each one of df checks AW (up to 3000) (previously an UL get rid of) 
        NOW, also check UL for Blazar colors!
        So if fin_ul_3d has any (and not already in all_3_index), add the k value to that! (aka make the 2nd round (UL) not inside the 
        previous h one, so don't have to check if k is already in all_3_index. 
        
        Instead, we add to 
        '''
       #Here the above, non-UL code ends. 
        
       
    #Each source will have it's own list, as p is going through t_ul, values that have an upper limit (no longer eliminating)
    #Recall, t_ul was part of df, but made a separate data table
    if k in all_3_index:
        continue #skip the rest of code since already has a source in all 3 dimensions
    h0_ul = []
    h1_ul = []
    h2_ul = []
    #Go through all the sources in each strip one at a time    
    for p in range(len(t_ul)):#len(t_ul) #since separate data set now, index is unique from original
        #p = 349
        for h in range(3):
            if h == 0: #W1-W2 vs W2-W3
                #Left side of bll, if W3 is UL, W2-W3 is lower limit, so if on left of strip, could be in
                #better to just say if right of right side you eliminate
                
                x1bl = 2.01
                y1bl = 0.37
                x2bl = 3.30
                y2bl = 1.17
                x3bl = 2.59
                y3bl = 1.20
                x4bl = 1.52
                y4bl = 0.51
                
                x1f = 2.90
                y1f = 0.85
                x2f = 3.81
                y2f = 1.17
                x3f = 3.29
                y3f = 1.67
                x4f = 2.29
                y4f = 1.08
                #for setting up the line, we will do a first approximation of just bll left line, 
                #and not bother with fsrq left line. FSRQ will be with W1W2
                
                #set n_ out here?
                #pointing(x,y, x1,y1,x2,y2)
                n_12bl = pointing(t_ul['w2w3'][p],t_ul['w1w2'][p],x1bl, y1bl, x2bl, y2bl)
                
                n_23bl = pointing(t_ul['w2w3'][p],t_ul['w1w2'][p],x2bl, y2bl, x3bl, y3bl)
                
                n_41bl = pointing(t_ul['w2w3'][p],t_ul['w1w2'][p],x1bl, y1bl, x4bl, y4bl)
                
                n_34bl = pointing(t_ul['w2w3'][p],t_ul['w1w2'][p], x4bl,y4bl, x3bl, y3bl) #flipped order to match my initial intuition on its direction
                
                n_23f = pointing(t_ul['w2w3'][p],t_ul['w1w2'][p],x2f, y2f, x3f, y3f)
                
                n_34f = pointing(t_ul['w2w3'][p],t_ul['w1w2'][p],x4f, y4f, x3f, y3f) #flipped order so matched my incorrect sign haha
                
                n_12f = pointing(t_ul['w2w3'][p],t_ul['w1w2'][p],x1f,y1f,x2f,y2f)
                
                '''
                n_12bl #1 -2 BLL 
                left: -
                right: +
                
                n_41bl #4-1 BLL
                under: -
                above: +
                
                n_23f #2-3 FSRQ
                under: -
                above: +
                
                n34bl: #3-4 BLL
                left: -
                right: +
                
                
                Currently written so that if at the lower limit in one dim could put it in the blazar strip,
                it is included. (hence, the other is actually at the lower limit)
                
                
                
                
                
                KEY: LOWER LIMIT = SMALLEST VALUE of the color
                     UPPER LIMIT = LARGEST VALUE of the color
                     UL = smallest number (smaller number is brighter)
                     LL = largest number (larger number is dimmer)
                
                If W1 is UL, thenLOWER limit on W1-W2  (smallest value, can get bigger) -> color value can be bigger
                
                If W2 is UL for W1-W2 then UPPER limit (largest value calculatd, can get smaller) -> color value can be smaller
                
                If W2 is UL for W2-W3 then LOWER limit (smallest mag value calculated, can get bigger) -> color value can be bigger
                
                If W3 is UL, for W3-W4 then LOWER limit (smallest mag value calculated, can get bigger) -> color value can be bigger
                
                If W3 is UL then W2-W3 is UPPER limit (largest magnitude of mag value calculated, can get smaller) -> color value can be smaller
                
                If W4 is UL, then UPPER limit (largest magnitude of mag value calculated, can get smaller) -> color value can be smaller
                
                
                
                For all of these, we approximate some borders with using either the FSRQ or BLL edge, but not both. 
                '''
                #h=0
                #With regards to W1-W2 vs. W2-W3 (y axis values)
                if (t_ul['W1_lim'][p] == 1):
                    #pass #will remain upper limit
                    #Since only doing 1D, split into an or conditional with the BLL and FSRQ portion
                    #of the plot - way easier, just do for each...
                    if ((t_ul['w2w3'][p]<x2f and t_ul['w2w3'][p]>x2bl and n_34f>0 and n_23f<0) or (t_ul['w2w3'][p]>x4bl and t_ul['w2w3'][p]<x2bl and n_34bl>0 and n_23bl<0)):
                        h0_ul.append(p)
                if t_ul['W2_lim'][p] == 1:
                    #affects W1-W2 -> LL 
                    #add due to Y
                    if (t_ul['w2w3'][p]<x2f and t_ul['w2w3'][p]>x4bl and n_12f<0): #approximate with n12f for both bll and fsrq
                        #print(n_)
                        h0_ul.append(p)
                    #W2-W3 (X axis)    
                    elif (t_ul['w1w2'][p]<y3f and t_ul['w1w2'][p]>y1bl and n_12bl<0 and n_23f<0):#approximate with n12bl line (anything to left of it)
                        h0_ul.append(p)
                
                
                if t_ul['W3_lim'][p] ==1:
                    #Add due to X
                    #n_ = pointing(t_ul['w2w3'][p],t_ul['w1w2'][p],x1bl_right, y1bl_right, x2bl_right, y2bl_right)
                    #left of bll right - will be ok for W2-W3
                    if (n_12bl>0 and t_ul['w1w2'][p]>y1bl and t_ul['w1w2'][p]<y3f and n_41bl>0): #have to check bottom left corner 
                        #n_23 = pointing(t_ul['w1w2'][p],t_ul['w2w3'][p],x1bl, y1bl, x2bl, y2bl)
                        #print('poingint from left BLL' + str(n_23))
                        h0_ul.append(p)
                               
                if t_ul['W4_lim'][p] == 1:
                    pass #not in this dimension
                    
                
            elif h ==1: #W2-W3 vs. W3-W4
                x1bl = 2.20
                y1bl = 1.65
                x2bl = 2.72
                y2bl = 2.57
                x3bl = 2.29
                y3bl = 3.30
                x4bl = 1.20
                y4bl = 1.96
                
                x1f = 2.25
                y1f = 2.22
                x2f = 3.04
                y2f = 3.05
                x3f = 2.67
                y3f = 3.70
                x4f = 1.68
                y4f = 2.85
                #relation to line
                n_12bl = pointing(t_ul['w3w4'][p],t_ul['w2w3'][p],x1bl, y1bl, x2bl, y2bl)
                
                n_41bl = pointing(t_ul['w3w4'][p],t_ul['w2w3'][p],x1bl, y1bl, x4bl, y4bl) # (-) if under, (+) if above
                
                n_34bl = pointing(t_ul['w3w4'][p],t_ul['w2w3'][p],x4bl, y4bl, x3bl, y3bl)
                
                n_23f = pointing(t_ul['w3w4'][p],t_ul['w2w3'][p],x2f, y2f, x3f, y3f)
                
                
                
                if t_ul['W1_lim'][p] == 1:
                    pass #not in this dimension
                
                if t_ul['W2_lim'][p] == 1:
                    if (t_ul['w3w4'][p]>x4bl and t_ul['w3w4'][p]<x2f and n_34bl>0 and n_23f<0):
                        h1_ul.append(p)
                if t_ul['W3_lim'][p] == 1:
                    #Y 
                    if (n_12bl < 0 and t_ul['w3w4'][p]>x4bl  and t_ul['w3w4'][p]<x2f and n_41bl>0): #t_ul['w2w3'][p]>y3f
                        #print(n_)
                        h1_ul.append(p)
                        #print(p)
                    #X    
                    elif (t_ul['w2w3'][p]>y1bl and t_ul['w2w3'][p]<y3f and n_12bl<0 and n_23f<0): #elif, since if already in it, don't need it again
                        h1_ul.append(p)
                        #print('pizza')
                    #elif (n_23f>0) and (n_)
                    
                    
                if t_ul['W4_lim'][p] == 1: #X
                    if (n_12bl>0 and t_ul['w2w3'][p]>y1bl and t_ul['w2w3'][p]<y3f and n_41bl>0): 
                        #n_23 = pointing(t_ul['w1w2'][p],t_ul['w2w3'][p],x1bl, y1bl, x2bl, y2bl)
                        #print('poingint from left BLL' + str(n_23))
                        h1_ul.append(p)
                        #print(p)
                
            elif h ==2: #W1-W2 vs. W3-W4
                x1bl = 2.05
                y1bl = 0.33
                x2bl = 2.83
                y2bl = 1.07
                x3bl = 2.28
                y3bl = 1.21
                x4bl = 1.20
                y4bl = 0.73
                
                x1f = 2.48
                y1f = 0.78
                x2f = 3.05
                y2f = 1.17
                x3f = 2.55
                y3f = 1.50
                x4f = 1.72
                y4f = 1.12
                #Point Lying
                n_12bl = pointing(t_ul['w3w4'][p],t_ul['w1w2'][p],x1bl, y1bl, x2bl, y2bl)
                
                n_41bl = pointing(t_ul['w3w4'][p],t_ul['w1w2'][p],x1bl, y1bl, x4bl, y4bl)
                
                n_23f = pointing(t_ul['w3w4'][p],t_ul['w1w2'][p],x2f, y2f, x3f, y3f)
                
                n_12f = pointing(t_ul['w3w4'][p],t_ul['w1w2'][p],x1f, y1f, x2f, y2f)
                
                n_34f = pointing(t_ul['w3w4'][p],t_ul['w1w2'][p],x4f, y4f, x3f, y3f)
            
                if (t_ul['W1_lim'][p] == 1):
                    if ((t_ul['w3w4'][p]>x4bl and t_ul['w3w4'][p]<x2f and n_34f>0  and n_41bl <0)): 
                        h2_ul.append(p)
                        
                if t_ul['W2_lim'][p] == 1:
                    #
                    #add due to Y
                    if (t_ul['w3w4'][p]<x2f and t_ul['w3w4'][p]>x4bl and n_12f<0 and n_41bl>0):
                        h2_ul.append(p)
                    
                if t_ul['W3_lim'][p] == 1:
                    if (t_ul['w1w2'][p]>y1bl and t_ul['w1w2'][p]<y3f and n_12f<0 and n_23f<0):
                        h2_ul.append()
                if t_ul['W4_lim'][p] == 1: #X
                    if (t_ul['w1w2'][p]>y1bl and t_ul['w1w2'][p]<y3f and n_34f >0 and n_41bl>0):
                        #(n_12bl<0) and (t_ul['w1w2'][p]>y1bl) and (t_ul['w1w2'][p]<y3f) and (n_23f<0): 
                        #n_23 = pointing(t_ul['w1w2'][p],t_ul['w2w3'][p],x1bl, y1bl, x2bl, y2bl)
                        #print('poingint from left BLL' + str(n_23))
                        h2_ul.append(p)
        
        
        #This is where we are still in p but not in h loop
        #Report all dimensions that have UL that could fall in the regions
        
        h_01 = set(h0_ul) & set(h1_ul)
        h_012 = set(h_01) & set(h2_ul)
        fin_ul_3 = list(h_012)
        if len(fin_ul_3)>0:
            all_3_ul.append(k) #if find it earlier in the range of p, then will just keep adding it #k is the index for the source detected, that is from the pandas dataframe we input
            blz_assoc_list.append(df['assoc'][k])
            break

        
    #combine lists
    # A_B = set(A) & set(B)
    # ABC = set(A_B) & set(C)
    # final_in3 = list(ABC)
    # print('Number of srcs in all 3 strips: ' + str(len(final_in3)) + '\n')
    # #Running list with indices of sources that are in all 3 blazar strips
    # if len(final_in3) > 0:
    #     all_3_index.append(k)


#Note: we can also possibly have new blazar candidates if we look at those that are in all 3 but are not XMed -> this will be for when pulling an r_95
#Then we don't have to do that last part by hand as we do now.... aka just do NOT in all_3_index
#print(all_3_index)




blazar_srcs_AW_3d = all_3_index + all_3_ul
blazar_srcs_AW_3d.sort()

print(blazar_srcs_AW_3d)
#Values in df that are NOT in blazar strip 3D
#l1 = df

df_list = list(df.index)
df_not_blazars = list(set(df_list).difference(blazar_srcs_AW_3d))


#Now, we will get rid of any remaining flare that has an 'assoc' that has already been determined by a different flaring
#Event to be a blazar by AW colors. 

nodup_assoc = list(set(list(blz_assoc_list))) #won't remove them now but plan would be to remove those that have previous assoc with
#already ruled out blazars. -> part of this may be associating, BUT ACTUALLY UNIQUE
#just because in this list doesn't mean any more survived!



#I based on what is slipping thru, it is sources just on the edge of the r_95. 

'''
2nd Catalog
'assoc' is based on 99% stat err + r_sys - solely positional. 
hence, if already in the 99% area and been ruled out, will take out. 

Oddly, missing 1 of the novae - likely due to being in 3FGL g-ray uncert. .... Need to check still
'''
#nodup_assoc is the list of all the associated things that came up with AW colors like a blazar for some detection -> remove those 
#we know are not blazars. 
#ANOTHER TEST -> Novae catalog and where fall in blazar strip. 
#Remove things we know are NOT blazars

#removing things that have blazar colors, that the association IS KNOWN TO BE GALACTIC!!!
#feel free to add or let me know of other sources
for each in nodup_assoc:
    if each == 'None':
        nodup_assoc.remove(each)
    elif each == 'LS I+61 303':
        nodup_assoc.remove(each)
    elif each == 'IC 443':
        nodup_assoc.remove(each)
    elif each == 'V407 Cyg':
        nodup_assoc.remove(each)
    elif each == 'Cyg X-3':
        nodup_assoc.remove(each)
    elif each == 'HESS J1303-631':
        nodup_assoc.remove(each)
    elif each == 'PSR B1259-63':
        nodup_assoc.remove(each)
    elif each == 'PSR J2032+4127':
        nodup_assoc.remove(each)
    elif each == 'PSR J0534+2200':
        nodup_assoc.remove(each)
    elif each == 'PSR J1826-1256':
        nodup_assoc.remove(each)
    elif each == 'NVSS J085238-312331':
        nodup_assoc.remove(each)
    elif each == 'PMN J1802-3940':
        nodup_assoc.remove(each)
    
    
'''
#If not in the list you will have to run each of these one at a time to confirm (could check to see if there is anything in this list and then run)
nodup_assoc.remove('None')
nodup_assoc.remove('LS I+61 303')
nodup_assoc.remove('IC 443')
nodup_assoc.remove('V407 Cyg') #means was previously caught - we should do a AW blazar strip test with novae as well. 
nodup_assoc.remove('Cyg X-3')
nodup_assoc.remove('LS I+61 303')
nodup_assoc.remove('HESS J1303-631')
nodup_assoc.remove('PSR B1259-63')
nodup_assoc.remove('PSR J2032+4127')
nodup_assoc.remove('LAT PSR J2032+4127')
nodup_assoc.remove('PSR J0534+2200')
nodup_assoc.remove('PSR J0248+6021')
nodup_assoc.remove('LAT PSR J1826-1256')
#nodup_assoc.remove('PSR J1826-1256')
nodup_assoc.remove('NVSS J085238-312331')
#nodup_assoc.remove('PMN J1802-3940)
#nodup_assoc.remove('')
'''




#are we interested in putting back in ones that failed and are known Galactic?

'''
Notorious Blazars
PMN86 J1802-394
B2 2114+33
'''



#Want the rows that ARE NOT in nodup_assoc
#Just say "if assoc by FAVA with known blazar at end, we remove" (could have done at beginning but )
no_prev_blz_candidates = df.loc[df_not_blazars]

for q in range(len(nodup_assoc)):
    no_prev_blz_candidates = no_prev_blz_candidates[no_prev_blz_candidates['assoc'].isin([nodup_assoc[q]]) == False]
    
no_prev_blz_candidates_fin = no_prev_blz_candidates[no_prev_blz_candidates['assoc'].isin(['PMN J1802-3940']) == False]
no_prev_blz_candidates_fin.reset_index(inplace=True)

#Uncomment below to get the flaring Gal. Non blazar candidates
#no_prev_blz_candidates_fin.to_csv('')



'''
Should probably have this output the AllWISE sources that do fall in all 3
'''
print('AW sources that Fall in all 3 blazar strips' + '\n')
print('Not using upper limits')

print(v_result[in3]['AllWISE'])
# we can also save this to a csv file
#uncomment below
aw_srcs_with_blz_colors = v_result[in3]
#aw_src_with_blz_colors.to_csv('AW_srcs_with_blz_colors.csv')
print('Using UL')
try:
    print(t_ul[fin_ul_3]['AllWISE'])
except NameError:
    print('No UL AW srcs')


'''
Final Notes: 
    this is fine since we have FAVA to take advantage of them using their own r_99 + r_sys (0.1deg) to do the associations. 
    Will need to discuss this down the road, since we are doing intervention to clean it up.
    
Next: tests on knonw blazars (BZCAT) and known Gal. srcs (in 4FGL) and Novae

Side project of interest: Crab flaring - number of cases and when. Crab vs. the Pulsar
'''


#and then make sure to remove None from this list.
#then remove any remaining that are also associated with those that are blazars 

#print("Index Values in df NOT in blazar_srcs_AW_3d:",df_not_blazars)
#Then, we add those that fall in all 3 to the all_3_index. We are then interested in all the sources that are NOT in all_3_index

'''#Uncomment this after
#candidates = df.loc[df_not_blazars]
candidates.reset_index(inplace=True)
candidates.to_csv('/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/historic_735_th_5p5/FAV_cand4.csv')
'''




'''
#W2-W3 vs. W3-W4
#P1-P2
x1_2bll = [2.20, 2.72]
y1_2bll = [1.65, 2.57]

x1_2fsq = [2.25, 3.04]
y1_2fsq = [2.22, 3.05]


#P2-P3
x2_3bll = [2.72, 2.29]
y2_3bll = [2.57,3.30]

x2_3fsq = [3.04,2.67]
y2_3fsq = [3.05,3.70]

#P3-P4
x3_4bll = [2.29, 1.20]
y3_4bll = [3.30, 1.96]

x3_4fsq = [2.67,1.68]
y3_4fsq = [3.70,2.85]

#P4-P1
x4_1bll = [1.20, 2.20]
y4_1bll = [1.96, 1.65]

x4_1fsq = [1.68,2.25]
y4_1fsq = [2.85,2.22]

plt.plot(x1_2bll, y1_2bll, c='b' , linestyle='-')#, zorder = 1)
plt.plot(x2_3bll, y2_3bll, c = 'b',linestyle= '-')#, zorder = 1)
plt.plot(x3_4bll, y3_4bll, c = 'b', linestyle='-')#, zorder = 1)
plt.plot(x4_1bll, y4_1bll, c = 'b', linestyle='-')#, zorder = 1)


plt.plot(x1_2fsq, y1_2fsq, c= 'r', linestyle='-')
plt.plot(x2_3fsq, y2_3fsq, c= 'r', linestyle='-')
plt.plot(x3_4fsq, y3_4fsq, c= 'r', linestyle='-')
plt.plot(x4_1fsq, y4_1fsq, c= 'r', linestyle='-')
plt.scatter(v_result['w2w3'],v_result['w1w2'])
plt.errorbar(v_result['w2w3'],v_result['w1w2'],yerr=v_result['w1w2_err'],xerr=v_result['w2w3_err'])






#W1-W2 vs. W3-W4
#P1-P2
x1_2bll = [2.05, 2.83]
y1_2bll = [0.33, 1.07]

x1_2fsq = [2.48, 3.05]
y1_2fsq = [0.78, 1.17]


#P2-P3
x2_3bll = [2.83, 2.28]
y2_3bll = [1.07,1.21]

x2_3fsq = [3.05,2.55]
y2_3fsq = [1.17,1.5]

#P3-P4
x3_4bll = [2.28, 1.20]
y3_4bll = [1.21, 0.73]

x3_4fsq = [2.55,1.72]
y3_4fsq = [1.50,1.12]

#P4-P1
x4_1bll = [1.20, 2.05]
y4_1bll = [0.73, 0.33]

x4_1fsq = [1.72,2.48]
y4_1fsq = [1.12,0.78]

plt.plot(x1_2bll, y1_2bll, c='b' , linestyle='-')#, zorder = 1)
plt.plot(x2_3bll, y2_3bll, c = 'b',linestyle= '-')#, zorder = 1)
plt.plot(x3_4bll, y3_4bll, c = 'b', linestyle='-')#, zorder = 1)
plt.plot(x4_1bll, y4_1bll, c = 'b', linestyle='-')#, zorder = 1)


plt.plot(x1_2fsq, y1_2fsq, c= 'r', linestyle='-')
plt.plot(x2_3fsq, y2_3fsq, c= 'r', linestyle='-')
plt.plot(x3_4fsq, y3_4fsq, c= 'r', linestyle='-')
plt.plot(x4_1fsq, y4_1fsq, c= 'r', linestyle='-')





'''
'''
#FAVA th 5.5 results-> 536 srcs for newly adjusted cuts up until AW check
fav_th5p5_gal = [0, 1, 2, 3, 4, 6, 7, 9, 10, 12, 13, 14, 15, 18, 19, 20, 21, 22, 24, 25, 26, 27, 28, 31, 32, 33, 34, 36, 37, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 62, 63, 64, 66, 67, 68, 70, 72, 73, 75, 76, 77, 78, 79, 81, 82, 83, 84, 86, 87, 88, 89, 92, 93, 94, 96, 97, 99, 101, 102, 103, 104, 105, 108, 109, 110, 111, 112, 114, 115, 116, 117, 119, 120, 121, 122, 123, 125, 126, 127, 128, 130, 131, 133, 134, 136, 137, 138, 142, 144, 145, 146, 147, 148, 149, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 165, 166, 167, 168, 169, 170, 171, 172, 174, 175, 177, 178, 179, 180, 181, 182, 183, 185, 186, 187, 188, 190, 193, 194, 195, 196, 197, 198, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 216, 217, 218, 219, 220, 222, 223, 224, 225, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 241, 242, 243, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 257, 259, 260, 262, 264, 267, 268, 269, 270, 272, 274, 275, 276, 277, 278, 279, 280, 281, 283, 284, 285, 286, 287, 288, 289, 291, 292, 293, 294, 295, 296, 297, 298, 300, 301, 302, 303, 304, 305, 307, 308, 309, 310, 311, 312, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 330, 331, 333, 334, 335, 336, 337, 338, 340, 341, 342, 344, 345, 346, 347, 350, 352, 353, 354, 355, 357, 358, 359, 360, 361, 363, 365, 366, 369, 373, 375, 376, 377, 378, 379, 381, 384, 386, 387, 388, 389, 391, 392, 393, 394, 395, 397, 398, 399, 402, 403, 404, 407, 408, 409, 412, 413, 418, 419, 423, 424, 425, 427, 429, 430, 431, 432, 434, 435, 438, 440, 441, 444, 445, 446, 447, 450, 452, 456, 457, 461, 463, 470, 472, 473, 474, 477, 478, 479, 480, 482, 493, 498, 510, 513, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 529, 530, 531, 532, 533, 534, 535]
    #XM this list with 


#Ok so can do union-intersection to get what is not in intersection (since union will include everything)
#or shortcut is set(a) ^ set(b)

Not_blazar_colors =  set(df.index) ^ set(fav_th5p5_gal)


candidates_1 = df.loc[Not_blazar_colors]

candidates_1.reset_index(inplace=True)
#candidates_1.to_csv('/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/historic_735_th_5p5/FAV_cand1.csv')

no_assoc_fav = candidates_1.loc[(candidates_1['assoc'] == 'None') & (candidates_1['fglassoc'] == 'None')]
no_assoc_fav.reset_index(inplace=True)
#no_assoc_fav.to_csv('/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/historic_735_th_5p5/no_assoc_cand1.csv')
'''
