#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 14:08:20 2023

@author: joffresd
"""

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

'''
Idea to catch blazars: 
    I think we are catching them, but running into issues with just outside our own bubbles. 
    So make list of 'assoc' values that correspond with these blazars, and then remove any remaining that have that assoc still. 
'''


def pointing (x,y,x1,y1,x2,y2):
    n_ = (x-x1)*(y2-y1) + (y-y1)*(-x2+x1)
    return n_

#Known 4FGL Blazar Sources
#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/4fgl_blazars.csv")


# known 4FGL Galactic sources
#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/historic_735_th_5p5/gal_flaring_notinbzcat_notblzin4fgl_combinedsrcs_noCrab1.csv")
#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/historic_735_th_5p5/pre_Cand_no_bzcat_4fgl3fglblzs_noCrab_2.csv")

#Gal test (instead of XM to closest one, will just do all of them) from 4FGL
#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/gal_known_4fgl_dr3_eC_removed.csv")

#Novae
#df = pd.read_csv('/Users/joffresd/Documents/Clemson/Research/FAVA/novae_xm_AW.csv')

df=pd.read_csv('/Users/joffresd/Documents/Clemson/Research/FAVA/gcvs_nova_srcs.csv')

#Crab only
#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/historic_735_th_5p5/Candidates_6/crab_and_pulsar_flares.csv")

#Threshold=5 result
#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/historic_735_th_5/final_flaring_noCrab.csv")

#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/historic_735_th_5/Useful_cross_matches/no_bz_4fgl_3fgl_blzs_noCrab.csv")

wknum = 758
#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/Weekly_Analysis/" +str(wknum)+ "/fava_cuts_" +str(wknum) + ".csv")
#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/Weekly_Analysis/TC_xmes_week_" + str(wknum) +".csv")

#From 4FGL unclassified
#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/historic_735_th_5/Assoc_post_proposal/unk.csv")


#final candidates, check with more AW sources
#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/historic_735_th_5/Paper_Cand_v2_10k/FAV_cand2.csv")


#FINAL FAVA TEST WITH UL
#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/historic_735_th_5/Paper_Cand_v2_10k/FAV_cand2.csv")

#Dr. Hartmann Test - 10<b<20 (rest of criteria)
#df = pd.read_csv("/Users/joffresd/Documents/Clemson/Research/FAVA/pyFAVA/historic_735_th_5/tests/DH2_postxm.csv")

all_3_index = []
all_3_ul = []
blz_assoc_list = []





for k in range(len(df)): #len(df['flareID']) # df['RAJ2000'][k]
    if k == 196:
        continue
    if k ==187:
        continue
    #k= 187
    #k = 21
    #k=0
    #First do just as one and then convert to loop through all of the sources
    ra =  df['RA'][k]#df['best_ra'][k] #from 4FGL #df['RAJ2000'][k]
    dec =  df['DEC'][k] #df['best_dec'][k] #from 4FGL #df['DEJ2000'][k]
    search_radius = 0.0008333333333333334#df['Conf_95_SemiMajor'][k]#df['best_r95'][k] #from 4FGL #df['Conf_95_SemiMajor'][k]
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
    '''
    5/23/23 Note
    IMPORTANT! Need to update code - I did UL with larger number as UL - WRONG
    For magnitudes, an UL will be the SMALLEST value
    '''
    
    
    
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
        #blz_assoc_list.append(df['assoc'][k])
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
                        h2_ul.append(p)
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
            #blz_assoc_list.append(df['assoc'][k])
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

#nodup_assoc = list(set(list(blz_assoc_list))) #won't remove them now but plan would be to remove those that have previous assoc with
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
'''
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
#Testing in parts
Blazars
1. 0-500 : 
blz1 = [0, 1, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 36, 37, 38, 39, 40, 41, 42, 44, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 73, 74, 75, 76, 77, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 249, 250, 251, 252, 253, 254, 255, 257, 258, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 278, 279, 280, 281, 282, 283, 284, 285, 286, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 305, 306, 307, 308, 309, 310, 311, 312, 314, 315, 316, 317, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 348, 349, 350, 351, 352, 354, 355, 356, 357, 358, 359, 362, 364, 366, 368, 369, 371, 373, 374, 375, 378, 379, 382, 383, 384, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 402, 403, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 417, 418, 419, 421, 423, 424, 425, 426, 427, 428, 429, 430, 432, 435, 436, 437, 438, 439, 440, 442, 443, 445, 446, 448, 449, 450, 451, 454, 457, 461, 462, 463, 464, 465, 467, 468, 469, 470, 471, 473, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 486, 487, 490, 492, 494, 495, 496, 498, 500]
    
2. 501-1200
blz2 = [501, 502, 503, 504, 505, 506, 508, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 625, 626, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 697, 698, 699, 700, 701, 702, 703, 704, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 726, 727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 747, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 819, 820, 821, 822, 823, 824, 825, 826, 828, 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 840, 842, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 865, 866, 867, 868, 869, 871, 872, 873, 874, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890, 891, 892, 893, 894, 895, 896, 897, 898, 899, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 931, 932, 933, 934, 935, 937, 938, 939, 940, 941, 942, 943, 944, 945, 946, 948, 950, 951, 952, 953, 954, 955, 956, 957, 958, 959, 960, 961, 962, 963, 964, 965, 966, 967, 968, 969, 970, 971, 972, 973, 974, 975, 976, 978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998, 999, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1177, 1178, 1179, 1180, 1181, 1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189, 1190, 1191, 1192, 1193, 1194, 1195, 1196, 1197, 1198, 1199]
    
3. 1200-2000
blz3 = [1200, 1201, 1202, 1203, 1204, 1205, 1206, 1208, 1209, 1210, 1211, 1212, 1213, 1214, 1215, 1216, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1225, 1226, 1227, 1228, 1229, 1230, 1231, 1232, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1241, 1242, 1243, 1244, 1245, 1246, 1247, 1248, 1249, 1250, 1251, 1252, 1253, 1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1262, 1263, 1264, 1265, 1266, 1268, 1269, 1270, 1271, 1272, 1273, 1274, 1275, 1276, 1277, 1278, 1279, 1280, 1281, 1282, 1283, 1284, 1285, 1286, 1287, 1288, 1289, 1290, 1291, 1292, 1293, 1294, 1295, 1296, 1297, 1298, 1299, 1300, 1301, 1302, 1303, 1304, 1305, 1306, 1307, 1308, 1309, 1311, 1312, 1313, 1314, 1315, 1316, 1317, 1318, 1319, 1320, 1321, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331, 1332, 1333, 1334, 1335, 1336, 1337, 1338, 1339, 1340, 1341, 1342, 1343, 1344, 1345, 1346, 1347, 1348, 1349, 1350, 1351, 1352, 1353, 1354, 1355, 1356, 1357, 1358, 1359, 1360, 1361, 1362, 1363, 1364, 1365, 1367, 1368, 1369, 1370, 1371, 1372, 1373, 1374, 1375, 1376, 1377, 1378, 1379, 1380, 1381, 1382, 1383, 1384, 1385, 1386, 1387, 1388, 1389, 1390, 1391, 1392, 1393, 1394, 1395, 1396, 1397, 1398, 1399, 1400, 1401, 1402, 1403, 1404, 1405, 1406, 1407, 1408, 1409, 1411, 1412, 1413, 1414, 1415, 1416, 1417, 1418, 1419, 1420, 1421, 1422, 1423, 1424, 1425, 1426, 1427, 1428, 1429, 1430, 1431, 1432, 1433, 1434, 1435, 1436, 1437, 1438, 1439, 1440, 1441, 1442, 1443, 1444, 1445, 1446, 1447, 1448, 1449, 1450, 1451, 1452, 1453, 1454, 1455, 1456, 1457, 1458, 1459, 1460, 1461, 1462, 1463, 1464, 1465, 1466, 1467, 1468, 1469, 1470, 1471, 1472, 1473, 1474, 1475, 1477, 1478, 1479, 1480, 1481, 1482, 1483, 1484, 1485, 1486, 1487, 1488, 1489, 1490, 1491, 1492, 1493, 1494, 1495, 1496, 1497, 1498, 1499, 1500, 1501, 1502, 1503, 1504, 1505, 1506, 1507, 1508, 1509, 1510, 1511, 1513, 1514, 1516, 1517, 1518, 1519, 1520, 1521, 1522, 1523, 1524, 1525, 1526, 1527, 1528, 1529, 1530, 1531, 1532, 1533, 1534, 1535, 1536, 1537, 1538, 1539, 1540, 1541, 1542, 1543, 1544, 1545, 1546, 1547, 1548, 1549, 1550, 1551, 1552, 1553, 1554, 1555, 1556, 1557, 1558, 1559, 1560, 1561, 1562, 1563, 1564, 1565, 1566, 1567, 1568, 1569, 1570, 1571, 1572, 1573, 1575, 1576, 1577, 1578, 1579, 1580, 1581, 1582, 1583, 1584, 1585, 1587, 1588, 1589, 1590, 1591, 1592, 1593, 1594, 1595, 1597, 1598, 1599, 1600, 1601, 1602, 1603, 1604, 1605, 1606, 1607, 1608, 1609, 1610, 1611, 1612, 1613, 1614, 1615, 1616, 1617, 1618, 1619, 1620, 1621, 1622, 1623, 1624, 1625, 1626, 1627, 1628, 1629, 1630, 1631, 1632, 1633, 1634, 1635, 1636, 1637, 1638, 1639, 1640, 1641, 1642, 1643, 1644, 1645, 1646, 1647, 1648, 1649, 1650, 1651, 1652, 1653, 1654, 1655, 1656, 1657, 1658, 1659, 1660, 1661, 1662, 1663, 1664, 1665, 1666, 1667, 1668, 1669, 1670, 1672, 1673, 1674, 1675, 1676, 1678, 1679, 1680, 1681, 1682, 1683, 1684, 1685, 1686, 1687, 1688, 1689, 1690, 1691, 1692, 1693, 1694, 1695, 1696, 1697, 1698, 1699, 1700, 1701, 1702, 1703, 1704, 1705, 1706, 1707, 1708, 1709, 1710, 1712, 1713, 1714, 1715, 1716, 1717, 1718, 1719, 1720, 1721, 1722, 1723, 1724, 1725, 1726, 1728, 1729, 1730, 1732, 1733, 1734, 1736, 1737, 1738, 1739, 1740, 1741, 1743, 1744, 1745, 1746, 1747, 1748, 1749, 1750, 1751, 1752, 1753, 1754, 1755, 1756, 1757, 1758, 1760, 1761, 1762, 1763, 1764, 1765, 1766, 1767, 1768, 1769, 1770, 1771, 1772, 1773, 1774, 1775, 1776, 1778, 1779, 1780, 1782, 1783, 1784, 1785, 1787, 1788, 1789, 1791, 1792, 1793, 1794, 1795, 1796, 1797, 1798, 1800, 1801, 1802, 1803, 1805, 1806, 1807, 1808, 1809, 1810, 1811, 1812, 1813, 1814, 1815, 1816, 1817, 1818, 1820, 1821, 1822, 1823, 1824, 1826, 1827, 1828, 1829, 1830, 1831, 1832, 1833, 1834, 1835, 1836, 1837, 1838, 1839, 1841, 1842, 1843, 1844, 1845, 1846, 1847, 1848, 1849, 1850, 1851, 1852, 1853, 1854, 1855, 1856, 1857, 1858, 1859, 1860, 1861, 1862, 1863, 1864, 1865, 1866, 1868, 1869, 1870, 1871, 1872, 1873, 1874, 1876, 1877, 1878, 1879, 1880, 1881, 1882, 1883, 1884, 1885, 1886, 1887, 1888, 1889, 1890, 1891, 1892, 1893, 1894, 1895, 1896, 1897, 1898, 1900, 1901, 1902, 1903, 1904, 1905, 1906, 1907, 1908, 1909, 1910, 1911, 1912, 1913, 1914, 1915, 1916, 1917, 1918, 1919, 1920, 1921, 1922, 1923, 1924, 1925, 1926, 1927, 1928, 1930, 1931, 1932, 1934, 1935, 1936, 1937, 1938, 1939, 1940, 1941, 1942, 1943, 1944, 1945, 1946, 1947, 1948, 1949, 1950, 1951, 1952, 1954, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999]
    
4. 2000-3000
blz4 = [2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029, 2030, 2031, 2032, 2033, 2034, 2035, 2036, 2037, 2038, 2039, 2040, 2041, 2042, 2043, 2044, 2045, 2046, 2047, 2048, 2049, 2050, 2051, 2052, 2053, 2054, 2055, 2056, 2057, 2058, 2059, 2060, 2061, 2062, 2063, 2064, 2065, 2066, 2067, 2068, 2069, 2070, 2071, 2072, 2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2083, 2084, 2085, 2087, 2088, 2089, 2090, 2091, 2092, 2093, 2094, 2095, 2096, 2097, 2098, 2099, 2100, 2101, 2102, 2103, 2104, 2105, 2106, 2107, 2108, 2109, 2110, 2111, 2112, 2114, 2116, 2117, 2118, 2119, 2120, 2121, 2122, 2123, 2124, 2125, 2126, 2127, 2128, 2129, 2130, 2131, 2132, 2133, 2134, 2135, 2136, 2137, 2138, 2139, 2140, 2141, 2142, 2143, 2144, 2145, 2146, 2147, 2148, 2149, 2150, 2151, 2152, 2153, 2154, 2155, 2156, 2157, 2158, 2159, 2160, 2161, 2162, 2163, 2164, 2165, 2166, 2167, 2168, 2169, 2170, 2171, 2172, 2173, 2174, 2175, 2176, 2177, 2178, 2179, 2180, 2181, 2182, 2183, 2184, 2185, 2186, 2187, 2188, 2189, 2190, 2191, 2192, 2193, 2194, 2195, 2196, 2197, 2198, 2199, 2200, 2201, 2202, 2203, 2204, 2205, 2206, 2207, 2208, 2209, 2210, 2211, 2212, 2213, 2214, 2215, 2216, 2217, 2218, 2219, 2220, 2221, 2222, 2223, 2224, 2225, 2226, 2227, 2228, 2229, 2230, 2231, 2232, 2233, 2234, 2235, 2236, 2237, 2238, 2239, 2240, 2241, 2242, 2243, 2244, 2245, 2246, 2247, 2248, 2249, 2250, 2251, 2252, 2253, 2254, 2255, 2256, 2257, 2258, 2259, 2260, 2261, 2262, 2263, 2264, 2265, 2266, 2267, 2268, 2269, 2270, 2271, 2272, 2273, 2274, 2275, 2276, 2277, 2278, 2279, 2280, 2281, 2282, 2283, 2284, 2285, 2286, 2287, 2288, 2289, 2290, 2291, 2292, 2293, 2294, 2295, 2296, 2297, 2298, 2299, 2300, 2301, 2302, 2303, 2304, 2305, 2306, 2307, 2308, 2309, 2310, 2311, 2312, 2313, 2314, 2315, 2316, 2317, 2318, 2319, 2320, 2321, 2322, 2323, 2324, 2325, 2326, 2327, 2328, 2329, 2330, 2331, 2332, 2333, 2334, 2335, 2336, 2337, 2338, 2339, 2340, 2341, 2342, 2343, 2344, 2345, 2346, 2347, 2348, 2349, 2350, 2351, 2352, 2353, 2354, 2355, 2356, 2357, 2358, 2359, 2360, 2361, 2362, 2363, 2364, 2365, 2366, 2367, 2368, 2369, 2370, 2371, 2372, 2373, 2374, 2375, 2376, 2377, 2378, 2379, 2380, 2381, 2382, 2383, 2384, 2385, 2386, 2387, 2388, 2389, 2390, 2391, 2392, 2393, 2394, 2395, 2396, 2397, 2398, 2399, 2400, 2401, 2402, 2403, 2404, 2405, 2407, 2408, 2409, 2411, 2412, 2413, 2414, 2415, 2416, 2417, 2418, 2419, 2420, 2421, 2422, 2423, 2424, 2425, 2426, 2427, 2428, 2429, 2430, 2431, 2432, 2433, 2434, 2435, 2436, 2437, 2438, 2439, 2440, 2441, 2442, 2443, 2444, 2445, 2447, 2448, 2449, 2450, 2451, 2452, 2453, 2454, 2455, 2456, 2457, 2458, 2459, 2460, 2461, 2462, 2463, 2464, 2465, 2466, 2467, 2468, 2469, 2470, 2471, 2472, 2473, 2474, 2475, 2476, 2477, 2478, 2479, 2480, 2481, 2482, 2483, 2484, 2485, 2486, 2487, 2488, 2489, 2490, 2491, 2492, 2493, 2494, 2495, 2496, 2497, 2498, 2499, 2500, 2501, 2502, 2503, 2504, 2505, 2506, 2507, 2508, 2509, 2510, 2511, 2512, 2513, 2514, 2515, 2516, 2517, 2518, 2519, 2520, 2521, 2522, 2523, 2524, 2525, 2526, 2527, 2528, 2529, 2530, 2531, 2532, 2533, 2534, 2535, 2536, 2537, 2538, 2539, 2540, 2541, 2542, 2543, 2544, 2545, 2546, 2547, 2548, 2549, 2550, 2551, 2552, 2553, 2554, 2555, 2556, 2557, 2558, 2559, 2560, 2561, 2562, 2563, 2564, 2565, 2566, 2567, 2568, 2569, 2570, 2571, 2572, 2573, 2574, 2575, 2576, 2577, 2578, 2579, 2580, 2581, 2582, 2583, 2584, 2585, 2586, 2587, 2588, 2589, 2590, 2591, 2592, 2593, 2594, 2595, 2596, 2597, 2598, 2599, 2600, 2601, 2602, 2603, 2604, 2605, 2606, 2607, 2608, 2609, 2610, 2611, 2612, 2613, 2614, 2615, 2616, 2617, 2618, 2619, 2620, 2621, 2622, 2623, 2624, 2625, 2626, 2627, 2628, 2629, 2630, 2631, 2632, 2633, 2634, 2635, 2636, 2637, 2638, 2639, 2640, 2641, 2642, 2643, 2644, 2645, 2646, 2648, 2649, 2650, 2651, 2652, 2653, 2654, 2655, 2656, 2657, 2658, 2659, 2660, 2661, 2662, 2663, 2664, 2665, 2666, 2667, 2668, 2669, 2670, 2671, 2672, 2673, 2675, 2676, 2677, 2678, 2679, 2680, 2681, 2682, 2683, 2684, 2685, 2686, 2687, 2688, 2689, 2690, 2691, 2692, 2693, 2694, 2695, 2696, 2697, 2698, 2699, 2700, 2702, 2703, 2704, 2705, 2706, 2709, 2710, 2711, 2712, 2713, 2714, 2715, 2716, 2717, 2718, 2719, 2720, 2721, 2722, 2723, 2724, 2725, 2726, 2727, 2728, 2729, 2730, 2731, 2732, 2733, 2734, 2735, 2736, 2737, 2738, 2739, 2740, 2741, 2742, 2743, 2744, 2745, 2746, 2747, 2748, 2749, 2750, 2751, 2752, 2753, 2754, 2755, 2756, 2757, 2758, 2759, 2760, 2761, 2762, 2763, 2764, 2765, 2766, 2767, 2768, 2769, 2770, 2771, 2772, 2773, 2774, 2775, 2776, 2777, 2778, 2779, 2780, 2781, 2782, 2783, 2784, 2785, 2786, 2787, 2788, 2789, 2790, 2791, 2792, 2793, 2794, 2795, 2796, 2797, 2798, 2799, 2800, 2801, 2802, 2803, 2804, 2805, 2806, 2807, 2808, 2809, 2810, 2811, 2812, 2813, 2814, 2815, 2816, 2817, 2818, 2819, 2820, 2821, 2822, 2823, 2824, 2825, 2826, 2827, 2828, 2829, 2830, 2831, 2832, 2833, 2834, 2835, 2836, 2837, 2838, 2839, 2840, 2841, 2842, 2843, 2844, 2845, 2846, 2847, 2848, 2849, 2850, 2851, 2852, 2853, 2854, 2856, 2857, 2858, 2859, 2860, 2861, 2862, 2863, 2864, 2865, 2866, 2867, 2868, 2869, 2870, 2871, 2872, 2873, 2874, 2875, 2876, 2877, 2878, 2879, 2880, 2881, 2882, 2883, 2884, 2885, 2886, 2887, 2888, 2889, 2890, 2891, 2892, 2893, 2894, 2895, 2896, 2897, 2898, 2899, 2900, 2901, 2902, 2903, 2904, 2905, 2906, 2907, 2908, 2909, 2910, 2911, 2912, 2913, 2914, 2915, 2916, 2917, 2918, 2919, 2920, 2921, 2922, 2923, 2924, 2925, 2926, 2927, 2928, 2929, 2930, 2931, 2932, 2933, 2934, 2935, 2936, 2937, 2938, 2939, 2940, 2941, 2942, 2943, 2944, 2945, 2946, 2947, 2948, 2949, 2950, 2951, 2952, 2953, 2954, 2955, 2956, 2958, 2959, 2960, 2961, 2962, 2963, 2964, 2965, 2966, 2967, 2968, 2969, 2970, 2971, 2972, 2974, 2975, 2976, 2977, 2978, 2979, 2980, 2981, 2982, 2983, 2984, 2985, 2986, 2987, 2988, 2989, 2990, 2991, 2992, 2993, 2994, 2995, 2996, 2997, 2998, 2999]

5. 3000-end
blz5 =  [3000, 3001, 3002, 3003, 3004, 3005, 3006, 3007, 3008, 3009, 3010, 3011, 3012, 3013, 3014, 3015, 3016, 3017, 3018, 3019, 3020, 3021, 3022, 3023, 3024, 3025, 3026, 3027, 3028, 3029, 3030, 3031, 3032, 3033, 3034, 3035, 3036, 3038, 3039, 3040, 3041, 3042, 3043, 3044, 3045, 3046, 3047, 3048, 3050, 3051, 3052, 3053, 3054, 3055, 3056, 3057, 3058, 3059, 3060, 3061, 3062, 3063, 3064, 3065, 3066, 3067, 3068, 3069, 3070, 3071, 3072, 3073, 3074, 3075, 3076, 3077, 3079, 3080, 3081, 3082, 3083, 3084, 3085, 3086, 3087, 3088, 3089, 3090, 3091, 3092, 3093, 3095, 3096, 3097, 3098, 3099, 3100, 3101, 3102, 3103, 3104, 3105, 3106, 3107, 3108, 3109, 3110, 3112, 3113, 3114, 3115, 3116, 3117, 3118, 3119, 3120, 3121, 3122, 3123, 3124, 3125, 3126, 3127, 3128, 3129, 3130, 3132, 3133, 3134, 3135, 3136, 3137, 3138, 3139, 3140, 3141, 3142, 3144, 3145, 3146, 3147, 3148, 3149, 3150, 3151, 3152, 3153, 3154, 3155, 3156, 3157, 3158, 3159, 3160, 3161, 3162, 3163, 3164, 3165, 3166, 3167, 3168, 3169, 3170, 3171, 3172, 3173, 3174, 3175, 3176, 3177, 3178, 3179, 3181, 3182, 3183, 3185, 3187, 3188, 3189, 3190, 3191, 3192, 3193, 3194, 3195, 3196, 3197, 3198, 3199, 3200, 3201, 3202, 3203, 3204, 3205, 3206, 3207, 3208, 3209, 3210, 3211, 3212, 3213, 3214, 3215, 3216, 3217, 3218, 3219, 3220, 3221, 3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229, 3230, 3231, 3232, 3233, 3234, 3235, 3236, 3237, 3238, 3239, 3240, 3241, 3242, 3243, 3244, 3245, 3246, 3247, 3248, 3249, 3250, 3251, 3252, 3253, 3254, 3255, 3256, 3257, 3258, 3259, 3260, 3261, 3262, 3263, 3264, 3265, 3266, 3267, 3268, 3269, 3270, 3271, 3272, 3273, 3274, 3275, 3276, 3277, 3278, 3279, 3280, 3281, 3282, 3283, 3284, 3285, 3286, 3287, 3288, 3289, 3290, 3291, 3292, 3293, 3294, 3295, 3296, 3297, 3298, 3299, 3300, 3301, 3302, 3303, 3304, 3305, 3306, 3307, 3308, 3309, 3310, 3311, 3312, 3313, 3314, 3315, 3316, 3317, 3318, 3319, 3320, 3321, 3322, 3323, 3324, 3325, 3326, 3327, 3328, 3329, 3330, 3331, 3332, 3333, 3334, 3335, 3336, 3337, 3338, 3339, 3340, 3341, 3342, 3343, 3344, 3345, 3346, 3347, 3348, 3349, 3350, 3351, 3352, 3353, 3354, 3355, 3356, 3357, 3358, 3359, 3360, 3361, 3362, 3363, 3364, 3365, 3366, 3367, 3368, 3369, 3370, 3371, 3372, 3373, 3374, 3375, 3376, 3377, 3378, 3379, 3380, 3381, 3382, 3383, 3384, 3385, 3386, 3387, 3388, 3389, 3390, 3391, 3392, 3393, 3394, 3395, 3396, 3397, 3398, 3399, 3400, 3401, 3402, 3403, 3404, 3405, 3406, 3407, 3408, 3409, 3410, 3411, 3412, 3413, 3414, 3415, 3416, 3417, 3418, 3419, 3420, 3421, 3422, 3423, 3424, 3425, 3426, 3427, 3428, 3429, 3430, 3431, 3432, 3433, 3434, 3435, 3436, 3437, 3438, 3439, 3440, 3441, 3442, 3443, 3444, 3445, 3446, 3448, 3449, 3450, 3451, 3452, 3453, 3454, 3455, 3456, 3457, 3458, 3459, 3460, 3461, 3462, 3463, 3464, 3465, 3466, 3467, 3468, 3469, 3470, 3471, 3472, 3473, 3474, 3475, 3476, 3477, 3478, 3479, 3480, 3481, 3482, 3483, 3484, 3485, 3486, 3487, 3488, 3489, 3490, 3491, 3492, 3493, 3494, 3495, 3496, 3497, 3498, 3499, 3500, 3501, 3502, 3503, 3504, 3505, 3506, 3507, 3508, 3509, 3510, 3511, 3512, 3513, 3514, 3515, 3516, 3517, 3518, 3520, 3521, 3522, 3523, 3524, 3525, 3526, 3527, 3528, 3529, 3530, 3531, 3532, 3533, 3534, 3535, 3536, 3538, 3539, 3540, 3541, 3542, 3543, 3544, 3545, 3546, 3547, 3548, 3549, 3550, 3551, 3552, 3553, 3554, 3555, 3556, 3557, 3558, 3559, 3560, 3561, 3562, 3563, 3564, 3565, 3566, 3567, 3568, 3569, 3570, 3571, 3572, 3573, 3574, 3575, 3576, 3577, 3578, 3579, 3580, 3581, 3582, 3583, 3584, 3585, 3586, 3587, 3588, 3589, 3590, 3591, 3592, 3593, 3594, 3595, 3596, 3597, 3598, 3599, 3600, 3601, 3602, 3603, 3604, 3605, 3606, 3607, 3608, 3609, 3610, 3611, 3612, 3613, 3614, 3615, 3616, 3617, 3618, 3619, 3620, 3621, 3622, 3623, 3624, 3625, 3626, 3627, 3628, 3629, 3630, 3631, 3632, 3633, 3634, 3635, 3636, 3637, 3638, 3639, 3640, 3641, 3642, 3643, 3644, 3645, 3646, 3647, 3648, 3649, 3650, 3651, 3652, 3653, 3654, 3655, 3656, 3657, 3658, 3659, 3660, 3661, 3662, 3663, 3664, 3665, 3666, 3667, 3668, 3669, 3670, 3671, 3672, 3673, 3674, 3675, 3676, 3677, 3678, 3679, 3680, 3681, 3682, 3683, 3684, 3685, 3686, 3687, 3688, 3689, 3690, 3691, 3692, 3693, 3694, 3695, 3696, 3697, 3698, 3699, 3700, 3701, 3702, 3703, 3704, 3705, 3706, 3707, 3708, 3709, 3710, 3711, 3712, 3713, 3714, 3715, 3716, 3717, 3718, 3719, 3720, 3721, 3722, 3723, 3724, 3725, 3726, 3727, 3728, 3729, 3730, 3731, 3732, 3733, 3734, 3735, 3736, 3737, 3738, 3739, 3740, 3741, 3742]

not_blzs = df.loc[~df.index.isin(blz)]
'''





'''
#Testing need to break into parts
Galactic
1. 0-300:  [2,
 5,
 7,
 8,
 9,
 11,
 14,
 15,
 17,
 18,
 19,
 20,
 21,
 22,
 24,
 25,
 26,
 30,
 31,
 32,
 33,
 34,
 35,
 36,
 37,
 38,
 40,
 41,
 42,
 44,
 45,
 46,
 47,
 48,
 49,
 51,
 54,
 56,
 57,
 58,
 59,
 60,
 64,
 65,
 66,
 67,
 69,
 70,
 71,
 72,
 73,
 74,
 76,
 78,
 81,
 84,
 89,
 94,
 99,
 100,
 107,
 109,
 111,
 113,
 115,
 116,
 117,
 118,
 120,
 121,
 124,
 126,
 128,
 129,
 130,
 131,
 134,
 136,
 137,
 140,
 142,
 143,
 147,
 149,
 151,
 152,
 154,
 156,
 158,
 159,
 167,
 173,
 178,
 181,
 184,
 188,
 196,
 200,
 207,
 209,
 216,
 217,
 218,
 223,
 232,
 244,
 245,
 251,
 266]


2. 300-450
[309,
 310,
 312,
 313,
 314,
 316,
 317,
 318,
 321,
 322,
 323,
 325,
 330,
 334,
 340,
 341,
 343,
 344,
 346,
 349,
 350,
 351,
 354,
 355,
 356,
 357,
 358,
 359,
 360,
 361,
 362,
 363,
 364,
 365,
 366,
 367,
 368,
 369,
 370,
 372,
 373,
 374,
 376,
 377,
 379,
 380,
 381,
 385,
 386,
 387,
 388,
 391,
 392,
 393,
 399,
 403,
 404,
 405,
 409,
 412,
 418,
 420,
 421,
 425,
 426,
 427,
 428,
 429,
 430,
 431,
 433,
 436,
 437,
 438,
 440,
 443,
 444,
 445,
 446,
 447,
 448]


3. 450-539:
[450,
 451,
 452,
 453,
 454,
 455,
 456,
 457,
 458,
 459,
 460,
 461,
 462,
 464,
 468,
 470,
 474,
 480,
 483,
 485,
 486,
 487,
 488,
 489,
 492,
 495,
 496,
 497,
 498,
 499,
 500,
 501,
 504,
 507,
 510,
 512,
 513,
 516,
 517,
 518,
 520,
 521,
 522,
 523,
 524,
 526,
 527,
 529,
 530,
 531,
 532,
 533,
 534,
 537,
 538]


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