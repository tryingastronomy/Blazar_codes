#Code utlized to find sources with a likelihood analysis technique, using an iterative method to eliminate systematic errors caused
#by poorly fit sources

#Will more thoroughly comment closer to graduation.

#IMPORTS
from fermipy import utils
utils.init_matplotlib_backend()
from fermipy.gtanalysis import GTAnalysis
import argparse
import numpy as np
from fermipy.castro import CastroData
import astropy.io.fits as pyfits
import os,sys
from operator import itemgetter
import math
from astropy.io import fits
from astropy.table import Table
import pandas as pd
print("")
print("Starting simulation")
print("")

usage = "usage: %(prog)s [config file]"
description = "Run fermipy analysis chain."
parser = argparse.ArgumentParser(usage=usage,description=description)
parser.add_argument('--config', default = 'config.yaml')
parser.add_argument('--source', default = None)

args = parser.parse_args()
gta = GTAnalysis(args.config,logging={'verbosity' : 3},
            fileio={'workdir_regex' : '\.xml$|\.npy$'})

gta.setup()
#File we will load every time for the simulation
gta.write_roi('prep',make_plots=False,save_model_map=False)
################## Define Functions ###################

# convert the fits catalog to dataframe
def catlog_to_df(fits_file):
    hdul = fits.open(fits_file)
    table = Table(hdul[1].data)
    names = [name for name in table.colnames if len(table[name].shape) <= 1]  # pandas can't covert multi-dimensional table to dataframe
    df = table[names].to_pandas()
    hdul.close()
    return df


### FUNCTION TO RE-OPTIMIZE SPECTRA AFTER REFINDING ###
#WHEN MANUALLY SET THE INDEX YOU NEED TO ALREADY INCLUDE THE SCALE FACTOR

def reoptsrc(names_of_srcs_list):
    #extended_4fgl = extended_list

    remove_list=[]
    SI_list = []
    SIerr_list = []
    indices=np.linspace(-0.5,-4.5,9)

    for each in names_of_srcs_list:
        source_name = each
        gta.free_sources(free=False)
        gta.free_source("galdiff")
        gta.free_source("isodiff")
        gta.free_source("isodiff",free=False,pars=["BreakValue","Index1","Index2"])
        gta.free_source(source_name)
        fit = gta.fit()
        # Write fit results:
        status = fit['fit_status']
        quality = fit['fit_quality']

        if status != 0 or quality != 3:
            conv = 0
        else:
            conv = 1
        # Change starting index value if source does not converge
        # Increase bounds of index value [-2,5]
        for ind in indices:
            if conv == 0:
                print('Testing with starting index', ind)
                gta.get_free_source_params(source_name)
                gta.set_parameter(source_name,par='Index',value=ind,scale=-1.0,bounds=[-5.,5.])
                fit = gta.fit()
                status  = fit['fit_status']
                quality = fit['fit_quality']
                if status == 0 and quality == 3:
                    conv+=1
                else:
                    continue
            else:
                break

        if conv==0:
            gta.set_parameter(source_name,par='Index',value=-2.0,scale=-1.0,bounds=[-5.,5.])
            gta.get_free_source_params(source_name)
            gta.free_source(source_name,free=False,pars=["Index"])
            gta.fit()
            remove_list.append(source_name)

        this_result = "quality: " + str(quality) + " status: " + str(status) + "\n"
        #I can't figure out how to access information in gta - so 
        gta.write_roi('refit')
        c = np.load('refit.npy',allow_pickle=True).flat[0]
        #just refitting the one src so we just want the 
        new_ind = c['sources'][source_name]['spectral_pars']['Index']['value']
        new_inderr = c['sources'][source_name]['spectral_pars']['Index']['error']
        SI_list.append(new_ind)
        SIerr_list.append(new_inderr)

    #also need to record and report the error of the spectral index valuei
    #include in output and then can collect these
    return(remove_list, SI_list, SIerr_list) 




### GLOBAL FIT FUNCTION ###
def global_fit(names_of_srcs_list,extended_list):
    # Free norm for all sources together, but only sources from the current TS iteration:
    #Should really only be 1 src if a refit of a previously detected source - but could have more...
    # Also refit Index1 of the isotropic
    gta.free_sources(free=False)
    #gta.free_sources(pars='norm',free=True,exclude=exclude_list) # instead of doing this we could just free all sources in the names of srcs list
    for every in names_of_srcs_list:
        gta.free_source(every, pars = ['norm'])
    gta.free_source("galdiff")
    gta.free_source("isodiff")
    gta.free_source("isodiff",free=False,pars=["BreakValue", "Index2"])
    #Don't free the Sun and moon parameters
    gta.free_source("diffuse00", free=False)
    gta.free_source("diffuse01", free=False)
    fit = gta.fit()

    # Free norm for all sources together aside 4FGL extended:
    # Also refit Index1 of the isotropic
    gta.free_sources(free=False)
    gta.free_sources(pars='norm',free=True, exclude = extended_4fgl)
    gta.free_source("galdiff")
    gta.free_source("isodiff")
    gta.free_source("isodiff",free=False,pars=["BreakValue", "Index2"])
    #Don't free the Sun and moon parameters
    gta.free_source("diffuse00", free=False)
    gta.free_source("diffuse01", free=False)
    fitglobal = gta.fit()
    return fitglobal #names_of_srcs_list, source_list # -> no return, just have the fit occur





######################################################

# Define main source list:
src_list = []
extended_4fgl = []


gta.print_model()
#Load in prefit
#gta.load_roi('sim_start_fit') #model saved with all 4FGL srcs
#simulate ROI
#Delete Sources-> then find them!
lat_4fgl = "/zfs/astrohe/sjoffre/Fermi/Catalogs/gll_psc_v32.fit"
hdu = fits.open(lat_4fgl)
data = hdu[1].data
cut_list = []
#wait might have problem that required .xml file

for each in data["Source_Name"]:
    if "e" not in each:
        cut_list.append(each)
        gta.delete_source(each)
    else:
        src_list.append(each)
        #extended_4fgl.append(each) -> potential problem from before
# Delete all 4fgl source except extended sources:
'''
hdu = fits.open(lat_4fgl)
data = hdu[1].data
cut_list = []
for each in data["Source_Name"]:
    if "e" not in each:
        cut_list.append(each)
'''
#gta.delete_sources(names=cut_list)


#Write XML file
#I can't remember why I had to use the XML file... - saving issue with new Extended sources
gta.write_xml('prefit_roi')

# Change spectral model of isodiff to broken power law:
gta.set_source_spectrum(name="isodiff", spectrum_type="BrokenPowerLaw")
gta.set_parameter(name="isodiff",par="Index1",scale=1,value=0,bounds=[-15,15])
gta.set_parameter(name="isodiff",par="Index2",scale=1,value=0,bounds=[-3,3])
gta.set_parameter(name="isodiff",par="Prefactor",scale=1,value=1,bounds=[0.1,5])
gta.set_parameter(name="isodiff",par="BreakValue",scale=1,value=30,bounds=[20,100])

# Initial fit for diffuse sources:
gta.free_sources(free=False)
gta.free_source("galdiff")
gta.free_source("isodiff")
gta.free_source("isodiff",free=False,pars=["BreakValue", "Index2"])
gta.fit()
gta.print_model()

#gta.write_roi('sim_start_fit',make_plots=False,save_model_map=False)
#gta.write_fits('fit_0_test')
gta.write_xml('fit_0')
#gta.tsmap('TS_map_start', make_plots=True)

#since we have deleted out the non-extended sources, we grab the remaining 4FGL srcs
#which will have to be 4FGL extended sources
sources_from_xml = []
with open('fit_0_00.xml') as g:
    for line in g:
        if line.find("source name=") !=-1:
            print(line)
            word_count = 0
            for j in range(len(line)):
                if word_count==1:
                    sources_from_xml.append(word)
                    break
                elif line[j] == "\"":
                    p=j+1
                    word=""
                    word_count = 1
                    while(line[p]!= "\""):
                            word = word + line[p]
                            p += 1
#Save the extended srcs, including to running src_list
ext_names = []
#extended_4fgl = []
for each in sources_from_xml:
    if '4FGL' in each:
        src_list.append(each)
        extended_4fgl.append(each)
print(extended_4fgl)
#--- INITIAL ITERATION STEPS TO FIND ROUGH POSITION ---#

gta.free_sources(free=False)
gta.free_source("galdiff")
gta.free_source("isodiff")
gta.free_source("isodiff",free=False,pars=["BreakValue", "Index2"])
gta.fit()
gta.print_model()

#NOTE: although labeled all 5 - will be TSsqrt4!!!
sqrt_TS = [8,6,4]
min_sep_ang = [0.5,0.5,0.5]
f2  = open("running_loglike.txt","w")
if os.path.isfile('ellipse_file.txt')==False:
    files_for_ellipse = open('ellipse_file.txt', 'w+')
    files_for_ellipse.write("#Source_name glon_ell[deg] glat_ell[deg] \n")
    files_for_ellipse.close()
else:
    files_for_ellipse = open('ellipse_file.txt', 'a')
    files_for_ellipse.write("#Source_name glon_ell[deg] glat_ell[deg] \n")
    files_for_ellipse.close()

#ITERATE Through all 3 TS levels
#Do a simple individual fit for each source - in-depth fit after refining position/spectrum of brightest TS
for i in range(0, len(sqrt_TS)):
    f  = open("TS_sqrt" + str(sqrt_TS[i])+"_initial_convergence.txt","w") #this is probably all screwed up later on
    #gta.free_sources(free=False)
    model = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
    sdict = gta.find_sources(model=model, sqrt_ts_threshold=sqrt_TS[i],min_separation=min_sep_ang[i], max_iter=10,sources_per_iter=20,tsmap_fitter='tsmap')
    #SAVE NAMES FOR EACH ITERATION
    if i == 0:
        name_listsrcfind8   = np.array([sdict['sources'][j]['name'] for j in range(0,len(sdict['sources']))])
        ts_listsrcfind8 = np.array([sdict['sources'][j]['ts'] for j in range(0,len(sdict['sources']))])
        G_listsrcfind8 = np.array([sdict['sources'][j]['spectral_pars']['Index']['value'] for j in range(0,len(sdict['sources']))])
        for each8 in name_listsrcfind8:
            gta.free_sources(free=False)
            gta.free_source("galdiff")
            gta.free_source("isodiff")
            gta.free_source("isodiff",free=False,pars=["BreakValue","Index1","Index2"])
            gta.free_source(each8)
            fit = gta.fit()
        print('name_listsrcfind8')
        print(name_listsrcfind8)
        #WRITE LOGLIKE OUTSIDE OF FOR LOOP
        '''
        ll = fit['loglike']
        dll = fit['dloglike']
        f2.write('TSsqrt8')
        f2.write("\n")
        f2.write(str(ll))
        f2.write("\n")
        f2.write(str(dll))
        f2.write("\n")
        '''

    elif i == 1:
        name_listsrcfind6   = np.array([sdict['sources'][j]['name'] for j in range(0,len(sdict['sources']))])
        ts_listsrcfind6 = np.array([sdict['sources'][j]['ts'] for j in range(0,len(sdict['sources']))])
        G_listsrcfind6 = np.array([sdict['sources'][j]['spectral_pars']['Index']['value'] for j in range(0,len(sdict['sources']))])
        for each6 in name_listsrcfind6:
            gta.free_sources(free=False)
            gta.free_source("galdiff")
            gta.free_source("isodiff")
            gta.free_source("isodiff",free=False,pars=["BreakValue","Index1","Index2"])
            gta.free_source(each6)
            fit = gta.fit()

        print('name_listsrcfind6')
        print(name_listsrcfind6)
        '''
        ll = fit['loglike']
        dll = fit['dloglike']
        f2.write('TSsqrt6')
        f2.write("\n")
        f2.write(str(ll))
        f2.write("\n")
        f2.write(str(dll))
        f2.write("\n")
        '''

    elif i ==2:
        name_listsrcfind5   = np.array([sdict['sources'][j]['name'] for j in range(0,len(sdict['sources']))])
        ts_listsrcfind5 = np.array([sdict['sources'][j]['ts'] for j in range(0,len(sdict['sources']))])
        G_listsrcfind5 = np.array([sdict['sources'][j]['spectral_pars']['Index']['value'] for j in range(0,len(sdict['sources']))])
        for each5 in name_listsrcfind5:
            gta.free_sources(free=False)
            gta.free_source("galdiff")
            gta.free_source("isodiff")
            gta.free_source("isodiff",free=False,pars=["BreakValue","Index1","Index2"])
            gta.free_source(each5)
            fit = gta.fit()
        print('name_listsrcfind5')
        print(name_listsrcfind5)
        '''
        ll = fit['loglike']
        dll = fit['dloglike']
        f2.write('TSsqrt5')
        f2.write("\n")
        f2.write(str(ll))
        f2.write("\n")
        f2.write(str(dll))
        f2.write("\n")
        '''
    '''
    status = fit['fit_status']
    quality = fit['fit_quality']
    write_list = [quality,status]
    f.write("\n")
    f.write('TSsqrt'+str(i))
    f.write("\n")
    f.write(str(write_list))
    f.write("\n")
    '''
    f.close()
    ###FIRST SRC FIND SAVED FILE###
    gta.write_roi('after_TS_sqrt'+str(sqrt_TS[i]), make_plots=False, save_model_map=False)
    #gta.tsmap('TS_map_Final_after_TS_sqrt'+str(sqrt_TS[i]), make_plots=True)
    gta.print_model()


### GLOBAL FIT HERE ### (no need to call function if only extended sources are in src_list
# Free norm for all sources together, but only sources from the current TS iteration:
# Also refit Index1 of the isotropic
gta.free_sources(free=False)
gta.free_sources(pars='norm',free=True,exclude=src_list) # at this point in time the only sources in src_list are the extended sources
gta.free_source("galdiff")
gta.free_source("isodiff")
gta.free_source("isodiff",free=False,pars=["BreakValue", "Index2"])
#Don't free the Sun and moon parameters
gta.free_source("diffuse00", free=False)
gta.free_source("diffuse01", free=False)
fit = gta.fit() ### -> HERE IS THE GLOBAL FIT after initial finding all sources
ll = fit['loglike']
dll = fit['dloglike']
f2.write('GF1')
f2.write("\n")
f2.write(str(ll))
f2.write("\n")
f2.write(str(dll))
f2.write("\n")

gta.write_roi('global_fit_1', make_plots=False, save_model_map=False)

# Write fit results:
status = fit['fit_status']
quality = fit['fit_quality']
#this_result = "quality: " + str(quality) + " status: " + str(status) + "\n"
write_list = [quality,status]
g  = open("global_1_convergence.txt","w")
g.write("\n")
g.write(str(write_list))
g.write("\n")
g.close()


#--PREPARE FOR RELOCALIZATION--#
#delete sources that have TS<25 from model

#-Update TS values from global fit-#
ts25_name=[]
ts25_index=[]
ts25_ts=[]
ts36_name=[]
ts36_index=[]
ts36_ts=[]
ts64_name=[]
ts64_index=[]
ts64_ts=[]

for u in range(len(gta.roi.sources)-4): #4 for the 4 diffuse components
    if gta.roi[gta.roi.sources[u].name]['Source_Name'] in name_listsrcfind5:
        #
        ts25_name.append(gta.roi[gta.roi.sources[u].name]['Source_Name'])
        ts25_ts.append(gta.roi[gta.roi.sources[u].name]['ts'])
        ts25_index.append(gta.roi[gta.roi.sources[u].name]['spectral_pars']['Index']['value'])
    elif gta.roi[gta.roi.sources[u].name]['Source_Name'] in name_listsrcfind6:
        ts36_name.append(gta.roi[gta.roi.sources[u].name]['Source_Name'])
        ts36_ts.append(gta.roi[gta.roi.sources[u].name]['ts'])
        ts36_index.append(gta.roi[gta.roi.sources[u].name]['spectral_pars']['Index']['value'])
    elif gta.roi[gta.roi.sources[u].name]['Source_Name'] in name_listsrcfind8:
        ts64_name.append(gta.roi[gta.roi.sources[u].name]['Source_Name'])
        ts64_ts.append(gta.roi[gta.roi.sources[u].name]['ts'])
        ts64_index.append(gta.roi[gta.roi.sources[u].name]['spectral_pars']['Index']['value'])



#to make clean - write function - for now don't mess with it. 
TS5dict = {'Name' : ts25_name, 'TS' : ts25_ts, 'Index' : ts25_index}
dfts5 = pd.DataFrame(data=TS5dict)
sorted_dfts5 = dfts5.sort_values(by='TS', ascending=False)
under_ts25 = sorted_dfts5[sorted_dfts5['TS']<25].reset_index(drop=True) #these we delete from the model
sorted_dfts5 = sorted_dfts5[sorted_dfts5['TS']>25].reset_index(drop=True) #these we call for the loops below
for x in range(len(under_ts25)):
    gta.delete_source(under_ts25['Name'][x])


#can also be from TSsqrt6 srcs so go through those as well
TS6dict = {'Name' : ts36_name, 'TS' : ts36_ts, 'Index' : ts36_index}
dfts6 = pd.DataFrame(data=TS6dict)
sorted_dfts61 = dfts6.sort_values(by='TS', ascending=False)
under_ts25_6 = sorted_dfts61[sorted_dfts61['TS']<25].reset_index(drop=True)
sorted_dfts6 = sorted_dfts61[sorted_dfts61['TS']>25].reset_index(drop=True)
#print('sorted_dfts61')
#print(sorted_dfts61)
#print('sorted_dfts6')
#print(sorted_dfts6['Name'])
#print('under_ts25_6')
#print(under_ts25_6['Name'])
if len(under_ts25_6['Name'])>0:
    for y in range(len(under_ts25_6)):
        #oh wtf, since index of the src is 1, does not exist...
        gta.delete_source(under_ts25_6['Name'][y])


TS8dict = {'Name' : ts64_name, 'TS' : ts64_ts, 'Index' : ts64_index}
dfts8 = pd.DataFrame(data=TS8dict)
sorted_dfts81 = dfts8.sort_values(by='TS', ascending=False)
under_ts25_8 = sorted_dfts81[sorted_dfts81['TS']<25].reset_index(drop=True) #these we delete from the model
sorted_dfts8 = sorted_dfts81[sorted_dfts81['TS']>25].reset_index(drop=True) #these we call for the loops below
for z in range(len(under_ts25_8)):
    gta.delete_source(under_ts25_8['Name'][z])

## --- RELOCALIZATION AND RE-OPTIMIZATION --- ##
#make sure to capture spectral index error values
#Update - model should use the src spectral index. 
f1  = open("running_convergence.txt","w")
remove_list8 = []
src_list8 = []
ts_list8 = []
G_list8 = []
Gerr_list8 = []
for n in range(len(sorted_dfts8)): #Adapted to be in descending order name_listsrcfind8:
    source_name = sorted_dfts8['Name'][n]
    #since not indexing grab the PL Index value from the dataframe with each
    #new_ind = sorted_dfts8.loc[sorted_dfts8['Name']==source_name].values[0][2]
    new_ind = sorted_dfts8['Index'][n]
    model = {'Index' : new_ind, 'SpatialModel' : 'PointSource'}
    #Delete only 1 at a time, starting with the brightest
    gta.delete_source(source_name)
    print('DELETE DETECTED SOURCE')
    #same as previous fit freeing
    gta.free_sources(free=False)
    gta.free_source("galdiff")
    gta.free_source("isodiff")
    gta.free_source("isodiff",free=False,pars=["BreakValue","Index1","Index2"])
    #model = {'Index' : 2.0, 'SpatialModel' : 'PointSource'} 
    ####SDICT_8####
    #change to only find 1 src with 1 iteration
    sdict_8 = gta.find_sources(model=model, sqrt_ts_threshold=8,min_separation=min_sep_ang[i], max_iter=1,sources_per_iter=1,tsmap_fitter='tsmap')
    #add warning if find more than 1#

    #NEED TO ALSO RECORD NAMES#
    name_listsrcfind8  = np.array([sdict_8['sources'][j]['name'] for j in range(0,len(sdict_8['sources']))])
    ts_listsrcfind8 = np.array([sdict_8['sources'][j]['ts'] for j in range(0,len(sdict_8['sources']))])
    G_listsrcfind8 = np.array([sdict_8['sources'][j]['spectral_pars']['Index']['value'] for j in range(0,len(sdict_8['sources']))])
    Gerr_listsrcfind8 = np.array([sdict_8['sources'][j]['spectral_pars']['Index']['error'] for j in range(0,len(sdict_8['sources']))])
    '''
    if len(name_listsrcfind8)>1:
        print('MORE THAN ONE SRC DETECTED FOR A SINGLE DELETION')
    for w in range(len(name_listsrcfind8)):
        src_list8.append(name_listsrcfind8[w])
        ts_list8.append(ts_listsrcfind8[w])
        G_list8.append(G_listsrcfind8[w])
        Gerr_list8.append(Gerr_listsrcfind8[w])
    '''
    #DO SOMETHING SIMILAR FOR RA/DEC POSR95 AS I DO FOR NAMES AND TS AS OPPOSED TO DOING THE LOOPS BELOW
    #REOPTSRC returns the remove list -> but we are only finding one at a time so you need to store them
    remove_8_1 = reoptsrc(name_listsrcfind8) # 0 is remove list, 1 is SI, 2 is SIerr
    if len(remove_8_1[0]) >0:
        remove_list8 = remove_list8 + remove_8_1[0]
    #remove_list8 = remove_list8 + remove_8_1[0] #incrementally add each one that fails to converge after repotimixing

    if len(name_listsrcfind8)>1:
        print('MORE THAN ONE SRC DETECTED FOR A SINGLE DELETION')
    #Ok, so add to the lists that we will make dataframes
    #which we will save so we have the index error values
    for w in range(len(name_listsrcfind8)):
        src_list8.append(name_listsrcfind8[w])
        ts_list8.append(ts_listsrcfind8[w])
        G_list8.append(remove_8_1[1][0])
        Gerr_list8.append(remove_8_1[2][0])
        

#GLOBAL FIT# -> out of TSsqrt8 for loop, need to input all sources detected above
gf = global_fit(src_list8,extended_4fgl)

gta.write_roi('post_8refit', make_plots=False, save_model_map=False)
#gta.tsmap('post_8refitTSmap', make_plots=True)


status = gf['fit_status']
quality = gf['fit_quality']
write_list = [quality,status]
f1.write("\n")
f1.write('GF 8 ')
f1.write("\n")
f1.write(str(write_list))
f1.write("\n")

ll = gf['loglike']
dll = gf['dloglike']
f2.write('GF 8 ')
f2.write("\n")
f2.write(str(ll))
f2.write("\n")
f2.write(str(dll))
f2.write("\n")
#now after global fits, we remove the sources that did not converge - hopefully better capture in next TS step
#src_list8 = src_list8.tolist() # -> no longer an array so don't need
'''
for each in src_list8:
    if each in remove_list8:
        src_list8.remove(each)
src_list = src_list + src_list8
'''

## FIND SOURCES IN TSsqrt6 ##
remove_list6=[]
src_list6=[]
ts_list6=[]
G_list6 = []
Gerr_list6 = []

if len(sorted_dfts6)>0:
    #Refit inidividual sources if sources still survive TS>25
    for m in range(len(sorted_dfts6)):
        source_name = sorted_dfts6['Name'][m]
        #new_ind = sorted_dfts6.loc[sorted_dfts6['Name']==source_name].values[0][2]
        new_ind = sorted_dfts6['Index'][m]
        model = {'Index' : new_ind, 'SpatialModel' : 'PointSource'}
        gta.delete_source(source_name)
        print('DELETE DETECTED SOURCE')
        gta.free_sources(free=False)
        gta.free_source("galdiff")
        gta.free_source("isodiff")
        gta.free_source("isodiff",free=False,pars=["BreakValue","Index1","Index2"])
        sdict_6 = gta.find_sources(model=model, sqrt_ts_threshold=6,min_separation=min_sep_ang[i], max_iter=1,sources_per_iter=1,tsmap_fitter='tsmap')
        # Add warning if find more than 1#
        name_listsrcfind6 = np.array([sdict_6['sources'][j]['name'] for j in range(0,len(sdict_6['sources']))])
        ts_listsrcfind6 = np.array([sdict_6['sources'][j]['ts'] for j in range(0,len(sdict_6['sources']))])
        G_listsrcfind6 = np.array([sdict_6['sources'][j]['spectral_pars']['Index']['value'] for j in range(0,len(sdict_6['sources']))])
        Gerr_listsrcfind6 = np.array([sdict_6['sources'][j]['spectral_pars']['Index']['error'] for j in range(0,len(sdict_6['sources']))])
        #INCLUDE OTHERS and save so we can see if improving...
        remove_6_1 = reoptsrc(name_listsrcfind6)
        if len(remove_6_1[0])>0:
            remove_list6.append(remove_6_1[0])
        #remove_list6.append(remove_6_1[0])

        if len(name_listsrcfind6)>1:
            print('MORE THAN ONE SRC DETECTED FOR A SINGLE DELETION')
        for w in range(len(name_listsrcfind6)):
            src_list6.append(name_listsrcfind6[w])
            ts_list6.append(ts_listsrcfind6[w])
            G_list6.append(remove_6_1[1][w])
            Gerr_list6.append(remove_6_1[2][w]) #since tuple need index

else:
    ### if no srcs in first fit are at TS>25, find, and re-optimize
    model = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
    sdict = gta.find_sources(model=model, sqrt_ts_threshold=6,min_separation=0.5, max_iter=10,sources_per_iter=20,tsmap_fitter='tsmap')
    
    #This time don't delete out TS<25 srcs
    #organize by TS and reoptimize
    name_listsrcfind6   = np.array([sdict['sources'][j]['name'] for j in range(0,len(sdict['sources']))])
    ts_listsrcfind6 = np.array([sdict['sources'][j]['ts'] for j in range(0,len(sdict['sources']))])
    G_listsrcfind6 = np.array([sdict['sources'][j]['spectral_pars']['Index']['value'] for j in range(0,len(sdict['sources']))])
    G_errlistsrcfind6 = np.array([sdict['sources'][j]['spectral_pars']['Index']['error'] for j in range(0,len(sdict['sources']))])
    for each6 in name_listsrcfind6:
        gta.free_sources(free=False)
        gta.free_source("galdiff")
        gta.free_source("isodiff")
        gta.free_source("isodiff",free=False,pars=["BreakValue","Index1","Index2"])
        gta.free_source(each6)
        fit = gta.fit()
    
    #don't remove TS<25 srcs
    #hopefully my relocalization method works ok for all new discoveries
    #don't bother with re-ordering from highest TS
    remove_6_1 = reoptsrc(name_listsrcfind6)
    if len(remove_6_1[0])>0:
        remove_list6.append(remove_6_1[0])

    for w in range(len(name_listsrcfind6)):
        src_list6.append(name_listsrcfind6[w])
        ts_list6.append(ts_listsrcfind6[w])
        G_list6.append(remove_6_1[1][w])
        Gerr_list6.append(remove_6_1[2][w])
    
    TS6dict = {'Name' : src_list6, 'TS' : ts_list6, 'Index' : G_list6, 'Indexerr' : Gerr_list6}
    dfts6 = pd.DataFrame(data=TS6dict)
    sorted_dfts6 = dfts6.sort_values(by='TS', ascending=False)


#GLOBAL FIT#
gf = global_fit(src_list6, extended_4fgl)
status = gf['fit_status']
quality = gf['fit_quality']
write_list = [quality,status]
f1.write("\n")
f1.write('GF 6 ')
f1.write("\n")
f1.write(str(write_list))
f1.write("\n")

ll = gf['loglike']
dll = gf['dloglike']
f2.write('GF 6')
f2.write("\n")
f2.write(str(ll))
f2.write("\n")
f2.write(str(dll))
f2.write("\n")

gta.write_roi('post_6refit', make_plots=False, save_model_map=False)
#gta.tsmap('post_6refitTSmap', make_plots=True)

### FIND SOURCES TS SQRT4 ###
remove_list5 = []
src_list5 = []
ts_list5 = []
G_list5 = []
Gerr_list5 = []
if len(sorted_dfts5)>0:
    #IF SOURCES IN FIRST ROUND SURVIVE TS>25
    for o in range(len(sorted_dfts5)):
        source_name = sorted_dfts5['Name'][o]
        new_ind = sorted_dfts5['Index'][o]
        #new_ind = sorted_dfts5.loc[sorted_dfts5['Name']==source_name].values[0][2]
        model = {'Index' : new_ind, 'SpatialModel' : 'PointSource'}
        gta.delete_source(source_name)
        print('DELETE DETECTED SOURCE')
        gta.free_sources(free=False)
        gta.free_source("galdiff")
        gta.free_source("isodiff")
        gta.free_source("isodiff",free=False,pars=["BreakValue","Index1","Index2"])
        #Maybe do something where find TSsqrt5 again, and then go through and do extra excesses, but for this let's go back down to TSsqrt4 (for the refit) and have the 10 iterations and 20 srcs per iteration since the issue seemed to be TSsqrt4 initially. 
        #what we also can do is do TSsqrt 8 and 6 alone -> refit those, and then do TSsqrt5 and TSsqrt4 -> refit those, and then do it all one last time
        #Although say 5, deciding to do 4 (no excess)

        sdict_5 = gta.find_sources(model=model, sqrt_ts_threshold=4,min_separation=min_sep_ang[i], max_iter=1,sources_per_iter=1,tsmap_fitter='tsmap')
        # Add warning if find more than 1#
        name_listsrcfind5 = np.array([sdict_5['sources'][j]['name'] for j in range(0,len(sdict_5['sources']))])
        ts_listsrcfind5 = np.array([sdict_5['sources'][j]['ts'] for j in range(0,len(sdict_5['sources']))])
        G_listsrcfind5 = np.array([sdict_5['sources'][j]['spectral_pars']['Index']['value'] for j in range(0,len(sdict_5['sources']))])
        Gerr_listsrcfind5 = np.array([sdict_5['sources'][j]['spectral_pars']['Index']['error'] for j in range(0,len(sdict_5['sources']))])
        
        remove_5_1 = reoptsrc(name_listsrcfind5)
        if len(remove_5_1[0])>0:
            remove_list5.append(remove_5_1)
        #BECAUSE WE ARE DOING 1 src at a time, but can find more than 1
        if len(name_listsrcfind5)>0:
            print('MORE THAN ONE SRC DETECTED FOR A SINGLE DELETION')
        for w in range(len(name_listsrcfind5)):
            src_list5.append(name_listsrcfind5[w])
            ts_list5.append(ts_listsrcfind5[w])
            G_list5.append(remove_5_1[1][w]) #the wth src will still be w in order of detected - instead of just 0, can be more
            Gerr_list5.append(remove_5_1[2][w])
else:
    model = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
    sdict = gta.find_sources(model=model, sqrt_ts_threshold=4,min_separation=0.5, max_iter=10,sources_per_iter=20,tsmap_fitter='tsmap')

    #This time don't delete out TS<25 srcs
    #organize by TS and reoptimize
    name_listsrcfind5   = np.array([sdict['sources'][j]['name'] for j in range(0,len(sdict['sources']))])
    ts_listsrcfind5 = np.array([sdict['sources'][j]['ts'] for j in range(0,len(sdict['sources']))])
    G_listsrcfind5 = np.array([sdict['sources'][j]['spectral_pars']['Index']['value'] for j in range(0,len(sdict['sources']))])
    Gerr_listsrcfind5 = np.array([sdict['sources'][j]['spectral_pars']['Index']['error'] for j in range(0,len(sdict['sources']))])
    for each5 in name_listsrcfind5:
        gta.free_sources(free=False)
        gta.free_source("galdiff")
        gta.free_source("isodiff")
        gta.free_source("isodiff",free=False,pars=["BreakValue","Index1","Index2"])
        gta.free_source(each5)
        fit = gta.fit()

    #don't remove TS<25 srcs
    #hopefully my relocalization method works ok for all new discoveries
    #don't bother with re-ordering from highest TS
    remove_5_1 = reoptsrc(name_listsrcfind5)
    if len(remove_5_1[0])>0:
        remove_list5.append(remove_5_1[0])

    for w in range(len(name_listsrcfind5)):
        src_list5.append(name_listsrcfind5[w])
        ts_list5.append(ts_listsrcfind5[w])
        G_list5.append(remove_5_1[1][w])
        Gerr_list5.append(remove_5_1[2][w])

    TS5dict = {'Name' : src_list5, 'TS' : ts_list5, 'Index' : G_list5, 'Indexerr' : Gerr_list5}
    dfts5 = pd.DataFrame(data=TS5dict)
    sorted_dfts5 = dfts5.sort_values(by='TS', ascending=False)

##WHAT TO DO IF DETECT MORE THAN ONE SRC??? -> will be same as the w value


# Final combination so we can get the SI and SIerr values
name_tot = src_list8+src_list6+src_list5
ts_tot = ts_list8 + ts_list6 + ts_list5
G_tot = G_list8 + G_list6 + G_list5
Gerr_tot = Gerr_list8 + Gerr_list6 + Gerr_list5

TS_final_dict = {'Name' : name_tot, 'TS' : ts_tot , 'Index' : G_tot, 'Indexerr' : Gerr_tot}
dfts_final = pd.DataFrame(data=TS_final_dict)
dfts_final.to_csv('spectral_values.csv')





#GLOBAL FIT#
gf = global_fit(src_list5, extended_4fgl)
status = gf['fit_status']
quality = gf['fit_quality']
write_list = [quality,status]
f1.write("\n")
f1.write('GF 5 ')
f1.write("\n")
f1.write(str(write_list))
f1.write("\n")

ll = gf['loglike']
dll = gf['dloglike']
f2.write('GF 5')
f2.write("\n")
f2.write(str(ll))
f2.write("\n")
f2.write(str(dll))
f2.write("\n")
'''
for each in src_list5:
    if each in remove_list5:
        src_list5.remove(each)
src_list = src_list + src_list5
'''
### GLOBAL FIT HERE ### -> all srcs not just recent ones (not calling the function since will just exclude extended
# Free norm for all sources together, but only sources from the current TS iteration:
# Also refit Index1 of the isotropic
gta.free_sources(free=False)
gta.free_sources(pars='norm',free=True,exclude=extended_4fgl)
gta.free_source("galdiff")
gta.free_source("isodiff")
gta.free_source("isodiff",free=False,pars=["BreakValue", "Index2"])
#Don't free the Sun and moon parameters
gta.free_source("diffuse00", free=False)
gta.free_source("diffuse01", free=False)
gta.fit() ### -> HERE IS THE GLOBAL FIT after initial finding all sources
gta.write_roi('global_fit_2', make_plots=False, save_model_map=False)
#gta.tsmap('global_fit_2', make_plots=True)


status = gf['fit_status']
quality = gf['fit_quality']
write_list = [quality,status]
f1.write("\n")
f1.write('GF Last ')
f1.write("\n")
f1.write(str(write_list))
f1.write("\n")

ll = gf['loglike']
dll = gf['dloglike']
f2.write('GF last')
f2.write("\n")
f2.write(str(ll))
f2.write("\n")
f2.write(str(dll))
f2.write("\n")



###Still need to implement something for refitting the sources that did not converge ###

#it appears we need to also record this potentially new spectral index and spectral index with error values
# solution - data frame of the values? 

srcs_to_refit = remove_list5 + remove_list6 + remove_list8

print(srcs_to_refit)

#if any of them are tuples, convert to lists
for q in range(len(srcs_to_refit)):
    if type(srcs_to_refit[q]) == tuple:
        srcs_to_refit[q] = list(srcs_to_refit[q])

#now we have list of lists
for r in range(len(srcs_to_refit)):
    if type(srcs_to_refit[r]) == list:
        srcs_to_refit[r] = srcs_to_refit[r][0]

#now make a list of strings
#srcs_to_refit = list(map(''.join,srcs_to_refit))

#since last step can have multiple sources need to flatten the list
#srcs_to_refit = sum(srcs_to_refit,[])
print('Sources to be refit')
#print(srcs_to_refit)
remove_list_end = []
src_names_rem =[]
si_refit = []
sierr_refit = []
#for j in range(len(gta.roi.sources[:])):


#Need to distinguish between list and list of lists - i can't find an easy way to do this
#if list of lists when call for len(srcs_to_refit[0]) unlikely to be 15 srcs (although could be
#fuck - roi43 has mix of list of a list and a single src...
fixed_refit_list = []
for each in srcs_to_refit:
    if len(each) == 15:
        fixed_refit_list.append(each)
    else:
        for every in each:
            fixed_refit_list.append(every)
    #len_src_refit = len(srcs_to_refit[0])

print(fixed_refit_list)
#now set to original name so I don't have to change the rest o fthe code 
srcs_to_refit = fixed_refit_list

indices=np.linspace(-0.5,-4.5,9)
for eachredo in srcs_to_refit:
    #names = gta.roi.sources[j].name
    gta.free_sources(free=False)
    gta.free_source("galdiff")
    gta.free_source("isodiff")
    gta.free_source("isodiff",free=False,pars=["BreakValue","Index1","Index2"])
    gta.free_source(eachredo)

    fit = gta.fit()

    status = fit['fit_status']
    quality = fit['fit_quality']


    if status != 0 or quality != 3:
        conv = 0
    else:
        conv = 1

    # Change starting index value if source does not converge
    # Increase bounds of index value [-5,5]


    for ind in indices:
        print('Testing with starting index', ind)
        if conv == 0:
            gta.get_free_source_params(eachredo)
            gta.set_parameter(eachredo,par='Index',value=ind,scale=-1.0,bounds=[-5.,5.])

            fit = gta.fit()

            status  = fit['fit_status']
            quality = fit['fit_quality']

            if status == 0 and quality == 3:
                conv+=1
            else:
                continue
        else:
            break
    if conv==0:
        gta.set_parameter(eachredo,par='Index',value=-2.0,scale=-1.0,bounds=[-5.,5.])
        gta.get_free_source_params(eachredo)
        gta.free_source(eachredo,free=False,pars=["Index"])
        gta.fit()
        remove_list_end.append(eachredo)

    gta.write_roi('refit_nonconvg',make_plots=False,save_model_map=False)
    rnc = np.load('refit_nonconvg.npy',allow_pickle=True).flat[0]
    new_ind = rnc['sources'][eachredo]['spectral_pars']['Index']['value']
    new_inderr = rnc['sources'][eachredo]['spectral_pars']['Index']['error']
    si_refit.append(new_ind)
    sierr_refit.append(new_inderr)

refit_dict = {'Name' : srcs_to_refit, 'Index' : si_refit, 'Indexerr' : sierr_refit}
dfts_final = pd.DataFrame(data=refit_dict)
dfts_final.to_csv('nonconvg_spectral_refit.csv')

    
#### Count Number ##### 
f = open('non_convergent.txt', 'w+')

how_many_end = len(remove_list_end)

if how_many_end!=0:
    for name in remove_list_end:
        f.write(name+"\n")
f.close()

### FINAL GLOBAL FIT ###
if how_many_end !=0:
    gta.free_sources(free=False)
    gta.free_sources(pars='norm',free=True,exclude=extended_4fgl) # at this point in time the only sources in src_list are the extended sources
    gta.free_source("galdiff")
    gta.free_source("isodiff")
    gta.free_source("isodiff",free=False,pars=["BreakValue", "Index2"])
    #Don't free the Sun and moon parameters
    gta.free_source("diffuse00", free=False)
    gta.free_source("diffuse01", free=False)
    fit = gta.fit() ### -> HERE IS THE GLOBAL FIT after initial finding all sources

gta.write_roi('global_post_refit',make_plots=False,save_model_map=False)
gta.residmap('SR_global_post_refit', make_plots=True)
gta.print_model()





####### Final cuts on TS #######
# Make final cuts on sources
delete_list = []
for s in range(0,len(gta.roi.sources)):

        source_name = gta.roi.sources[s].name
        ts = gta.roi[source_name]['ts']

        if ts < 9:
                delete_list.append(source_name)

for each in delete_list:
        gta.delete_source(each)

# Write final model:
gta.print_model()
gta.write_roi("Final_model",make_plots=False,save_model_map=False)


## Final cut on edge cuts ## 
######################################
# Make edge cuts

delete_list = []
for s in range(0,len(gta.roi.sources)):

        source_name = gta.roi.sources[s].name
        offset_edge = gta.roi[source_name]['offset_roi_edge']

        if offset_edge>-2.5:
            if 'e' in source_name:
                continue
            else:
               delete_list.append(source_name)

for each in delete_list:
        gta.delete_source(each)

# Write final model with edge cuts:
gta.print_model()
gta.write_roi("Final_model_with_2.5deg_edge_cut",make_plots=False,save_model_map=False)
gta.print_params()









f1.close()
f2.close()
