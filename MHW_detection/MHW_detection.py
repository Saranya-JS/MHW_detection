#### MHW-detection code 
#### Inspired by Eric Oliver's codes and also follows the MHW framework of Hobday et al. (2016)  
# Written by Panini Dasgupta and Saranya JS 


import numpy as np
from datetime import date
import pandas as pd




def season_oliver(s,L=365,no_harmonics=3)
       """
       This function is adopted from Eric Oliver code
       input :
       s            : Time series
       L            : No of days in a year
       no_harmonics : Number of harmonics.
       
       Return :
       smooth       : smoothed climatology
       anomaly      : anomaly after removing climatology.
       
       
       """
    
    n = len(s)
    time = np.mat(np.arange(1,n+1)/(1.*L))

    #   set up mean and harmonics to fit data

    P = np.mat(np.ones((n,1)))
    K = no_harmonics
    
    for k in range(1, K+1):
        P = np.concatenate((P, np.cos(k*2*np.pi*time.T)), 1)
        P = np.concatenate((P, np.sin(k*2*np.pi*time.T)), 1)

    #   Remove seasonal cycle by harmonic regression
    beta = (np.linalg.inv(P.T*P)*P.T)*s
    smooth = P*beta
    anomaly = s - perc

    return smooth,anomaly



def select_seg_(i,len_1,len_0):
      """
       This function is adopted from Hobday et al. (2016)
       input 
       i         :   (Temperature >= Threshold)*1
       len_1     :   5
       len_0     :   2
       
       Return 
       ind1      :   where Temperature > Threshold
       ind2      :   where Temperature < Threshold
       """
   
    ind1 = []
    ind2 = []
    z = []
   
    t2 = i[1:]
    t3 = i[0:-1];
    z =  t2-t3;
   
    ind1 = np.where(z == 1)[0] +1;
    ind2 = np.where(z == -1 ) [0]+1 ;
    if (ind1[0] > ind2[0]):
        ind2 = ind2[1:]
    if (ind1[-1]>ind2[-1]):
        ind1 = ind1[:-1]
    
    print(len(ind1))
    
    ind_1 = np.where((ind2-ind1)< len_1)[0]
    
    print(len(ind_1))
    
    ind1 = np.delete(ind1,ind_1)
    ind2 = np.delete(ind2,ind_1)
    k = ind1[1:] - ind2[:-1]
    w = np.where(k<=len_0)[0]
    while len(w)>0:
        for j in w:
                ind2[j] = ind2[j+1]
        ind1 = np.delete(ind1,j+1)
        ind2 = np.delete(ind2,j+1)
        k = ind1[1:] - ind2[:-1]
        w = np.where(k<=len_0)[0]
     
   
    return ind1,ind2


def mhw_param(df,ind1,ind2):
"""""
        input 
        df     :   Dataframe of temperature anomaly
        ind1   :
        ind2   :
        Return
        mhw_df :  Dataframe with MHW statistics
"""""  
      
    mhw_mean= []
    mhw_max= []
    mhw_cum= []
    mhw_duration= []
    no_days = []
    max_date = []
    ## reading sst values ###
    
    ssta = df.values.squeeze()
    time = pd.Series(df.index)
    sd = time.iloc[ind1] #### starting dates
    ed = time.iloc[ind2-1] ### ending dates 
    
    ### prepare the dataframe
    mhw_df= pd.DataFrame()
    
    for i in range(len(ind1)):
        no_days.append(ind2[i]-ind1[i])
        mhw_mean.append(np.mean(ssta[ind1[i]:ind2[i]]))
        mhw_max.append(np.max(ssta[ind1[i]:ind2[i]]))
        mhw_cum.append(np.sum(ssta[ind1[i]:ind2[i]]))
        max_date.append(np.where(df==np.max(ssta[ind1[i]:ind2[i]]))[0].tolist())
    
    max_date = np.array(max_date).squeeze()
    md = time.iloc[max_date] 
    sd1 = sd.to_frame().reset_index(drop =True)
    ed1 = ed.to_frame().reset_index(drop =True)
    md1 = md.to_frame().reset_index(drop =True)
    
    
    nod   = pd.DataFrame(no_days)
    mean1 = pd.DataFrame(mhw_mean)
    max1  = pd.DataFrame(mhw_max)
    cum1  = pd.DataFrame(mhw_cum)
    
    mhw_df = pd. concat([sd1,ed1,md1,nod,mean1,max1,cum1],axis=1,ignore_index=1)
    mhw_df.columns = ['start_date','end_date','Max_date','No of days','mean_intensity','max_intensity','cumulative_intensity']
    return mhw_df








def MHW_Area_calculation(df,DS,l):
"""""
  This function give us the area of MHW during the peak dates
        input 
        df     :   Dataframe of temperature, anomaly, threshold and climatology with MHW peak dates
        DS     :   Dataset of the region with latitude,longitude and temperature 
        l      :   list of MHW peak dates
        Return
        MHW_area_df:  Dataframe of MHW area
"""""  


    no_grid = np.zeros((len(df)))

    for i in range(len(df)):

            threshold = df.iloc[i,0]
            time = df.iloc[i].name
            temp = DS.sst.sel(time=time)
            grid= np.where((temp> threshold))[0]
            no_grid[i] = len(grid)
    grid = pd.DataFrame(no_grid, columns=['no_grid'])
    MHW_area_df = grid.set_index([l]) 
    MHW_area_df["resolution"] = "756.25"  #0.25 degree resolution
    MHW_area_df["area"] = ""
    MHW_area_df["area"] = MHW_area_df.no_grid * MHW_area_df.resolution.astype(float)

    return MHW_area_df
