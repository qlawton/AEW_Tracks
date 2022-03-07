## Run Boostrapped Significance

def bootstrap_compare(array1, array2, n_comps, thresh_perc = 95, axis = 0):
    import numpy as np
    import xarray as xr
    import numpy.ma as ma
    '''This function runs a bootstrapped statistical significance test and returns 
    the results. Bootstrapping is done with replacement.'''
    array1_mini_comp = []
    array2_mini_comp = []
    thresh = (100 - thresh_perc)/2/100.
    stat_thresh = int(round(n_comps*thresh)) #Get our threshold cutoff
    
    #First, loop over for the first array (array 1)
    for ni in range(n_comps):
        comp_sample = array1.shape[axis] #get number of samples
        index = np.random.choice(array1.shape[axis], comp_sample, replace = True) #Random index samples, with replacement
        try:
            if axis == 0:
                temp_array = array1[index,:,:]
            elif axis == 1:
                temp_array = array1[:,index,:]
            elif axis == 2:
                temp_array = array1[:,:, index]
            else:
                raise Exception('This function only supports axis sorting between 0 and 2')
        except: 
            raise Exception('Error pulling out random indices')
        temp_final_mean = np.nanmean(temp_array, axis = axis)
        
        if array1_mini_comp == []:
            array1_mini_comp = temp_final_mean
        else:
            array1_mini_comp = np.dstack((array1_mini_comp, temp_final_mean))
    #Next, loop over second array (array 2)
    for ni in range(n_comps):
        comp_sample = array2.shape[axis] #get number of samples
        index = np.random.choice(array2.shape[axis], comp_sample, replace = True) #Random index samples, with replacement
        try:
            if axis == 0:
                temp_array = array2[index,:,:]
            elif axis == 1:
                temp_array = array2[:,index,:]
            elif axis == 2:
                temp_array = array2[:,:,index]
            else:
                raise Exception('This function only supports axis sorting between 0 and 2')
        except: 
            raise Exception('Error pulling out random indices')
        temp_final_mean = np.nanmean(temp_array, axis = axis)
        
        if array2_mini_comp == []:
            array2_mini_comp = temp_final_mean
        else:
            array2_mini_comp = np.dstack((array2_mini_comp, temp_final_mean))
            
            
    #Finally, spit out final arrays
    array1_mini_comp = np.array(array1_mini_comp)
    array2_mini_comp = np.array(array2_mini_comp)
    sorted_array1 = np.sort(array1_mini_comp, axis = 2)
    sorted_array2 = np.sort(array2_mini_comp, axis = 2)
    
    #Now, define our significance arrays
    sig_1 = sorted_array1[:,:,stat_thresh] > sorted_array2[:,:,-stat_thresh] #One version of significance testing
    sig_2 = sorted_array1[:,:,-stat_thresh]< sorted_array2[:,:,stat_thresh] #Another version
    
    sig_total = np.zeros(sig_1.shape)
    sig_total[sig_1 == 1] = 1
    sig_total[sig_2 == 1] = 1
    
    diff_mask = ma.masked_where(sig_total == 0, (array1.mean(axis = 0) - array2.mean(axis = 0))) #Create a difference array with non-sig things masked out
    diff_filter = ma.masked_where(sig_total == 1, (array1.mean(axis = 0) - array2.mean(axis = 0))) #Create a difference array with non-sig things masked out
    
    return diff_mask, diff_filter, sig_total, sorted_array1, sorted_array2