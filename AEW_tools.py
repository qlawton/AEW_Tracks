#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 12:29:29 2022

@author: Quinton Lawton, University of Miami 
"""

import numpy as np

def get_dist_meters(lon, lat):
        earth_circ = 6371*2*np.pi*1000 #earth's circumference in meters
        lat_met = earth_circ/360 #get the number of meters in a degree latitude (ignoring "bulge")
        lat_dist = np.gradient(lat, axis=0)*lat_met
        lon_dist = np.gradient(lon, axis=1)*np.cos(np.deg2rad(lat))*lat_met
        return lon_dist, lat_dist 


def rad_mask(i, j, dx, dy, radius, res):
        #start = tm.time()
        
        earth_circ = 6371*2*np.pi*1000 #earth's circumference in meters
        lat_met = earth_circ/360 #get the number of meters in a degree latitude (ignoring "bulge")
        lat_km_lat = lat_met/1000/1 #Updated to fix high latitude errors...
        lat_km_lon = lat_met/1000/2 #Updated to work up to 60N...
        buffer = int(np.ceil(radius/lat_km_lat/res))
        buffer_j = int(np.ceil(radius/lat_km_lon/res))
        boolean_array = np.zeros(np.shape(dx), dtype=bool)
        #print(buffer)
        #i_array = np.zeros(np.shape(dy))
        #j_array = np.zeros(np.shape(dx))
        dy = -dy

        # Before doing this, we recognize something -- there is a max radius in the x and y direction that will
        # be computed. To save on computational time, we need to slice out the extraneous points that we already
        # know from this NOT to be in the circle. How? We know the absolute highest distance between gridspaces
        # will be the maximum value of latitude on the elliptical earth... ~ 111 km. So slice out roughly a 112km "box" 
        # plus one gridbox for buffer. 
        i_st = (i-(buffer+1))
        i_end = (i+(buffer+1))
        j_st = (j-(buffer_j+1))
        j_end = (j+(buffer_j+1))
        #print(i_st, i_end)
        #print(j_st, j_end)

        new_i = i-i_st
        new_j = j-j_st

        dy_slc = dy[i_st:i_end,j_st:j_end]
        dx_slc = dx[i_st:i_end,j_st:j_end]
        
        #print(np.shape(dy_slc), np.shape(dy))

            
        i_array_sub = np.zeros(np.shape(dy_slc))
        j_array_sub = np.zeros(np.shape(dx_slc))

        i_array_sub[new_i,:] = 0
        i_array_sub[(new_i+1):,:] = np.add.accumulate(dy_slc[(new_i+1):,:])
        i_array_sub[(new_i-1)::-1,:] = np.add.accumulate(dy_slc[(new_i-1)::-1,:])

        j_array_sub[:,new_j] = 0
        j_array_sub[:,(new_j+1):] = np.add.accumulate(dx_slc[:,(new_j+1):], axis=1)
        j_array_sub[:,(new_j-1)::-1] = np.add.accumulate(dx_slc[:,(new_j-1)::-1], axis=1)

        radial_array = (np.sqrt(np.square(i_array_sub)+np.square(j_array_sub))/1000) < radius # in km
        boolean_array[i_st:i_end,j_st:j_end] = radial_array
        #end = tm.time()
        #print(end - start)
        return boolean_array
    
def lon_lat_dxdy(coord1, coord2):
    import geopy.distance
    ''' Get the dx, dy, and total distance between two latlon coordinates. Output in kilometers.'''
    dxcoord1 = (coord1[0], coord1[1])
    dxcoord2 =  (coord2[0], coord1[1]) #Not a mistake, want to only look at dx

    dycoord1 = (coord1[0], coord1[1])
    dycoord2 = (coord1[0], coord2[1]) #Not a mistake, want to only look at dy

    dx_abs = geopy.distance.distance(dxcoord1, dxcoord2).km
    dy_abs = geopy.distance.distance(dycoord1, dycoord2).km
    dtot = geopy.distance.distance(coord1, coord2).km
    
    # Now, determine which direction these changes have occured
    if coord2[0]<coord1[0]:
        dx = -dx_abs
    else:
        dx = dx_abs
    if coord2[1]<coord1[1]:
        dy = -dy_abs
    else:
        dy = dy_abs
    
    return dx, dy
def stm_motion(init_idx, lon_list, lat_list, dt_list):
    '''Calculate the motion of system in each direction'''

    #Pull dx, dy values
    coord1 = (lat_list[init_idx], lon_list[init_idx])
    coord2 = (lat_list[init_idx+1], lon_list[init_idx+1])
    #print(coord1, coord2)
    dx, dy = lon_lat_dxdy(coord1, coord2)
    #print(coord1, import geopy.distancecoord2)

    curr_time = datetime.datetime.strptime(dt_list[init_idx+1], "%Y%m%d%H")
    prev_time = datetime.datetime.strptime(dt_list[init_idx], "%Y%m%d%H")
    timedelta = curr_time - prev_time #time delta 
    #print(timedelta)
    ### ------------------
    timedelta_sec = timedelta.seconds #Convert to float in seconds
    ### ------------------

    vel_x = (dx/timedelta_sec)*1000 #Take dx (in km), divide by dt (seconds), and then convert to m/s 
    vel_y = (dy/timedelta_sec)*1000 #Take dy (in km), divide by dt (seconds), and then convert to m/s

    return vel_y, vel_x

def stm_motion_centered(init_idx, lon_list, lat_list, dt_list):
    '''Calculate the motion of system in each direction'''
    
    if init_idx != 0 and init_idx != len(lon_list)-1: #So basically, if we can take a centered average
        #Pull dx, dy values
        coord1 = (lat_list[init_idx-1], lon_list[init_idx-1])
        coord2 = (lat_list[init_idx+1], lon_list[init_idx+1])
        #print(coord1, coord2)
        dx, dy = lon_lat_dxdy(coord1, coord2)
        #print(coord1, import geopy.distancecoord2)

        curr_time = datetime.datetime.strptime(dt_list[init_idx+1], "%Y%m%d%H")
        prev_time = datetime.datetime.strptime(dt_list[init_idx-1], "%Y%m%d%H")
        timedelta = curr_time - prev_time #time delta 
        #print(timedelta)
        ### ------------------
        timedelta_sec = timedelta.seconds #Convert to float in seconds
        ### ------------------

        vel_x = (dx/timedelta_sec)*1000 #Take dx (in km), divide by dt (seconds), and then convert to m/s 
        vel_y = (dy/timedelta_sec)*1000 #Take dy (in km), divide by dt (seconds), and then convert to m/s
        
    elif init_idx == 0: #If we are starting out at index 0, take a forward speed
        #Pull dx, dy values
        coord1 = (lat_list[init_idx], lon_list[init_idx])
        coord2 = (lat_list[init_idx+1], lon_list[init_idx+1])
        #print(coord1, coord2)
        dx, dy = lon_lat_dxdy(coord1, coord2)
        #print(coord1, import geopy.distancecoord2)

        curr_time = datetime.datetime.strptime(dt_list[init_idx+1], "%Y%m%d%H")
        prev_time = datetime.datetime.strptime(dt_list[init_idx], "%Y%m%d%H")
        timedelta = curr_time - prev_time #time delta 
        #print(timedelta)
        ### ------------------
        timedelta_sec = timedelta.seconds #Convert to float in seconds
        ### ------------------

        vel_x = (dx/timedelta_sec)*1000 #Take dx (in km), divide by dt (seconds), and then convert to m/s 
        vel_y = (dy/timedelta_sec)*1000 #Take dy (in km), divide by dt (seconds), and then convert to m/s
        
    elif init_idx == len(lon_list)-1: #If we are at end, use just a previous point to define
        #Pull dx, dy values
        coord1 = (lat_list[init_idx-1], lon_list[init_idx-1])
        coord2 = (lat_list[init_idx], lon_list[init_idx])
        #print(coord1, coord2)
        dx, dy = lon_lat_dxdy(coord1, coord2)
        #print(coord1, import geopy.distancecoord2)

        curr_time = datetime.datetime.strptime(dt_list[init_idx], "%Y%m%d%H")
        prev_time = datetime.datetime.strptime(dt_list[init_idx-1], "%Y%m%d%H")
        timedelta = curr_time - prev_time #time delta 
        #print(timedelta)
        ### ------------------
        timedelta_sec = timedelta.seconds #Convert to float in seconds
        ### ------------------

        vel_x = (dx/timedelta_sec)*1000 #Take dx (in km), divide by dt (seconds), and then convert to m/s 
        vel_y = (dy/timedelta_sec)*1000 #Take dy (in km), divide by dt (seconds), and then convert to m/s

    return vel_y, vel_x

def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        #print(idx)
        return idx, array[idx]
def choose_bin(bin_bounds, lon_in, lat_in, lat_min, lat_max):
    if lon_in<bin_bounds[0] or lon_in>bin_bounds[-1]:
        #Outside range, don't include 
        return np.NaN
    if lat_in > lat_max or lat_in <lat_min:
        #Outside latitude bounds
        return np.NaN
    
    closest_i = find_nearest(bin_bounds, lon_in)[0]
    
    if bin_bounds[closest_i] <= lon_in:
        out_i = closest_i
        
    elif bin_bounds[closest_i]> lon_in:
        out_i = closest_i -1
    
    return out_i
        