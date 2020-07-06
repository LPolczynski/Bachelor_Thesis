from math import pi
import shutil

def water_diffusion(w_diff,pi_t,pi_t_d):
    dV_exchange = - w_diff*(pi_t - pi_t_d)
    dV_exchange_d = + w_diff*(pi_t - pi_t_d)
    return dV_exchange,dV_exchange_d


def osmolyte_diffusion(c_diff, c_i, c_i_d, r, r_d):
    dc_i_diff = - c_diff*(c_i - c_i_d)/(4/3*pi*r**3)
    dc_i_diff_d = + c_diff*(c_i - c_i_d)/(4/3*pi*r_d**3)
    return dc_i_diff, dc_i_diff_d

def get_growth_species(spec,times,t,i):
    return spec[i] + (spec[i+1]-spec[i])*(t-times[i])/(times[i+1]-times[i])

def get_time_index(times,t,i_start = None):
    right_index = False
    if i_start:
        i = i_start
    else:
        i = 0
    while not right_index:
        if times[i] <= t and t <= times[i+1]:
            right_index = True
        else:
            i += 1
    return i

def save_sim_file(simfile, save_dir):
    src = simfile
    dst = save_dir
    shutil.copy(src, dst)