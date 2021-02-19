import cvxpy as cp
import numpy as np
#import environment10
#import subproblem22
#import subproblem21
#import subproblem1
import scipy.io as sio
import random

import random
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from matplotlib import cm
#import compute_rate_assist
import os
import time
import os.path as osp
import cvxpy as cp
import numpy as np
import matplotlib
import math
import random
from scipy.stats import rice
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from matplotlib import cm
#import main10_5
import scipy.io as sio
import os
import time
import os.path as osp

import simulation_environment
import basic_problem

def state_update_func(env):
    I = env["I"]
    M = env["M"]
    h_ul = env["h_ul"]
    bandwidth_max = env["bandwidth_max"]
    tx_power_m_max = env["tx_power_m_max"]
    comp_bs = env["comp_bs"]
    comp_m = env["comp_m"]
    task = env["task"]
    tau = env["tau"]
    data = env["data"]
    M_list = env["M_list"]
    time_scale = env["time_scale"]

    rate_m_itr = bandwidth_max * np.log2(1 + tx_power_m_max * h_ul)
    average_rate = bandwidth_max * np.log2(1 + tx_power_m_max * h_ul) / M
    t = 0
    time_interval = 1/time_scale
    state_h_ul = []

    state_alpha = []
    state_beta = []
    state_omega = []
    state_rate = []

    state_local_remain_data = []
    state_bs_remain_data = []
    state_local_remain_data[t] = max(state_local_remain_data[t-1]*state_alpha[t-1]-state_rate[t-1]*time_interval,0)+ max(state_local_remain_data[t-1]*(1-state_alpha[t-1])-comp_m*time_interval/tau,0)

    state_bs_remain_data[t] = max(state_bs_remain_data[t-1]-comp_bs*state_beta[t-1]*time_interval,0)+min(state_local_remain_data[t-1]*state_alpha[t-1], state_rate[t-1]*time_interval)


    time_scale_num = 1000
    local_state = []
    bs_state = []
    local_state.append()
    state_update = {"h_ul": h_ul, "I": I, "M": M, "tx_power_m_max": tx_power_m_max, "bandwidth_max": bandwidth_max,
       "comp_m": comp_m, "comp_bs": comp_bs,
       "task": task, "data": data, "M_list": M_list, "time_scale": time_scale,"state_local_remain_data":state_local_remain_data,"state_bs_remain_data":state_bs_remain_data}

    return state_update