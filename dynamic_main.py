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

import multi_time_env
import multi_time_problem
import copy

now = time.strftime("%Y-%m-%d-%H_%M_%S", time.localtime(time.time()))
seed = 1
np.random.seed ( seed )
random.seed ( seed )
main_env = simulation_environment.one_cell(now, seed)
print(main_env["I"])
I = main_env["I"]
M = main_env["M"]
h_ul = main_env["h_ul"]
bandwidth_max = main_env["bandwidth_max"]
tx_power_m_max = main_env["tx_power_m_max"]
comp_bs = main_env["comp_bs"]
comp_m = main_env["comp_m"]
task = main_env["task"]
tau = main_env["tau"]
data = main_env["data"]
M_list = main_env["M_list"]
time_scale = main_env["time_scale"]
time_interval = 1/time_scale
T = 100 #total circle

#initialize the iteration
main_alpha_v = 0.5*np.ones([main_env["M"],main_env["I"]])
harvest_time_circle = np.random.randint(30,50,[main_env["M"],main_env["I"]])


result = []

one_cell_AoI = []
single_user_AoI = []
single_user_circle = []
one_cell_circle = []
single_user_change_circle = []
single_user_circle_start_time = []
for i in range(0,T):
    single_user_circle_start_time.append(np.zeros([M, I]))


one_cell_circle_start_time = np.zeros([T, I])
one_cell_circle_flag  = np.zeros([T,I])
for i in range (0,I):
    one_cell_circle_flag[0:T,i] = np.linspace(1,101,100)


#single_user_AoI.append( np.random.randint(5,10,[main_env["M"],main_env["I"]]))
single_user_AoI.append( np.zeros([main_env["M"],main_env["I"]]))
one_cell_AoI.append(np.zeros(I))
one_cell_AoI[0] = np.max(single_user_AoI[0],axis=0)
single_user_circle.append(np.ones([M,I]))
one_cell_circle.append(np.ones(I))
single_user_change_circle.append(np.zeros([M,I]))
single_user_change_circle.append(np.zeros([M,I]))
single_circle_now = np.ones([M,I])
one_cell_circle_now = np.ones(I)

state_alpha = []
state_beta = []
state_omega = []
state_rate = []
#intialize AoI and state
state_local_remain_data = []
state_bs_remain_data = []



# update state
state_local_remain_data.append(data)
state_bs_remain_data.append(np.zeros([M,I]))


# solve optimization problem
main_alpha_v = 0.1*np.ones([main_env["M"],main_env["I"]])
main_beta_v = 0.1*np.ones([main_env["M"],main_env["I"]])
main_omega_v = 0.1*np.ones([main_env["M"],main_env["I"]])

result_store = []

result_store.append(basic_problem.solve_problem_func(env=main_env, alpha_parameter=main_alpha_v,beta_parameter=main_beta_v,omega_parameter=main_omega_v,single_user_AoI_parameter=single_user_AoI[0],one_cell_AoI_parameter=one_cell_AoI[0]))
print(result_store[0])
flag = 1
result_value = []
alpha_itr = []
beta_itr = []
omega_itr = []
result_value.append(result_store[0]["problem6.value"])
alpha_itr.append(result_store[0]["alpha.value"])
beta_itr.append(result_store[0]["beta.value"])
omega_itr.append(result_store[0]["omega.value"])

