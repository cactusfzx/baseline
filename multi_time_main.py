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

#initialize the iteration
main_alpha_v = 0.5*np.ones([main_env["M"],main_env["I"]])
harvest_time_circle = np.random.randint(30,50,[main_env["M"],main_env["I"]])


result = []

one_cell_AoI = []
single_user_AoI = []

single_user_AoI.append( np.random.randint(5,10,[main_env["M"],main_env["I"]]))
one_cell_AoI.append(max(single_user_AoI[0]))

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
main_alpha_v = 0.5*np.ones([main_env["M"],main_env["I"]])
main_beta_v = 0.1*np.ones([main_env["M"],main_env["I"]])
main_omega_v = 0.1*np.ones([main_env["M"],main_env["I"]])

result_store = []

result_store.append(basic_problem.solve_problem_func(env=main_env, alpha_parameter=main_alpha_v,beta_parameter=main_beta_v,omega_parameter=main_omega_v))
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
for itr in range(1,100):
    flag = flag + 1
    result_store.append(
        basic_problem.solve_problem_func(env=main_env, alpha_parameter=alpha_itr[itr-1], beta_parameter=beta_itr[itr-1],
                                         omega_parameter=beta_itr[itr-1]))
    result_value.append(result_store[itr]["problem6.value"])

    if result_value[itr] == 0:
        result_store[itr] = result_store[itr-1]


    alpha_itr.append(result_store[itr]["alpha.value"])
    beta_itr.append(result_store[itr]["beta.value"])
    omega_itr.append(result_store[itr]["omega.value"])
    if abs(result_value[itr]-result_value[itr-1])/result_value[itr-1] <1e-7:
        #update variable
        state_alpha.append(alpha_itr[itr])
        state_beta.append(beta_itr[itr])
        state_omega.append((omega_itr[itr]))
        state_rate.append(state_omega[0] * bandwidth_max * np.log2(1 + tx_power_m_max * h_ul))
        break
#t is the scheduling cirlcle one t = 10ms

for t in range(1,100):

    temp1 = np.zeros([M, I])
    temp2 = np.zeros([M, I])
    single_user_AoI.append(np.zeros([M,I]))
    one_cell_AoI.append(0)
    #update state
    for m in range(0,M):
        for i in range(0,I):

            temp1[m,i] = max(state_local_remain_data[t-1][m,i]*state_alpha[t-1][m,i]-state_rate[t-1][m,i]*time_interval,0)+ max(state_local_remain_data[t-1][m,i]*(1-state_alpha[t-1][m,i])-comp_m[m,i]*time_interval/tau[m,i],0)
            temp2[m,i] = max(state_bs_remain_data[t-1][m,i]-comp_bs*state_beta[t-1][m,i]*time_interval/tau[m,i],0)+min(state_local_remain_data[t-1][m,i]*state_alpha[t-1][m,i], state_rate[t-1][m,i]*time_interval)

            if temp1[m,i] <=0 and temp2[m,i]<=0:
                single_user_AoI[t][m,i] = 0
            else:
                single_user_AoI[t][m, i] = single_user_AoI[t-1][m, i] + 1

            if np.all(single_user_AoI) == 0:
                one_cell_AoI[t] = 0
            else:
                one_cell_AoI[t] = one_cell_AoI[t-1] + 1

    state_local_remain_data.append(temp1)
    state_bs_remain_data.append(temp2)

    multi_time_para = multi_time_env.one_cell_multi_time(seed=seed, time=t, main_env=main_env,
                                                        state=[state_local_remain_data[t],state_bs_remain_data[t]])
    #solve optimization problem
    result_store = []
    result_store.append(
        multi_time_problem.solve_multi_time_problem_func(time=t,multi_time_env=multi_time_para, alpha_parameter=main_alpha_v, beta_parameter=main_beta_v,
                                         omega_parameter=main_omega_v))
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
    for itr in range(1, 100):
        flag = flag + 1
        result_store.append(
            multi_time_problem.solve_multi_time_problem_func(time=t,multi_time_env=multi_time_para, alpha_parameter=alpha_itr[itr - 1],
                                             beta_parameter=beta_itr[itr - 1],
                                             omega_parameter=beta_itr[itr - 1]))
        result_value.append(result_store[itr]["problem6.value"])

        if result_value[itr] == 0:
            result_store[itr] = result_store[itr - 1]

        alpha_itr.append(result_store[itr]["alpha.value"])
        beta_itr.append(result_store[itr]["beta.value"])
        omega_itr.append(result_store[itr]["omega.value"])
        if result_value[itr - 1]==0:
            state_alpha.append(alpha_itr[itr])
            state_beta.append(beta_itr[itr])
            state_omega.append((omega_itr[itr]))
            state_rate.append(state_omega[t] * bandwidth_max * np.log2(1 + tx_power_m_max * h_ul))
            break

        else:

            if abs(result_value[itr] - result_value[itr - 1]) / result_value[itr - 1] < 1e-7:
                # update variable
                state_alpha.append(alpha_itr[itr])
                state_beta.append(beta_itr[itr])
                state_omega.append((omega_itr[itr]))
                state_rate.append(state_omega[t] * bandwidth_max * np.log2(1 + tx_power_m_max * h_ul))
                break


#iteration = np.arange(1,100,[100])
'''
print(result_store)
plt.subplot(449)
plt.title("iteration result_value")
plt.plot(np.arange(1,flag+1), result_value[0:flag], '-^')
'''
# store the result to matlab file
os.mkdir ( now )
save_fn00 = 'Single_user_AoI.mat'
save_array = single_user_AoI
sio.savemat ( osp.join ( now, save_fn00 ), {'Single_user_AoI': save_array} )  # 和上面的一样，存在了array变量的第一行

# construct the single user AoI matrix into the vector and draw a 3D figure


plt.show()

