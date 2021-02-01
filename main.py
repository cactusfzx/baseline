# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
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
global seed
import simulation_environment
import basic_problem

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
now = time.strftime("%Y-%m-%d-%H_%M_%S", time.localtime(time.time()))
seed = 1
main_env = simulation_environment.one_cell(now, seed)
print(main_env["I"])

#initialize the iteration
main_alpha_v = 0.1*np.ones([main_env["M"],main_env["I"]])

result = []
result.append(basic_problem.solve_problem_func(env=main_env, alpha_parameter=main_alpha_v)["problem6.value"])
print(result[0])
flag = 1
alpha_itr = []
alpha_itr.append(basic_problem.solve_problem_func(env=main_env, alpha_parameter=main_alpha_v)["alpha.value"])
for itr in range(1,100):
    flag = flag + 1
    result.append(basic_problem.solve_problem_func(env=main_env, alpha_parameter=alpha_itr[itr-1])["problem6.value"])
    if result[itr] ==0:
        result[itr] = result[itr-1]
    alpha_itr.append(basic_problem.solve_problem_func(env=main_env, alpha_parameter=alpha_itr[itr-1])["alpha.value"])
    if abs(result[itr]-result[itr-1])/result[itr-1] <1e-7:
        break

#iteration = np.arange(1,100,[100])
print(result)
plt.subplot(449)
plt.title("iteration result")
plt.plot(np.arange(1,flag+1), result[0:flag], '-^')
plt.show()


