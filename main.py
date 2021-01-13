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
env = simulation_environment.one_cell(now, seed)
print(env["I"])

basic_problem.solve_problem_func(env)
