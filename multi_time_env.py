import cvxpy as cp
import numpy as np
import matplotlib
import math
import random
from scipy.stats import rice
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from matplotlib import cm
# import main
import scipy.io as sio
import os
import time
import os.path as osp


def one_cell_multi_time(seed, time, main_env, state):



    np.random.seed ( seed )
    random.seed ( seed )
    interference = 1
    psi = 4
    rou = 0
    time_scale = 100 #1000ms=1 100ms=10 10ms=100 1ms=1000
    scenario = [0, 1]
    # ---------------simulation environment---------------------
    cell_num = [1, 3, 7]
    scenario = 'uma'  # uma umi free
    cell_radius = [100]
    if scenario == 'uma':
        ISD = 500
    elif scenario == 'umi':
        ISD = 200
    elif scenario == 'free':
        ISD = 500

    # ISD = 500 #base station distance
    cell_radius = ISD / 2
    task_num = [2, 5, 10, 15, 20, 25, 30, 50]
    # ue_num = [2, 3, 4, 5, 6, 7]
    I = cell_num[0]
    # M = task_num[0]
    M = task_num[2]
    # N = ue_num

    # ----------computing environment---------#
    comp_m = 1e9 * np.random.uniform ( 0.5, 1, (M, I) )
    # comp_sc = 1e9 * 2*np.ones([M,I])
    # comp_bs = 1e9 * np.random.uniform ( 2,10 , I )
    comp_bs = 6* 1e9
    # comp_ue = 1e9 * np.random.uniform ( 2, 4, (N, I) )
    #data = 1e3 * np.random.uniform ( 300, 500, (M, I) )
    tau = main_env["tau"]
    local_remain_data = state[0]
    bs_remain_data = state[1]
    local_task = np.multiply ( tau, local_remain_data )
    bs_task = np.multiply(tau, bs_remain_data)
    M_list = np.arange ( 1, M+1, 1 )












    # ----------cell location----------#
    height_bs = 25  # uma
    height_bs = 10  # umi
    cell_locx = 0
    cell_locy = 0

    # ----------user distribution----------#
    m_theta = np.random.uniform ( 0, 2 * math.pi, (M, 1) )
    m_x_theta = np.cos ( m_theta )
    m_y_theta = np.sin ( m_theta )
    distance_m2D = np.random.uniform ( 35, ISD / 2, (M, 1) )
    height_m = np.random.uniform ( 1.5, 2, (M, 1) )

    # ----------distance between user and the base station----------#

    cell_m_locx = cell_locx+distance_m2D * m_x_theta
    cell_m_locy = cell_locy+distance_m2D * m_y_theta
    distance_m3D = np.sqrt ( np.square ( distance_m2D )+np.square ( height_bs-height_m ) )

    # ----------communication enviroment----------#
    tx_power_bs_max = math.pow ( 10, 49 / 10 )  # downlink power 49dBm
    tx_power_m_max = 1e-3*math.pow ( 10, 28 / 10 )  # uplink power 28dBm
    tx_power_m_max = math.pow(10, 28 / 10)  # uplink power 28dBm
    #bandwidth_max = 20 * 1e6  # bandwidth 20MHz
    carrier = 180  # KHz
    sigma2 = -174  # dBm/Hz
    sigma2_power = -100  # dBm
    bandwidth_max = 180*1e3*M

    reference = 1
    factor = math.pow ( 10, -100 / 10 )

    # ----------pathloss setting 3GPP TR 38.901 V16.1.0 (2019-12) page 27--------#
    if scenario == 'uma':
        # PL = 32.4 + 20log10(fc)(fc=GHZ) + 30log10(d_3D)(d_3D meter) + sigma^2=7.8 uma
        PL_m = 32.4+ 20 * np.log10 ( 2.5 )+ 30 * np.log10 ( distance_m3D )+np.random.normal ( loc=0.0, scale=7.8,
                                                                                            size=(M, 1) )  # dB
        h_ul = np.power(10, -1 * PL_m / 10) / factor
    elif scenario == 'umi':
        # PL = 35.3log10(d_3D)(d_3D meter)+ 22.4+21.3log10(fc)(fc=GHZ)  -0.3*(h_ut-1.5)+ sigma^2=7.82 umai
        PL_m = 35.3 * np.log10 ( distance_m3D )+22.4+21.3 * np.log10 ( 2.5 )-0.3 * (height_m-1.5)+np.random.normal (
            loc=0.0, scale=7.82, size=(M, 1) )
        h_ul = np.power(10, -1 * PL_m / 10) / factor
    elif scenario == 'free':
        PL_m = 10*np.log10(1*np.square(3*1e8/2e9)/np.square(4*3.14*distance_m3D))
        h_ul = np.square(3*1e8/2e9)/np.square(4*3.14*distance_m3D)/factor



    rate_m = bandwidth_max / M * np.log2 ( 1+tx_power_m_max * h_ul )
    print("distance2d",distance_m2D)
    print("normal", np.random.normal ( loc=0.0, scale=7.8,size=(M, 1) ))
    print("hul", h_ul)
    print ( rate_m / 1e6 )
    print ( M_list )
    #print ( data / rate_m )

    '''
    y = data / rate_m
    print ( y )
    plt.subplot ( 331 )
    plt.title ( "user_rate" )
    plt.bar ( x=M_list, height=(rate_m / 1e6).reshape ( M ), width=1 )

    plt.subplot ( 332 )
    plt.title ( "full_transmit_delay" )
    plt.bar ( x=M_list, height=y.reshape ( M ), width=1 )

    plt.subplot ( 333 )
    plt.title ( "local_computing_delay" )
    plt.bar ( x=M_list, height=(task / comp_m).reshape ( M ), width=1 )

    plt.subplot ( 334 )
    plt.title ( "offloading_computing_delay" )
    plt.bar ( x=[0], height=(np.sum ( task ) / comp_bs).reshape ( 1 ), width=1 )
    '''
    #plt.show ()

    one_cell_env = {"h_ul": h_ul, "I": I, "M": M, "tx_power_m_max": tx_power_m_max, "bandwidth_max":bandwidth_max, "comp_m":comp_m, "comp_bs":comp_bs,
                   "local_task":local_task, "bs_task":bs_task,"local_remain_data":local_remain_data,"bs_remain_data":bs_remain_data, "M_list":M_list,"time_scale":time_scale, "tau":tau}
    return one_cell_env
