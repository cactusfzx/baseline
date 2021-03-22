import cvxpy as cp
import numpy as np
import matplotlib.pyplot as plt
from cvxpy import SolverError

# import simulation_environment
import matplotlib
import random
#import compute_rate_assist


# random.seed(0)

def solve_multi_time_problem_func(time,multi_time_env,alpha_parameter,beta_parameter,omega_parameter,single_user_AoI_parameter,one_cell_AoI_parameter,single_user_circle_parameter):
    # ----------parameters----------
    I = multi_time_env["I"]
    M = multi_time_env["M"]
    h_ul = multi_time_env["h_ul"]
    bandwidth_max = multi_time_env["bandwidth_max"]
    tx_power_m_max = multi_time_env["tx_power_m_max"]
    comp_bs = multi_time_env["comp_bs"]
    comp_m = multi_time_env["comp_m"]

    bs_remain_data = multi_time_env["bs_remain_data"]
    local_remain_data = multi_time_env["local_remain_data"]
    local_task= multi_time_env["local_task"]
    bs_task = multi_time_env["bs_task"]


    M_list = multi_time_env["M_list"]
    time_scale = multi_time_env["time_scale"]
    rate_m_itr = bandwidth_max * np.log2 ( 1+tx_power_m_max * h_ul )
    average_rate = bandwidth_max * np.log2 ( 1+tx_power_m_max * h_ul )/M
    alpha_v = alpha_parameter
    beta_v = beta_parameter
    omega_v = omega_parameter
    single_user_AoI = single_user_AoI_parameter
    onne_cell_AoI = one_cell_AoI_parameter
    single_user_circle = single_user_circle_parameter

    for j in range(0,I):
        print("before",j)
        single_user_AoI[0:M,j] = single_user_AoI[0:M,j]*(1+np.max(single_user_circle[0:M,j])-single_user_circle[0:M,j])
        print("after",j)
    # ----------itration variable parameter----------
    # this part is defined to use the iteration variable in the subproblem

    # ----------variables----------
    # define the variables
    # bandwidth allocation factor w
    omega = cp.Variable([M,I])
    # offloading factor alpha
    alpha = cp.Variable([M,I])
    # bs sever allocation factor beta
    beta = cp.Variable([M,I])
    eta = cp.Variable([M,I])
    mu = cp.Variable([M,I])

    # max delay T
    T = cp.Variable()

    # additional variables
    t1 = cp.Variable()
    t2 = cp.Variable()
    t3 = cp.Variable()




    # ----------problem formulation----------

    # --------objective function--------
    objective_func = T

    objective1 = cp.Minimize ( objective_func )

    # --------constraints--------
    # offloading constraints
    c1 = [alpha >= 0, alpha <= 1]
    # bandwidth allocation constraints
    c4_5 = [cp.sum(omega) <= 1, omega >= (1e-7)]
    # bs server allocation constraints
    c2_3 = [cp.sum(beta) <= 1, beta >= (1e-7)]

    c6 = [t1>=single_user_AoI]
    c6 = c6 + [cp.multiply ( 1-alpha, local_task )/comp_m*time_scale+single_user_AoI<=t1]
    #c7b=[-cp.log(t3)-cp.log(beta)+np.log2(alpha_v)+cp.multiply(1/alpha_v/np.log(2),(alpha-alpha_v))+np.log2((local_task+bs_task)/comp_bs*time_scale)<=0]
    #c8b=[-cp.log(t2)-cp.log(omega)+np.log2(alpha_v)+cp.multiply(1/alpha_v/np.log(2),(alpha-alpha_v))+np.log2(local_remain_data/rate_m_itr*time_scale)<=0]
    #rewrite the c7b and c8b in elementwise method
    c7b = [t3>=1/2*single_user_AoI]
    c8b = [t2>=1/2*single_user_AoI]
    for i in range(0,I):
        for m in range(0,M):
            if local_task[m,i]+bs_task[m,i] >= (1e-5):
                c7b= c7b + [-cp.log(t3-1/2*single_user_AoI[m,i])-cp.log(beta[m,i])+np.log2(alpha_v[m,i])+cp.multiply(1/alpha_v[m,i]/np.log(2),(alpha[m,i]-alpha_v[m,i]))+np.log2((local_task[m,i]+bs_task[m,i])/comp_bs*time_scale)<=0]

            if local_remain_data[m,i] >=(1e-5):
                c8b = c8b + [-cp.log(t2-1/2*single_user_AoI[m,i])-cp.log(omega[m,i])+np.log2(alpha_v[m,i])+cp.multiply(1/alpha_v[m,i]/np.log(2),(alpha[m,i]-alpha_v[m,i]))+np.log2(local_remain_data[m,i]/rate_m_itr[m,i]*time_scale)<=0]



    #c8 = [cp.multiply(eta,(local_task+bs_task)/comp_bs)*time_scale<=t3]
    #c7 = [cp.multiply ( mu, local_remain_data / rate_m_itr)*time_scale<=t2]

    R1a2 = []
    for i in range(0,M):
        for j in range(0,I):

            R1a2 =R1a2 + [cp.norm(cp.vstack([2*(alpha[i,j]),eta[i,j]-beta[i,j]]),2)<=eta[i,j]+beta[i,j]]
            #assert  R1a2[0].is_dqcp()

    R1a3 = [alpha<=1e7*beta+(1e-7)*eta-1,alpha<=eta]

    R2a3 = [alpha<=1e7*omega+(1e-7)*mu-1,alpha<=mu]

    R2a4 = [eta<=1e7,mu<=1e7]


    # local computing constraints R1
    #R1 = [cp.ceil(cp.multiply(alpha, task)/comp_m) <= T]
    R1 = [cp.ceil ( cp.multiply ( 1-alpha, local_task ) / comp_m ) -T <= 0]
    #R1 = [cp.ceil(alpha) - 10 <= 0]
    #assert R1[0].is_dqcp ()
    # offloading constraints R2
    R2a = cp.ceil(cp.multiply(cp.multiply(alpha,(local_task+bs_task)/comp_bs),cp.inv_pos(beta)))
    R2b = cp.ceil(cp.multiply(cp.multiply(alpha,local_remain_data/bandwidth_max*np.log2 ( 1+tx_power_m_max * h_ul )),cp.inv_pos(omega)))
    R2 = [R2a+R2b <= T]
    R2a1 = cp.ceil((cp.multiply(eta,(local_task+bs_task)/comp_bs)))

    objective2 = cp.Minimize(1/2*cp.max(cp.ceil ( cp.multiply ( alpha, local_task ) / comp_m )))
    objective2 = cp.Minimize(cp.max(cp.ceil(cp.multiply(1-alpha, (local_task+bs_task)/beta/comp_bs))))


    obj3_part1 = 1/2*cp.max(cp.ceil(cp.multiply(1-alpha,local_task/comp_m)))
    obj3_part1 = 1 / 2 * cp.ceil ( cp.max ( cp.multiply ( 1-alpha, local_task / comp_m ) ) )
    obj3_part2 = 1/2*cp.max(cp.ceil(cp.multiply(eta,(local_task+bs_task)/comp_bs))+cp.ceil(cp.multiply(mu,local_remain_data/rate_m_itr)))
    obj3_part2 = 1 / 2 * cp.ceil (cp.max( cp.multiply ( eta, (local_task+bs_task) / comp_bs ) ))+1/2*cp.ceil ( cp.max(cp.multiply ( mu, local_remain_data / rate_m_itr ) ) )

    obj3_part3 = cp.abs(obj3_part1-obj3_part2)
    objective3 = cp.Minimize(obj3_part1+obj3_part2+obj3_part3)


    #assert obj3_part1.is_dqcp ()
    #assert obj3_part1.is_dcp()
    #assert (cp.max(cp.ceil(cp.multiply(eta,local_task/comp_bs)))).is_dqcp()
    #assert (cp.max(cp.ceil(cp.multiply(eta,local_task/comp_bs)))).is_dcp()
    #assert (cp.ceil(cp.multiply(mu,local_remain_data/rate_m_itr))).is_dqcp
    #assert (cp.ceil(alpha)+cp.ceil(beta)).is_dqcp()
    #assert obj3_part2.is_dqcp ()
    #assert obj3_part3.is_dqcp ()
    #assert objective3.is_dqcp ()


    t = cp.Variable([M,I])
    objective2 = cp.Minimize ( cp.max ( cp.ceil ( cp.multiply ( t, (local_task+bs_task) / comp_bs ) ) ) )
    c4 = [(1-alpha)/beta <= t]

    objective6 = cp.Minimize(1/2*cp.ceil(t1+t2+t3+cp.abs(t1-t2-t3))+3/2)

    rho = 1
    upsilon = 1
    varsigma = 1e-27
    obj7_func1 = 1/2*cp.ceil(t1+t2+t3+cp.abs(t1-t2-t3))+3/2
    obj7_func1 = 1 / 2 * (t1 + t2 + t3 + cp.abs(t1 - t2 - t3))
    obj7_func2 = cp.sum(cp.multiply(1-alpha,varsigma*local_task*np.square(comp_m)))+cp.sum(tx_power_m_max*cp.multiply ( mu, local_remain_data / rate_m_itr))
    objective7 = cp.Minimize(rho*obj7_func1+upsilon*obj7_func2)





    try:
        # ----------probalem solve and results----------
        # problem6 = cp.Problem(objective6, c1+c2_3+c4_5+c6+c7+c8+R1a3+R2a3+R2a4)
        problem6 = cp.Problem(objective6, c1 + c2_3 + c4_5 + c6 + c7b + c8b)
        problem6.solve(qcp=True, verbose=True,solver=cp.ECOS)
        print()
        print("problem1 solve: ", problem6.value)
        print("alpha.value", alpha.value)
        print("beta.value", beta.value)
        print("omega.value", omega.value)

        np.ceil(local_remain_data / rate_m_itr * M * time_scale)
        # ----------data collection and depict-----------

        # plt.clf()
        # plt.subplot(441)
        # plt.title("user_rate in M/s")
        # # plt.bar(x=M_list, height=(average_rate/1e6).reshape(M), width=1)
        # plt.plot(M_list, omega.value * rate_m_itr / 1e6, '-*', color='b', label="optimized rate")
        # plt.plot(M_list, average_rate / 1e6, '-o', color='r', label="average bandwidth allocation  rate")
        #
        # plt.legend()
        #
        # plt.subplot(442)
        # plt.title("local_remain_data transmitting delay")
        # plt.plot(M_list, local_remain_data / average_rate * time_scale, '-*', color='b',
        #          label="full date transmitting delay with average bandwidth allocation")
        # plt.plot(M_list, (alpha.value * local_remain_data) / (omega.value * rate_m_itr) * time_scale, '-o', color='r',
        #          label="optimized transmitting delay")
        # # plt.bar(x=M_list, height=np.ceil(local_remain_data / average_rate*time_scale).reshape(M), width=1)
        # plt.legend(fontsize='xx-small')
        #
        # plt.subplot(443)
        # plt.title("local_computing_delay")
        # bar_width = 0.3
        # plt.bar(x=M_list, height=np.ceil(time_scale * local_task / comp_m).reshape(M), width=bar_width,
        #         label='full local computing delay')
        # plt.bar(x=M_list + bar_width, height=np.ceil(time_scale * (1 - alpha.value) * local_task / comp_m).reshape(M),
        #         width=bar_width, label='remain local computing delay')
        #
        # plt.plot(M_list, problem6.value * np.ones(M), 'o', color='m', label="optimized delay")
        # plt.legend(fontsize='xx-small')  # 显示图例，即label
        # plt.xticks(x=M_list + bar_width / 2)  # 显示x坐标轴的标签,即tick_label,调整位置，使其落在两个直方图中间位置
        #
        # plt.subplot(444)
        # plt.title("offloading_computing_delay")
        #
        # bar_width = 0.3  # 设置柱状图的宽度
        # plt.bar(x=M_list, height=(np.ceil((local_task+bs_task) / (comp_bs / M)) * time_scale).reshape(M), width=bar_width,
        #         label='average full offloading computing delay')
        # plt.bar(x=M_list + bar_width,
        #         height=(np.ceil(local_task * alpha.value / (beta.value * comp_bs) * time_scale)).reshape(M), width=bar_width,
        #         label='optimized offloading computing delay')
        #
        # # 绘制并列柱状图
        #
        # plt.legend()  # 显示图例，即label
        # plt.xticks(x=M_list + bar_width / 2)  # 显示x坐标轴的标签,即tick_label,调整位置，使其落在两个直方图中间位置
        #
        # plt.subplot(445)
        # plt.title("optimized local computing_delay gain")
        # plt.plot(M_list, (problem6.value - np.ceil(local_task / comp_m * time_scale)).reshape(M), 'x', color='r')
        # plt.plot(M_list, problem6.value * np.ones(M), 'o', color='m', label="optimized delay")
        # plt.plot(M_list, t1.value * np.ones(M), 'v', color='g', label="local computing delay t1")
        # plt.plot(M_list, t3.value * np.ones(M), '^', color='k', label="offloading computing delay t3")
        # plt.plot(M_list, t2.value * np.ones(M), '*', color='b', label="transmitting delay t2")
        # plt.legend(fontsize='xx-small')
        #
        # plt.subplot(446)
        # plt.title("optimized alpha beta")
        # plt.plot(M_list, alpha.value, '-v', label='alpha')
        # plt.plot(M_list, beta.value, '-x', label='beta')
        # plt.plot(M_list, omega.value, '^', label='omega')
        # plt.legend()
        #
        # plt.subplot(447)
        # plt.title("optimized omega")
        # plt.plot(M_list, omega.value, '-^')
        #
        # plt.subplot(448)
        # plt.title("optimized beta")
        # plt.plot(M_list, beta.value, '-^')

        '''
        plt.subplot(449)
        plt.title("optimized eta")
        plt.plot(M_list, eta.value,'-^')

        plt.subplot(4,4,10)
        plt.title("optimized mu")
        plt.plot(M_list, mu.value,'-^')

        plt.subplot(4,4,11)
        plt.title("recalculated omega")
        plt.plot(M_list, alpha.value/eta.value,'-^')

        plt.subplot(4,4,12)
        plt.title("recalculated beta")
        plt.plot(M_list, alpha.value/mu.value,'-^')
        '''

        print("func1.value", obj7_func1.value)
        print("local computing energy", cp.sum(cp.multiply(1 - alpha, varsigma * local_task * np.square(comp_m))).value)
        print("offloading energy", cp.sum(tx_power_m_max * cp.multiply(mu, local_remain_data / rate_m_itr)).value)
        print("func2.value", obj7_func2.value)



        problem_result = {"problem6.value": problem6.value, "alpha.value": alpha.value,"beta.value":beta.value,"omega.value":omega.value}

        return problem_result


    except SolverError:
        problem6_result = 0


        problem_result = {"problem6.value": problem6_result, "alpha.value": alpha_v,"beta.value":beta_v,"omega.value":omega_v}

        return problem_result

        pass

    #plt.show()



