import cvxpy as cp
import numpy as np
import matplotlib.pyplot as plt
from cvxpy import SolverError

# import simulation_environment
import matplotlib
import random
#import compute_rate_assist


# random.seed(0)

def solve_problem_func(env):
    # ----------parameters----------
    I = env["I"]
    M = env["M"]
    h_ul = env["h_ul"]
    bandwidth_max = env["bandwidth_max"]
    tx_power_m_max = env["tx_power_m_max"]
    comp_bs = env["comp_bs"]
    comp_m = env["comp_m"]
    task = env["task"]
    data = env["data"]
    M_list = env["M_list"]
    time_scale = env["time_scale"]
    rate_m_itr = bandwidth_max * np.log2 ( 1+tx_power_m_max * h_ul )

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
    c1 = [alpha >= 1e-5, alpha <= 1]
    # bandwidth allocation constraints
    c4_5 = [cp.sum(omega) <= 1, omega >= 1e-5]
    # bs server allocation constraints
    c2_3 = [cp.sum(beta) <= 1, beta >= 1e-5]

    c6 = [cp.multiply ( 1-alpha, task )/comp_m*time_scale<=t1]
    c7 = [cp.multiply(eta,task/comp_bs)*time_scale<=t2]
    c8 = [cp.multiply ( mu, data / rate_m_itr)*time_scale<=t3]
    R1a2 = []
    for i in range(0,M):
        for j in range(0,I):

            R1a2 =R1a2 + [cp.norm(cp.vstack([2*(alpha[i,j]),eta[i,j]-beta[i,j]]),2)<=eta[i,j]+beta[i,j]]
            #assert  R1a2[0].is_dqcp()

    R1a3 = [alpha>=eta+beta-1, alpha<=beta,alpha<=eta]

    R2a3 = [alpha>=omega+mu-1,alpha<=omega,alpha<=mu]


    # local computing constraints R1
    #R1 = [cp.ceil(cp.multiply(alpha, task)/comp_m) <= T]
    R1 = [cp.ceil ( cp.multiply ( 1-alpha, task ) / comp_m ) -T <= 0]
    #R1 = [cp.ceil(alpha) - 10 <= 0]
    #assert R1[0].is_dqcp ()
    # offloading constraints R2
    R2a = cp.ceil(cp.multiply(cp.multiply(alpha,task/comp_bs),cp.inv_pos(beta)))
    R2b = cp.ceil(cp.multiply(cp.multiply(alpha,data/bandwidth_max*np.log2 ( 1+tx_power_m_max * h_ul )),cp.inv_pos(omega)))
    R2 = [R2a+R2b <= T]
    R2a1 = cp.ceil((cp.multiply(eta,task/comp_bs)))

    objective2 = cp.Minimize(1/2*cp.max(cp.ceil ( cp.multiply ( alpha, task ) / comp_m )))
    objective2 = cp.Minimize(cp.max(cp.ceil(cp.multiply(1-alpha, task/beta/comp_bs))))


    obj3_part1 = 1/2*cp.max(cp.ceil(cp.multiply(1-alpha,task/comp_m)))
    obj3_part1 = 1 / 2 * cp.ceil ( cp.max ( cp.multiply ( 1-alpha, task / comp_m ) ) )
    obj3_part2 = 1/2*cp.max(cp.ceil(cp.multiply(eta,task/comp_bs))+cp.ceil(cp.multiply(mu,data/rate_m_itr)))
    obj3_part2 = 1 / 2 * cp.ceil (cp.max( cp.multiply ( eta, task / comp_bs ) ))+1/2*cp.ceil ( cp.max(cp.multiply ( mu, data / rate_m_itr ) ) )

    obj3_part3 = cp.abs(obj3_part1-obj3_part2)
    objective3 = cp.Minimize(obj3_part1+obj3_part2+obj3_part3)


    #assert obj3_part1.is_dqcp ()
    #assert obj3_part1.is_dcp()
    #assert (cp.max(cp.ceil(cp.multiply(eta,task/comp_bs)))).is_dqcp()
    #assert (cp.max(cp.ceil(cp.multiply(eta,task/comp_bs)))).is_dcp()
    #assert (cp.ceil(cp.multiply(mu,data/rate_m_itr))).is_dqcp
    #assert (cp.ceil(alpha)+cp.ceil(beta)).is_dqcp()
    #assert obj3_part2.is_dqcp ()
    #assert obj3_part3.is_dqcp ()
    #assert objective3.is_dqcp ()


    t = cp.Variable([M,I])
    objective2 = cp.Minimize ( cp.max ( cp.ceil ( cp.multiply ( t, task / comp_bs ) ) ) )
    c4 = [(1-alpha)/beta <= t]

    objective6 = cp.Minimize(1/2*cp.ceil(t1+t2+t3+cp.abs(t1-t2-t3))+3/2)

    rho = 1
    upsilon = 1
    varsigma = 1e-27
    obj7_func1 = 1/2*cp.ceil(t1+t2+t3+cp.abs(t1-t2-t3))+3/2
    obj7_func1 = 1 / 2 * (t1 + t2 + t3 + cp.abs(t1 - t2 - t3))
    obj7_func2 = cp.sum(cp.multiply(1-alpha,varsigma*task*np.square(comp_m)))+cp.sum(tx_power_m_max*cp.multiply ( mu, data / rate_m_itr))
    objective7 = cp.Minimize(rho*obj7_func1+upsilon*obj7_func2)



    # ----------probalem solve and results----------
    problem6 = cp.Problem(objective6, c1+c2_3+c4_5+c6+c7+c8+R1a3+R2a3)
    problem6.solve(qcp=True, verbose=True)
    print()
    print("problem1 solve: ", problem6.value)
    print("alpha.value", alpha.value)
    print("beta.value", beta.value)
    print("omega.value", omega.value)


    y = np.ceil(data / rate_m_itr/M*10)
    # ----------data collection and depict-----------
    plt.subplot(331)
    plt.title("user_rate")
    plt.bar(x=M_list, height=(rate_m_itr/1e6/M).reshape(M), width=1)

    plt.subplot(332)
    plt.title("average bandwidth allocation full_transmit_delay")
    plt.bar(x=M_list, height=y.reshape(M), width=1)

    plt.subplot(333)
    plt.title("local_computing_delay")
    plt.bar(x=M_list, height=np.ceil(time_scale*task / comp_m).reshape(M), width=1)

    plt.subplot(334)
    plt.title("offloading_computing_delay")
    plt.bar(x=M_list, height=(np.ceil(task *M/ comp_bs)*time_scale).reshape(M), width=1)

    plt.subplot(335)
    plt.title("optimized local computing_delay gain")
    plt.plot(M_list, (problem6.value-np.ceil(task / comp_m*time_scale)).reshape(M),color='r')
    plt.plot(M_list, t1.value*np.ones(M),color='g')
    plt.plot(M_list, t2.value * np.ones(M),color='k')
    plt.plot(M_list, t3.value * np.ones(M),color='b')

    plt.subplot(336)
    plt.title("optimized alpha beta")
    plt.plot(M_list, alpha.value ,'-v')
    plt.plot(M_list, beta.value,'-x')

    plt.subplot(337)
    plt.title("optimized omega")
    plt.plot(M_list, omega.value,'-^')

    print("func1.value", obj7_func1.value)
    print("local computing energy", cp.sum(cp.multiply(1-alpha,varsigma*task*np.square(comp_m))).value)
    print("offloading energy", cp.sum(tx_power_m_max*cp.multiply ( mu, data / rate_m_itr)).value)
    print("func2.value", obj7_func2.value)




    plt.show()


