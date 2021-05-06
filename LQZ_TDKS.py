#!/usr/bin/python3.7
import numpy as np
from numpy import linalg as LA
import os
from scipy import sparse
from scipy.sparse import linalg
from scipy.sparse.linalg import splu
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
#from lqz-plot import td_plot

class constant:
    # M = num of sites
    M = 50
    # so total 1D length is l = M * ∆x
    l = 0.001
    # ∆M = length of two near site: x_i+1 - x_i
    delta_x = l / (M - 1)
    # NOTICE: less propergate, higher accurate, not time smaller interval
    num_t = 100
    delta_t = 0.000001
    # t is time
    t_tot = (num_t - 1) * delta_t
    A = 1
    omega = 0.5

def fsys():
# get current working dir
    work_dir = os.getcwd()
    dir_tdevec = work_dir + '/TDEVEC'
    if (os.path.isdir(dir_tdevec)):
        print("dir TDEVEC exist")
    else:
# all strings are True
        print("really biult TDEVEC directory and remove old one? True or 1 means do it, no input means don't")
        buildornot = bool(input())
        if (buildornot):
            #print(buildornot)
        #    os.rmdir(dir_tdevec)
            os.mkdir(dir_tdevec)
        else:
            return dir_tdevec, buildornot
    return dir_tdevec, buildornot

# define 1d equidistant grid to store real lattice points.
def gen_grid(M, l):
    x = np.linspace(0, l, M, endpoint = True)
    #print(x)
    return x

# define single-partical Hamiltonian.
#Kin = (phi(j+1) - 2 * phi(j) + phi(j-1))/(delta_x * delta_x)
def gen_kin(M, delta_x):
    #Kin = np.zeros([ M, M ])
    Kin = 2 * np.identity(M)
    for i in range(0, M-1):
        Kin[ i ][ i+1 ] = -1.0
        Kin[ i+1 ][ i ] = -1.0
    K = 0.5 * Kin / delta_x**2
    return K
# Potential Part
# use harmonic ocsillator
def gen_pot(x):
    vx = 0.5 * x**2
    V = np.diag(vx)
    return V 

# static Vex
def gen_vex0(M):
    vex = 0.3 * np.identity(M)
    return vex

# TD-External potential Vex (M)
# v(x_j, t)=A*f(t)*x_j*sin(ωt)
# f(t) is the envelope function
def gen_vex(x, M, t):
    ft = np.exp(-(t**2)/2)
    Vex = constant.A * ft * x * np.sin(constant.omega * t)


# solve GS eigen problem and normalized
def h_solve(x, delta_x, M, H):
    #print(H)
# GS Hamiltonian
    #H0 = H0_k + H0_v
# eigh solve diag, eig solve generous 
    evalue, estate = LA.eigh(H)
    #print(estate)
# Normalized eigenstate
    wfnorm = 0.0
    for i in range(M):
        wfnorm = wfnorm + delta_x * np.vdot(estate[ i ], estate[ i ])
    #print(wfnorm)
# linalg.norm returns the sqrt(a^2+b^2+...)=normal
    #normal = LA.norm(estate)
    #est_norm = estate/normal
    est_norm = estate/np.sqrt(wfnorm)
    #print(normal, np.sqrt(wfnorm))
    return evalue, est_norm

# define Crank–Nicholson algorithm
# e^(-iĤΔτ)≈(1-iĤ∆τ/2)/(1+iĤΔτ/2)
# (1+i*Ĥ(τ_(j+½))∆τ/2)ψ(τ_(j+1))=(1-i*Ĥ(τ_(j+½))∆τ/2)ψ(τ_j)
def Crank_Nicholson(delta_t, psi_j0, Ht_j0, M):
    Ht_1 = 1 - 1j * Ht_j0 * delta_t * 0.5
    Ht_2 = 1 + 1j * Ht_j0 * delta_t * 0.5
    sumpsi=0.0
    PSI = np.zeros([ M, M ], complex)
# method1: solve linear function
  #  for i in range(M):
  #      B = Ht_1.dot(psi_j0[ i ])
  #      PSI[ i ] = LA.solve(Ht_2, B)
  #      #PSI[ i ] = sparse.linalg.spsolve(Ht_2, B)
  #      sumpsi=sumpsi+constant.delta_x * np.vdot(PSI[ i ], PSI[ i ])
        #sumpsi = sumpsi + constant.delta_x * np.vdot(np.matmul(np.linalg.inv(Ht_2), B),np.matmul(np.linalg.inv(Ht_2), B))
#
# method2: use ψ_j+1 = UA^-1*UB * ψ_j
    U = np.matmul(np.linalg.inv(Ht_2), Ht_1)
    #PSI = U.dot(psi_j0)
    PSI = np.matmul(U, psi_j0)
    for i in range(M):
        sumpsi=sumpsi+constant.delta_x * np.vdot(PSI[ i ], PSI[ i ])

    print(sumpsi)
    #b = Ht_1 * phi_j0
    #print(Ht_2, b)
    #phi_t = LA.solve(Ht_2, b)
    #phi_t = LA.tensorsolve(Ht_2, b)
    return PSI

# Aψ_j+1=Bψ_j
def new_CN(delta_t, psi_j0, Ht_j0, M):
    Ht_1 = 1 - 1j * Ht_j0 * delta_t * 0.5
    Ht_2 = 1 + 1j * Ht_j0 * delta_t * 0.5
    U1 = sparse.csr_matrix(Ht_2)   # Here's the initialization of the sparse matrix.
    U2 = sparse.csr_matrix(Ht_1)
    #PSI = np.zeros((M, num_t),complex) # M times num_t array to store all solutions
    #PSI[:,0] = psi_j0 # psi(x,0)
    #LU = sparse.linalg.splu(U1) # compute LU-decomposition of U1
    #for t in range(num_t): # loop over time-steps
    #    B = U2.dot(PSI[:, t])
    #    PSI[:, t + 1] = LU.solve(B) # solve system of equations for each time step
    PSI = np.zeros([ M, M ], complex)
    LU = sparse.linalg.splu(U1) # compute LU-decomposition of U1
    sumpsi=0.0
    for i in range(M):
        B = U2.dot(psi_j0[ i ])
        PSI[ i ] = LU.solve(B) # solve system of equations for each time step
        #sumpsi=sumpsi+constant.delta_x*np.vdot(PSI[ i ], PSI[ i ])
    #print(sumpsi)
    return PSI




# iterate calculate C-N
#def iterate_cal():

def check_norm(phi, M):
# norm of wf at each time steps, and check it conserved.
    wfnorm = 0.0
    for i in range(M):
        wfnorm = wfnorm + constant.delta_x * np.vdot(phi[ i ], phi[ i ])
    with open('checknorm.dat', 'a+') as f_check1:
        f_check1.write("%.8f    %.8f\n" %(np.real(wfnorm), np.imag(wfnorm)))
    #print(wfnorm)

# check under static Vex, ψ(t)=ψ(0)*exp(-iεt)
def check_tdwf(phi0, eva0, t_tot, phit, M, x):
    print("print the check_tdwf or not, if not don't input anything?")
    write_check_tdwf = input()
    for i in range(M):
        wf_check = phi0[ :, i ] * np.exp(-1j * eva0[ i ] * t_tot)
        wf_tddft = phit[ :, i ]
        f, ax = plt.subplots()
        ax.plot(x, wf_check, color='red')
        ax.plot(x, wf_tddft, color='blue')
        plt.savefig("wf_check"+str(i)+".png")
        #if (os.path.isfile("wf_check.dat")):
        #    print("wf_check.dat exist, delet?")
        #    wff_del = input()
        if (write_check_tdwf):
            with open("wf_check_iet.dat", 'a+') as f_check2:
                f_check2.write("phi_%d\n" %(i))
                for item1 in wf_check:
                    f_check2.write("%.8f    %.8f\n" %(np.real(item1), np.imag(item1)))
                f_check2.write("\n")
            with open("wf_check_cn.dat", 'a+') as f_check3:
                f_check3.write("phi_%d\n" %(i))
                for item2 in wf_tddft:
                    f_check3.write("%.8f    %.8f\n" %(np.real(item2), np.imag(item2)))
                f_check3.write("\n")


def main():
# File system
    dir_tdevec, buildornot = fsys()
    x = gen_grid(constant.M, constant.l)
    H0_k = gen_kin(constant.M, constant.delta_x)
    H0_v = gen_pot(x)
    H0 = H0_k + H0_v
    eva0, est0 = h_solve(x, constant.delta_x, constant.M, H0)
# Initialize 1st step est, and phi
    Vex = gen_vex0(constant.M)
    Ht = H0_k + H0_v + Vex
    eva_t, est_t = h_solve(x, constant.delta_x, constant.M, Ht)
    print(Ht)
    print(eva_t)
    phi_t = est0.astype(complex)
    #phi = np.zeros([ len(est0[ 0 ]), len(est0[ 1 ]) ])
# instantaneous eigenstate INST_EIGST
    INST_EIGST = np.zeros((constant.M, constant.M, constant.num_t),complex)
    INST_EIGST[ :, :, 0 ] = phi_t
    #phi_t = Crank_Nicholson(constant.delta_t, phi_0, Ht, constant.M)
    #INST_EIGST[ :, :, 1 ] = phi_t
    print("\n")
# Time iteration
    for t in range(1, constant.num_t, 1):
        INST_EIGST[ :, :, t ] = Crank_Nicholson(constant.delta_t, phi_t, Ht, constant.M)
        #INST_EIGST[ :, :, t ] = new_CN(constant.delta_t, phi_t, Ht, constant.M)
        phi_t = INST_EIGST[ :, :, t ]
        #check_norm(phi_t, constant.M)
        #check_tdwf(est_t, eva_t, constant.t_tot, phi_t, constant.M)
    check_tdwf(INST_EIGST[ :, :, 0 ], eva_t, constant.t_tot, INST_EIGST[ :, :, constant.num_t - 1 ], constant.M, x)
# Ploting wf
  #  for ti in range(constant.num_t):
  #      for mi in range(constant.M):
  #          #td_plot.wf_plot(x, INST_EIGST[ :, :, 9 ], ax, f)
  #          f, ax = plt.subplots()
  #          ax.plot(x, INST_EIGST[ :, mi, ti ], color='blue')
  #  #        ax.plot(x, v
  #          plt.savefig("wf"+"_s"+str(mi)+"_t"+str(ti)+".png")
    if (buildornot):
        for t in range(constant.num_t):
            td_file = dir_tdevec + '/TDEVEC_step'+str(t)+'.dat'
            np.savetxt(td_file, phi_t, fmt='%.8f')

if __name__ == "__main__":
    main()
