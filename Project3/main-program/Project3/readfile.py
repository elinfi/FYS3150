# This Python file uses the following encoding: utf-8
import numpy as np
import matplotlib.pyplot as plt

def read_a(filename):
    infile = open(filename, 'r')
    n_length = 7
    res_ = np.zeros(n_length)
    lambd_ = np.zeros(n_length)
    n_ = np.zeros(n_length)
    time_ = np.zeros(n_length)
    rel_err_ = np.zeros(n_length)
    i = 0
    anal = np.ones(n_length)*((5*np.pi**2)/(16**2))
    for line in infile:
        if line[0] == '-':
            print(n_, res_)
            plt.figure(1)
            #plt.plot(n_, res_ , label="lambda = %2.3f" %lambd_[2])
            #plt.legend()
            #plt.figure(2)
            plt.loglog(time_, rel_err_, label="lambda = %2.3f" %lambd_[2])
            plt.legend()

            i=0
        else:
            res, n, lambd, time, rel_err = line.split()
            #print(res, n, lambd, time)
            res = float(res)
            lambd = float(lambd)
            rel_err = float(rel_err)
            n = int(n)
            time = float(time)
            res_[i] = float(res)
            lambd_[i] = float(lambd)
            n_[i] = int(n)
            time_[i] = time
            rel_err_[i] = abs(rel_err)



            i +=1
    #plt.title("Results for varying N for Gauss-Legendre", fontsize = 20)
    #plt.xlabel("number of steps N", fontsize = 16)
    #plt.ylabel("result of integral", fontsize = 16)
    #plt.plot(n_, anal, '-', label='analytical')
    #plt.legend()
    #plt.show()
    plt.title("Relative error for varying N and lambda" , fontsize = 20)
    plt.xlabel("time [s]", fontsize = 16)
    plt.ylabel("relative error", fontsize = 16)
    plt.legend()
    plt.show()

#read_a("a_long.txt")


def read_b(filename):
    infile = open(filename, 'r')
    n_length = 7
    res_ = np.zeros(n_length)
    lambd_ = np.zeros(n_length)
    n_ = np.zeros(n_length)
    time_ = np.zeros(n_length)
    rel_err_ = np.zeros(n_length)
    i = 0
    n=0
    anal = np.ones(n_length)*((5*np.pi**2)/(16**2))
    for line in infile:
        res, n, time, err = line.split()
        print(res, n)
        res = float(res)
        n = int(n)
        time = float(time)
        err = float(err)
        res_[i] = float(res)
        n_[i] = int(n)
        time_[i] = float(time)
        rel_err_[i] = abs(float(err))

        i +=1

    plt.title("relativ error for different calculation times for laguerre")
    plt.ylabel("relative error", fontsize = 16)
    plt.xlabel("time [s]", fontsize = 16)
    print(rel_err_, time_)
    plt.loglog(time_, rel_err_, 'o')
    #plt.legend()
    plt.show()

#read_b("b.txt")

def read_c(filename):
    infile = open(filename, 'r')
    n_length = 3553
    infile.readline()
    res_ = np.zeros(n_length)
    lambd_ = np.zeros(n_length)
    n_ = np.zeros(n_length)
    variance_ = np.zeros(n_length)
    i = 0

    anal = np.ones(n_length)*((5*np.pi**2)/(16**2))
    for line in infile:
        if line[0] == '-':
            print(res_)
            plt.subplot(2,1,1)
            plt.plot(n_[:-1], res_[:-1] , label="lambda = %2.3f" %lambd_[2])
            plt.legend()

            i=0
        else:
            print(line)
            res, lambd, n, variance = line.split()
            #print(res, lambd, n)
            res = float(res)
            lambd = float(lambd)
            n = int(n)
            variance = float(variance)
            res_[i] = float(res)
            lambd_[i] = float(lambd)
            n_[i] = int(n)
            variance_[i] = variance




            i +=1
    plt.subplot(2,1,1)
    plt.plot(n_, anal, '-', label='analytical')
    plt.title("Results for varying N for simple monte carlo", fontsize = 20)
    plt.xlabel("number of samplings N", fontsize = 16)
    plt.ylabel("result of integral", fontsize = 16)
    #plt.axis([0, 5000000, 0.17, 0.46])
    plt.legend()
    plt.subplot(2,1,2)
    #plt.show()
    plt.plot(n_[:-1], variance_[:-1], '-r')
    plt.xlabel("variance")
    plt.ylabel("number of samplings N")
    plt.legend()
    plt.show()

    plt.show()
#read_c("c_extreme.txt")

def read_c(filename):
    infile = open(filename, 'r')
    n_length = int(infile.readline().strip())
    res_ = np.zeros(n_length)
    lambd_ = np.zeros(n_length)
    n_ = np.zeros(n_length)
    variance_ = np.zeros(n_length)
    i = 0

    anal = np.ones(n_length)*((5*np.pi**2)/(16**2))
    for line in infile:
        if line[0] == '-':
            plt.subplot(2,1,1)
            plt.plot(n_[3:-1], res_[3:-1] , label="lambda = %2.3f" %lambd_[2])
            plt.legend()

            i=0
        else:
            res, n, variance = line.split()
            print(res, n)
            res = float(res)

            n = int(n)
            variance = float(variance)
            res_[i] = float(res)

            n_[i] = int(n)
            variance_[i] = variance




            i +=1
    plt.subplot(2,1,1)
    plt.plot(n_, anal, '-', label='analytical')
    plt.title("Result of the improved Monte Carlo", fontsize = 20)
    plt.xlabel("number of samplings N", fontsize = 16)
    plt.ylabel("result of integral", fontsize = 16)
    #plt.axis([0, 5000000, 0.17, 0.46])
    plt.legend()
    plt.subplot(2,1,2)
    #plt.show()
    plt.plot(n_[1:-1], variance_[1:-1], '-r')
    plt.ylabel("variance")
    plt.xlabel("number of samplings N")
    plt.legend()
    plt.show()

    plt.show()
read_c("boogaloo_500000.txt")



