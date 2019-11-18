# This Python file uses the following encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt


def readfile(filename):
        infile = open(filename, 'r')
        data = infile.readlines()
        data_ = data[0].split()
        n = len(data[0])
        E_exp = []
        E_variance = []
        M_exp = []
        M_variance = []
        acceptance = []
        n = []

        for i, elem in enumerate(data):
            obj = elem.split()

            E_exp.append(float(obj[0]))
            E_variance.append(float(obj[1]))
            M_exp.append(float(obj[2]))
            M_variance.append(float(obj[3]))
            acceptance.append(float(obj[4]))
            n.append(int(obj[5]))

        return E_exp, E_variance, M_exp, M_variance, acceptance, n

#oppgave b

"""
E_exp, E_variance, M_exp, M_variance, acceptance , n= readfile('20-1-false_final.txt')


print(E_variance)
plt.plot(n, E_exp, label='false')
plt.xlabel('number of Montecarlo Cycles')
plt.ylabel('expectance energy')


#Hva er det som er interresant å plotte og se på her??


E_exp, E_variance, M_exp, M_variance, M_abs , n= readfile('20-1-true_final.txt')

print(np.sqrt(E_variance))
plt.plot(n, E_exp, label= 'true')
plt.legend()
plt.show()
"""

#PLOTTING THE NUMBER OF ACCEPTED CHANGES AGAINST MONTE CARLO CYCLES. NOTE, EACH CYCLE HAS A NEW MATRIX,
#MEANING SOME RANDOM CONFIGURATIONS MAY NEED LESS CHANGES TO COME CLOSE TO TERMINAL STATE
#THIS SHOULD CAUSE THE FLUCTUATIONS IN THE RANDOM STATE
#THE LOWER TEMPERATURES ALSO HAVE FEWER STATES THEY CAN OCCUPY (LIKE THREE OR SOMETHING), THIS MEANS
#WE EXPECT TO SEE LESS MOVEMENT IN THEESE
E_exp1_false, E_variance1_false, M_exp1_false, M_variance1_false, acceptance1_false , n1_false= readfile('20-1-false_final.txt')
plt.plot(n1_false, acceptance1_false, label='ordered')


E_exp1_true, E_variance1_true, M_exp1_true, M_variance1_true, acceptance1_true , n1_true= readfile('20-1-true_final.txt')
plt.plot(n1_false, acceptance1_true, label='random')
plt.xlabel('number of Montecarlo cycles')
plt.ylabel('accepted configuration changes')
plt.title('configuration changes for random and ordered matrices.T=1')
plt.legend()
plt.show()

E_exp24_false, E_variance24_false, M_exp24_false, M_variance24_false, acceptance24_false , n24_false= readfile('20-24-false_final.txt')
plt.plot(n24_false, acceptance24_false, label='ordered')

E_exp24_true, E_variance24_true, M_exp24_true, M_variance24_true, acceptance24_true , n24_true= readfile('20-24-true_final.txt')
plt.plot(n24_false, acceptance24_true, label='random')
plt.xlabel('number of Montecarlo cycles')
plt.ylabel('accepted configuration changes')
plt.title('configuration changes for random and ordered matrices.T=24')
plt.legend()
plt.grid='on'
plt.show()


#oppgave c
def readfile_energies(filename):
        infile = open(filename, 'r')
        data = infile.readlines()
        energies = np.zeros(len(data))
        acceptance = np.zeros(len(data))
        for i, elem in enumerate(data):
            data_ = elem.split()


            energies[i] = float(data_[0].strip())
            acceptance[i] = float(data_[1].strip())
        return energies, acceptance



data24, acceptance24 = readfile_energies('20-2.400000-1-ENERGIES_acceptance.txt')
data1, acceptance1 = readfile_energies('20-1.000000-1-ENERGIES_acceptance.txt')
a = np.hstack((data1   [5000:])) #slicing data to remove the nasty early values
print(a)
plt.hist(a, bins=5, density='true')
plt.xlabel('energy [E/J]')
plt.ylabel('probability')
plt.title('Histogram of energies for T=1.0')

mean = np.mean(data1)
sigma = np.std(data1)
print(sigma)
print(mean)

x24 = np.arange(0, len(acceptance24))
x1 = np.arange(0, len(acceptance1))

#plt.plot(np.linspace(mean-sigma,mean+sigma,2), np.linspace(0.0048, 0.0048,2),markersize = 20,marker ="|",color = "r",label = "Standard deviation")
plt.show()


#LAGRE ACCEPTED ENDRINGER FOR HVER ENESTE MC CYCLE, OG PLOTTE
#MOT ANTALL N SOM I OPPGAVE C??


#oppgave C 2

readfile_energies('20-1.000000-0-ENERGIES_acceptance.txt')


#Hva er det vi ønsker å tolke ut av de to forskjellige ordnet/random matrise plottene?
"""
Vi skal se om ordnet matrise er nærmere terminalisert enn en tilfeldig matrise.
Eksempelvis er en 2x2 matrise som er ordnet allerede terminalisert.
Går en ordnet 20x20 matrise fortere til terminalisering?
"""

#Hvilken energi er det de ber oss plotte i oppgave d? og hva er intensjonen med sansynligheten P(E)??
"""
Vi ønsker å se hvor sansynlig en energi state er. Lagre alle energi (IKKE EXPECTED ENERGY) og skriv disse til fil.
plott disse så i et histogram
"""
#Hvorfor regner vi ut De før vi flipper??
"""
Man kan bare se på et atom og regne ut hva energien til denne er hvis den er flippet, uten å flippe den. Forskjellen
er jo bare et fortegn, som man lett kan legge inn i utregningen av De.
"""

#SKrive hver energi for hver eneste MC syklus, eller bare skrive ut average energy for hver N-te antall sykluser?
