import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.special import comb
from scipy.linalg import expm
def mybin(x,N):
    S = bin(x).replace('0b','')
    S = S.rjust(N,'0')
    return S



N = 12
al_even = 0
al_odd = 0
even_basis_state_2 = []
even_basis_state_10 = []

odd_basis_state_2 = []
odd_basis_state_10 = []

J = 1
delta= 3
mu = 0

E_even = []
E_odd = []
G = []
# define the basis
for X in range(2**N):
    S = np.binary_repr(X, width=N)
    if S.count('1')%2==0:
        even_basis_state_2.append(S)
        even_basis_state_10.append(X)
        al_even += 1

    else:
        odd_basis_state_2.append(S)
        odd_basis_state_10.append(X)
        al_odd += 1

H_even = np.zeros((al_even, al_even))
H_odd = np.zeros((al_odd, al_odd))

even_basis_state_10 = np.array(even_basis_state_10)
odd_basis_state_10 = np.array(odd_basis_state_10)

even_basis_state_2 = np.array(even_basis_state_2)
odd_basis_state_2 = np.array(odd_basis_state_2)
Chi = []

#define the hamiltonian
for ka in range(al_even):
    for position in range(int(N-1)):
        if even_basis_state_2[ka][position]+even_basis_state_2[ka][position + 1] == '00':
            demo1 = even_basis_state_10[ka] ^ (2 ** (N - 1 - position) + 2 ** (N - 2 - position))
            x = np.where(even_basis_state_10 == demo1)
            H_even[x, ka] = -delta
            H_even[ka, x] = -delta

    for position in range(int(N-1)):
        if int(even_basis_state_2[ka][position]) + int(even_basis_state_2[ka][position + 1]) ==1:
            demo1 = even_basis_state_10[ka] ^ (2 ** (N - 1 - position) + 2 ** (N - position - 2))
            x = np.where(even_basis_state_10 == demo1)
            H_even[x, ka] = -J
            H_even[ka, x] = -J
    if even_basis_state_2[ka][0]+even_basis_state_2[ka][-1] == '00':
        demo1 = even_basis_state_10[ka] ^ (2 ** (N-1) + 2 ** (0))
        x = np.where(even_basis_state_10 == demo1)
        H_even[x, ka] = delta
        H_even[ka, x] = delta

    if int(even_basis_state_2[ka][0]) + int(even_basis_state_2[ka][-1]) == 1:
        demo1 = even_basis_state_10[ka] ^ (2 ** (N-1) + 2 ** (0))
        x = np.where(even_basis_state_10 == demo1)
        H_even[x, ka] = J
        H_even[ka, x] = J
    S = even_basis_state_2[ka]
    n = S.count('1')
    H_even[ka,ka] = 2*mu*n-mu*N

for ka in range(al_odd):
    for position in range(int(N - 1)):
        if odd_basis_state_2[ka][position] + odd_basis_state_2[ka][position + 1] == '00':
            demo1 = odd_basis_state_10[ka] ^ (2 ** (N - 1 - position) + 2 ** (N - 2 - position))
            x = np.where(odd_basis_state_10 == demo1)
            H_odd[x, ka] = -delta
            H_odd[ka, x] = -delta

    for position in range(int(N - 1)):
        if int(odd_basis_state_2[ka][position]) + int(odd_basis_state_2[ka][position + 1]) == 1:
            demo1 = odd_basis_state_10[ka] ^ (2 ** (N - 1 - position) + 2 ** (N - position - 2))
            x = np.where(odd_basis_state_10 == demo1)
            H_odd[x, ka] = -J
            H_odd[ka, x] = -J

    if odd_basis_state_2[ka][0] + odd_basis_state_2[ka][-1] == '00':
        demo1 = odd_basis_state_10[ka] ^ (2 ** (N - 1) + 2 ** (0))
        x = np.where(odd_basis_state_10 == demo1)
        H_odd[x, ka] = -delta
        H_odd[ka, x] = -delta

    if int(odd_basis_state_2[ka][0]) + int(odd_basis_state_2[ka][-1]) == 1:
        demo1 = odd_basis_state_10[ka] ^ (2 ** (N - 1) + 2 ** (0))
        x = np.where(odd_basis_state_10 == demo1)
        H_odd[x, ka] = -J
        H_odd[ka, x] = -J
    S = odd_basis_state_2[ka]
    n = S.count('1')
    H_odd[ka,ka] = 2*mu*n-mu*N

# find the eigenvalue
e_even, v_even = np.linalg.eig(H_even)
e_odd, v_odd = np.linalg.eig(H_odd)
e_energy = np.hstack((e_even,e_odd))

e_even_range = np.argsort(np.real(e_even))
e_odd_range = np.argsort(np.real(e_odd))
e_range = np.argsort(np.real(e_energy))

e_even = np.round(e_even[e_even_range],4)
e_odd = np.round(e_odd[e_odd_range],4)
e_energy = np.round(e_energy[e_range],4)
count = 0
d = 0
for n in range(int(N)):
    count = count + 1
    d = d + np.sqrt(1+mu**2-2*mu*np.cos(2*np.pi*n/N))

H_even = -H_even
