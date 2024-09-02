# Molten salt fluid properties module
# Written by CdS
# 14 March 2024
 
# Import modules
from numpy import array
 
# Function MScp
def MoltenSalt(T,output):
 
    if T < 496.:
        print('T below Min, T adjusted to 496K')
        T = 496.
        check = 1
    elif T > 873.:
        print('T above Max, T adjusted to 873K')
        T = 873.
        check = 2
    else:
        check = 0
 
    c_p = 1396 + 0.172*T
    c_v = 2724.69
    k = 0.45    
    mu = 0.1800696829 - 5.005107238e-4*T + 4.697206e-7*T**2 - 1.474e-10*T**3
    rho = 2263.628 - 0.636*T
 
    solution = array([c_p, c_v, k, mu, rho, check])
 
    if output == 'All':
        return solution
    elif output == 'c_p':
        return solution[0]
    elif output == 'c_v':
        return solution[1]
    elif output == 'k':
        return solution[2]
    elif output == 'mu':
        return solution[3]
    elif output == 'rho':
        return solution[4]
   
# Test
print(MoltenSalt(500+273.15,'c_p'))
 