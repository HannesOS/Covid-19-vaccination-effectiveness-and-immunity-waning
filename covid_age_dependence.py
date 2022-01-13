import numpy as np

# no_immunity_loss = NIL
# 1_year_immunity = 1Y
# 3_year immunity = 3Y

age_group_names = np.array(['0 - 4', '5 - 9', '10 - 14', '15 - 19', '20 - 24', '25 - 29',
                        '30 - 34', '35 - 39', '40 - 44', '45 - 49', '50 - 54', '55 - 59',
                         '60 - 64', '65 - 69', '70 - 74', '75+'])

#Age dependet parameters for each scenario in the form param_array_scenario = [param_agegroup0, param_agegroup1,..., param_agegroup15]
theta_array_NIL = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.03, 0.03, 0.15, 0.15, 0.28, 0.28])
pi_array_NIL = np.array([0.04, 0.00, 0.02, 0.00, 0.03, 0.07, 0.07, 0.18, 0.11, 0.15, 0.16, 0.13, 0.12, 0.13, 0.08, 0.03])
tau_array_NIL = np.array([0.00, 0.00, 0.00, 0.00, 0.01, 0.00, 0.00, 0.00, 0.02, 0.03, 0.01, 0.06, 0.07, 0.08, 0.18, 0.33])
mu_array_NIL = np.array([0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.11, 0.07, 0.04, 0.16, 0.14, 0.24, 0.24, 0.34, 0.53])
rho_array_NIL = np.array([0.17, 0.17, 0.12, 0.12, 0.14, 0.14, 0.12, 0.12, 0.10, 0.10, 0.09, 0.09, 0.08, 0.08, 0.07, 0.07])
omega_array_NIL = np.array([0.03, 0.03, 0.07, 0.07, 0.05, 0.05, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04])

theta_array_1Y = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.10, 0.10, 0.22, 0.22])
pi_array_1Y = np.array([0.05, 0.00, 0.03, 0.01, 0.03, 0.07, 0.04, 0.11, 0.12, 0.17, 0.20, 0.16, 0.15, 0.16, 0.10, 0.03])
tau_array_1Y = np.array([0.00, 0.00, 0.00, 0.00, 0.01, 0.00, 0.00, 0.00, 0.02, 0.03, 0.02, 0.06, 0.08, 0.09, 0.18, 0.34])
mu_array_1Y = np.array([0.09, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.11, 0.07, 0.04, 0.16, 0.14, 0.23, 0.23, 0.33, 0.51])
rho_array_1Y = np.array([0.17, 0.17, 0.12, 0.12, 0.14, 0.14, 0.12, 0.12, 0.10, 0.10, 0.09, 0.09, 0.08, 0.08, 0.07, 0.07])
omega_array_1Y = np.array([0.03, 0.03, 0.07, 0.07, 0.05, 0.05, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04])

theta_array_3Y = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.15, 0.15, 0.27, 0.27])
pi_array_3Y = np.array([0.03, 0.00, 0.02, 0.00, 0.02, 0.05, 0.08, 0.21, 0.16, 0.21, 0.15, 0.12, 0.11, 0.12, 0.08, 0.03])
tau_array_3Y = np.array([0.00, 0.00, 0.00, 0.00, 0.01, 0.00, 0.00, 0.00, 0.02, 0.03, 0.01, 0.06, 0.07, 0.08, 0.17, 0.32])
mu_array_3Y = np.array([0.09, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.11, 0.07, 0.04, 0.16, 0.14, 0.23, 0.23, 0.33, 0.51])
rho_array_3Y = np.array([0.17, 0.17, 0.12, 0.12, 0.14, 0.14, 0.12, 0.12, 0.10, 0.10, 0.09, 0.09, 0.08, 0.08, 0.07, 0.07])
omega_array_3Y = np.array([0.03, 0.03, 0.07, 0.07, 0.05, 0.05, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04])


#Age dependet initial conditions for each compartment in the form comparment0_array = [compartment0_agegroup0, compartment0_agegroup1,.., compartment0_agegroup15]
S0_array = np.array([386989.00, 402879.00, 454523.00, 490150.00, 471175.00, 468934.00, 465720.00, 552733.00, 672346.00, 676536.00, 637053.00, 633486.00, 586943.00, 541284.00, 508977.00, 1014532.00])
E0_array = np.array([1144.00, 1631.00, 2014.00, 2436.00, 2906.00, 2764.00, 2798.00, 3263.00, 3910.00, 4130.00, 3732.00, 3535.00, 3092.00, 2295.00, 1924.00, 6077.00])
IA0_array = np.array([811.00, 1337.00, 1472.00, 1794.00, 2030.00, 2018.00, 1987.00, 2372.00, 2769.00, 2925.00, 2649.00, 2339.00, 2004.00, 1519.00, 1170.00, 3433.00])
IS0_array = np.array([1012.00, 1668.00, 1836.00, 2238.00, 2532.00, 2517.00, 2478.00, 2958.00, 3453.00, 3648.00, 3304.00, 2917.00, 2499.00, 1895.00, 1459.00, 4281.00])
H0_array = np.array([8.00, 9.00, 4.00, 5.00, 20.00, 20.00, 49.00, 59.00, 128.00, 129.00, 274.00, 272.00, 509.00, 468.00, 1226.00, 2473.00])
ICU0_array = np.array([0.00, 0.00, 0.00, 0.00, 3.00, 2.00, 6.00, 7.00, 31.00, 32.00, 80.00, 79.00, 143.00, 132.00, 81.00, 164.00])
M0_array = np.array([1.00, 0.00, 0.00, 2.00, 0.00, 8.00, 11.00, 13.00, 30.00, 63.00, 95.00, 184.00, 361.00, 570.00, 844.00, 8602.00])
R0_array = np.array([46237.00, 48319.00, 45091.00, 48697.00, 71778.00, 71417.00, 93545.00, 111017.00, 101557.00, 102270.00, 97991.00, 97329.00, 81211.00, 74749.00, 33910.00, 68359.00])


def get_age_group_names():
    return age_group_names

def get_age_dependent_parameter_arrays(parameter=None, waning='NIL'):
    if(waning == 'NIL'):
        if parameter == None:
            return theta_array_NIL, pi_array_NIL, tau_array_NIL, mu_array_NIL, rho_array_NIL, omega_array_NIL
        elif parameter == 'theta':
            return theta_array_NIL
        elif parameter == 'pi':
            return pi_array_NIL
        elif parameter == 'tau':
            return tau_array_NIL
        elif parameter == 'mu':
            return mu_array_NIL
        elif parameter == 'rho':
            return rho_array_NIL
        elif parameter == 'omega':
            return omega_array_NIL
        else:
            exit(f'Parameter {parameter} does not exist.')
    elif waning == '1Y':
        if parameter == None:
            return theta_array_1Y, pi_array_1Y, tau_array_1Y, mu_array_1Y, rho_array_1Y, omega_array_1Y
        elif parameter == 'theta':
            return theta_array_1Y
        elif parameter == 'pi':
            return pi_array_1Y
        elif parameter == 'tau':
            return tau_array_1Y
        elif parameter == 'mu':
            return mu_array_1Y
        elif parameter == 'rho':
            return rho_array_1Y
        elif parameter == 'omega':
            return omega_array_1Y
        else:
            exit(f'Parameter {parameter} does not exist.')
    elif waning == '3Y':
        if parameter == None:
            return theta_array_3Y, pi_array_3Y, tau_array_3Y, mu_array_3Y, rho_array_3Y, omega_array_3Y
        elif parameter == 'theta':
            return theta_array_3Y
        elif parameter == 'pi':
            return pi_array_3Y
        elif parameter == 'tau':
            return tau_array_3Y
        elif parameter == 'mu':
            return mu_array_3Y
        elif parameter == 'rho':
            return rho_array_3Y
        elif parameter == 'omega':
            return omega_array_3Y
        else:
            exit(f'Parameter {parameter} does not exist.')
    else:
        exit('Please choose between no immunity loss (NIL), 1 year immunity (1Y) and 3 year immunity (3Y)')

def get_age_dependent_paramters_of_age_group(age_group=0, waning='NIL'):
    return [el[age_group] for el in get_age_dependent_parameter_arrays(waning=waning)]
        
def print_age_dependent_parameter_table(min_agegroup=0, max_agegroup=len(age_group_names), waning='NIL'):
    table = []
    for age in range(min_agegroup, max_agegroup):
        row = get_age_dependent_paramters_of_age_group(age, waning)
        table.append(row)
        print(row)
    return table

def get_age_dependent_initial_conditions(compartment=None):
    S_v0_array = np.zeros_like(S0_array)
    E_v0_array = np.zeros_like(E0_array)
    IA_v0_array = np.zeros_like(IA0_array)
    IS_v0_array = np.zeros_like(IS0_array)
    H_v0_array = np.zeros_like(H0_array)
    ICU_v0_array = np.zeros_like(ICU0_array)
    M_v0_array = np.zeros_like(M0_array)
    R_v0_array = np.zeros_like(R0_array)
    if compartment == None:
        return (S0_array, E0_array, IA0_array, IS0_array, H0_array, ICU0_array, M0_array, R0_array,
                S_v0_array, E_v0_array, IA_v0_array, IS_v0_array, H_v0_array, ICU_v0_array, M_v0_array, R_v0_array) 
    elif compartment == 'S':
        return S0_array
    elif compartment == 'E':
        return E0_array
    elif compartment == 'IA':
        return IA0_array
    elif compartment == 'IS':
        return IS0_array
    elif compartment == 'H':
        return H0_array
    elif compartment == 'ICU':
        return ICU0_array
    elif comparment == 'M':
        return M0_array
    elif compartment == 'R':
        return R0_array
    elif compartment == 'N':
        return N
    else:
        exit('Compartment does not exist.')

def get_age_dependent_initial_conditions_of_age_group(age_group=0):
    return [el[age_group] for el in get_age_dependent_initial_conditions()]
