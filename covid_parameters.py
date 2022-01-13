import numpy as np
import pandas as pd
from covid_age_dependence import *
from covid_interventions import *


beta = 0.1088   # transmission probability
#betas adjusted for R0 = 1, 1.5, 2, ..., 9.5, 10
betas_adjusted_for_R0 = np.array([0.0272, 0.0408, 0.0544, 0.0679, 0.0815,
                                     0.0951, 0.1088, 0.1223, 0.1359, 0.1495,
                                      0.1631, 0.1767, 0.1902, 0.2038, 0.2174,
                                       0.2310, 0.2446, 0.2582, 0.2718])

beta_R0_32 = 0.0870

epsilon = 0.2632    # latent period in days^-1


r_s = 0.2941    # infectious period (symptomatic) in days^−1
r_a = 0.2941    # infectious period (assymptomatic) in days^−1


theta_NIL = get_age_dependent_parameter_arrays('theta', waning='NIL')
theta_1Y = get_age_dependent_parameter_arrays('theta', waning='1Y')
theta_3Y = get_age_dependent_parameter_arrays('theta', waning='3Y')


p = 0.5550  # proportion of symptomatic cases


rho_NIL = get_age_dependent_parameter_arrays('rho', waning='NIL')
rho_1Y = get_age_dependent_parameter_arrays('rho', waning='1Y')
rho_3Y = get_age_dependent_parameter_arrays('rho', waning='3Y')


pi_NIL = get_age_dependent_parameter_arrays('pi', waning='NIL')
pi_1Y = get_age_dependent_parameter_arrays('pi', waning='1Y')
pi_3Y = get_age_dependent_parameter_arrays('pi', waning='3Y')


tau_NIL = get_age_dependent_parameter_arrays('tau', waning='NIL')
tau_1Y = get_age_dependent_parameter_arrays('tau', waning='1Y')
tau_3Y = get_age_dependent_parameter_arrays('tau', waning='3Y')


omega_NIL = get_age_dependent_parameter_arrays('omega', waning='NIL')
omega_1Y = get_age_dependent_parameter_arrays('omega', waning='1Y')
omega_3Y = get_age_dependent_parameter_arrays('omega', waning='3Y')


mu_NIL = get_age_dependent_parameter_arrays('mu', waning='NIL')
mu_1Y = get_age_dependent_parameter_arrays('mu', waning='1Y')
mu_3Y = get_age_dependent_parameter_arrays('mu', waning='3Y')


alpha_A = 0.5000    # asymptomatic reduction in transmission
alpha_S = 1


#### Intervention parameters

ypsilon = vaccination_rate_for_all_age_groups()
ypsilon_5_to_11_yo_vaccinated = vaccination_rate_for_all_age_groups(vaccination_plan_for_5_to_11_yo=[221, 278], vaccination_coverage_for_5_to_11_yo=0.8)

e = 0.7000  # reduction of susceptibility upon exposure


f = 0.8333  # reduction of the risk of symptomatic disease upon acquisition


li_and_ti_NIL = get_contact_reductions(waning='NIL')
li_and_ti_1Y = get_contact_reductions(waning='1Y')
li_and_ti_3Y = get_contact_reductions(waning='3Y')



#### Exploratory parameters


gamma_NIL = 0   # rate of loss of natural immunity (No-immunity-Loss)
gamma_1Y = 1/365    # rate of loss of natural immunity (1Y immunity)
gamma_3Y = 1/1095   # rate of loss of natural immunity (3Y immunity)


gamma_v_NIL = 0     # rate of loss of immunity conferred by vaccination (No-immunity-Loss)
gamma_v_1Y = 1/365     # rate of loss of immunity conferred by vaccination (1Y immunity)
gamma_v_3Y = 1/1095     # rate of loss of immunity conferred by vaccination (3Y immunity)


#####
NPI_CONTACTS_REDUCTION = NPI_CONTACTS_REDUCTION
#####



def load_parameters(vaccinate_5_to_11_yo=False):
    if vaccinate_5_to_11_yo:
        yps = ypsilon_5_to_11_yo_vaccinated
    else:
        yps = ypsilon
    return (beta, epsilon, r_s, r_a, theta_NIL, theta_1Y, theta_3Y, p, rho_NIL, rho_1Y, rho_3Y,
            pi_NIL, pi_1Y, pi_3Y, tau_NIL, tau_1Y, tau_3Y, omega_NIL, omega_1Y, omega_3Y,
            mu_NIL, mu_1Y, mu_3Y, alpha_S, alpha_A, yps, e, f, li_and_ti_NIL, li_and_ti_1Y, li_and_ti_3Y,
            gamma_NIL, gamma_1Y, gamma_3Y, gamma_v_NIL, gamma_v_1Y, gamma_v_3Y, NPI_CONTACTS_REDUCTION)

def get_waning_dependent_gamma(waning):
    if waning == 'NIL':
        return gamma_NIL, gamma_v_NIL
    elif waning == '1Y':
        return gamma_1Y, gamma_v_1Y
    elif waning == '3Y':
        return gamma_v_3Y, gamma_v_3Y
    else:
        exit('Please choose between no immunity loss (NIL), 1 year immunity (1Y) and 3 year immunity (3Y)')

def get_waning_dependent_li_and_ti(waning):
    if waning == 'NIL':
        return li_and_ti_NIL 
    elif waning == '1Y':
        return li_and_ti_1Y 
    elif waning == '3Y':
        return li_and_ti_3Y 
    else:
        exit('Please choose between no immunity loss (NIL), 1 year immunity (1Y) and 3 year immunity (3Y)')

def load_waning_dependent_parameters(waning='NIL'):
    theta, pi, tau, mu, rho, omega = get_age_dependent_parameter_arrays(waning=waning)
    gamma, gamma_v = get_waning_dependent_gamma(waning=waning)
    li_and_ti = get_waning_dependent_li_and_ti(waning=waning)
    return  theta, pi, tau, mu, rho, omega, gamma, gamma_v, li_and_ti
    

def load_initial_conditions():
    return get_age_dependent_initial_conditions()

def vaccination_rate_from_age(vaccinations_rates, age=2):
    for i in range(len(vaccinations_rates)):
        if i < age:
            vaccinations_rates[i] = 0
    return vaccinations_rates

def load_contact_matrix():
    df = pd.read_excel('data\\contacts_age_groups_portugal.xlsx', header=None)
    df = df.to_numpy()
    return df

def lambda_of_age_group(age_group, beta, alpha_S, contact_matrix, I_S_j, I_S_v_j, N_j, alpha_A, I_A_j, I_A_v_j):
    contacts_of_age_group = contact_matrix[age_group]
    sum = 0
    for j in range(len(contact_matrix)):
        sum += alpha_S * contacts_of_age_group[j] * (I_S_j + I_S_v_j) / N_j + alpha_A * contacts_of_age_group[j] * (I_A_j + I_A_v_j) / N_j
    return beta * sum

def all_lambdas(beta, alpha_S, contact_matrix, I_S, I_S_v, N, alpha_A, I_A, I_A_v):
    result = []
    for i in range(len(contact_matrix)):
        contacts_of_age_group = contact_matrix[i]
        sum = 0
        for j in range(len(contact_matrix)):
            sum += alpha_S * contacts_of_age_group[j] * (I_S[j] + I_S_v[j]) / N[j] + alpha_A * contacts_of_age_group[j] * (I_A[j] + I_A_v[j]) / N[j]
        result.append(beta * sum)
    return np.array(result)

def get_beta_for_R0(R0):
    if R0 == 3.2:
        return beta_R0_32
    elif R0 % 0.5 == 0 and R0 >= 1 and R0 <= 10:
        index = int(R0 * 2 - 2)
        return betas_adjusted_for_R0[index]
    else:
        exit(f'R0 {R0} not supported')
