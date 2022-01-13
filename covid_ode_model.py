import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
from covid_age_dependence import age_group_names
from covid_parameters import load_parameters, load_waning_dependent_parameters, load_contact_matrix, get_beta_for_R0, load_initial_conditions, vaccination_rate_from_age, all_lambdas
from covid_reproduction import R_het
from covid_interventions import get_vaccination_coverage_of_age_group, get_contact_reductions, get_vaccination_date_of_age_group
from covid_visualisations import *

def aggregate_compartments(compartments):
    comp_sums = []
    for compartment in compartments:
        comp_sums.append(np.sum(compartment, axis=0))
    return np.array(comp_sums)

def aggregate_vaccinated_and_unvaccinated(compartments):
    split_index = int(len(compartments)/2)
    unvaccinated = compartments[0:split_index]
    vaccinated = compartments[split_index:len(compartments)]
    aggs = []
    for i in range(len(unvaccinated)):
        agg = unvaccinated[i] + vaccinated[i]
        aggs.append(agg)
    return np.array(aggs)

(beta, epsilon, r_s, r_a, theta_NIL, theta_1Y, theta_3Y, p, rho_NIL, rho_1Y, rho_3Y,
            pi_NIL, pi_1Y, pi_3Y, tau_NIL, tau_1Y, tau_3Y, omega_NIL, omega_1Y, omega_3Y,
            mu_NIL, mu_1Y, mu_3Y, alpha_S, alpha_A, ypsilon, e, f, li_and_ti_NIL, li_and_ti_1Y, li_and_ti_3Y,
            gamma_NIL, gamma_1Y, gamma_3Y, gamma_v_NIL, gamma_v_1Y, gamma_v_3Y, NPI_CONTACTS_REDUCTION) = load_parameters(vaccinate_5_to_11_yo=False)

contact_matrix = load_contact_matrix()
vaccination_plan = get_vaccination_date_of_age_group()

init = load_initial_conditions()
N = np.sum(init, axis=0)



####################################### FIGURE 1 ###############################################
def figure1():
    coverage = get_vaccination_coverage_of_age_group()
    R0s = [1.5, 2.5, 3.5, 4.5, 5.5]
    NPIs = [False, True]
    for R0 in R0s:
        for NPI in NPIs:
            Rt_plots(R0, alpha_A, alpha_S, contact_matrix, N, epsilon, p, r_a, r_s, NPI=NPI)
################################################################################################

def model(t, x, params, contact_matrix):
    (beta, epsilon, r_s, r_a, theta_NIL, theta_1Y, theta_3Y, p, rho_NIL, rho_1Y, rho_3Y,
            pi_NIL, pi_1Y, pi_3Y, tau_NIL, tau_1Y, tau_3Y, omega_NIL, omega_1Y, omega_3Y,
            mu_NIL, mu_1Y, mu_3Y, alpha_S, alpha_A, ypsilon, e, f, li_and_ti_NIL, li_and_ti_1Y, li_and_ti_3Y,
            gamma_NIL, gamma_1Y, gamma_3Y, gamma_v_NIL, gamma_v_1Y, gamma_v_3Y, NPI_CONTACTS_REDUCTION) = params


    S, E, IA, IS, H, ICU, M, R, S_v, E_v, IA_v, IS_v, H_v, ICU_v, M_v, R_v = x
    lambda_ = all_lambdas(beta, alpha_S, contact_matrix, IS, IS_v, N, alpha_A, IA, IA_v)
    theta, pi, tau, mu, rho, omega, gamma, gamma_v, li_and_ti = load_waning_dependent_parameters(waning=waning)

    dS_dt = -lambda_ * S - ypsilon * S + gamma * R + gamma_v * (R_v + S_v)
    dE_dt = lambda_ * S - epsilon * E
    dIA_dt = (1 - p) * epsilon * E - r_a * IA
    dIS_dt = p * epsilon * E - r_s * IS
    dH_dt = theta * r_s * IS  - rho * H
    dICU_dt = pi * rho * H - omega * ICU
    dM_dt = mu * omega * ICU + tau * rho * H
    dR_dt = (1 - omega) * r_s * IS + r_a * IA + (1 - pi - tau) * rho * H + (1 - mu) * omega * ICU - ypsilon * R - gamma * R

    dS_v_dt = -(1 - e) * lambda_ * S_v + ypsilon * S - gamma_v * S_v
    dE_v_dt = (1 - e) * lambda_ * S_v - epsilon * E_v
    dIA_v_dt = (1 - (1 - f) * p) * epsilon * E_v - r_a * IA_v
    dIS_v_dt = (1 - f) * p * epsilon * E_v - r_s * IS_v
    dH_v_dt = theta * r_s * IS_v - rho * H_v
    dICU_v_dt = pi * rho * H_v - omega * ICU_v
    dM_v_dt = mu * omega * ICU_v + tau * rho * H_v
    dR_v_dt = (1 - omega) * r_s * IS_v + r_a * IA_v + (1 - pi - tau) * rho * H_v + (1 - mu) * omega * ICU_v + ypsilon * R - gamma_v * R_v

    return np.array([dS_dt, dE_dt, dIA_dt, dIS_dt, dH_dt, dICU_dt, dM_dt, dR_dt,
                     dS_v_dt, dE_v_dt, dIA_v_dt, dIS_v_dt, dH_v_dt, dICU_v_dt, dM_v_dt, dR_v_dt])

def solve(x0, params, contact_matrix, vaccination_plan, NPI, waning, amount_of_classes = 16, amount_of_age_groups = 16, bounds=[0, 600], dt=0.01):
    (beta, epsilon, r_s, r_a, theta_NIL, theta_1Y, theta_3Y, p, rho_NIL, rho_1Y, rho_3Y,
            pi_NIL, pi_1Y, pi_3Y, tau_NIL, tau_1Y, tau_3Y, omega_NIL, omega_1Y, omega_3Y,
            mu_NIL, mu_1Y, mu_3Y, alpha_S, alpha_A, ypsilon, e, f, li_and_ti_NIL, li_and_ti_1Y, li_and_ti_3Y,
            gamma_NIL, gamma_1Y, gamma_3Y, gamma_v_NIL, gamma_v_1Y, gamma_v_3Y, NPI_CONTACTS_REDUCTION) = params

    ts = np.arange(bounds[0], bounds[1], dt)
    solved = np.zeros((amount_of_classes, amount_of_age_groups, ts.size))
    contact_reductions = get_contact_reductions(waning=waning)
    current_ypsilon = np.copy(ypsilon)



    for i in range(len(solved)):
        for j in range(len(solved[i])):
            solved[i][j][0] = x0[i][j]

    for i in range(1, ts.size):
        t = (i-1) * dt
        if i % 2000 == 0:
            print(f'Solving t={t}...')
        y = solved[:, :, i-1]
        current_contacts = update_contact_matrix(t, contact_matrix, contact_reductions)
        current_ypsilon = update_vaccination(t, current_ypsilon, vaccination_plan)        
        params = (beta, epsilon, r_s, r_a, theta_NIL, theta_1Y, theta_3Y, p, rho_NIL, rho_1Y, rho_3Y,
            pi_NIL, pi_1Y, pi_3Y, tau_NIL, tau_1Y, tau_3Y, omega_NIL, omega_1Y, omega_3Y,
            mu_NIL, mu_1Y, mu_3Y, alpha_S, alpha_A, current_ypsilon, e, f, li_and_ti_NIL, li_and_ti_1Y, li_and_ti_3Y,
            gamma_NIL, gamma_1Y, gamma_3Y, gamma_v_NIL, gamma_v_1Y, gamma_v_3Y, NPI_CONTACTS_REDUCTION)

        new_values = model(t, y, params, current_contacts)
        for j in range(len(solved)):
            for k in range(len(solved[j])):
                solved[j][k][i] = y[j][k] + new_values[j][k] * dt 

    return ts, solved

def update_contact_matrix(t, contact_matrix, contact_reductions, NPI=True):
    new_c = np.copy(contact_matrix)
    if NPI:
        new_c = new_c * (1 - NPI_CONTACTS_REDUCTION)

    for reduction in contact_reductions:
        t_start, t_end, strength = reduction
        if t >= t_start and t < t_end:
            new_c = new_c * (1 - strength)
            break
    return new_c

def update_vaccination(t, ypsilon, vaccination_plan):
    new_yps = np.copy(ypsilon)
    for i in range(len(vaccination_plan)):
        vaccination_start, vaccination_end = vaccination_plan[i]
        if vaccination_start != None and vaccination_end != None:
            if t < vaccination_start or t > vaccination_end:
                new_yps[i] = 0
    return new_yps



x0 = init
params = load_parameters() 
contact_matrix = load_contact_matrix()
NPI = False
wanings = ['NIL']


for waning in wanings:
    ts, solved = solve(x0, params, contact_matrix, vaccination_plan=vaccination_plan, NPI=NPI, waning=waning)

    (S, E, IA, IS, H, ICU, M, R) = aggregate_vaccinated_and_unvaccinated(aggregate_compartments(solved))
    plt.plot(ts, H)
plt.ylabel('Hospital occupancy')
plt.legend(['1Y waning', '3Y waning', 'no immunity waning'])
plt.show()

#aggregate_vaccinated_and_unvaccinated((S, E, IA, IS, H, ICU, M, R, S_v, E_v, IA_v, IS_v, H_v, ICU_v, M_v, R_v)

#(S, E, IA, IS, H, ICU, M, R, S_v, E_v, IA_v, IS_v, H_v, ICU_v, M_v, R_v) = aggregate_compartments(solved)


#plt.plot(ts, S, label='S')
#plt.show()
#plt.plot(ts, IA, label='IA')
#plt.plot(ts, IS, label='IS')
#plt.show()
#plt.plot(ts, M, label='M')
#plt.plot(ts, R, label='R')
#plt.show()
#plt.legend()

