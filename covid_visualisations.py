import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
from covid_age_dependence import age_group_names
from covid_parameters import get_beta_for_R0
from covid_reproduction import R_het
from covid_interventions import get_vaccination_coverage_of_age_group


def Rt_plots(R0, alpha_A, alpha_S, contact_matrix, N, epsilon, p, r_a, r_s, NPI=False):
    if NPI:
        contact_matrix = contact_matrix * 0.6

    coverage = get_vaccination_coverage_of_age_group()
    beta = get_beta_for_R0(R0)


    ages = np.arange(1, 16, 1)
    es = np.arange(0.5, 1, 0.03571428571)

    f = 0.5
    result = np.empty((len(es), len(ages)))
    for i, e in enumerate(es):
        coverage_ = np.copy(coverage)
        for j, age in enumerate(ages):
            for k in range(len(coverage)):
                if k < age:
                    coverage_[k] = 0 
            S_v_dash = N * coverage_
            S_dash = N - S_v_dash
            Rt = R_het(beta, alpha_A, alpha_S, contact_matrix, S_dash, S_v_dash, e, N, epsilon, p, r_a, r_s, f)
            result[i][j] = Rt
            print(f'e = {e}. Vaccination start: {age_group_names[age]}.  R_t = ', Rt)

    cmap = colors.ListedColormap(['palegreen', 'limegreen', 'darkgreen', 'orange', 'darkorange', 'tomato', 'indianred', 'brown', 'darkred'])
    bounds=[0, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 5]
    norm = colors.BoundaryNorm(bounds, cmap.N)


    fig, ax = plt.subplots(1,1)
    fig.set_dpi(180)
    img = ax.imshow(result.T, origin='upper', extent=[0, 1, 0, 1], cmap=cmap, norm=norm)
    ax.set_xlabel('Vaccine effectiveness')
    ax.set_ylabel('Vaccine strategy (youngest vaccinated age group)')
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_xticklabels(['50%', '60%', '70%', '80%', '90%', '100%'])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_yticklabels(['>75', '>50', '>45', '>30', '>20','>5',])
    if NPI:
        ax.set_title(f'$R_0$ = {R0} with NPIs')
    else:
        ax.set_title(f'$R_0$ = {R0}')
    plt.show()