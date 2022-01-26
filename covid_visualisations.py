import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
from covid_age_dependence import age_group_names
from covid_parameters import get_beta_for_R0
from covid_reproduction import R_het
from covid_interventions import get_vaccination_coverage_of_age_group
from covid_ode_model import solve


def Rt_plots(R0, alpha_A, alpha_S, contact_matrix, N, epsilon, p, r_a, r_s, NPI=False, legend=True):
    if NPI:
        contact_matrix = contact_matrix * 0.53

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
    bounds=[0, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 7]
    norm = colors.BoundaryNorm(bounds, cmap.N)


    fig, ax = plt.subplots(1,1)
    fig.set_dpi(180)
    img = ax.imshow(result.T, origin='upper', extent=[0, 1, 0, 1], cmap=cmap, norm=norm)
    ax.set_xlabel('Vaccine effectiveness')
    ax.set_ylabel('Vaccine strategy (youngest vaccinated age group)')
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_xticklabels(['50%', '60%', '70%', '80%', '90%', '100%'])
    ax.set_yticks([0.5/15, 2.5/15, 4.5/15, 6.5/15, 8.5/15, 10.5/15, 12.5/15, 14.5/15])
    ax.set_yticklabels(['>75', '>65', '>55', '>45', '>35', '>25', '>15','>5',])
    if NPI:
        ax.set_title(f'$R_0$ = {R0} with NPIs')
    else:
        ax.set_title(f'$R_0$ = {R0}')

    if legend:
        legend_dict = { '$<0.5$': 'palegreen', '$0.5 - 0.75$' : 'limegreen', '$0.75 - 1.0$' : 'darkgreen' ,
             '$1.0 - 1.5$': 'orange', '$1.5 - 2.0$':'darkorange' , '$2.0 - 2.5$':'tomato', '$2.5 - 3.0$':'indianred', '$3.0-  3.5$':'brown', '$>3.5$' : 'darkred'}

        patchList = []
        for key in legend_dict:
            data_key = mpatches.Patch(color=legend_dict[key], label=key)
            patchList.append(data_key)
        ax.legend(handles=patchList, loc='best', bbox_to_anchor=(0.85, 0.5, 0.5, 0.5), frameon=False, title='$R_t$')

    fig.savefig(f'figures\\{R0}_{NPI}.png')
    plt.show()

def time_series(x0, params, contact_matrix, NPI, wanings, max_x=715):

    xticks = np.arange(10, max_x, 70)
    print(xticks)
    xlabels = ['feb 2021', 'april 2021', 'june 2021', 'aug 2021', 'oct 2021', 'dec 2021',
            'feb 2022', 'april 2022', 'june 2022', 'aug 2022', 'oct 2022']
    fig, ax = plt.subplots()
    fig.set_dpi(180)
    fig.set_size_inches(10,4)
    colors = ['mediumblue', 'darkred', 'forestgreen']
    for i, waning in enumerate(wanings):

        ts, solved = solve(x0, params, contact_matrix, vaccination_plan=vaccination_plan, NPI=NPI, waning=waning, bounds=[0, max_x])

        (S, E, IA, IS, H, ICU, M, R) = aggregate_vaccinated_and_unvaccinated(aggregate_compartments(solved))
        ax.plot(ts, ICU, linewidth=1.6, color=colors[i])
    ax.set_ylabel('Intensive care unit occupancy')


    ax.set_xticks(xticks)
    ax.set_xlim(0,max_x)
    ax.hlines(280, 0, max_x, linestyles='dashed', color='black')
    ax.set_xticklabels(xlabels, fontsize=6)
    ax.set_ylim(-100, 2420)

    plt.legend(['1 year immunity', '3 year immunity', 'no immunity loss'], frameon=False)
    plt.show()