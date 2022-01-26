import numpy as np

NPI_CONTACTS_REDUCTION = 0.47

#Contact reductions in the form [[t_{from}, t_{till}, contact reduction]]
contact_reductions_NIL = np.array([[0, 32, 0.8527], [32, 139, 0.6657], [139, 191, 0.6260], [191, 248, 0.6097]])
contact_reductions_1Y = np.array([[0, 43, 0.8230], [43, 158, 0.6070], [158, 187, 0.6898], [187, 248, 0.5955]])
contact_reductions_3Y = np.array([[0, 39, 0.8295], [39, 160, 0.6663], [160, 203, 0.6100], [203, 248, 0.6652]])

#The vaccination plan for each age-group (see covid_age_dependence.py) in the form 
#vaccination_plan = [[t_{from}_agegroup0, t_{till}_agegroup0], [t_{from}_agegroup1, t_{till}_agegroup1],..., [t_{from}_agegroup15, t_{till}_agegroup15]]] 
vaccination_plan = np.array([[None, None], [None, None], [221, 278], [221, 278], [204, 264], [204, 264], [190, 264],
                             [190, 264], [173, 248], [173, 248], [159, 234], [159, 234], [16, 217], [16, 217], [16, 217], [16, 217]]) 

#vaccination_plan = np.array([[None, None], [None, None], [221, 221+14], [221, 221+14], [204, 204+14], [204, 204+14], [190, 190+14],
#                             [190, 190+14], [173, 173+14], [173, 173+14], [159, 159+14], [159, 159+14], [16, 16+14], [16, 16+14], [16, 16+14], [16, 16+14]]) 


vaccination_plan_as_date = np.array([[None, None], [None, None], ['01.09.2021', '29.10.2021'], ['01.09.2021', '29.10.2021'], ['15.08.2021', '15.10.2021'],
                                    ['15.08.2021', '15.10.2021'], ['01.08.2021', '15.10.2021'], ['01.08.2021', '15.10.2021'], ['15.07.2021', '29.09.2021'],
                                    ['15.07.2021', '29.09.2021'], ['01.07.2021', '15.09.2021'], ['01.07.2021', '15.09.2021'], ['08.02.2021', '29.08.2021'],
                                    ['08.02.2021', '29.08.2021'], ['08.02.2021', '29.08.2021'], ['08.02.2021', '29.08.2021']])

#We assume that roughly that in the class of 10 to 14 year olds, all ages are represented with the same proportion
#hense vaccination coverage_{10-14} = (0 + 0 + 0.8 + 0.8)/4 = 0.4
vaccination_coverage = np.array([0, 0, 0.4, 0.8, 0.85, 0.85, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.99, 0.99, 0.99, 0.99])

def get_contact_reductions(waning='NIL'):
    if waning == 'NIL':
        return contact_reductions_NIL
    elif waning == '1Y':
        return contact_reductions_1Y
    elif waning == '3Y':
        return contact_reductions_3Y
    else:
        exit('Please choose between no immunity loss (NIL), 1 year immunity (1Y) and 3 year immunity (3Y)')


def get_vaccination_date_of_age_group(age_group=None, vaccination_plan_for_5_to_11_yo=None, as_date=False):
    if as_date:
        plan = vaccination_plan_as_date
    else:
        plan = vaccination_plan
    if vaccination_plan_for_5_to_11_yo != None:
        plan[1] = vaccination_plan_for_5_to_11_yo
    if age_group == None:
        return plan
    else:
        return plan[age_group]

def get_vaccination_coverage_of_age_group(age_group=None, vaccination_coverage_for_5_to_11_yo=None):
    cov = vaccination_coverage
    if vaccination_coverage_for_5_to_11_yo != None:
        cov[1] = vaccination_coverage_for_5_to_11_yo
        cov[2] = (vaccination_coverage_for_5_to_11_yo * 2 + 0.8 * 2) / 4
    if age_group != None:
        cov = cov[age_group]
    return cov

        
def vaccination_rate_of_age_group(age_group,  vaccination_plan_for_5_to_11_yo=None, vaccination_coverage_for_5_to_11_yo=None):
    P_i = get_vaccination_coverage_of_age_group(age_group, vaccination_coverage_for_5_to_11_yo=vaccination_coverage_for_5_to_11_yo)
    t_start, t_end = get_vaccination_date_of_age_group(age_group, vaccination_plan_for_5_to_11_yo=vaccination_plan_for_5_to_11_yo) 
    if t_start == None or t_end == None:
        return 0
    T_i = t_end - t_start + 1
    return -np.log((1 - P_i) / T_i)

def vaccination_rate_for_all_age_groups(vaccination_plan_for_5_to_11_yo=None, vaccination_coverage_for_5_to_11_yo=None):
    rates = []
    for age in range(len(vaccination_coverage)):
        rates.append(vaccination_rate_of_age_group(age, vaccination_plan_for_5_to_11_yo=vaccination_plan_for_5_to_11_yo,
                                                        vaccination_coverage_for_5_to_11_yo=vaccination_coverage_for_5_to_11_yo ))
    return np.array(rates)

def vaccine_efficacy(e, f):
    return 1 - (1 - e) * (1 - f)
