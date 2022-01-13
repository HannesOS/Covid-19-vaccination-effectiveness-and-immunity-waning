import numpy as np
from scipy.linalg import block_diag
from scipy.optimize import curve_fit
np.set_printoptions(threshold=np.inf)


def new_infections_matrix(beta, alpha_A, alpha_S, contact_matrix, S_dash, S_v_dash, e, N):
    matrix = []
    blocks = []
    for i in range(len(contact_matrix)):
        row = np.array([])
        for j in range(len(contact_matrix[i])):
            block = F_i_j(i, j, beta, alpha_A, alpha_S, contact_matrix, S_dash, S_v_dash, e, N)
            if len(row) < 1:
                row = block
            else:
                row = np.concatenate((row, block), axis=1)
        if len(matrix) < 1:
            matrix = row
        else:
            matrix = np.concatenate((matrix, row), axis=0)
    return np.array(matrix)


def F_i_j(i, j, beta, alpha_A, alpha_S, contact_matrix, S_dash, S_v_dash, e, N):
    matrix = []
    matrix.append(np.array([0,
            beta * alpha_A * contact_matrix[i][j] * S_dash[i] / N[j],
            beta * alpha_S * contact_matrix[i][j] * S_dash[i] / N[j],
            0,
            beta * alpha_A * contact_matrix[i][j] * S_v_dash[i] / N[j],
            beta * alpha_S * contact_matrix[i][j] * S_v_dash[i] / N[j]]))
    matrix.append(np.array([0, 0, 0, 0, 0, 0]))
    matrix.append(np.array([0, 0, 0, 0, 0, 0]))
    matrix.append(np.array([0,
            beta * alpha_A * contact_matrix[i][j] * (1 - e) * S_dash[i] / N[j],
            beta * alpha_S * contact_matrix[i][j] * (1 - e) * S_dash[i] / N[j],
            0,
            beta * alpha_A * contact_matrix[i][j] * (1 - e) * S_v_dash[i] / N[j],
            beta * alpha_S * contact_matrix[i][j] * (1 - e) * S_v_dash[i] / N[j]]))
    matrix.append(np.array([0, 0, 0, 0, 0, 0]))
    matrix.append(np.array([0, 0, 0, 0, 0, 0]))
    return np.array(matrix)


def compartment_transition_matrix(epsilon, p, r_a, r_s, f, contact_matrix):
    blocks = []
    for i in range(len(contact_matrix)):
        block = V_i(i, epsilon, p, r_a, r_s, f)
        blocks.append(block)
    matrix = block_diag(*blocks)
    return matrix


def V_i(i, epsilon, p, r_a, r_s, f):
    matrix = []
    matrix.append(np.array([epsilon, 0, 0, 0, 0, 0]))
    matrix.append(np.array([-(1 - p) * epsilon, r_a, 0, 0, 0, 0]))
    matrix.append(np.array([-p * epsilon, 0, r_s, 0, 0, 0]))
    matrix.append(np.array([0, 0, 0, epsilon, 0, 0]))
    matrix.append(np.array([0, 0, 0, -(1 - (1 - f) * p) * epsilon, r_a, 0]))
    matrix.append(np.array([0, 0, 0, -(1 - f) * p * epsilon, 0, r_s]))
    return np.array(matrix)


def next_generation_matrix(beta, alpha_A, alpha_S, contact_matrix, S_dash, S_v_dash, e, N, epsilon, p, r_a, r_s, f):
    F = new_infections_matrix(beta, alpha_A, alpha_S, contact_matrix, S_dash, S_v_dash, e, N)
    V = compartment_transition_matrix(epsilon, p, r_a, r_s, f, contact_matrix)
    return np.matmul(F, np.linalg.inv(V))


def R_het(beta, alpha_A, alpha_S, contact_matrix, S_dash, S_v_dash, e, N, epsilon, p, r_a, r_s, f):
    eigenvalues = np.linalg.eigvals(next_generation_matrix(beta, alpha_A, alpha_S, contact_matrix, S_dash, S_v_dash, e, N, epsilon, p, r_a, r_s, f))
    dominant = np.max(eigenvalues)
    return dominant.real


def adjust_beta_to_fit_R0(R0, alpha_A, alpha_S, contact_matrix, e, N, epsilon, p, r_a, r_s, f):
    S_dash = N
    S_v_dash = np.zeros_like(N)

    max_iteration = 250000
    best_beta = 0
    best_error = 1e10
    betas = np.linspace(0,1,max_iteration)
    for beta in betas:
        Rt = R_het(beta, alpha_A, alpha_S, contact_matrix, S_dash, S_v_dash, e, N, epsilon, p, r_a, r_s, f)
        error = np.abs(R0 - Rt)
        if error < best_error:
            best_beta = beta
            best_error = error
        if error > best_error:
            break
    return best_beta

'''
beta = 0.1088 
alpha_A = 1.1
alpha_S = 0.5
contact_matrix = [[2, 1], [1, 2]] 
S_dash = [10, 20]
S_v_dash = [5, 10]
e = 0.6
N = [150, 200]
epsilon = 0.1
p = 0.5
r_a = 0.1
r_s = 0.2
f = 0.5
'''

#print(R_het(beta, alpha_A, alpha_S, contact_matrix, S_dash, S_v_dash, e, N, epsilon, p, r_a, r_s, f))