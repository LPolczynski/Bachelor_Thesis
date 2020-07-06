from para import *



def ODE_V_G1(t, y, para):
    # parameters
    c_e = para[0]
    R = para[1]
    T = para[2]
    Lp = para[3]
    d = para[4]
    phi = para[5]
    pi_tc = para[6]
    nu = para[7]
    k_uptake = para[8]
    k_maintenance = para[9]
    E = para[10]

    phi_d = para[11]
    E_d = para[12]
    c_diff = para[13]
    w_diff = para[14]

    k_A = para[15]
    k_D = para[16]
    k_b = para[17]
    k_cost = para[18]

    # integration variables
    V_os = y[0]
    V_b = y[1]
    r = y[2]
    R_ref = y[3]
    c_i = y[4]
    pi_t = y[5]

    B_A = y[6]
    B_R = y[7]
    B = y[8]

    V_os_d = y[9]
    V_b_d = y[10]
    r_d = y[11]
    R_ref_d = y[12]
    c_i_d = y[13]
    pi_t_d = y[14]

    B_A_d = y[15]
    B_R_d = y[16]
    B_d = y[17]

    # dependent variables
    A = 4 * pi * r ** 2
    pi_i = c_i * R * T
    pi_e = c_e * R * T
    plastic_expansion = max(pi_t - pi_tc, 0)

    # derivatives
    dV_os = - Lp * A * (pi_t + (pi_e - pi_i))
    dV_b = (k_D / (k_b - k_D)) * dV_os
    dr = 1 / 3 * (3 / 4 / pi * k_b / (k_b - k_D)) ** (1 / 3) * dV_os / V_os ** (2 / 3)
    dR_ref = phi * R_ref * r / (2 * d) * plastic_expansion
    dc_i = 3 * k_uptake / r - k_maintenance - 3 * dr / r * (c_i + k_D * k_cost)
    dpi_t = E * 2 * d / (1 - nu) * (dr / r ** 2 - dR_ref / (R_ref * r)) - dr / r * pi_t

    dB_A = k_A * 8 * pi * R_ref * dR_ref
    dB_R = k_D * k_b / (k_b - k_D) * dV_os - dB_A
    dB = dB_A + dB_R

    dV_os_d = 0
    dV_b_d = 0
    dr_d = 0
    dR_ref_d = 0
    dc_i_d = 0
    dpi_t_d = 0

    dB_A_d = 0
    dB_R_d = 0
    dB_d = 0

    f = [dV_os, dV_b, dr, dR_ref, dc_i, dpi_t, dB_A, dB_R, dB, dV_os_d,
         dV_b_d, dr_d, dR_ref_d, dc_i_d, dpi_t_d, dB_A_d, dB_R_d, dB_d]
    return f


def ODE_V_G2(t, y, para):
    # parameters
    c_e = para[0]
    R = para[1]
    T = para[2]
    Lp = para[3]
    d = para[4]
    phi = para[5]
    pi_tc = para[6]
    nu = para[7]
    k_uptake = para[8]
    k_maintenance = para[9]
    E = para[10]

    phi_d = para[11]
    E_d = para[12]
    c_diff = para[13]
    w_diff = para[14]

    k_A = para[15]
    k_D = para[16]
    k_b = para[17]
    k_cost = para[18]

    # integration variables
    V_os = y[0]
    V_b = y[1]
    r = y[2]
    R_ref = y[3]
    c_i = y[4]
    pi_t = y[5]

    B_A = y[6]
    B_R = y[7]
    B = y[8]

    V_os_d = y[9]
    V_b_d = y[10]
    r_d = y[11]
    R_ref_d = y[12]
    c_i_d = y[13]
    pi_t_d = y[14]

    B_A_d = y[15]
    B_R_d = y[16]
    B_d = y[17]

    # dependent variables
    A = 4 * pi * r ** 2
    pi_i = c_i * R * T
    pi_e = c_e * R * T
    plastic_expansion = max(pi_t - pi_tc, 0)

    A_d = 4 * pi * r_d ** 2
    pi_i_d = c_i_d * R * T
    plastic_expansion_d = max(pi_t_d - pi_tc, 0)

    # derivatives
    if coupling:    #coupling variable is stored in the parameter file
        dV_ex = water_diffusion(w_diff, pi_t, pi_t_d)
        dV_exchange = dV_ex[0]
        dV_exchange_d = dV_ex[1]

        dc_diff = osmolyte_diffusion(c_diff, c_i, c_i_d, r, r_d)
        dc_i_diff = dc_diff[0]
        dc_i_diff_d = dc_diff[1]

        dV_os = - Lp * (pi_t + (pi_e - pi_i)) + dV_exchange
        dV_b = (k_D / (k_b - k_D)) * dV_os
        dr = 1 / 3 * (3 / 4 / pi * k_b / (k_b - k_D)) ** (1 / 3) * dV_os / V_os ** (2 / 3)
        dR_ref = phi * R_ref * r / (2 * d) * plastic_expansion
        dc_i = 3 * k_uptake / r - k_maintenance - 3 * dr / r * (c_i + k_D * k_cost) + dc_i_diff
        dpi_t = E * 2 * d / (1 - nu) * (dr / r ** 2 - dR_ref / (R_ref * r)) - dr / r * pi_t

        dB_A = 0
        dB_R = 0
        dB = 0

        dV_os_d = - Lp * (pi_t_d + (pi_e - pi_i_d)) + dV_exchange_d
        dV_b_d = (k_D / (k_b - k_D)) * dV_os_d
        dr_d = 1 / 3 * (3 / 4 / pi * k_b / (k_b - k_D)) ** (1 / 3) * dV_os_d / V_os_d ** (2 / 3)
        dR_ref_d = phi_d * R_ref_d * r_d / (2 * d) * plastic_expansion_d
        dc_i_d = 3 * k_uptake / r_d - k_maintenance - 3 * dr_d / r_d * (c_i_d + k_D * k_cost) + dc_i_diff_d
        dpi_t_d = E_d * 2 * d / (1 - nu) * (dr_d / r_d ** 2 - dR_ref_d / (R_ref_d * r_d)) - dr_d / r_d * pi_t_d

        dB_A_d = k_A * 8 * pi * R_ref_d * dR_ref_d
        dB_R_d = k_D * k_b / (k_b - k_D) * dV_os_d - dB_A_d
        dB_d = dB_A_d + dB_R_d
    else:
        dV_os = 0
        dV_b = 0
        dr = 0
        dR_ref = 0
        dc_i = 0
        dpi_t = 0

        dB_A = 0
        dB_R = 0
        dB = 0

        dV_os_d = - Lp * A_d * (pi_t_d + (pi_e - pi_i_d))
        dV_b_d = (k_D / (k_b - k_D)) * dV_os_d
        dr_d = 1 / 3 * (3 / 4 / pi * k_b / (k_b - k_D)) ** (1 / 3) * dV_os_d / V_os_d ** (2 / 3)
        dR_ref_d = phi_d * R_ref_d * r_d / (2 * d) * plastic_expansion_d
        dc_i_d = 3 * k_uptake / r_d - k_maintenance - 3 * dr_d / r_d * (c_i_d + k_D * k_cost)
        dpi_t_d = E_d * 2 * d / (1 - nu) * (dr_d / r_d ** 2 - dR_ref_d / (R_ref_d * r_d)) - dr_d / r_d * pi_t_d

        dB_A_d = k_A * 8 * pi * R_ref_d * dR_ref_d
        dB_R_d = k_D * k_b / (k_b - k_D) * dV_os_d - dB_A_d
        dB_d = dB_A_d + dB_R_d

    f = [dV_os, dV_b, dr, dR_ref, dc_i, dpi_t, dB_A, dB_R, dB, dV_os_d,
         dV_b_d, dr_d, dR_ref_d, dc_i_d, dpi_t_d, dB_A_d, dB_R_d, dB_d]
    return f