from utility_f import *




def ODEs(X, t, p):
    # parameters
    prod_Cln = p[0]
    prod_Clb = p[1]
    deg_mCLN = p[2]
    deg_Cln = p[3]
    deg_mCLB = p[4]
    deg_Clb = p[5]

    # get radius/surface
    list_t = p[-1][0]
    list_r_m = p[-1][1]
    list_r_d = p[-1][2]
    t_start = p[-1][3]
    list_B_R_m = p[-1][4]
    list_B_R_d = p[-1][5]
    i_start = p[-1][6]

    t_now = t_start + t
    i_now = get_time_index(list_t, t_now, i_start)
    r_m = get_growth_species(list_r_m, list_t, t_now, i_now)
    r_d = get_growth_species(list_r_d, list_t, t_now, i_now)
    V_m = r_m ** 3 * 4 / 3 * pi
    V_d = r_d ** 3 * 4 / 3 * pi
    B_R_m = get_growth_species(list_B_R_m, list_t, t_now, i_now)
    B_R_d = get_growth_species(list_B_R_d, list_t, t_now, i_now)

    # get cyclins
    mCLN = X[0]
    Cln = X[1]
    mCLB = X[2]
    Clb = X[3]

    dmCLN = -deg_mCLN * mCLN
    dCln = prod_Cln * mCLN * (B_R_m + B_R_d) / (V_m + V_d) - deg_Cln * Cln
    dmCLB = -deg_mCLB * mCLB
    dClb = prod_Clb * mCLB * B_R_d / (V_d + V_m) - deg_Clb * Clb
    return [dmCLN, dCln, dmCLB, dClb]