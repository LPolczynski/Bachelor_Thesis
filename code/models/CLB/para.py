from math import pi


#parameters and initial values#
###############################
#constants
coupling = False
threshold = {'Cln':30,'Clb':6}
dur = 15

#parameter
#cyclins
prod_Cln_G1 = 1   #fl/(pg*min)
prod_Clb_G1 = 0   #fl/(pg*min)
deg_mCLN_G1 = .05
deg_Cln_G1 = 0.1
deg_mCLB_G1 = .5
deg_Clb_G1 = 0.5


prod_Cln_G2 = 0   #fl/(pg*min)
prod_Clb_G2 = 1   #fl/(pg*min)
deg_mCLN_G2 = .5
deg_Cln_G2 = 0.5
deg_mCLB_G2 = .05
deg_Clb_G2 = 0.1

g1_parameters = [prod_Cln_G1,prod_Clb_G1,deg_mCLN_G1,deg_Cln_G1,deg_mCLB_G1,deg_Clb_G1,'res_growth'] # is for storage of temp_sol
g2_parameters = [prod_Cln_G2,prod_Clb_G2,deg_mCLN_G2,deg_Cln_G2,deg_mCLB_G2,deg_Clb_G2,'res_growth']

#growth
ce = 240  #mM
R = 8.314  #J/(mol*K)
T = 303  #K
Lp = 1.19*10**(-6)  #mu m/(s*Pa)
d = 0.115  #mu m
phi = 10**(-4)   #1/(s*Pa)
pi_tc = 0.2*10**6  #Pa
nu = 0.5  #arb.unit
growth = 20
k_uptake = 0.35 *growth #mmol/((mu m)**2*s)
k_maintenance = 0.3 *growth  #mmol/((mu m)**3*s)
modulus_adjustment = (1-nu**2)**(-1)
E3d = 2.58*10**6  #Pa
E = modulus_adjustment * E3d

k_D = 0.5 # g/ml --> g/cm**3 --> pg/(mu m)**3 (https://www.merckmillipore.com/DE/de/product/Yeast-extract,MDA_CHEM-111926?bd=1
                                              #millipore Merck, Hefeextrakt, density)
k_b = 1.1 # g/ml --> g/cm**3 --> pg/(mu m)**3 (source:W.Baldwin,H.Kubitchek;1984;Journal of Bacteriology)
k_A = d*k_b  # pg/(mu m)**2
k_cost = 0  #mmol/pg

if coupling == False:
    phi_d = phi
    E_d = E
    c_diff = 0   # (mu m)**3/min
    w_diff = 0   # (mu m)**3/(min*Pa)
elif coupling == True:
    phi_d = phi * 100
    E_d =  E *1.28
    c_diff = 1.   # (mu m)**3/mins
    w_diff =  1.  # (mu m)**3/(min*Pa)

growth_parameters = [ce,R,T,Lp,d,phi,pi_tc,nu,k_uptake,k_maintenance,
                     E,phi_d,E_d,c_diff,w_diff,k_A,k_D,k_b,k_cost]


#initial values
#cyclins
init_species_cyclins = {'mCLN':0., 'Cln':0., 'mCLB':0., 'Clb':0.}   #number of molecules

#growth
#mother
V_os = 10 #(mu m)**3
V_b = k_D/(k_b-k_D)*V_os   #(mu m)**3
r = (3/4/pi*k_b/(k_b-k_D)*V_os)**(1/3)   #mu m
c_i = 319.17   #mM
pi_t = 0.2*10**(6)   #Pa
R_ref = r/(1 + (1 - nu) * (pi_t * r) / (E * 2 * d))   #mu m
B_A = k_A * 4*pi*R_ref**2   #pg
B_R = k_D*k_b/(k_b-k_D)*V_os - B_A   #pg
if B_R < 0:
    print('Invalid start volume, B_R = ', B_R, '!')
B = B_A + B_R   #pg

#daughter
V_os_d = .7   #(mu m)**3
V_b_d = k_D/(k_b-k_D)*V_os_d   #(mu m)**3
r_d = (3/4/pi*k_b/(k_b-k_D)*V_os_d)**(1/3)   #mu m
c_i_d = 319.17   #mM
pi_t_d = 0.2*10**(6)   #Pa
R_ref_d = r_d/(1 + (1 - nu) * (pi_t_d * r_d) / (E_d * 2 * d))   #mu m
B_A_d = k_A * 4*pi*R_ref_d**2   #pg
B_R_d = k_D*k_b/(k_b-k_D)*V_os_d - B_A_d   #pg
if B_R_d < 0:
    print('Invalid start volume, B_R_d = ', B_R_d, '!')
B_d = B_A_d + B_R_d  #pg

init_species_growth_m = {'V_os':V_os, 'V_b':V_b, 'r':r, 'R_ref':R_ref, 'c_i':c_i, 'pi_t':pi_t,
                         'B_A':B_A, 'B_R':B_R, 'B':B}
init_species_growth_d = {'V_os':V_os_d, 'V_b':V_b_d, 'r':r_d, 'R_ref':R_ref_d, 'c_i':c_i_d, 'pi_t':pi_t_d,
                        'B_A':B_A_d, 'B_R':B_R_d, 'B':B_d}
