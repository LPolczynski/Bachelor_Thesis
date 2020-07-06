import numpy as np
from scipy import integrate
import copy

from utility_f import *
from ODE_growth import *
from ODE_cyc import *
from para import *



class cell:
    def __init__(self, species_cyclins, species_growth_m, sim_length, t):
        self.lifetime = 0
        self.time_of_birth = t
        self.sim_length = sim_length - t
        if species_cyclins == 'init_species_cyclins':
            self.species_cyclins = self.initialise_species(init_species_cyclins, self.sim_length + 1)
        if species_growth_m == 'init_species_growth':
            self.species_growth_m = self.initialise_species(init_species_growth_m, self.sim_length + 1)
        else:
            self.species_cyclins = self.initialise_species(species_cyclins, self.sim_length + 1)
            self.species_growth_m = self.initialise_species(species_growth_m, self.sim_length + 1)
        self.species_growth_d = self.initialise_species(init_species_growth_d,
                                                        self.sim_length + 1)  # daughter species for volume growth are always init_species_growth_d at cell initialization
        self.phase = 0
        self.para_cyclins = copy.deepcopy(g1_parameters)
        self.para_growth = copy.deepcopy(growth_parameters)
        self.k_mRNA = 1  # scaling of mRNA burst size with cell size
        self.gen = 0
        self.t_in_g1 = 0
        self.times_in_g1 = []
        self.t_in_s = 0
        self.t_in_G2 = 0
        self.times_in_G2 = []
        self.t_in_M = 0
        self.temp_sol = None  # storage space for solution of precalculated volume model
        self.time_conversion = 1  # conversion ratio between time of cyclin and growth model
        self.pre_dur = dur + 5
        self.i_start = 0  # index for faster search of right time index before each cyclin integration step
        self.threshold = threshold

    def simulate(self):
        daughters = []
        self.integration_growth()
        for step in range(self.sim_length):
            self.lifetime += 1
            self.transcribe(self.phase)
            self.integration_cyclins(1)
            daughter = self.update_phase(1)
            if daughter:
                daughters.append(daughter)
        return self, daughters

    def integration_cyclins(self, t_step):
        '''- retrieves the results of growth simulation that are needed for this timestep
           - integrates cyclin ODEs
           - stores cyclin results
           - stores approximated results of growth simulation at timepoints of cyclin integration'''

        self.para_cyclins[-1] = self.prepare_integration_cyclins(t_step)  # get radius,metabolic biomass, timepoints
        # and corresponding index from
        # growth simulation

        X0 = [val[self.lifetime - 1] for val in self.species_cyclins.values()]
        t_points = np.arange(0, t_step + 1, 1)
        X = integrate.odeint(ODEs, X0, t_points, (self.para_cyclins,))  # ,hmax=0.1
        for i, spec in enumerate(self.species_cyclins):
            self.species_cyclins[spec][self.lifetime] = X[1, i]

        res_growth = self.get_growth_results(t_step)
        spec_nr = len(self.species_growth_m)
        for i, spec in enumerate(self.species_growth_m):
            self.species_growth_m[spec][self.lifetime] = res_growth[i]
            self.species_growth_d[spec][self.lifetime] = res_growth[spec_nr + i]

    def integration_growth(self):
        '''- retrieve results of last growth simulation at new starting point
           - integrate long enough timespan to ensure the next cell phase is finished before end
           - save results in intermediate store self.temp_sol
           - adjust timepoints of growth simulation with conversion factor'''

        self.i_start = 0  # reset start index before each growth simulation
        X0 = self.prepare_integration_growth()  # [val[self.lifetime] for val in self.species_growth_m.values()] + [val[self.lifetime] for val in self.species_growth_d.values()]
        if self.phase == 0:
            if len(self.times_in_g1) > 2:
                dur_sim = (self.times_in_g1[-1] + self.times_in_g1[-2] + self.times_in_g1[-3]) * 2
            else:
                dur_sim = self.pre_dur
            self.temp_sol = integrate.solve_ivp(ODE_V_G1, [0, dur_sim], X0, method='BDF', args=[self.para_growth],
                                                dense_output=True, max_step=1.)
        if self.phase == 1:
            if len(self.times_in_G2) > 2:
                dur_sim = (self.times_in_G2[-1] + self.times_in_G2[-2] + self.times_in_G2[-3]) * 2
            else:
                dur_sim = self.pre_dur
            self.temp_sol = integrate.solve_ivp(ODE_V_G2, [0, dur_sim], X0, method='BDF', args=[self.para_growth],
                                                dense_output=True, max_step=1.)

        self.temp_sol.t = self.temp_sol.t * self.time_conversion  # adjust time scales
        # plot_growth_species(self.temp_sol)

    def prepare_integration_cyclins(self, t_step):
        t_start = self.get_t_start()
        t_end = t_start + t_step
        self.i_start = get_time_index(self.temp_sol.t, t_start, self.i_start)
        res_t = self.temp_sol.t
        res_r_m = self.temp_sol.y[2]  # total radius of mother
        res_r_d = self.temp_sol.y[11]  # total radius of daughter
        res_B_R = self.temp_sol.y[7]  # metabolic biomass of mother
        res_B_R_d = self.temp_sol.y[16]  # metabolic biomass of daughter
        return [res_t, res_r_m, res_r_d, t_start, res_B_R, res_B_R_d, self.i_start]

    def prepare_integration_growth(self):
        initial_values = []
        for spec in self.species_growth_m:
            initial_values.append(self.species_growth_m[spec][self.lifetime])
        if self.phase == 0:
            for spec in self.species_growth_d:
                initial_values.append(self.species_growth_d[spec][self.lifetime])
        if self.phase == 1:
            for spec in self.species_growth_d:
                if (spec == 'c_i' or spec == 'pi_t') and coupling == True:
                    initial_values.append(self.species_growth_m[spec][self.lifetime])
                else:
                    initial_values.append(self.species_growth_d[spec][self.lifetime])
        return initial_values

    def initialise_species(self, species, my_sim_length):
        my_species = {}
        for k in species:
            my_species[k] = np.zeros(my_sim_length)
            my_species[k][0] = copy.deepcopy(species[k])
        return my_species

    def update_phase(self, t_step):
        daughter = None
        if self.phase == 0:
            if self.species_cyclins['Cln'][self.lifetime] >= self.threshold['Cln']:
                print('im in G2 now!')
                print(self.lifetime)
                self.phase = 1
                self.times_in_g1.append(self.t_in_g1)
                self.para_cyclins = copy.deepcopy(g2_parameters)
                self.t_in_g1 = 0
                self.integration_growth()
            else:
                self.t_in_g1 += t_step
        elif self.phase == 1:
            if self.t_in_s >= 25:
                self.phase = 2
                self.t_in_s = 0
            else:
                self.t_in_s += t_step
        elif self.phase == 2:
            if self.species_cyclins['Clb'][self.lifetime] >= self.threshold['Clb']:
                self.phase = 3
                self.times_in_G2.append(self.t_in_G2)
                self.t_in_G2 = 0
            else:
                self.t_in_G2 += t_step
        elif self.phase == 3:
            if self.t_in_M >= 5:
                daughter = self.divide()
                print(self.lifetime)
            else:
                self.t_in_M += t_step
        return daughter

    def divide(self):
        print('division')
        daughter_species_cyclins = copy.deepcopy(init_species_cyclins)  # get structure of dictionary
        daughter_species_growth = copy.deepcopy(init_species_growth_m)  # to store results
        V_m = 4 / 3 * pi * self.species_growth_m['r'][self.lifetime] ** 3
        V_d = 4 / 3 * pi * self.species_growth_d['r'][self.lifetime] ** 3
        for k in daughter_species_cyclins:  # fill dictionary with results
            daughter_species_cyclins[k] = self.species_cyclins[k][self.lifetime] * V_d / (V_m + V_d)
            self.species_cyclins[k][self.lifetime] = self.species_cyclins[k][self.lifetime] - daughter_species_cyclins[
                k]
        for k in daughter_species_growth:
            daughter_species_growth[k] = self.species_growth_d[k][self.lifetime]
            self.species_growth_d[k][self.lifetime] = init_species_growth_d[k]
        self.para_cyclins = g1_parameters
        self.phase = 0
        self.integration_growth()
        return cell(daughter_species_cyclins, daughter_species_growth, self.time_of_birth + self.sim_length,
                    self.time_of_birth + self.lifetime)

    def transcribe(self, phase):
        samp = np.random.randint(1, 11, 1)
        if phase == 0:
            if samp[0] < 5:
                self.species_cyclins['mCLN'][self.lifetime - 1] += 1
        elif phase == 2:
            if samp[0] < 5:
                self.species_cyclins['mCLB'][self.lifetime - 1] += 1

    def get_t_start(self):
        if self.phase == 0:
            t_start = self.t_in_g1
        elif self.phase == 1:
            t_start = self.t_in_s
        elif self.phase == 2:
            t_start = 25 + self.t_in_G2
        elif self.phase == 3:
            t_start = 25 + self.times_in_G2[-1] + self.t_in_M
        return t_start

    def get_growth_results(self, t_step):
        t_start = self.get_t_start()
        t_end = t_start + t_step
        i = get_time_index(self.temp_sol.t, t_end)
        res_growth = []
        for j in range(len(self.temp_sol.y)):
            res_growth.append(get_growth_species(self.temp_sol.y[j], self.temp_sol.t, t_end, i))
        return res_growth
