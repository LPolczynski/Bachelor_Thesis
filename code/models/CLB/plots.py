import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pylab as p
import pickle

from cellclass import *


def plot_pop_volume_life(pop, xlab, ylab, save_name, int_1, int_2, save_dir):
    font = {'size': 14}
    matplotlib.rc('font', **font)

    logit = True
    f1 = p.figure()
    tmp = []

    X = 'r'
    # zahl=2000
    # tmp2=zeros((zahl,zahl),dtype=float)

    zahl = max([len(cell.species_growth_m[X]) for cell in pop])  # "zahl" is the simulation time of the oldest cell
    nr_cells = len(pop)
    tmp2 = np.zeros((nr_cells, zahl), dtype=float)  # 2D array with one array of length(zahl) for every cell

    tmp_mn = []
    tmp_sd = []
    i = 1
    for k1, cell in enumerate(pop):
        if np.any(cell.species_growth_m[X]):
            # assemble results of species X for cell number k1
            for k2, m, d in zip(range(zahl - len(cell.species_growth_m[X]), zahl), cell.species_growth_m[X],
                                cell.species_growth_d[X]):
                tmp2[k1, k2] = 4 / 3 * pi * (m ** 3 + d ** 3)
            # plot in range between int_1 and int_2, check if log scale should be used
            if not logit and len(cell.species_growth_m[X]) >= int_1 and i <= int_2:
                p.plot(range(zahl - len(cell.species_growth_m[X]), zahl), cell.species_growth_m[X], alpha=0.1)
                i += 1
            elif len(cell.species_growth_m[X]) >= int_1 and i <= int_2:
                p.plot(range(zahl - len(cell.species_growth_m[X]), zahl), np.log10(cell.species_growth_m[X]),
                       alpha=0.1)  # log-scale
                # p.xlim(400,600)
                i += 1
    # for every timestep calculate mean and standard deviation of all cells living at that time
    for m in range(zahl):
        if not logit:
            tmp_list = [n for n in tmp2[:, m] if n != 0]
        else:
            tmp_list = [np.log10(n) for n in tmp2[:, m] if n != 0]  # log-scale
        tmp_mn.append(np.mean(tmp_list))
        tmp_sd.append(np.std(tmp_list))
    tmp_mn = np.array(tmp_mn)
    tmp_sd = np.array(tmp_sd)
    print(tmp_mn[-1], tmp_sd[-1])

    p.text(50, 1, '%s size %s var' % (round(tmp_mn[-1], 2), round(tmp_sd[-1], 2)))

    x_values = [i * 1 for i in range(zahl - len(tmp_mn), zahl)]

    p1, = p.plot(x_values, tmp_mn, color='#25597c', lw=3)
    p.fill_between(x_values, tmp_mn - tmp_sd, tmp_mn + tmp_sd, color='#82bcd3', alpha=1)

    p.plot(x_values, tmp_mn - tmp_sd, color='black', lw=0.3)  # plot black line along standart deviation lower end
    p.plot(x_values, tmp_mn + tmp_sd, color='black', lw=0.3)  # plot black line along standart deviation upper end

    r = p.Rectangle((0, 0), .2, .2, facecolor='#82bcd3', alpha=1, edgecolor='black',
                    lw=0.3)  # creates rectangle patch for legend use.

    p.legend([p1, r], ('Mean', 'Standard deviation'), loc='lower right', shadow=True)

    p.yticks([np.log10(q) for q in np.linspace(10, 250, 7)], np.linspace(10, 250, 7))

    p.xlabel(xlab)
    p.ylabel(ylab)
    # p.title('linear')
    f1.savefig(save_dir + save_name)




def plot_species(cell, save_dir, cell_nr, plot_conc=False):
    spec_nr = int(len(cell.species_cyclins) / 2) + len(cell.species_growth_m) + 2
    fig, ax = plt.subplots(spec_nr, 2, figsize=(32, 64))
    i = 0
    h = 0
    for spec in cell.species_cyclins:
        life_length = len(cell.species_cyclins[spec])
        t = np.linspace(0, life_length - 1, life_length)
        if plot_conc:
            ax[i][h].plot(t,
                          cell.species_cyclins[spec] / (cell.species_growth_m['V_os'] + cell.species_growth_m['V_b'] +
                                                        cell.species_growth_d['V_os'] + cell.species_growth_d['V_b']))
        else:
            ax[i][h].plot(t, cell.species_cyclins[spec])
        ax[i][h].set_xlabel('t')
        ax[i][h].legend([spec], shadow=True)
        if h == 0:
            h = 1
        elif h == 1:
            h = 0
            i += 1
    j = 0
    for spec in cell.species_growth_m:
        life_length = len(cell.species_growth_m[spec])
        t = np.linspace(0, life_length - 1, life_length)
        ax[i + j][0].plot(t, cell.species_growth_m[spec])
        ax[i + j][0].set_xlabel('t')
        ax[i + j][0].legend([spec], shadow=True)
        j += 1
    k = 0
    for spec in cell.species_growth_d:
        life_length = len(cell.species_growth_d[spec])
        t = np.linspace(0, life_length - 1, life_length)
        ax[i + k][1].plot(t, cell.species_growth_d[spec])
        ax[i + k][1].set_xlabel('t')
        ax[i + k][1].legend([spec], shadow=True)
        k += 1
    ax[i + k][0].plot(t, cell.species_growth_m['B_R'] / (4 / 3 * pi * cell.species_growth_m['r'] ** 3))
    ax[i + k][0].set_xlabel('t')
    ax[i + k][0].legend(['B_R/V'], shadow=True)
    ax[i + k][1].plot(t, cell.species_growth_d['B_R'] / (4 / 3 * pi * cell.species_growth_d['r'] ** 3))
    ax[i + k][1].set_xlabel('t')
    ax[i + k][1].legend(['B_R/V'], shadow=True)
    ax[i + k + 1][0].plot(t, (4 / 3 * pi * cell.species_growth_m['r'] ** 3))
    ax[i + k + 1][0].set_xlabel('t')
    ax[i + k + 1][0].legend(['V'], shadow=True)
    ax[i + k + 1][1].plot(t, (4 / 3 * pi * cell.species_growth_d['r'] ** 3))
    ax[i + k + 1][1].set_xlabel('t')
    ax[i + k + 1][1].legend(['V'], shadow=True)
    fig.savefig(save_dir + 'cell_%d' % cell_nr)
    plt.close(fig)
    print('figure save')


def plot_G1(pop_dir):
    pop = pickle.load(open(pop_dir + 'pop' + '.p', 'rb'))
    max_gen_nr = max([len(cell.times_in_g1) for cell in pop])
    mean_all_gen = []
    std_all_gen = []
    for gen in range(max_gen_nr):
        gen_array = [cell.times_in_g1[gen] for cell in pop if len(cell.times_in_g1) > gen]
        mean_all_gen.append(np.mean(gen_array))
        std_all_gen.append(np.std(gen_array))

    X = np.arange(max_gen_nr)
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.75,0.75])
    ax.bar(X + 0.00, mean_all_gen, color='b', width=0.25, label='Time in G1 mean')
    ax.bar(X + 0.25, std_all_gen, color='g', width=0.25, label='Time in G1 std')
    ax.set_xlabel('genealogical age')
    ax.set_ylabel('time (min)')
    ax.legend()
    fig.savefig(pop_dir + 'TG1')
    #plt.show()


def plot_G2(pop_dir):
    pop = pickle.load(open(pop_dir + 'pop' + '.p', 'rb'))
    max_gen_nr = max([len(cell.times_in_G2) for cell in pop])
    mean_all_gen = []
    std_all_gen = []
    for gen in range(max_gen_nr):
        gen_array = [cell.times_in_G2[gen] for cell in pop if len(cell.times_in_G2) > gen]
        mean_all_gen.append(np.mean(gen_array))
        std_all_gen.append(np.std(gen_array))

    X = np.arange(max_gen_nr)
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.75,0.75])
    ax.bar(X + 0.00, mean_all_gen, color='b', width=0.25, label='Time in G2 mean')
    ax.bar(X + 0.25, std_all_gen, color='g', width=0.25, label='Time in G2 std')
    ax.set_xlabel('genealogical age')
    ax.set_ylabel('time (min)')
    ax.legend()
    fig.savefig(pop_dir + 'TG2')
    #plt.show()


def plot_G1G2_md(dire):
    colors = ['b', 'g', 'r', 'c', 'm', 'y']
    fig = plt.figure()
    i = 0
    for p_dir in dire:
        pop = pickle.load(open(p_dir + 'pop' + '.p', 'rb'))
        TG1d = np.mean([cell.times_in_g1[0] for cell in pop if cell.times_in_g1])
        TG2d = np.mean([cell.times_in_G2[0] for cell in pop if cell.times_in_G2])
        t1 = []
        t2 = []
        for cell in pop:
            if len(cell.times_in_g1) > 1:
                t1.append(np.mean(cell.times_in_g1[1:]))
            if len(cell.times_in_G2) > 1:
                t2.append(np.mean(cell.times_in_G2[1:]))
        TG1 = np.mean(t1)
        TG2 = np.mean(t2)
        plt.plot(TG2, TG1, '%s' % colors[i] + 's', label='mother')
        plt.plot(TG2d, TG1d, '%s' % colors[i] + '^', label='daughter')
        i += 1
    plt.xlabel('G2')
    plt.ylabel('G1')
    plt.legend()
    fig.savefig('./G1G2_md')
    #plt.show()

def plot_birthsize_G1(pop_dir):
    pop = pickle.load(open(pop_dir + 'pop' + '.p', 'rb'))
    TG1 = []
    Vbirth = []
    for cell in pop:
        if cell.times_in_g1:
            TG1.append(cell.times_in_g1[0])
            Vbirth.append(cell.species_growth_m['r'][0] ** 3 * 4 * pi / 3)
    fig = plt.figure()
    plt.plot(Vbirth, TG1, 'bo')
    plt.xlabel('V at birth(fL)')
    plt.ylabel('Duration of first G1 phase(min)')
    fig.savefig(pop_dir + 'birthsize_G1')
    #plt.show()


def plot_deltaVG1_birthsize(pop_dir):
    pop = pickle.load(open(pop_dir + 'pop' + '.p', 'rb'))
    VG1 = []
    Vbirth = []
    for cell in pop:
        if cell.times_in_g1:
            Vbirth.append(cell.species_growth_m['r'][0] ** 3 * 4 * pi / 3)
            TG1 = cell.times_in_g1[0]
            VG1.append(
                cell.species_growth_m['r'][TG1] ** 3 * 4 * pi / 3 - cell.species_growth_m['r'][0] ** 3 * 4 * pi / 3)
    fig = plt.figure()
    plt.plot(Vbirth, VG1, 'bo')
    plt.xlabel('V at birth(fL)')
    plt.ylabel('Delta V in first G1')
    fig.savefig(pop_dir + 'birthsize_VG1')
    #plt.show()


def plot_Vdiv_VG1(pop_dir):
    pop = pickle.load(open(pop_dir + 'pop' + '.p', 'rb'))
    VG1 = []
    Vdiv = []

    TS = 25
    TM = 5
    for cell in pop:
        if cell.times_in_G2:
            pre_div_nr = len(cell.times_in_G2)
            t = sum(cell.times_in_g1[:pre_div_nr]) + sum(cell.times_in_G2[:pre_div_nr]) + (TS + TM) * pre_div_nr
            if t <= cell.lifetime:
                div_nr = pre_div_nr
            else:
                div_nr = pre_div_nr - 1
            if len(cell.times_in_g1) > div_nr:
                gen_nr = div_nr
            else:
                gen_nr = div_nr - 1
            for i in range(gen_nr):
                Tdiv = cell.times_in_g1[i] + TS + cell.times_in_G2[i] + TM
                Vdiv.append(cell.species_growth_m['r'][Tdiv] ** 3 * 4 * pi / 3)
                TG1 = cell.times_in_g1[i + 1]
                VG1.append(cell.species_growth_m['r'][Tdiv + TG1] ** 3 * 4 * pi / 3 - cell.species_growth_m['r'][
                    Tdiv] ** 3 * 4 * pi / 3)
    fig = plt.figure()
    plt.plot(Vdiv, VG1, 'bo')
    plt.xlabel('V at division(fL)')
    plt.ylabel('Delta V in next G1')
    fig.savefig(pop_dir + 'Vdivision_VG1')
    #plt.show()


def plot_Vdiv_gen(pop_dir):
    pop = pickle.load(open(pop_dir + 'pop' + '.p', 'rb'))
    Gen = []
    Vdiv = []
    TS = 25
    TM = 5
    for cell in pop:
        if cell.times_in_G2:
            pre_div_nr = len(cell.times_in_G2)
            t = sum(cell.times_in_g1[:pre_div_nr]) + sum(cell.times_in_G2[:pre_div_nr]) + (TS + TM) * pre_div_nr
            if t <= cell.lifetime:
                gen_nr = pre_div_nr
            else:
                gen_nr = pre_div_nr - 1
            for i in range(gen_nr):
                Tdiv = cell.times_in_g1[i] + TS + cell.times_in_G2[i] + TM
                Vdiv.append(cell.species_growth_m['r'][Tdiv] ** 3 * 4 * pi / 3)
                Gen.append(i)
    fig = plt.figure()
    plt.plot(Gen, Vdiv, 'bo')
    plt.xlabel('genealogical age')
    plt.ylabel('V at division(fL)')
    fig.savefig(pop_dir + 'Vdivision_gen')
    #plt.show()


def plot_Vdiv_Vbirth(pop_dir):
    pop = pickle.load(open(pop_dir + 'pop' + '.p', 'rb'))
    Vbirth = []
    Vdiv = []
    TS = 25
    TM = 5
    for cell in pop:
        if cell.times_in_G2:
            pre_div_nr = len(cell.times_in_G2)
            t = sum(cell.times_in_g1[:pre_div_nr]) + sum(cell.times_in_G2[:pre_div_nr]) + (TS + TM) * pre_div_nr
            if t <= cell.lifetime:
                gen_nr = pre_div_nr
            else:
                gen_nr = pre_div_nr - 1
            for i in range(gen_nr):
                Tdiv = cell.times_in_g1[i] + TS + cell.times_in_G2[i] + TM
                Vdiv.append(cell.species_growth_m['r'][Tdiv] ** 3 * 4 * pi / 3)
                Vbirth.append(cell.species_growth_d['r'][Tdiv] ** 3 * 4 * pi / 3)
    fig = plt.figure()
    plt.plot(Vdiv, Vbirth, 'bo')
    plt.xlabel('V at division(fL)')
    plt.ylabel('V of new daughter(fL)')
    fig.savefig(pop_dir + 'Vdivision_gen')
    #plt.show()


def plot_size_distr(dire, growth_rates, labels):
    color = ['bo', 'go', 'ro', 'co', 'mo', 'yo']
    fig, ax = plt.subplots(1, 2)
    fig.tight_layout(pad=3.0)
    for i in range(len(dire)):
        dire[i] = [dire[i] + '%s' % str(rate) + '/' for rate in growth_rates]
        mn = []
        std = []
        mn_log = []
        std_log = []
        for p_dir in dire[i]:
            pop = pickle.load(open(p_dir + 'pop' + '.p', 'rb'))
            cell_nr = len(pop)
            cell_array = np.zeros(cell_nr)
            for k, cell in enumerate(pop):
                cell_array[k] = 4 * pi / 3 * (cell.species_growth_m['r'][-1] ** 3 + cell.species_growth_d['r'][-1] ** 3)
            mean_vol = np.mean(cell_array)
            std_vol = np.std(cell_array)

            cell_array_log = np.log10(cell_array)
            mean_vol_log = np.mean(cell_array_log)
            std_vol_log = np.std(cell_array_log)

            mn.append(mean_vol)
            std.append(std_vol)
            mn_log.append(mean_vol_log)
            std_log.append(std_vol_log)

        ax[0].plot(growth_rates, mn, color[i], label=labels[i])
        ax[1].plot(growth_rates, std_log, color[i], label=labels[i])
    ax[0].set_xlabel('growth factor')
    ax[1].set_xlabel('growth factor')
    ax[0].set_ylabel('mean cell volume (fL)')
    ax[1].set_ylabel('standard variation of log(volume)')
    ax[0].legend(loc='upper left')
    ax[1].legend(loc='upper right')
    fig.savefig('./size_distr')
    #plt.show()

