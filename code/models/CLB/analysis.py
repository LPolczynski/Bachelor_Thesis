from plots import *
from joblib import Parallel, delayed

def main():
    dire = ['./results/2020-07-05-growth-5/','./results/2020-07-05-growth-10/','./results/2020-07-05-growth-15/']
    plot_G1G2_md(dire)

    pop_dir = './results/2020-07-05-growth-15/'
    plot_birthsize_G1(pop_dir)

    pop_dir = './results/2020-07-05-growth-15/'
    plot_deltaVG1_birthsize(pop_dir)

    pop_dir = './results/2020-07-05-growth-15/'
    plot_Vdiv_VG1(pop_dir)

    pop_dir = './results/2020-07-05-growth-15/'
    plot_Vdiv_gen(pop_dir)

    pop_dir = './results/2020-07-05-growth-15/'
    plot_Vdiv_Vbirth(pop_dir)

    dire = ['../CLB/results/2020-07-05-growth-','../CLBnotLOC/results/2020-07-05-growth-','../noCLB/results/2020-07-05-growth-']
    growth = [5,10,15,20]
    labels=['Clb local','fixed G2', 'Clb not local']
    plot_size_distr(dire,growth,labels)

    pop_dir = './results/2020-07-05-growth-15/'
    plot_G1(pop_dir)
    plot_G2(pop_dir)

    plot_pop_volume_life(pop, 'Time (min)', 'Cell volume (fL)', 'lifecycle_V_' + filename_sim + '.pdf', 1, 10000, save_dir)
    Parallel(n_jobs=-1)(delayed(plot_species)(pop[cell_nr], save_dir, cell_nr) for cell_nr in range(min(len(pop), 10)))

if __name__=='__main__':
    main()