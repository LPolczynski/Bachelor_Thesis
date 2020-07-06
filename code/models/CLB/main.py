import os
from time import time
from datetime import datetime

from cellclass import *
from plots import *
from para import *
from simulation import *

def main():
    date = datetime.today().strftime('%Y-%m-%d')
    hour = str(datetime.now().time().hour)
    minute = str(datetime.now().time().minute)
    second = str(datetime.now().time().second)
    timepoint = date + '-' + hour + '-' + minute + '-' + second
    starttime = time()
    pop = simulation(2, dur)
    stoptime = time()
    print(len(pop))
    print(stoptime - starttime, ' s')

    save_dir = './results/%s/' % timepoint
    filename_sim = 'plot'
    os.makedirs(save_dir)
    simfile = 'para.py'
    save_sim_file(simfile, save_dir)

    pickle.dump(pop, open(save_dir + 'pop' + '.p', 'wb'))

if __name__=='__main__':
    main()