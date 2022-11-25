import glob
import sipm.util.submit as submit

dirs = glob.glob('/scratch/gpfs/GALBIATI/data/sipm/reflector_studies/2022-11-22/**/*/')

scheduler = submit.Scheduler(dirs=dirs)
scheduler.submit()