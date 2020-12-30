import numpy as np, os, run_analytic_skewspec as local
from quicksub import *

def create_runfile(f,i=0,Om=0.3,Mhnum=50):
    
    add('import numpy as np, run_analytic_skewspec as local',f,ini=True)
    add("for j, Mh in enumerate(local.params_array(.5,1.2,"+str(Mhnum)+")):",f)
    add("    local.compute_skewspec("+str(i)+",j,Om="+str(Om)+",As=Mh*(0.3/"+str(Om)+")**2*2e-9)",f)
        

def jobfile(tag,**kwargs):

    #set run file
    f_run = 'tmp_job_run_'+tag+'.py'
    create_runfile(f_run,Om,**kwargs)

    # set job file
    f_sub = 'tmp_job_sub_'+tag+'.sh'
    set_sbatch_params(f_sub,tag,mem='32G',t='0-04:00',email=True)
    add('source ~/.bashrc.ext',f_sub)
    add('py4so',f_sub)
    add('python '+f_run,f_sub)
    
    # submit
    #os.system('sbatch '+f_sub)
    os.system('sh '+f_sub)
    os.system('rm -rf '+f_run+' '+f_sub)


Ommin, Ommax, Omnum = .2, .4, 1
Oms = local.params_array(Ommin,Ommax,Omnum)

for i, Om in Oms:
    print(Om)
    jobfile('Om_'+str(Om),i=i,Om=Om,Mhnum=1)

