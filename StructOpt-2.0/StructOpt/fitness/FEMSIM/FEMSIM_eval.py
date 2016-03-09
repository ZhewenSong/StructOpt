import subprocess
import shlex
import sys
import json
import time
import numpy as np
import ase.io
import os
import logging
import math
import shutil
from mpi4py import MPI

class FEMSIM_eval(object):
    def __init__(self):
        self.args = self.read_inputs()

        #self.step_number = 0

        self.vk = np.multiply(self.args['thickness_scaling_factor'], self.vk)  # Multiply the experimental data by the thickness scaling factor


    def read_inputs(self):
        args = json.load(open('femsim_inp.json'))
        data = open(args['vk_data_filename']).readlines()
        data.pop(0)  # Comment line
        data = [line.strip().split()[:2] for line in data]
        data = [[float(line[0]), float(line[1])] for line in data]
        k, vk = zip(*data)
        # Set k and vk data for chi2 comparison
        self.k = np.array(k)
        self.vk = np.array(vk)
        return args


    def update_parameters(self, **kwargs):
        #self.step_number += 1

        #self.args['aberation_coef'] += 1

        for key, value in kwargs.items():
            self.args[key] = value

    def evaluate_fitness(self, Optimizer, individ):
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.Get_rank()

        if rank==0:
            femsimfiles = '{filename}-rank0/FEMSIMFiles'.format(filename=Optimizer.filename)
            try:
                os.mkdir(femsimfiles)
            except OSError:
                pass
            
            ntimes=int(math.ceil(1.*len(individ)/comm.Get_size()))
            nadd=int(ntimes*comm.Get_size()-len(individ))
            maplist=[[] for n in range(ntimes)]
            strt=0
            for i in range(len(maplist)):
                maplist[i]=[indi for indi in individ[strt:comm.Get_size()+strt]]
                strt+=comm.Get_size()
            for i in range(nadd):
                maplist[len(maplist)-1].append(None)
        else:
            ntimes=None
        ntimes = comm.bcast(ntimes,root=0)
        outs=[]
        for i in range(ntimes):
            if rank==0:
                one=maplist[i]
            else:
                one=None
            ind = comm.scatter(one,root=0)

            if ind == None:
                rank = MPI.COMM_WORLD.Get_rank()
                stro = 'Evaluated none individual on {0}\n'.format(rank)
                out = (None, stro)
            else:
                try:
                    os.mkdir('{filename}-rank0/FEMSIMFiles/rank-{r}'.format(filename=Optimizer.filename,r=rank))
                except OSError:
                    pass
                out = self.evaluate_indiv(Optimizer, ind, rank)

            outt = comm.gather(out,root=0)
            if rank == 0:
                outs.extend(outt)
        return outs

    def evaluate_indiv(self, Optimizer, individ, rank):

        logger = logging.getLogger(Optimizer.loggername)

        logger.info('Received individual HI = {0} for FEMSIM evaluation'.format(
            individ.history_index))
        
        cwd = os.getcwd()
        
        os.chdir('{filename}-rank0/FEMSIMFiles/rank-{r}'.format(filename=Optimizer.filename,r=rank))
        
        if not os.path.isfile(self.args['vk_data_filename']):
            shutil.copy(os.path.join(cwd,self.args['vk_data_filename']),self.args['vk_data_filename'])
        paramfilename = self.args['parameter_filename']
        paramfilename = paramfilename.split('.')
        paramfilename[-2] = '{head}_{i}'.format(head=paramfilename[-2], i=individ.history_index) 
        paramfilename = '.'.join(paramfilename)
        
        self.write_paramfile(paramfilename, Optimizer, individ, rank)

        base = 'indiv{i}'.format(i=individ.history_index)
        self.run_femsim(base=base, paramfilename=paramfilename)
        vk = self.get_vk_data(base)

        chisq = self.chi2(vk)
        logger.info('M:finish chi2 evaluation, chi2 = {0} @ rank = {1}'.format(chisq,rank))
        stro = 'Evaluated individual {0} @ rank = {1}\n'.format(individ.history_index,rank)
        
        os.chdir(cwd)
        return chisq, stro


    def write_paramfile(self, paramfilename, Optimizer, individ, rank):
        # Write structure file to disk so that the fortran femsim can read it in
        ase.io.write('structure_{i}.xyz'.format(i=individ.history_index), individ[0])
        # Replace comment line
        lines = open('structure_{i}.xyz'.format(i=individ.history_index)).readlines()
        lines[1] = "{} {} {}\n".format(self.args['xsize'], self.args['ysize'], self.args['zsize'])
        open('structure_{i}.xyz'.format(i=individ.history_index), 'w').write(''.join(lines))

        with open(paramfilename, 'w') as f:
            f.write('# Parameter file for generation {gen}, individual {i} @ rank {r}\n'.format(gen=Optimizer.generation, i=individ.history_index, r=rank))
            f.write('{}\n'.format('structure_{i}.xyz'.format(i=individ.history_index)))
            f.write('{}\n'.format(self.args['vk_data_filename']))
            f.write('{}\n'.format(self.args['Q']))
            f.write('{} {} {}\n'.format(self.args['nphi'], self.args['npsi'], self.args['ntheta']))
            f.write('{}\n'.format(self.args['thickness_scaling_factor']))



    def run_femsim(self, base, paramfilename):
        self.run_subproc('{femsim_command} {base} {paramfilename}'.format(femsim_command=os.getenv('FEMSIM_COMMAND'),base=base, paramfilename=paramfilename))


    def run_subproc(self, args):
        """ args should be the string that you would normally run from bash """
        #print("Running (via python): {0}".format(args))
        sargs = shlex.split(args)
        p = subprocess.Popen(sargs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = []
        for nextline in iter(p.stdout.readline, ""):
            sys.stdout.write(nextline)
            output.append(nextline)
            sys.stdout.flush()
        poutput = p.stdout.read()
        perr = p.stderr.read()
        preturncode = p.wait()
        if(preturncode != 0):
            print("{0} exit status: {1}".format(args, preturncode))
            print("{0} failed: {1}".format(args, perr))
        return ''.join(output)

    def get_vk_data(self, base):
        # Sleep until we can get the file
        # There may be in issue if the file is only partially written to when the open command gets run...
        # Let's hope that doesn't happen. If it does, I will convert this to a `if f in os.listdir()` command, followed by another short sleep.
        while True:
            try:
                data = open('vk_initial_{base}.txt'.format(base=base)).readlines()
                break # final? or initial?
            except IOError:
                time.sleep(5.0)
        #data.pop(0)  # I don't see the comment line in the vk output data
        data = [line.strip().split()[:2] for line in data]
        data = [[float(line[0]), float(line[1])] for line in data]
        vk = np.array([vk for k, vk in data])
        return vk


    def chi2(self, vk):
        return np.sum(((self.vk - vk) / self.vk)**2) / len(self.k)

