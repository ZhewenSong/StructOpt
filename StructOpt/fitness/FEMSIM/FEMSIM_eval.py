import subprocess
import shlex
import sys
import json
import time
import numpy as np
import ase.io


def run_subproc(args):
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


class FEMSIM_eval(object):
    def __init__(self):
        self.args = self.read_inputs()

        self.step_number = 0

        self.vk = np.multiply(self.args['thickness_scaling_factor'], self.vk)  # Multiply the experimental data by the thickness scaling factor


    def read_inputs(self):
        args = json.load(open('inputs.json'))
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
        self.step_number += 1

        #self.args['aberation_coef'] += 1

        for key, value in kwargs.items():
            self.args[key] = value


    def evaluate_fitness(self, individuals):
        chi2s = [None for _ in range(len(individuals))]
        for i, individual in enumerate(individuals):
            paramfilename = self.args['parameter_filename']
            paramfilename = paramfile.split('.')
            paramfilename[-2] = '{head}_{i}'.format(head=paramfilename[-2], i=i)
            paramfilename = '.'.join(paramfilename)

            self.write_paramfile(paramfilename, i, individual)

            base = 'step{i}'.format(i=i)
            self.run_femsim(base=base, paramfilename=paramfilename)
            vk = self.get_vk_data(base=base)

            chi2[i] = self.chi2(vk, range(individuals))
        return chi2s


    def write_paramfile(self, paramfilename, i, individual):
        # Write structure file to disk so that the fortran femsim can read it in
        ase.io.write('structure_{i}.xyz'.format(i=i), individual)
        # Replace comment line
        lines = open('structure_{i}.xyz'.format(i=i)).readlines()
        lines[1] = "{} {} {}".format(self.args['xsize'], self.args['ysize'], self.args['zsize'])
        open('structure_{i}.xyz'.format(i=i), 'w').write(''.join(lines))

        with open(paramfilename, 'w') as f:
            f.write('# Parameter file for step {step}, individual {i}\n'.format(step=self.step_number, i=i))
            f.write('{}\n'.format('structure_{i}.xyz'.format(i=i)))
            f.write('{}\n'.format(self.args['vk_data_filename']))
            f.write('{}\n'.format(self.args['Q']))
            f.write('{} {} {}\n'.format(self.args['nphi'], self.args['npsi'], self.args['ntheta']))
            f.write('{}\n'.format(self.args['thickness_scaling_factor']))


    def run_femsim(self, base, paramfilename):
        run_subproc('femsim-hrmc/femsim {base} {paramfilename}'.format(base=base, paramfilename=paramfilename))


    def get_vk_data(base):
        # Sleep until we can get the file
        # There may be in issue if the file is only partially written to when the open command gets run...
        # Let's hope that doesn't happen. If it does, I will convert this to a `if f in os.listdir()` command, followed by another short sleep.
        while True:
            try:
                data = open('vk_final_{base}.txt'.format(base=base)).readlines()
                break
            except IOError:
                time.sleep(5.0)
        data.pop(0)  # Comment line
        data = [line.strip().split()[:2] for line in data]
        data = [[float(line[0]), float(line[1])] for line in data]
        vk = np.array([vk for k, vk in data])
        return vk


    def chi2(self, vk):
        return np.sum(((self.vk - vk) / self.vk)**2) / len(self.k)

