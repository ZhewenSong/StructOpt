class FEM_eval(object):
    def __init__(self, ind):
        self.Optimizer, self.individ = ind
        #self.fitness = 0
        self.args = self.read_inputs()

        #self.step_number = 0

        self.vk = np.multiply(self.args['thickness_scaling_factor'], self.vk)  # Multiply the experimental data by the thickness scaling factor

        #self.structure = None  # This gets initialized for the first time when the optimizer calls update_parameters and passes in a new structure
        #self.structure = ase.io.read('Zr50Cu35Al15_start.xyz')  # TODO Delete

    def read_inputs(self):
        args = json.load(open('fem_inputs.json'))
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

        #for key, value in kwargs.items():
        #    self.args[key] = value
        return

    def evaluate_fitness(self):
        self.write_paramfile(individual)
        vk = self.run_femsim()
        chi2 = self.chi2(vk)
        stro = ''
        return chi2, stro

    def write_paramfile(self, individ):
        # Write structure file to disk so that the fortran femsim can read it in
        ase.io.xyz.write_xyz(open('structure.xyz', 'w'), self.individ[0], comment="{} {} {}".format(self.args['xsize'], self.args['ysize'], self.args['zsize']))

        with open(self.args['parameter_filename'], 'w') as f:
            f.write('# Paramfile for step {}\n'.format(self.step_number))
            f.write('{}\n'.format(self.args['model_filename']))
            f.write('{}\n'.format(self.args['vk_data_filename']))
            f.write('{}\n'.format(self.args['Q']))
            f.write('{} {} {}\n'.format(self.args['nphi'], self.args['npsi'], self.args['ntheta']))
            f.write('{}\n'.format(self.args['thickness_scaling_factor']))

    def run_femsim(self):
        #subproc.run('femsim step{i} {parameterfile}'.format(i=optimizer.step_number, parameterfile=args['parameter_filename']))
        #data = open('vk_final.txt').readlines()
        data = open(self.args['vk_data_filename']).readlines()  # TODO Delete the following lines and correctly run femsim instead
        data.pop(0)  # Comment line
        data = [line.strip().split()[:2] for line in data]
        data = [[float(line[0]), float(line[1])] for line in data]
        vk = np.array([vk for k, vk in data])
        return vk

    def chi2(self, vk):
        return np.sum(((self.vk - vk) / self.vk)**2) / len(self.k)
