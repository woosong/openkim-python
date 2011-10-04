
class descriptor:
    def __init__(self, name, units='fixed', model=False):
        if model:
            self.type = 'MODEL'
        else:
            self.type = 'TEST'
        self.name = name
        self.units = units
        self.particle_types = {}
        self.conventions = {}
        self.model_input = {}
        self.model_output = {}
        self.model_params = {}

    def build_species(self):
        str = 'SUPPORTED_ATOM/PARTICLES_TYPES:\n'
        for spec in self.particle_types:
            str += spec + '\tspec\t' + self.particle_types[spec].__str__() + '\n'
        return str

    def build_conventions(self):
        str = 'CONVENTIONS:\n'
        for conv in self.conventions:
            str += conv + '\tdummy\n'
        return str

    def build_model_input(self):
        str = 'MODEL_INPUT:\n'
        for name in self.model_input:
            data = self.model_input[name]
            str += name + '\t' + '\t'.join(data) + '\n'
        return str

    def build_model_output(self):
        str = 'MODEL_OUTPUT:\n'
        for name in self.model_output:
            data = self.model_output[name]
            str += name + '\t' + '\t'.join(data) + '\n'
        return str

    def build_model_params(self):
        if len(self.model_params) == 0:
            return ''
        str = 'MODEL_PARAMETERS:\n'
        for name in self.model_params:
            data = self.model_params[name]
            str += name + '\t' + '\t'.join(data) + '\n'
        return str

    def build_descriptor(self):
        str = self.type + '_NAME := ' + self.name + '\n'
        str += 'SystemOfUnitsFix := ' + self.units + '\n'
        str += self.build_species()
        str += self.build_conventions()
        str += self.build_model_input()
        str += self.build_model_output()
        str += self.build_model_params()
        return str


