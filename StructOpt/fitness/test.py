class test(object):
    def __init__(self, opt):
        self.opt = opt
        for k, v in self.opt.items():
            setattr(self, k, v)
