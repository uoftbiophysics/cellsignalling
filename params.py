from settings import GLOB_C, GLOB_K_ON, GLOB_K_OFF, GLOB_K_P


class Params:
    def __init__(self, c=GLOB_C, k_on=GLOB_K_ON, k_off=GLOB_K_OFF, k_p=GLOB_K_P):
        self.c = c
        self.k_on = k_on
        self.k_off = k_off
        self.k_p = k_p
        self.x = c * k_on / k_off
        self.pss = self.x / (1 + self.x)
        self.r = c * k_on + k_off

    def unpack(self):
        return self.c, self.k_on, self.k_off, self.k_p, self.x, self.pss, self.r


DEFAULT_PARAMS = Params()
