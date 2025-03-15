from abc import abstractmethod
from skewsmithnormalform import pir
from skewsmithnormalform import skew_functions as sf 

class sidesED(pir.PIR):
    @abstractmethod
    def norm(self):
        pass
    def __lt__(self, x):
        return self.norm() < x.norm()
    def __gt__(self, x):
        return self.norm() > x.norm()
    @abstractmethod
    def get_left_q(self, x):
        pass
    @abstractmethod
    def get_right_q(self, x):
        pass
    def __leftfloordiv__(self, x):
        return self.get_left_q(x)
    def __rightfloordiv__(self, x):
        return self.get_right_q(x)
    @abstractmethod
    def get_left_r(self, x):
        pass
    @abstractmethod
    def get_right_r(self, x):
        pass
    def __leftmod__(self, x):
        return self.get_left_r(x)
    def __rightmod__(self, x):
        return self.get_right_r(x)
    def extended_grcd(a, b):
        x0 = type(b).getOne()
        x1 = type(b).getZero()
        y0 = type(b).getZero()
        y1 = type(b).getOne()
        while b != type(b).getZero():
            tempa = a
            tempb = b
            q = tempa.get_right_q(tempb)
            b = tempa.get_right_r(tempb)
            a = tempb
            tempx0 = x0
            x0 = x1
            x1 = tempx0 - q.__rightmul__(x0)
            tempy0 = y0
            y0 = y1
            y1 = tempy0 - q.__rightmul__(y0)
        return [a, x0, y0]
    def extended_glcd(a, b):
        x0 = type(b).getOne()
        x1 = type(b).getZero()
        y0 = type(b).getZero()
        y1 = type(b).getOne()
        while b != type(b).getZero():
            tempa = a
            tempb = b
            q = tempa.get_left_q(tempb)
            b = tempa.get_left_r(tempb)
            a = tempb
            tempx0 = x0
            x0 = x1
            x1 = tempx0 - q.__leftmul__(x0)
            tempy0 = y0
            y0 = y1
            y1 = tempy0 - q.__leftmul__(y0)
        return [a, x0, y0]