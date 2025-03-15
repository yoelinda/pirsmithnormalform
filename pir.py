from abc import ABC, abstractmethod
class InvalidInitialContent(Exception):
    pass

class PIR(ABC):
    @abstractmethod
    def __str__(self):
        pass
    @abstractmethod
    def __eq__(self, x):
        pass
    @abstractmethod
    def __ne__(self, x):
        pass
    @abstractmethod
    def __neg__(self):
        pass
    @abstractmethod
    def __add__(self, x):
        pass
    @abstractmethod
    def __sub__(self, x):
        pass
    @abstractmethod
    def __leftmul__(self, x):
        pass
    @abstractmethod
    def __rightmul__(self, x):
        pass
    @abstractmethod
    def __leftfloordiv__(self, x):
        pass
    @abstractmethod
    def __rightfloordiv__(self, x):
        pass
    @abstractmethod
    def __leftmod__(self, x):
        pass
    @abstractmethod
    def __rightmod__(self, x):
        pass
    @staticmethod
    @abstractmethod
    def getZero():
        pass
    @staticmethod
    @abstractmethod
    def getOne():
        pass
    @abstractmethod
    def isUnit(self):
        pass
    def isrightUnitMultipleOf(self, x):
        if not (self.get_right_r(x)) == self.getZero():
            return False
        if not (self.get_right_q(x)).isUnit():
            return False
        return True
    def isleftUnitMultipleOf(self, x):
        if not (self.get_left_r(x)) == self.getZero():
            return False
        if not (self.get_left_q(x)).isUnit():
            return False
        return True
    @abstractmethod
    def extended_grcd(a, b):
        pass
    @abstractmethod
    def extended_glcd(a, b):
        pass