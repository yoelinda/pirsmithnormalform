from smithnormalform import ed, pir
import galois
from smithnormalform import skew_functions as sf 
field = galois.GF(4) 
class GFX(ed.sidesED):
    def __init__(self, content): 
        if isinstance(content, galois.Poly):
            self.a = content
        else:
            raise pir.InvalidInitialContent
    def __str__(self):
        return str(self.a)
    def __eq__(self, x): 
        return self.a == x.a
    def __ne__(self, x): 
        return self.a != x.a
    def __neg__(self):
        return GFX(-self.a)
    def __add__(self, x):
        return GFX(self.a + x.a)
    def __sub__(self, x):
        return GFX(self.a - x.a)
    def __leftmul__(self, x): 
        return GFX( sf.skewmultiplication(x.a, self.a)  )
    def __mul__(self, x): 
        return GFX(sf.skewmultiplication( self.a , x.a))
    def __rightmul__(self, x):
        return GFX(sf.skewmultiplication( self.a , x.a))
    def get_left_q(self, x):
        Q,R = sf.skew_Leuclideandivision(self.a , x.a)
        return GFX(Q)
    def get_right_q(self, x):
        Q,R = sf.skew_Reuclideandivision(self.a , x.a)
        return GFX(Q)
    def get_left_r(self, x):
        Q,R = sf.skew_Leuclideandivision(self.a , x.a)
        return GFX(R)
    def get_right_r(self, x):
        Q,R = sf.skew_Reuclideandivision(self.a , x.a)
        return GFX(R)
    def norm(self):
        return self.a(field(0)) 
    @staticmethod
    def getZero():
        return GFX(   galois.Poly.Zero(field)   )
    @staticmethod
    def getOne():
        return GFX(  galois.Poly.One(field) )
    def isUnit(self):
        number = sf.p**sf.m
        return ((int(self.a) < number) and int(self.a) >0) 
    def istotaldivisorof(self, x):
        return sf.istotaldivisorof(self.a, x.a)