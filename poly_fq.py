from smithnormalform import ed, pir
import galois
GF2 = galois.GF(2)
class GFX(ed.ED):
    def __init__(self, content):
        if isinstance(content, galois.Poly):
            self.a = content 
        else:
            raise pir.InvalidInitialContent
    def __str__(self):
        return str(self.a)
    def __hash__(self):
        return hash(("com.corbinmcneill.smithnormalform.z", self.a))
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
    def __mul__(self, x):
        return GFX(self.a * x.a)
    def get_q(self, x):
        return GFX(self.a // x.a)
    def get_r(self, x):
        return GFX(self.a % x.a)
    def norm(self):
        return self.a(GF2(0)) 
    @staticmethod
    def getZero():
        return GFX(   galois.Poly.Zero(GF2)   )
    @staticmethod
    def getOne():
        return GFX(  galois.Poly.One(GF2) )
    def isUnit(self):
        return (self.a == galois.Poly.One(GF2)) 