from math import*
import numpy as np

class DualNum():
    def __init__(self,val,dot):
        self.val = val
        self.dot = dot

        self.type = type(val)

    def checkType(self,other):
        if other.__class__!= DualNum and other.__class__ is not list:
            other = DualNum(other, 0)
        elif other.__class__ is list:
            dot =[]
            if(type(other[0]) is DualNum):
                for i in range(len(other)):
                    dot.append(other[i].dot)
                    other[i] = other[i].val

            else:
                dot =[0]*len(other)
            return DualNum(np.array(other),np.array(dot))
        return other

    def __str__(self):
        return "(" + str(self.val) + "," + str(self.dot) + ")"

    def __add__(self, other):
        other = self.checkType(other)
        val = self.val + other.val
        dot = self.dot + other.dot
        return DualNum(val, dot)

    def __sub__(self, other):
        other = self.checkType(other)
        val =self.val-other.val
        dot = self.dot-other.dot
        return DualNum(val, dot)

    def __mul__(self, other):
        other = self.checkType(other)

        val = np.dot(self.val,other.val)
        dot = np.dot(self.val,other.dot) + np.dot(self.dot,other.val)
        return DualNum(val, dot)

    def __truediv__(self, other):
        other = self.checkType(other)
        val = self.val / other.val
        dot = (np.dot(self.dot,other.val) )
        dot = (dot- np.dot(self.val,other.dot))
        dot = dot/np.power(other.val,2)
        return DualNum(val, dot)

    def __neg__(self):
        val = -self.val
        dot = -self.dot
        return DualNum(val, dot)


    def pow(self, power, modulo=None):

        val = np.power(self.val, power)
        if power<1 and self.val ==0:
            #print("here")
            dot = 0
        else:
            dot = power * np.power(self.val, power - 1) * self.dot

        return DualNum(val, dot)

    def __radd__(self, other):
        other = self.checkType(other)
        val = self.val + other.val
        dot = self.dot + other.dot
        return DualNum(val, dot)

    def __rsub__(self, other):
        other = self.checkType(other)
        val = other.val-self.val
        dot = other.dot-self.dot
        return DualNum(val, dot)

    def __round__(self, n=None):
        return round(self.val,n)

    def __rmul__(self, other):
        other = self.checkType(other)
        val = np.dot(self.val, other.val)
        dot = np.dot(self.val, other.dot) + np.dot(self.dot, other.val)
        return DualNum(val, dot)
    def __abs__(self):
        return np.abs(self.val)

    def __rdiv__(self, other):
        other = self.checkType(other)
        out = other/self
        return out
    def __rtruediv__(self, other):
        other = self.checkType(other)
        out = other / self
        return out

    def __lt__(self, other):
        compVal = other
        if(type(other).__name__=="DualNum"):
            compVal = other.val
        if(self.val<compVal):
            return True
        else:
            return False

    def __gt__(self, other):
        compVal = other

        if (type(other).__name__ == "DualNum"):
            compVal = other.val
        if (self.val > compVal):
            return True
        else:
            return False

    def __ge__(self,other):
        compVal = other

        if (type(other).__name__ == "DualNum"):
            compVal = other.val
        if (self.val >= compVal):
            return True
        else:
            return False

    def __pow__(self, power, modulo=None):
        val = np.power(self.val, power)

        if power<1 and self.val ==0:
            #print("here")
            dot = 0
        else:
            dot = power * np.power(self.val, power - 1) * self.dot

        return DualNum(val, dot)

    def __iter__(self):
        if type(self.val) is list or type(self.val)==type(np.array(0)):
            for i in range(len(self.val)):
                yield DualNum(self.val[i],self.dot[i])
        else:
            yield self
    def sin(self):
        val = sin(self.val)
        dot = cos(self.val)*self.dot
        return DualNum(val, dot)

    def cos(self):
        val = cos(self.val)
        dot = -1*sin(self.val)*self.dot
        return DualNum(val, dot)

    def tan(self):
        val = sin(self.val)/cos(self.val)
        dot = 1/pow(cos(self.val),2) *self.dot
        return DualNum(val, dot)


    def log(self):
        base = e
        val = np.log(self.val)
        dot = 1/(self.val * np.log(base))*self.dot
        return  DualNum(val, dot)

    def rec(self):
        val = 1/self.val
        dot = -1/pow(self.val,2)*self.dot
        return DualNum(val, dot)

    def arcsin(self):
        val = asin(self.val)
        dot = 1/sqrt(1-pow(self.val,2))*self.dot
        return DualNum(val, dot)

    def arccos(self):
        val = acos(self.val)
        dot = -1/sqrt(1-pow(self.val,2))*self.dot
        return DualNum(val, dot)

    def arctan(self):
        val = atan(self.val)
        dot = 1/(1+pow(self.val,2))*self.dot
        return DualNum(val, dot)

    def exp(self):
        val = np.exp(self.val)
        dot = np.exp(self.val)*self.dot
        return DualNum(val, dot)

