'''
 ____________________________________________________________________________________________________________________
|    __________________________________________________________________________________________________________     |
|    |                                   @author: Sadegh Shokatian                                             |    |
|    |                                Supervisor: Stephen G. Urquhart                                          |    |
|    |                                    department of Chemistry                                              |    |
|    |                                   University of Saskatchewan                                            |    |
|    |_________________________________________________________________________________________________________|    |
|___________________________________________________________________________________________________________________|
'''
import math # To do all math opertions
import re # To split a line in a text file
import os # To make folders
import shutil # To copy files with different options
import sys # To exit the program when its necessary
a = 7.4 # a unit cell
b = 4.93# b unit cell
c = 16.82/2.0# c unit cell
n = 2 # number of layers of crystal along x and y (right now its 2 layer along z)
alkane_angle = 48.8 # Angle of alkane chain with XY plane
path = os.path.abspath(os.path.join(('Orthorhombic.py')))
path2 = os.path.split(os.path.dirname(path))[1]
Nn = str(path2) # Title of jobs that will be submitted with deMon2k
ZLayer = 3
###

def precision(x,p): # To round all coordinate numbers to any significant figures that we define 
    """
    returns a string representation of x formatted with a precision of p

    """

    x = float(x)

    if x == 0.:
        return "0." + "0"*(p-1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)

    if n < math.pow(10, p - 1):
        e = e -1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= math.pow(10,p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return "".join(out)

###
class Atoms: # Returns all the atoms with their name and coordinate in the molecular structure
    def __init__(self,_name,_x,_y,_z):
        self.x = _x
        self.y = _y
        self.z = _z
        self.q = 1.0
        self.name = _name
    def __repr__(self):
        return "Name:" + str(self.name) + "X: " +str(precision(self.x,7))+ "Y: "+str(precision(self.y,7))+ "Z: " + str(precision(self.z,7)) +"\n"
    def __str__(self):
        return "Name:" + str(self.name) + "X: " +str(precision(self.x,7))+ "Y: "+str(precision(self.y,7))+ "Z: " + str(precision(self.z,7)) +"\n"
    def xangle(self,_name,_x,_y,_z):
        return math.acos((self.x)/(math.sqrt(((self.x)**2.0)+((self.y)**2)+((self.z)**2.0))))
    def yangle(self,_name,_x,_y,_z):
        return math.acos((self.y)/(math.sqrt(((self.x)**2.0)+((self.y)**2)+((self.z)**2.0))))
    def zangle(self,_name,_x,_y,_z):
        return math.acos((self.z)/(math.sqrt(((self.x)**2.0)+((self.y)**2)+((self.z)**2.0))))
    def dist(self,j):
        mySum = ((self.x - j.x)*(self.x - j.x)) + ((self.y - j.y)*(self.y - j.y)) + ((self.z - j.z)*(self.z - j.z))
        return math.sqrt(mySum)
###
def GetFirstLastOneToLastC(molecule): # Returns the coordinates of first carbon and last carbon in the molecule
    result = []
    
    carbons = []
    for i in molecule:
        if i.name == 'C':
            carbons.append([i,0])

    hydrogens = []
    for i in molecule:
        if i.name == 'H':
            hydrogens.append(i)
        
    for i in carbons:
        for j in hydrogens:
            if i[0].dist(j) < 1.20:
                i[1]= i[1]+1

    firstC = Atoms('C',0,0,0)
    lastC = Atoms('C',0,0,0)
    secondLastC = Atoms('C',0,0,0)
    
    temp1 = 0
    temp2 = 1
    for i in carbons:
        if i[1] == 3:
            if temp1 ==0:
                firstC = i[0]
                temp1 = temp1 + 1
            elif temp1 == 1:
                lastC = i[0]                                
        
             
            
    for i in carbons:
        if lastC.dist(i[0]) < 1.70 and lastC.dist(i[0]) > 1.1:
            secondLastC = i[0]
    
    result.append(firstC)
    result.append(lastC)
    result.append(secondLastC)
    
    return result
###    
def TransToOrigin(molecule,firstC): # transfr the whole molecule to the origin (0,0,0)
    trans = [-1.0 * float(firstC.x),-1.0 * float(firstC.y),-1.0 * float(firstC.z)]
    for i in molecule:
        i.x = i.x + trans[0]
        i.y = i.y + trans[1]
        i.z = i.z + trans[2]
    return molecule
###
def RotToXZPlane(molecule, vector): # rotates the molecule to put it in the XZ plane
    a = vector[0]
    b = vector[1]
    c = vector[2]
    l = math.sqrt(a**2 + b**2 + c**2)
    v = math.sqrt(a**2 + c**2)
    if l == 0:
        return molecule
    elif l != 0:
        if b ==0:
            return molecule
        else:
            for i in molecule:
                x0 = i.x
                y0 = i.y
                z0 = i.z
                q0 = i.q
                i.x = x0
                i.y = (z0*(b/l)) - (y0*(v/l)) 
                i.z = - 1.0 * (y0*(b/l)) - (z0*(v/l))
                i.q = q0
            return molecule
###
def RotToZAxis(molecule,vector): # rotates the molecule to put it on the Z axis
    a = vector[0]
    b = vector[1]
    c = vector[2]
    l = math.sqrt(a**2 + c**2)
    if l == 0:
        return molecule
    elif l != 0:
        if a == 0 and b == 0:
            return molecule
        else:
            for i in molecule:
                x0 = i.x
                y0 = i.y
                z0 = i.z
                q0 = i.q
                i.x = -1.0 * (x0 * (c/l)) - (z0 * (a/l))
                i.y = y0
                i.z = (x0 * (a/l)) - (z0 * (c/l))
                i.q = q0
            return molecule
###
def TransXtoZero(molecule,vector):# rotates the vector between first and second carbon so it will be in the zero angle with X axis
    a = vector[0]
    b = vector[1]
    c = vector[2]
    l = math.sqrt(a**2 + b**2)
    if l == 0:
        return molecule
    elif l != 0:
        if b == 0:
            return molecule
        else:
            for i in molecule:
                x0 = i.x
                y0 = i.y
                z0 = i.z
                q0 = i.q
                i.x = -1.0 * (x0 * (a/l)) - (y0 * (b/l))
                i.y = (x0 * (b/l)) - (y0 * (a/l))  
                i.z = z0
                i.q = q0
            return molecule
###

def RotToXAngle(molecule,teta):# rotates the vector between first and second carbon so it will be in the defined angle with X axis
    cost = math.cos(teta)
    sint = math.sin(teta)
    for i in molecule:
        x0 = i.x
        y0 = i.y
        z0 = i.z
        q0 = i.q
        i.x = (x0 * cost) + (y0 * sint)
        i.y = (y0 * cost) - (x0 * sint)
        i.z = z0
        i.q = q0
    return molecule

###
def FirstLayer(molecule):# Makes the first layer of the crystal
    result = []
    Z = ZLayer
    if Z == 1:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
    elif Z ==2:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
    elif Z == 3:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
    else:
        print "please choose a value between 1 and 3 for Z"
        sys.exit()
    return result
#
def SecondLayer(molecule):# Makes the second layer of the crystal
    result = []
    Z = ZLayer
    if Z == 1:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
    elif Z ==2:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            
    elif Z == 3:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
    else:
        print "please choose a value between 1 and 3 for Z"
        sys.exit()
    return result
#
def ThirdLayer(molecule):# Makes the third layer of the crystal
    result = []
    Z = ZLayer
    if Z == 1:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
    elif Z ==2:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            
    elif Z == 3:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            '''result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])'''
    elif Z == 0:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c * 2.0)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z- float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z- float(c * 2.0)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c * 2.0)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c * 2.0)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c * 2.0)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c * 2.0)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c * 2.0)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c * 2.0)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c * 2.0)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c * 2.0)),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b * 2.0)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b * 2.0)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
    else:
        print "please choose a value between 1 and 3 for Z"
        sys.exit()    
    return result
#
def FourthLayer(molecule):# Makes the forth layer of the crystal
    result = []
    Z = ZLayer
    if Z == 1:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            #
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
    elif Z ==2:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            
    elif Z == 3:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            
    else:
        print "please choose a value between 1 and 3 for Z"
        sys.exit()
    return result
###
def FifthLayer(molecule):# Makes the fifth layer of the crystal
    result = []
    Z = ZLayer
    if Z == 1:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            #
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            #
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
    elif Z ==2:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
    elif Z == 3:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
    else:
        print "please choose a value between 1 and 3 for Z"
        sys.exit()
    return result
###
def SixthLayer(molecule):# Makes the sixth layer of the crystal
    result = []
    Z = ZLayer
    if Z == 1:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            #
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            #
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            #
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(1.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(1.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(1.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(1.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(1.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(1.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(1.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(1.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
    elif Z ==2:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(1.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(1.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(1.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(1.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(1.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(1.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(1.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(1.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(1.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(1.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(1.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(1.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(1.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(1.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(1.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(1.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
    elif Z == 3:
        for m in molecule:
            xn = round((-1.0 * m.y),7)
            yn = round(m.x ,7)
            zn = round(m.z,7)
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z+ float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y),7)),(precision((m.z- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn+ float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn- float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn- float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z + float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x+float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            result.append([m.name,(precision((m.x-float(a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z - float(c)),7))])
            #
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn + float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn+float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float(b/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((3.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((5.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn+float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float(a/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            result.append([m.name,(precision((xn-float((3.0 * a)/2.0)),7)),(precision((yn-float((5.0 * b)/2.0)),7)),(precision((zn - float(c)),7))])
            #
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(1.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(1.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(1.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(1.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(1.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(1.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x+float(1.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(1.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(1.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(1.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(1.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(1.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(1.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(1.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x+float(1.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(1.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z +  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y-float(1.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(1.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x-float(3.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x-float(1.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x+float(1.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(2.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y+float(1.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(1.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(2.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x+float(3.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x+float(2.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x+float(1.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x-float(1.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
            result.append([m.name,(precision((m.x-float(2.0 * a)),7)),(precision((m.y-float(3.0 * b)),7)),(precision((m.z -  float(c)),7))])
    else:
        print "please choose a value between 1 and 3 for Z"
        sys.exit()
    return result
###                    

class ChDir:# Submits Q jobs on grex
    def __init__(self, newpath):
        self.savedpath = os.getcwd()
        os.chdir(newpath)
        os.system('qsub d2k-batch.run')
    def __finalize__(self):
        if str(i) == n_carbon:
            sys.exitfunc()
    def __del__(self):
        os.chdir(self.savedpath)

###
with open('Sources//Geometry.xyz','r+') as file: # opens the initial geometry file as an input    
    molecule = []
    for line in file:        
        parts = re.split("[\s]+", line)[:-1];
        atom = Atoms(parts[0],float(parts[1]),float(parts[2]),float(parts[3]))
        molecule.append(atom)

###
'''Following lines will find the first and the last carbon of molecule,
then a vector between these two carbon will be defined,
based on this vector, the whole molecule will be transfered to origin (0,0,0),
then, the whole molecule will be rotated to XZ plane, again, based on the new defined vector
 between the first and the last carbon,
 next step will be rotating the whole molecule to put it along Z axis,
 then, a new vector between the first and the second carbon will be defined,
 using this new vector, molecule will be rotated so the angle between this vector with X axis will be zero,
 final step is rotating the whole molecule, based on the vector defined in last step, to a defined angle of
 orthorhomic structure
 '''

findResult = GetFirstLastOneToLastC(molecule)
firstC = findResult[0]
lastC = findResult[1]
secondLastC = findResult[2]
vector = [(lastC.x - firstC.x),(lastC.y - firstC.y),(lastC.z - firstC.z)]

###

molecule = TransToOrigin(molecule,firstC)

###
findResult = GetFirstLastOneToLastC(molecule)
firstC = findResult[0]
lastC = findResult[1]
secondLastC = findResult[2]
vector = [(lastC.x - firstC.x),(lastC.y - firstC.y),(lastC.z - firstC.z)]

molecule = RotToXZPlane(molecule,vector)




###

findResult = GetFirstLastOneToLastC(molecule)
firstC = findResult[0]
lastC = findResult[1]
secondLastC = findResult[2]
vector = [(lastC.x - firstC.x),(lastC.y - firstC.y),(lastC.z - firstC.z)]
molecule = RotToZAxis(molecule,vector)

###

findResult = GetFirstLastOneToLastC(molecule)
lastC = findResult[1]
secondLastC = findResult[2]
vector = [(lastC.x - secondLastC.x),(lastC.y - secondLastC.y),(lastC.z - secondLastC.z)]
molecule = TransXtoZero(molecule,vector)

###
findResult = GetFirstLastOneToLastC(molecule)
lastC = findResult[1]
secondLastC = findResult[2]
vector = [(lastC.x - secondLastC.x),(lastC.y - secondLastC.y),(lastC.z - secondLastC.z)]
teta = math.radians(alkane_angle) 
molecule = RotToXAngle(molecule,teta)

###

'''In this part, depend on numbers of layers defined at the beginning of the program, layers of crystal will be made'''
 
if n == 1:
    crystal = FirstLayer(molecule)
elif n == 2:
    crystal = SecondLayer(molecule)
elif n == 3:
    crystal = ThirdLayer(molecule)
elif n == 4:
    crystal =  FourthLayer(molecule)
elif n == 5:
    crystal =  FifthLayer(molecule)
elif n == 6:
    crystal =  SixthLayer(molecule)
elif n ==0:
    crystal = '['      
else:
    print 'The limitation is six layers right know, I will calculate the default crystal with just one layer'
    crystal = FirstLayer(molecule)
   
final_crystall = str(crystal)

'''
In this step, all the unnecessary string will be removed from crystal and molecule file.
'''

with open('Sources//cry.xyz','w') as final:
    final.write(final_crystall)
final.close()
with open('Sources//cry.xyz','r+') as fin:
    with open('Sources//crys.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace('], ','\n'))
#         
fin.close()
fout.close()
with open('Sources//crys.xyz','r+') as fin:
    with open('Sources//cryst.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace('[',''))
fin.close()
fout.close()
with open('Sources//cryst.xyz','r+') as fin:
    with open('Sources//crysta.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace(',','    '))
fin.close()
fout.close()
with open('Sources//crysta.xyz','r+') as fin:
    with open('Sources//crystal.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace("'",''))
fin.close()
fout.close()
with open('Sources//crystal.xyz','r+') as fin:
    with open('Sources//crystall.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace("]]",''))
fin.close()
fout.close()
os.remove('Sources//crysta.xyz')
os.remove('Sources//cryst.xyz')
os.remove('Sources//crys.xyz')
os.remove('Sources//cry.xyz')
os.remove('Sources//crystal.xyz')

###


###
mol = str(molecule)
with open('Sources//mo.xyz','w') as final:
    final.write(mol)
final.close()
with open('Sources//mo.xyz','r+') as fin:
    with open('Sources//mol.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace('], ','\n'))
#         
fin.close()
fout.close()
with open('Sources//mol.xyz','r+') as fin:
    with open('Sources//mole.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace('[',''))
fin.close()
fout.close()
with open('Sources//mole.xyz','r+') as fin:
    with open('Sources//molec.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace(',','    '))
fin.close()
fout.close()
with open('Sources//molec.xyz','r+') as fin:
    with open('Sources//molecu.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace("'",''))
fin.close()
fout.close()
with open('Sources//molecu.xyz','r+') as fin:
    with open('Sources//molecul.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace("]]",''))
fin.close()
fout.close()
with open('Sources//molecul.xyz','r+') as fin:
    with open('Sources//molecule1.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace("]",''))
            
            
fin.close()
fout.close()

with open('Sources//molecule1.xyz','r+') as fin:
    with open('Sources//molecule2.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace('     Name:',''))
fin.close()
fout.close()
with open('Sources//molecule2.xyz','r+') as fin:
    with open('Sources//molecule3.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace('Name:','    '))
fin.close()
fout.close()
with open('Sources//molecule3.xyz','r+') as fin:
    with open('Sources//molecule4.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace("q: 1.0",''))
fin.close()
fout.close()
with open('Sources//molecule4.xyz','r+') as fin:
    with open('Sources//molecule5.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace("X:",'    '))
fin.close()
fout.close()
with open('Sources//molecule5.xyz','r+') as fin:
    with open('Sources//molecule6.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace("Y:",'    '))
fin.close()
fout.close()
with open('Sources//molecule6.xyz','r+') as fin:
    with open('Sources//molecule7.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace("Z:",'    '))
fin.close()
fout.close()
with open('Sources//molecule7.xyz','r+') as fin:
    with open('Sources//molecule.xyz','w') as fout:
        for line in fin:
            fout.write(line.replace("    H",'H'))
fin.close()
fout.close()
os.remove('Sources//mo.xyz')
os.remove('Sources//mol.xyz')
os.remove('Sources//mole.xyz')
os.remove('Sources//molec.xyz')
os.remove('Sources//molecu.xyz')
os.remove('Sources//molecul.xyz')
os.remove('Sources//molecule6.xyz')
os.remove('Sources//molecule5.xyz')
os.remove('Sources//molecule4.xyz')
os.remove('Sources//molecule1.xyz')
os.remove('Sources//molecule3.xyz')
os.remove('Sources//molecule2.xyz')
os.remove('Sources//molecule7.xyz')
###

###

###
'''In this step, number of carbons in the molecule will be calculated so based on this number we can make input files
'''

n_carbon = 0
c_counter = 1
with open('Sources//molecule.xyz','r') as file:
    for line in file:
        if line[0] == 'C':
            n_carbon = n_carbon + 1
n_carbon = n_carbon + 1
ff = 'Inputs'
if os.path.exists(ff):
    shutil.rmtree(ff)
    os.mkdir(ff)
else:
    os.mkdir(ff)
for i in range(n_carbon):
    if os.path.exists(ff+'//c'+str(i)):
        shutil.rmtree(ff+'//c'+str(i))
        os.mkdir(ff+'//c'+str(i))
    else:
        os.mkdir(ff+'//c'+str(i))

''' In this step, based on numbers of carbons in the molecule, input files will be made in their corresponding 
folder.
First, lines of input.txt file will be copied to new input files and then
the coordinates of crystal will be appended to each input file
in the second loop of this step changing of C to C1 will happen
'''
        
for c_counter in range(n_carbon):
    attempt = 0
    with open('Sources//molecule.xyz','r') as file:
        with open('Sources//input.txt','r') as infile:
            with open(ff+'//c'+str(c_counter)+'//c'+str(c_counter)+'.inp','w') as new:
                for line in file:
                    if line[0] == 'C':
                        shutil.copy2('Sources//input.txt', ff+'//c'+str(c_counter)+'//c'+str(c_counter)+'.inp')
              
for c_counter in range(n_carbon):
    attempt = 0
    with open('Sources//molecule.xyz','r') as file:
        with open('Sources//crystall.xyz','r') as cryfile:
            with open(ff+'//c'+str(c_counter)+'//c'+str(c_counter)+'.inp','a') as new:
                for line in file:
                    if line[0] == 'C':
                        attempt = attempt + 1
                        if attempt < c_counter:
                            new.write(line.replace('C','C2'))
                        elif attempt > c_counter:
                            new.write(line.replace('C','C2'))
                        elif attempt == c_counter:
                            new.write(line.replace('C','C1'))
                            c_counter = c_counter + 1
                            attempt = attempt + n_carbon
                    else:
                        new.write(line)

###

for c_counter in range(n_carbon):
    attempt = 0
    with open('Sources//crystall.xyz','r') as cryfile:
        with open(ff+'//c'+str(c_counter)+'//c'+str(c_counter)+'.inp','a') as new:
            for lines in cryfile:
                new.write(lines)

###
with open('Sources//crystall.xyz','a') as cryfile:
    with open('Sources//molecule.xyz','r') as file:
        cryfile.write('\n')
        for line in file:
            cryfile.write(line)
cryfile.close()
file.close()
#shutil.rmtree(ff+'//c0') 
shutil.move('Sources//crystall.xyz', ff+'//crystall.xyz')

###


'''This part is to make an output folder so all files will be copied here based on their type'''
'''oo = 'outputs'
if os.path.exists(oo):
    shutil.rmtree(oo)
    os.mkdir(oo)
    os.mkdir(oo +'//xry')
    os.mkdir(oo +'//out')
    os.mkdir(oo +'//rst')
    os.mkdir(oo +'//xas')
    os.mkdir(oo +'//inp')
    os.mkdir(oo +'//mol')
else:
    os.mkdir(oo)
    os.mkdir(oo +'//xry')
    os.mkdir(oo +'//out')
    os.mkdir(oo +'//rst')
    os.mkdir(oo +'//xas')
    os.mkdir(oo +'//inp')
    os.mkdir(oo +'//mol')'''

'''
In this step, all the necessary files for calculation, such as batch file and
 xray files will be copied to each folder.
 '''

Server = raw_input('Which server are you using? p (plato) or g(Grex)?')

if os.path.isfile('d2k-batch.run'):
    os.remove('d2k-batch.run')

if Server == 'p' or Server == 'P':
    shutil.copy2('Sources//Plato//d2k-batch_Plato.run','Sources//d2k-batch.run')
    shutil.copy2('Sources//Plato//xray2k.run', 'Sources//xray2k.run')
    
elif Server == 'g' or Server == 'G':
    shutil.copy2('Sources//Grex//d2k-batch_Grex.run','Sources//d2k-batch.run')
    shutil.copy2('Sources//Grex//xray2k.run', 'Sources//xray2k.run')
else:
    print 'Please indicate the right server!'   
for i in range(n_carbon):
    with open('Sources//d2k-batch.run','r') as file:
        with open(ff+'//c'+str(i)+'//d2k-batch.run','w') as batch:
            for line in file:
                batch.write(line.replace('TBA','c'+str(i)+'.inp'))
    shutil.copy2('Sources//Xray//xray.inp', ff+'//c'+str(i))
    shutil.copy2('Sources//xray2k.run', ff+'//c'+str(i))
    shutil.copy2('Sources//Xray//xray2k.x', ff+'//c'+str(i))
aa = os.curdir

os.remove('Sources//d2k-batch.run')
os.remove('Sources//molecule.xyz')
os.remove('Sources//xray2k.run')


'''
in this step, all the Q jobs will be submitted
if you want to run the program on a system other than westgrid,please activate this step

'''

'''for i in range(1,n_carbon):
    ChDir(ff+'//c'+str(i))
'''



print('Its all done...')
