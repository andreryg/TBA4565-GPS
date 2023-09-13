import numpy as np
from scipy.constants import speed_of_light
from itertools import compress

class GPS(object):
    t = 558000
    we = 7.2921151467 * 10**-5
    GM = 3.986005 * 10**14

    def __init__(self, p):
        self.toe = p[7]
        self.sqrta = p[6]
        self.e = p[4]
        self.M0 = p[2]
        self.w = p[13]
        self.i0 = p[11]
        self.lambda0 = p[9]
        self.deltan = p[1]
        self.i = p[15]
        self.omega = p[14]
        self.cuc = p[3]
        self.cus = p[5]
        self.crc = p[12]
        self.crs = p[0]
        self.cic = p[8]
        self.cis = p[10]
        self.P = p[16]
        self.dt = p[17]
        self.dion = p[18]
        self.dtrop = p[19]

    def disregard_corrections(self):
        self.deltan, self.i, self.omega, self.cuc, self.cus, self.crc, self.crs, self.cic, self.cis = [eval(i) for i in ["0"]*9]
        return True

    def _tk(self):
        tk = self.t-self.toe
        if tk > 302400:
            tk -= 604800
        elif tk < -302400:
            tk += 604800
        return tk

    def _lambdak(self):
        return self.lambda0 + (self.omega-self.we)*self._tk() - self.we*self.toe

    def _Mk(self):
        return self.M0 + (np.sqrt(self.GM/(self.sqrta**6)) + self.deltan) * self._tk()

    def _EK(self):
        Ek = self._Mk()
        for i in range(3):
            Ek = Ek + (self._Mk()-Ek+self.e*np.sin(Ek))/(1-self.e*np.cos(Ek))
        return Ek

    def _fk(self):
        return 2 * np.arctan(np.sqrt((1+self.e)/(1-self.e))*np.tan(self._EK()/2))

    def _ik(self):
        return self.i0 + self.i*self._tk() + self.cic*np.cos(2*(self.w+self._fk())) + self.cis*np.sin(2*(self.w+self._fk()))

    def _uk(self):
        return self.w + self._fk() + self.cuc*np.cos(2*(self.w+self._fk())) + self.cus*np.sin(2*(self.w+self._fk()))

    def _rk(self):
        #print((self.sqrta**2)*(1-self.e*np.cos(self._EK())) + self.crc*np.cos(self.w+self._fk())**2 + self.crs*np.sin(self.w+self._fk())**2)
        return (self.sqrta**2)*(1-self.e*np.cos(self._EK())) + self.crc*np.cos(2*(self.w+self._fk())) + self.crs*np.sin(2*(self.w+self._fk()))

    def _R3(self, theta):
        return np.array([[np.cos(theta), np.sin(theta), 0], [-np.sin(theta), np.cos(theta), 0], [0, 0, 1]])

    def _R1(self, theta):
        return np.array([[1, 0, 0], [0, np.cos(theta), np.sin(theta)], [0, -np.sin(theta), np.cos(theta)]])
    
    def calculateXYZ(self):
        self.X, self.Y, self.Z = self._R3(-self._lambdak())@self._R1(-self._ik())@self._R3(-self._uk())@np.array([[self._rk()], [0], [0]])
        self.X, self.Y, self.Z = float(self.X), float(self.Y), float(self.Z)
        return True
    
    def set_rho(self, Xr, Yr, Zr):
        self.rho = float(np.sqrt((self.X-Xr)**2 + (self.Y-Yr)**2 + (self.Z-Zr)**2))
        return True

def N(lat, a=6378137, b=6356752.3141):
    return a**2 / (np.sqrt(a**2*np.cos(lat)**2+b**2*np.sin(lat)**2))

def geodetic_iteration(x, y, z, e2):
    p = np.sqrt(x**2 + y**2)
    lat0 = np.arctan(z / p * (1-e2)**-1)
    while True:
        N0 = N(lat0)
        h = p / np.cos(lat0) - N0
        lat = np.arctan((z / p) * (1 - e2* N0/(N0+h))**-1)
        if abs(lat-lat0) <= 10**-8:
            h = p / np.cos(lat) - N0
            return np.rad2deg(lat), h
        else:
            lat0 = lat

def main():
    with open('ephemeries.txt') as f:
        lines = f.readlines()
    lines = lines[7:]
    for i in range(0, 7):
        lines[i] += ''.join(lines[i+1:i+8])
        del lines[i+1:i+8]
    relevant_values = [0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0]
    for i,v in enumerate(lines):
        v = list(compress(v.split(), relevant_values))
        lines[i] = [float(x.replace('D', 'E')) for x in v]
    #Task 1
    SV08 = GPS(lines[0]+[22550792.660, 0.00013345632, 3.344, 4.055])
    SV10 = GPS(lines[1]+[22612136.900, 0.00004615571, 2.947, 4.297])
    SV21 = GPS(lines[2]+[20754631.240, -0.00015182034, 2.505, 2.421])
    SV24 = GPS(lines[3]+[23974471.500, 0.00026587520, 3.644, 9.055])
    SV17 = GPS(lines[4]+[24380357.760, -0.00072144074, 6.786, 9.756])
    SV03 = GPS(lines[5]+[24444143.500, 0.00022187057, 4.807, 10.863])
    SN14 = GPS(lines[6]+[22891323.280, -0.00013020719, 4.598, 4.997])

    SV08.calculateXYZ()
    SV10.calculateXYZ()
    SV21.calculateXYZ()
    SV24.calculateXYZ()
    SV17.calculateXYZ()
    SV03.calculateXYZ()
    SN14.calculateXYZ()
    print(SV08.X, SV08.Y, SV08.Z)
    print(SV10.X, SV10.Y, SV10.Z)
    print("ANSWER", SV03.X, SV03.Y, SV03.Z)
    print(SV21.X, SV21.Y, SV21.Z)

    #Task 2
    """ SV06.disregard_corrections()
    SV06.calculateXYZ()
    SV10.disregard_corrections()
    SV10.calculateXYZ()
    SN16.disregard_corrections()
    SN16.calculateXYZ()
    SV21.disregard_corrections()
    SV21.calculateXYZ() 
    print(SV06.X, SV06.Y, SV06.Z)
    print(SV10.X, SV10.Y, SV10.Z)
    print(SN16.X, SN16.Y, SN16.Z)
    print(SV21.X, SV21.Y, SV21.Z) """


    #Task 3
    lat = np.deg2rad(63.2)
    long = np.deg2rad(10.2)
    a, b, h = 6378137, 6356752.3141, 400
    Xr, Yr, Zr = np.array([(N(lat)+h)*np.cos(lat)*np.cos(long), (N(lat)+h)*np.cos(lat)*np.sin(long), ((b**2/a**2)*N(lat)+h)*np.sin(lat)])

    print("XYZ: ",Xr, Yr, Zr)

    #Task 4
    for i in range(10):
        SV08.set_rho(Xr, Yr, Zr)
        SV10.set_rho(Xr, Yr, Zr)
        SV21.set_rho(Xr, Yr, Zr)
        SV24.set_rho(Xr, Yr, Zr)
        SV17.set_rho(Xr, Yr, Zr)
        SV03.set_rho(Xr, Yr, Zr)
        SN14.set_rho(Xr, Yr, Zr)
        
        A = np.array([
            [-(SV08.X-Xr)/SV08.rho, -(SV08.Y-Yr)/SV08.rho, -(SV08.Z-Zr)/SV08.rho, -speed_of_light],
            [-(SV10.X-Xr)/SV10.rho, -(SV10.Y-Yr)/SV10.rho, -(SV10.Z-Zr)/SV10.rho, -speed_of_light],
            [-(SV21.X-Xr)/SV21.rho, -(SV21.Y-Yr)/SV21.rho, -(SV21.Z-Zr)/SV21.rho, -speed_of_light],
            [-(SV24.X-Xr)/SV24.rho, -(SV24.Y-Yr)/SV24.rho, -(SV24.Z-Zr)/SV24.rho, -speed_of_light],
            [-(SV17.X-Xr)/SV17.rho, -(SV17.Y-Yr)/SV17.rho, -(SV17.Z-Zr)/SV17.rho, -speed_of_light],
            [-(SV03.X-Xr)/SV03.rho, -(SV03.Y-Yr)/SV03.rho, -(SV03.Z-Zr)/SV03.rho, -speed_of_light],
            [-(SN14.X-Xr)/SN14.rho, -(SN14.Y-Yr)/SN14.rho, -(SN14.Z-Zr)/SN14.rho, -speed_of_light]
        ])
        L = np.array([
            [SV08.P - SV08.rho - speed_of_light*SV08.dt - SV08.dion - SV08.dtrop],
            [SV10.P - SV10.rho - speed_of_light*SV10.dt - SV10.dion - SV10.dtrop],
            [SV21.P - SV21.rho - speed_of_light*SV21.dt - SV21.dion - SV21.dtrop],
            [SV24.P - SV24.rho - speed_of_light*SV24.dt - SV24.dion - SV24.dtrop],
            [SV17.P - SV17.rho - speed_of_light*SV17.dt - SV17.dion - SV17.dtrop],
            [SV03.P - SV03.rho - speed_of_light*SV03.dt - SV03.dion - SV03.dtrop],
            [SN14.P - SN14.rho - speed_of_light*SN14.dt - SN14.dion - SN14.dtrop]
        ])
        
        dXr, dYr, dZr, dTr = np.linalg.inv(A.T@A)@A.T@L
        Xr += float(dXr)
        Yr += float(dYr)
        Zr += float(dZr)

    print("XYZ: ", Xr, Yr, Zr)
    print("dTr", dTr)

    #Task 5
    Qx = np.linalg.inv(A.T@A)
    PDOP = np.sqrt(Qx[0][0] + Qx[1][1] + Qx[2][2])

    print(PDOP)

    #Task 6
    long = np.rad2deg(np.arctan(Yr/Xr))
    e = (a**2-b**2)/a**2
    lat, h = geodetic_iteration(Xr, Yr, Zr, e)

    print(lat, long, h)

    #Task 7

    print(dTr)


if __name__=="__main__":
    main()