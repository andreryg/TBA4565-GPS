import numpy as np
from scipy.constants import speed_of_light


class GPS(object):
    t = 129600
    we = 7.2921151467 * 10**-5
    GM = 3.986005 * 10**14

    def __init__(self, toe, sqrta, e, M0, w, i0, lambda0, deltan, i, omega, cuc, cus, crc, crs, cic, cis):
        self.toe = toe
        self.sqrta = sqrta
        self.e = e
        self.M0 = M0
        self.w = w
        self.i0 = i0
        self.lambda0 = lambda0
        self.deltan = deltan
        self.i = i
        self.omega = omega
        self.cuc = cuc
        self.cus = cus
        self.crc = crc
        self.crs = crs
        self.cic = cic
        self.cis = cis

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
    for i,v in enumerate(lines):
        print(i,v)
    #Task 1
    SV06 = GPS(1.295840*10**5, 5.153618*10**3, 5.747278*10**-3, -2.941505, -1.770838, 9.332837*10**-1, 2.123898, 5.243075*10**-9, -6.853856*10**-10, -8.116052*10**-9, -1.184642*10**-6, 7.672235*10**-6, 2.146562*10**2, -2.140625*10**1, 2.980232*10**-8, -1.117587*10**-8)
    SV10 = GPS(1.296*10**5, 5.153730*10**3, 7.258582*10**-3, 4.044839*10**-1, 4.344642*10**-1, 9.71311*10**-1, -2.006987, 4.442685*10**-9, 2.521533*10**-10, -8.495353*10**-9, 4.714354*10**-6, -1.825392*10**-7, 3.868750*10**2, 8.978125*10**1, 3.725290*10**-9, 8.940696*10**-8)
    SN16 = GPS(1.296*10**5, 5.153541*10**3, 3.506405*10**-3, 1.808249, -7.60081*10**-1, 9.624682*10**-1, 1.122991, 4.937348*10**-9, 2.367955*10**-10, -8.054621*10**-9, 9.49949*10**-7, 5.437061*10**-6, 2.709062*10**2, 1.515625*10**1, 6.332993*10**-8, -2.421438*10**-8)
    SV21 = GPS(1.29584*10**5, 5.153681*10**3, 1.179106*10**-2, 3.122437, -2.904128, 9.416507*10**-1, -3.042819, 4.445542*10**-9, -4.035882*10**-11, -7.757823*10**-9, 6.897374*10**-6, 1.069344*10**-5, 1.630625*10**2, 1.329375*10**2, -1.080334*10**-7, -8.009374*10**-8)




    SV06.calculateXYZ()
    SV10.calculateXYZ()
    SN16.calculateXYZ()
    SV21.calculateXYZ()
    print(SV06.X, SV06.Y, SV06.Z)
    print(SV10.X, SV10.Y, SV10.Z)
    print(SN16.X, SN16.Y, SN16.Z)
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
    lat = np.deg2rad(47.1)
    long = np.deg2rad(15.5)
    a, b, h = 6378137, 6356752.3141, 400
    Xr, Yr, Zr = np.array([(N(lat)+h)*np.cos(lat)*np.cos(long), (N(lat)+h)*np.cos(lat)*np.sin(long), ((b**2/a**2)*N(lat)+h)*np.sin(lat)])

    print(Xr, Yr, Zr)

    #Task 4
    PSV06, PSV10, PSN16, PSV21 = 20509078.908, 23568574.070, 23733776.587, 22106790.995
    for i in range(10):
        SV06.set_rho(Xr, Yr, Zr)
        SV10.set_rho(Xr, Yr, Zr)
        SN16.set_rho(Xr, Yr, Zr)
        SV21.set_rho(Xr, Yr, Zr)
        
        A = np.array([
            [-(SV06.X-Xr)/SV06.rho, -(SV06.Y-Yr)/SV06.rho, -(SV06.Z-Zr)/SV06.rho, -speed_of_light],
            [-(SV10.X-Xr)/SV10.rho, -(SV10.Y-Yr)/SV10.rho, -(SV10.Z-Zr)/SV10.rho, -speed_of_light],
            [-(SN16.X-Xr)/SN16.rho, -(SN16.Y-Yr)/SN16.rho, -(SN16.Z-Zr)/SN16.rho, -speed_of_light],
            [-(SV21.X-Xr)/SV21.rho, -(SV21.Y-Yr)/SV21.rho, -(SV21.Z-Zr)/SV21.rho, -speed_of_light],
        ])
        L = np.array([
            [PSV06 - SV06.rho],
            [PSV10 - SV10.rho],
            [PSN16 - SN16.rho],
            [PSV21 - SV21.rho]
        ])
        
        dXr, dYr, dZr, dTr = np.linalg.inv(A)@L
        Xr += float(dXr)
        Yr += float(dYr)
        Zr += float(dZr)

    print(Xr, Yr, Zr)
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