import numpy as np
from scipy.constants import speed_of_light
import math

class satellite(object):
    def __init__(self, coordinates_A, carrier_phase_A, coordinates_B, carrier_phase_B):
        self.X_A = coordinates_A[0]
        self.Y_A = coordinates_A[1]
        self.Z_A = coordinates_A[2]
        self.carrier_phase_A = carrier_phase_A
        self.X_B = coordinates_B[0]
        self.Y_B = coordinates_B[1]
        self.Z_B = coordinates_B[2]
        self.carrier_phase_B = carrier_phase_B
    
    def set_rho_A(self, Xr, Yr, Zr):
        self.rho_A = float(np.sqrt((self.X_A-Xr)**2 + (self.Y_A-Yr)**2 + (self.Z_A-Zr)**2))
        return True

    def set_rho_B(self, Xr, Yr, Zr):
        self.rho_B = float(np.sqrt((self.X_B-Xr)**2 + (self.Y_B-Yr)**2 + (self.Z_B-Zr)**2))
        return True

def ellipsoidal_to_cartesian(latitude, longitude, h, a=6378137, b=6356752.3141):
    latitude, longitude = np.deg2rad(latitude), np.deg2rad(longitude)
    N = a**2 / (np.sqrt(a**2*np.cos(latitude)**2+b**2*np.sin(latitude)**2))
    X = (N+h) * np.cos(latitude) * np.cos(longitude)
    Y = (N+h) * np.cos(latitude) * np.sin(longitude)
    Z = (((b**2 / a**2) * N + h) * np.sin(latitude))
    return X, Y, Z

def a(satellite, main_satellite, r):
    return [-(satellite.X_B - r[0])/satellite.rho_B + (main_satellite.X_B - r[0])/main_satellite.rho_B, -(satellite.Y_B - r[1])/satellite.rho_B + (main_satellite.Y_B - r[1])/main_satellite.rho_B, -(satellite.Z_B - r[2])/satellite.rho_B + (main_satellite.Z_B - r[2])/main_satellite.rho_B]

def phi(satellite, main_satellite):
    return satellite.carrier_phase_B - main_satellite.carrier_phase_B - satellite.carrier_phase_A + main_satellite.carrier_phase_A

def N_value(lat, a=6378137, b=6356752.3141):
    return a**2 / (np.sqrt(a**2*np.cos(lat)**2+b**2*np.sin(lat)**2))

def geodetic_iteration(x, y, z, e2):
    #print("x", x, "y", y, "z", z)
    p = np.sqrt(x**2 + y**2)
    lat0 = np.arctan(z / p * (1-e2)**-1)
    while True:
        N0 = N_value(lat0)
        h = p / np.cos(lat0) - N0
        lat = np.arctan((z / p) * (1 - e2* N0/(N0+h))**-1)
        if abs(lat-lat0) <= 10**-8:
            h = p / np.cos(lat) - N0
            return np.rad2deg(lat), h
        else:
            lat0 = lat

def main():
    A_latitude = -32.003884648
    A_longitude = 115.894802001
    A_h = 23.983
    B_latitude = -31.9
    B_longitude = 115.75
    B_h = 50

    #TASK 1
    A_X, A_Y, A_Z = ellipsoidal_to_cartesian(A_latitude, A_longitude, A_h)
    B_X, B_Y, B_Z = ellipsoidal_to_cartesian(B_latitude, B_longitude, B_h)
    print(A_X, A_Y, A_Z)
    print(B_X, B_Y, B_Z)

    P = 1/(10*0.002**2) * np.array([
        [4,-1,-1,-1,0,0,0,0],
        [-1,4,-1,-1,0,0,0,0],
        [-1,-1,4,-1,0,0,0,0],
        [-1,-1,-1,4,0,0,0,0],
        [0,0,0,0,4,-1,-1,-1],
        [0,0,0,0,-1,4,-1,-1],
        [0,0,0,0,-1,-1,4,-1],
        [0,0,0,0,-1,-1,-1,4]
    ])

    #TASK 2
    f = 1575.42 * 10**6
    wavelength = speed_of_light/f
    t1 = 172800
    t2 = 175020
    AB154_t1 = satellite([-26916298.03, -2738678.66, -11996368.48], 143588831.82, [-26916297.22, -2738678.33, -11996370.37], 144732675.77)
    AB155_t1 = satellite([7525690.37, 19506497.68, -20949608.30], 130653024.41, [7525691.76, 19506498.31, -20949607.21], 131796866.21)
    AB159_t1 = satellite([-13393419.83, 11968129.58, -23519577.49], 126717695.82, [-13393418.26, 11968130.27, -23519578.04], 127861541.27)
    AB174_t1 = satellite([-21266103.93, 19852330.46, 5502623.54], 135590276.69, [-21266103.64, 19852330.18, 5502625.68], 136734190.29)
    AB181_t1 = satellite([-25697148.88, 8478510.63, -11983680.13], 132277486.86, [-25697149.62, 8478511.10, -11983678.21], 133421353.12)
    AB154_t2 = satellite([-28824855.18, -3331684.88, -5833025.24], 148445614.13, [-28824855.42, -3331684.93, -5833024.04], 147792428.87)
    AB155_t2 = satellite([2847363.30, 17829892.74, -23451064.94], 129436683.11, [2847362.36, 17829892.49, -23451065.25], 128783480.33)
    AB159_t2 = satellite([-18112952.18, 10371378.82, -20977224.27], 129024672.38, [-18112953.03, 10371378.61, -20977223.64], 128371483.32)
    AB174_t2 = satellite([-21527125.49, 20297088.12, -1169497.61], 131540267.25, [-21527125.43, 20297088.12, -1169498.87], 130887130.25)
    AB181_t2 = satellite([-23104316.90, 6456482.66, -17326698.31], 133044764.25, [-23104316.37, 6456482.17, -17326699.19], 132391579.66)

    satellites = [AB154_t1, AB155_t1, AB159_t1, AB174_t1, AB181_t1, AB154_t2, AB155_t2, AB159_t2, AB174_t2, AB181_t2]

    for xdd in range(10):
        #print(B_X, B_Y, B_Z)
        for s in satellites:
            s.set_rho_A(A_X, A_Y, A_Z)
            s.set_rho_B(B_X, B_Y, B_Z)
            
        A = np.array([
            a(AB155_t1, AB154_t1, [B_X, B_Y, B_Z])+[wavelength, 0, 0, 0],
            a(AB159_t1, AB154_t1, [B_X, B_Y, B_Z])+[0, wavelength, 0, 0],
            a(AB174_t1, AB154_t1, [B_X, B_Y, B_Z])+[0, 0, wavelength, 0],
            a(AB181_t1, AB154_t1, [B_X, B_Y, B_Z])+[0, 0, 0, wavelength],
            a(AB155_t2, AB154_t2, [B_X, B_Y, B_Z])+[wavelength, 0, 0, 0],
            a(AB159_t2, AB154_t2, [B_X, B_Y, B_Z])+[0, wavelength, 0, 0],
            a(AB174_t2, AB154_t2, [B_X, B_Y, B_Z])+[0, 0, wavelength, 0],
            a(AB181_t2, AB154_t2, [B_X, B_Y, B_Z])+[0, 0, 0, wavelength]
        ])
        #print(A[0])
        deltaL = np.array([
            [wavelength * phi(AB155_t1, AB154_t1) - AB155_t1.rho_B + AB154_t1.rho_B + AB155_t1.rho_A - AB154_t1.rho_A],
            [wavelength * phi(AB159_t1, AB154_t1) - AB159_t1.rho_B + AB154_t1.rho_B + AB159_t1.rho_A - AB154_t1.rho_A],
            [wavelength * phi(AB174_t1, AB154_t1) - AB174_t1.rho_B + AB154_t1.rho_B + AB174_t1.rho_A - AB154_t1.rho_A],
            [wavelength * phi(AB181_t1, AB154_t1) - AB181_t1.rho_B + AB154_t1.rho_B + AB181_t1.rho_A - AB154_t1.rho_A],
            [wavelength * phi(AB155_t2, AB154_t2) - AB155_t2.rho_B + AB154_t2.rho_B + AB155_t2.rho_A - AB154_t2.rho_A],
            [wavelength * phi(AB159_t2, AB154_t2) - AB159_t2.rho_B + AB154_t2.rho_B + AB159_t2.rho_A - AB154_t2.rho_A],
            [wavelength * phi(AB174_t2, AB154_t2) - AB174_t2.rho_B + AB154_t2.rho_B + AB174_t2.rho_A - AB154_t2.rho_A],
            [wavelength * phi(AB181_t2, AB154_t2) - AB181_t2.rho_B + AB154_t2.rho_B + AB181_t2.rho_A - AB154_t2.rho_A]
        ])
        dB_X, dB_Y, dB_Z, dN1, dN2, dN3, dN4 = np.linalg.inv(A.T@P@A)@A.T@P@deltaL
        deltaX = np.array([
            [dB_X, dB_Y, dB_Z, dN1, dN2, dN3, dN4]
        ])
        #print(deltaX)
        #print(dB_X, dB_Y, dB_Z, dN1, dN2, dN3, dN4)
        #print(np.rad2deg(np.arctan(B_Y/B_X)))
        B_X += float(dB_X)
        B_Y += float(dB_Y)
        B_Z += float(dB_Z)
        
    print(A_X, B_X, A_Y, B_Y, A_Z, B_Z)
    print(dN1, dN2, dN3, dN4)
    #print(np.cov(A))
    

    Qx = np.linalg.inv(A.T@P@A)
    print(np.diag(Qx))
    f2 = lambda x: np.sqrt(x)
    dank = f2(np.diag(Qx))
    print("ddd", dank)
    

    """#FJERN
    a_=6378137
    b_=6356752.3141
    e = (a_**2-b_**2)/a_**2
    print("e", e)
    B_long = np.rad2deg(np.arctan(B_Y/B_X)+np.pi)
    B_lat, B_h = geodetic_iteration(B_X, B_Y, B_Z, e)
    print(B_lat, B_long, B_h)"""

    #TASK 3
    for xdd in range(10):
        #print(B_X, B_Y, B_Z)
        for s in satellites:
            s.set_rho_A(A_X, A_Y, A_Z)
            s.set_rho_B(B_X, B_Y, B_Z)

        A = np.array([
            a(AB155_t1, AB154_t1, [B_X, B_Y, B_Z]),
            a(AB159_t1, AB154_t1, [B_X, B_Y, B_Z]),
            a(AB174_t1, AB154_t1, [B_X, B_Y, B_Z]),
            a(AB181_t1, AB154_t1, [B_X, B_Y, B_Z]),
            a(AB155_t2, AB154_t2, [B_X, B_Y, B_Z]),
            a(AB159_t2, AB154_t2, [B_X, B_Y, B_Z]),
            a(AB174_t2, AB154_t2, [B_X, B_Y, B_Z]),
            a(AB181_t2, AB154_t2, [B_X, B_Y, B_Z])
        ])
        #print(A[0])
        deltaL = np.array([
            wavelength * phi(AB155_t1, AB154_t1) - AB155_t1.rho_B + AB154_t1.rho_B + AB155_t1.rho_A - AB154_t1.rho_A - wavelength * dN1[0],
            wavelength * phi(AB159_t1, AB154_t1) - AB159_t1.rho_B + AB154_t1.rho_B + AB159_t1.rho_A - AB154_t1.rho_A - wavelength * dN2[0],
            wavelength * phi(AB174_t1, AB154_t1) - AB174_t1.rho_B + AB154_t1.rho_B + AB174_t1.rho_A - AB154_t1.rho_A - wavelength * dN3[0],
            wavelength * phi(AB181_t1, AB154_t1) - AB181_t1.rho_B + AB154_t1.rho_B + AB181_t1.rho_A - AB154_t1.rho_A - wavelength * dN4[0],
            wavelength * phi(AB155_t2, AB154_t2) - AB155_t2.rho_B + AB154_t2.rho_B + AB155_t2.rho_A - AB154_t2.rho_A - wavelength * dN1[0],
            wavelength * phi(AB159_t2, AB154_t2) - AB159_t2.rho_B + AB154_t2.rho_B + AB159_t2.rho_A - AB154_t2.rho_A - wavelength * dN2[0],
            wavelength * phi(AB174_t2, AB154_t2) - AB174_t2.rho_B + AB154_t2.rho_B + AB174_t2.rho_A - AB154_t2.rho_A - wavelength * dN3[0],
            wavelength * phi(AB181_t2, AB154_t2) - AB181_t2.rho_B + AB154_t2.rho_B + AB181_t2.rho_A - AB154_t2.rho_A - wavelength * dN4[0]
        ])
        #print(AB154_t1.rho_B)
        dB_X, dB_Y, dB_Z = np.linalg.inv(A.T@P@A)@A.T@P@deltaL
        #print(dB_X, dB_Y, dB_Z)
        B_X += float(dB_X)
        B_Y += float(dB_Y)
        B_Z += float(dB_Z)
        """N1 += float(dN1)
        N2 += float(dN2)
        N3 += float(dN3)
        N4 += float(dN4)
        N5 += float(dN1)
        N6 += float(dN2)
        N7 += float(dN3)
        N8 += float(dN4)"""
    print(A_X, B_X, A_Y, B_Y, A_Z, B_Z)
    """x = np.array([
        [dX_B],
        [dY_B],
        [dZ_B],
        [AB155_t1.carrier_phase_B - AB154_t1.carrier_phase_B - AB155_t1.carrier_phase_A + AB154_t1.carrier_phase_A],
        [AB159_t1.carrier_phase_B - AB154_t1.carrier_phase_B - AB159_t1.carrier_phase_A + AB154_t1.carrier_phase_A],
        [AB174_t1.carrier_phase_B - AB154_t1.carrier_phase_B - AB174_t1.carrier_phase_A + AB154_t1.carrier_phase_A],
        [AB181_t1.carrier_phase_B - AB154_t1.carrier_phase_B - AB181_t1.carrier_phase_A + AB154_t1.carrier_phase_A]
    ]) """ 

    a_=6378137
    b_=6356752.3141
    e = (a_**2-b_**2)/a_**2
    B_long = np.rad2deg(np.arctan(B_Y/B_X)+np.pi)
    B_lat, B_h = geodetic_iteration(B_X, B_Y, B_Z, e)
    print(B_lat, B_long, B_h)

    Qx = np.linalg.inv(A.T@P@A)
    print(np.diag(Qx))
    f2 = lambda x: np.sqrt(x)
    dank2 = f2(np.diag(Qx))
    print(dank2)

    for xdd in range(10):
        #print(B_X, B_Y, B_Z)
        for s in satellites:
            s.set_rho_A(A_X, A_Y, A_Z)
            s.set_rho_B(B_X, B_Y, B_Z)

        A = np.array([
            a(AB155_t1, AB154_t1, [B_X, B_Y, B_Z]),
            a(AB159_t1, AB154_t1, [B_X, B_Y, B_Z]),
            a(AB174_t1, AB154_t1, [B_X, B_Y, B_Z]),
            a(AB181_t1, AB154_t1, [B_X, B_Y, B_Z]),
            a(AB155_t2, AB154_t2, [B_X, B_Y, B_Z]),
            a(AB159_t2, AB154_t2, [B_X, B_Y, B_Z]),
            a(AB174_t2, AB154_t2, [B_X, B_Y, B_Z]),
            a(AB181_t2, AB154_t2, [B_X, B_Y, B_Z])
        ])
        #print(A[0])
        deltaL = np.array([
            wavelength * phi(AB155_t1, AB154_t1) - AB155_t1.rho_B + AB154_t1.rho_B + AB155_t1.rho_A - AB154_t1.rho_A - wavelength * round(dN1[0]),
            wavelength * phi(AB159_t1, AB154_t1) - AB159_t1.rho_B + AB154_t1.rho_B + AB159_t1.rho_A - AB154_t1.rho_A - wavelength * round(dN2[0]),
            wavelength * phi(AB174_t1, AB154_t1) - AB174_t1.rho_B + AB154_t1.rho_B + AB174_t1.rho_A - AB154_t1.rho_A - wavelength * round(dN3[0]),
            wavelength * phi(AB181_t1, AB154_t1) - AB181_t1.rho_B + AB154_t1.rho_B + AB181_t1.rho_A - AB154_t1.rho_A - wavelength * round(dN4[0]),
            wavelength * phi(AB155_t2, AB154_t2) - AB155_t2.rho_B + AB154_t2.rho_B + AB155_t2.rho_A - AB154_t2.rho_A - wavelength * round(dN1[0]),
            wavelength * phi(AB159_t2, AB154_t2) - AB159_t2.rho_B + AB154_t2.rho_B + AB159_t2.rho_A - AB154_t2.rho_A - wavelength * round(dN2[0]),
            wavelength * phi(AB174_t2, AB154_t2) - AB174_t2.rho_B + AB154_t2.rho_B + AB174_t2.rho_A - AB154_t2.rho_A - wavelength * round(dN3[0]),
            wavelength * phi(AB181_t2, AB154_t2) - AB181_t2.rho_B + AB154_t2.rho_B + AB181_t2.rho_A - AB154_t2.rho_A - wavelength * round(dN4[0])
        ])
        #print(AB154_t1.rho_B)
        dB_X, dB_Y, dB_Z = np.linalg.inv(A.T@P@A)@A.T@P@deltaL
        #print(dB_X, dB_Y, dB_Z)
        B_X += float(dB_X)
        B_Y += float(dB_Y)
        B_Z += float(dB_Z)
        """N1 += float(dN1)
        N2 += float(dN2)
        N3 += float(dN3)
        N4 += float(dN4)
        N5 += float(dN1)
        N6 += float(dN2)
        N7 += float(dN3)
        N8 += float(dN4)"""
    print(A_X, B_X, A_Y, B_Y, A_Z, B_Z)
    """x = np.array([
        [dX_B],
        [dY_B],
        [dZ_B],
        [AB155_t1.carrier_phase_B - AB154_t1.carrier_phase_B - AB155_t1.carrier_phase_A + AB154_t1.carrier_phase_A],
        [AB159_t1.carrier_phase_B - AB154_t1.carrier_phase_B - AB159_t1.carrier_phase_A + AB154_t1.carrier_phase_A],
        [AB174_t1.carrier_phase_B - AB154_t1.carrier_phase_B - AB174_t1.carrier_phase_A + AB154_t1.carrier_phase_A],
        [AB181_t1.carrier_phase_B - AB154_t1.carrier_phase_B - AB181_t1.carrier_phase_A + AB154_t1.carrier_phase_A]
    ]) """ 

    a_=6378137
    b_=6356752.3141
    e = (a_**2-b_**2)/a_**2
    B_long = np.rad2deg(np.arctan(B_Y/B_X)+np.pi)
    B_lat, B_h = geodetic_iteration(B_X, B_Y, B_Z, e)
    print(B_lat, B_long, B_h)

    h_values = []
    other = []
    
    for i in range(math.ceil(dank[3])+1):
        for j in range(math.ceil(dank[4])+1):
            for k in range(math.ceil(dank[5])+1):
                for l in range(math.ceil(dank[6])+1):
                    for xdd in range(10):
                        #print(B_X, B_Y, B_Z)
                        for s in satellites:
                            s.set_rho_A(A_X, A_Y, A_Z)
                            s.set_rho_B(B_X, B_Y, B_Z)

                        A = np.array([
                            a(AB155_t1, AB154_t1, [B_X, B_Y, B_Z]),
                            a(AB159_t1, AB154_t1, [B_X, B_Y, B_Z]),
                            a(AB174_t1, AB154_t1, [B_X, B_Y, B_Z]),
                            a(AB181_t1, AB154_t1, [B_X, B_Y, B_Z]),
                            a(AB155_t2, AB154_t2, [B_X, B_Y, B_Z]),
                            a(AB159_t2, AB154_t2, [B_X, B_Y, B_Z]),
                            a(AB174_t2, AB154_t2, [B_X, B_Y, B_Z]),
                            a(AB181_t2, AB154_t2, [B_X, B_Y, B_Z])
                        ])
                        #print(A[0])
                        deltaL = np.array([
                            wavelength * phi(AB155_t1, AB154_t1) - AB155_t1.rho_B + AB154_t1.rho_B + AB155_t1.rho_A - AB154_t1.rho_A - wavelength * round(dN1[0])+i,
                            wavelength * phi(AB159_t1, AB154_t1) - AB159_t1.rho_B + AB154_t1.rho_B + AB159_t1.rho_A - AB154_t1.rho_A - wavelength * round(dN2[0])+j,
                            wavelength * phi(AB174_t1, AB154_t1) - AB174_t1.rho_B + AB154_t1.rho_B + AB174_t1.rho_A - AB154_t1.rho_A - wavelength * round(dN3[0])+k,
                            wavelength * phi(AB181_t1, AB154_t1) - AB181_t1.rho_B + AB154_t1.rho_B + AB181_t1.rho_A - AB154_t1.rho_A - wavelength * round(dN4[0])+l,
                            wavelength * phi(AB155_t2, AB154_t2) - AB155_t2.rho_B + AB154_t2.rho_B + AB155_t2.rho_A - AB154_t2.rho_A - wavelength * round(dN1[0])+i,
                            wavelength * phi(AB159_t2, AB154_t2) - AB159_t2.rho_B + AB154_t2.rho_B + AB159_t2.rho_A - AB154_t2.rho_A - wavelength * round(dN2[0])+j,
                            wavelength * phi(AB174_t2, AB154_t2) - AB174_t2.rho_B + AB154_t2.rho_B + AB174_t2.rho_A - AB154_t2.rho_A - wavelength * round(dN3[0])+k,
                            wavelength * phi(AB181_t2, AB154_t2) - AB181_t2.rho_B + AB154_t2.rho_B + AB181_t2.rho_A - AB154_t2.rho_A - wavelength * round(dN4[0])+l
                        ])
                        #print(AB154_t1.rho_B)
                        dB_X, dB_Y, dB_Z = np.linalg.inv(A.T@P@A)@A.T@P@deltaL
                        #print(dB_X, dB_Y, dB_Z)
                        B_X += float(dB_X)
                        B_Y += float(dB_Y)
                        B_Z += float(dB_Z)
                        
                    print("m", B_X, B_Y, B_Z)
                    

                    a_=6378137
                    b_=6356752.3141
                    e = (a_**2-b_**2)/a_**2
                    B_long = np.rad2deg(np.arctan(B_Y/B_X)+np.pi)
                    B_lat, B_h = geodetic_iteration(B_X, B_Y, B_Z, e)
                    #print(B_lat, B_long, B_h)
                    h_values.append(B_h)
                    other.append([round(dN1[0])+i,round(dN2[0])+j,round(dN3[0])+k,round(dN4[0])+l])

                    for xdd in range(10):
                        #print(B_X, B_Y, B_Z)
                        for s in satellites:
                            s.set_rho_A(A_X, A_Y, A_Z)
                            s.set_rho_B(B_X, B_Y, B_Z)

                        A = np.array([
                            a(AB155_t1, AB154_t1, [B_X, B_Y, B_Z]),
                            a(AB159_t1, AB154_t1, [B_X, B_Y, B_Z]),
                            a(AB174_t1, AB154_t1, [B_X, B_Y, B_Z]),
                            a(AB181_t1, AB154_t1, [B_X, B_Y, B_Z]),
                            a(AB155_t2, AB154_t2, [B_X, B_Y, B_Z]),
                            a(AB159_t2, AB154_t2, [B_X, B_Y, B_Z]),
                            a(AB174_t2, AB154_t2, [B_X, B_Y, B_Z]),
                            a(AB181_t2, AB154_t2, [B_X, B_Y, B_Z])
                        ])
                        #print(A[0])
                        deltaL = np.array([
                            wavelength * phi(AB155_t1, AB154_t1) - AB155_t1.rho_B + AB154_t1.rho_B + AB155_t1.rho_A - AB154_t1.rho_A - wavelength * round(dN1[0])-i,
                            wavelength * phi(AB159_t1, AB154_t1) - AB159_t1.rho_B + AB154_t1.rho_B + AB159_t1.rho_A - AB154_t1.rho_A - wavelength * round(dN2[0])-j,
                            wavelength * phi(AB174_t1, AB154_t1) - AB174_t1.rho_B + AB154_t1.rho_B + AB174_t1.rho_A - AB154_t1.rho_A - wavelength * round(dN3[0])-k,
                            wavelength * phi(AB181_t1, AB154_t1) - AB181_t1.rho_B + AB154_t1.rho_B + AB181_t1.rho_A - AB154_t1.rho_A - wavelength * round(dN4[0])-l,
                            wavelength * phi(AB155_t2, AB154_t2) - AB155_t2.rho_B + AB154_t2.rho_B + AB155_t2.rho_A - AB154_t2.rho_A - wavelength * round(dN1[0])-i,
                            wavelength * phi(AB159_t2, AB154_t2) - AB159_t2.rho_B + AB154_t2.rho_B + AB159_t2.rho_A - AB154_t2.rho_A - wavelength * round(dN2[0])-j,
                            wavelength * phi(AB174_t2, AB154_t2) - AB174_t2.rho_B + AB154_t2.rho_B + AB174_t2.rho_A - AB154_t2.rho_A - wavelength * round(dN3[0])-k,
                            wavelength * phi(AB181_t2, AB154_t2) - AB181_t2.rho_B + AB154_t2.rho_B + AB181_t2.rho_A - AB154_t2.rho_A - wavelength * round(dN4[0])-l
                        ])
                        #print(AB154_t1.rho_B)
                        dB_X, dB_Y, dB_Z = np.linalg.inv(A.T@P@A)@A.T@P@deltaL
                        #print(dB_X, dB_Y, dB_Z)
                        B_X += float(dB_X)
                        B_Y += float(dB_Y)
                        B_Z += float(dB_Z)
                    
                    #print(A_X, B_X, A_Y, B_Y, A_Z, B_Z)
                    

                    a_=6378137
                    b_=6356752.3141
                    e = (a_**2-b_**2)/a_**2
                    B_long = np.rad2deg(np.arctan(B_Y/B_X)+np.pi)
                    B_lat, B_h = geodetic_iteration(B_X, B_Y, B_Z, e)
                    #print(B_lat, B_long, B_h)
                    h_values.append(B_h)
                    other.append([round(dN1[0])-i,round(dN2[0])-j,round(dN3[0])-k,round(dN4[0])-l])

    print(np.array(h_values) - 23.787)
    print(h_values)
    
    

if __name__=="__main__":
    main()