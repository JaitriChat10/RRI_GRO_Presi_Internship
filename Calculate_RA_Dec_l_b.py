import numpy as np
from datetime import datetime, timezone

#Function for Julian Day
def Julian_Day(year, month, day, hour, minute, sec):
    a = (14 - month) // 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    
    jdn = day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
    JD = jdn + (hour - 12) / 24 + minute / 1440 + sec / 86400
    return JD

#Function for Local Sidereal Time
def LST(east_long, dt):
    jd = Julian_Day(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    jdt = int(jd + 0.5) - 0.5

    d = jdt - 2451545
    T = d / 36525
   
    gmst = 280.46061837 + 360.98564736629*(jd - 2451545.0) + 0.000387933*(T**2) - (T**3)/38710000           #Greenwich Mean Sidereal Time
    gmst = gmst % 360

    lmst = (gmst + east_long) 
    return lmst

#Function for conversion from Horizontal to Equatorial System
def Horizontal_Equatorial(a, A, lat, lst_deg):
    a = np.radians(a)           #Altitude
    A = np.radians(A)           #Azimuth
    lat = np.radians(lat)       #Latitude
    lst = np.radians(lst_deg)   #LST

    # Declination
    sin_dec = np.sin(a) * np.sin(lat) + np.cos(a) * np.cos(lat) * np.cos(A)
    dec = np.asin(sin_dec)

    # Hour Angle
    cos_h = (np.sin(a) - np.sin(lat)*np.sin(dec)) / (np.cos(lat)*np.cos(dec))
    h = np.acos(cos_h)

    # Adjust HA depending on azimuth
    if np.sin(A) > 0:
        h = 2*np.pi - h

    # Right Ascension
    ra = lst - h

    ra = ra % (2*np.pi)

    return np.degrees(ra), np.degrees(dec)

#Function for conversion from Equatorial to Galactic Coordinate System
def Equatorial_Galactic(ra_deg, dec_deg):
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)

    # Galactic North Pole (J2000)
    ra_gp = np.radians(192.85948)
    dec_gp = np.radians(27.12825)
    l_N = np.radians(123.0)

    sin_b = (np.sin(dec)*np.sin(dec_gp)) + (np.cos(dec)*np.cos(dec_gp)*np.cos(ra - ra_gp))
    b = np.arcsin(sin_b)

    sin_ln_l = (np.cos(dec)*np.sin(ra - ra_gp))/np.cos(b)
    l = l_N - np.arcsin(sin_ln_l)
    return np.degrees(l), np.degrees(b)
    
    
#Interactive cell for entering required details
if __name__ == "__main__":
    print("Enter observation site and antenna orientation details:")

    lat = float(input("Latitude (in degrees, +N): "))
    lon = float(input("Longitude (in degrees, +E): "))
    alt = float(input("Altitude (in degrees): "))
    az = float(input("Azimuth (in degrees, 0=N, 90=E): "))

    date_str = input("Enter observation time (YYYY-MM-DD HH:MM:SS): ")
    dt = datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S")
    dt = dt.astimezone(timezone.utc)

    # Calculate LST
    lst = LST(lon, dt) 

    # Convert to RA-Dec
    ra, dec = Horizontal_Equatorial(alt, az, lat, lst)

    # Convert to Galactic
    glon, glat = Equatorial_Galactic(ra, dec)

    print("RESULTS:")
    print(f"Right Ascension (in degree): {ra:.2f}")
    print(f"Declination (in degree): {dec:.2f}")
    print(f"Galactic Longitude (in degree): {glon:.2f}")
    print(f"Galactic Latitude (in degree): {glat:.2f}")
