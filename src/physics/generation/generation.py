
import numpy as np
import pandas as pd
from scipy.integrate import quad
import matplotlib.pyplot as plt



def temperature_integral_numerical(T_start, T_end, E_a, R, num_points=100):
    """
    Calculate the temperature integral for non-isothermal kinetics.
    Integral: I(E, T) = ∫[T_start, T_end] exp(-E/(RT)) dT
    """
    integrand = lambda T, E_a, R: np.exp(-E_a / (R * T))
    integral_result = quad(integrand, T_start, T_end, args=(E_a, R))[0]
    return integral_result

def calculate_conversion_non_isothermal(S_i0, A, E_a, T_intervals, HRs):
    """
    Calculate total conversion based on provided parameters for non-isothermal kinetics.
    """
    R = 1.987e-3  # Gas constant, kcal/(mol·K)
    conversions = []
    temperatures = []
    current_fracs = S_i0
    df = pd.DataFrame(columns=['Myr', 'Temperature Interval', 'Current Conversion'])  # Initialize df here

    for j in range(len(T_intervals)):
        T_start, T_end = T_intervals[j]
        T_start, T_end = T_start + 273.15, T_end + 273.15  # Convert to Kelvin
        Hr = HRs[j]

        for i in range(len(S_i0)):
            I_E_T = temperature_integral_numerical(T_start, T_end, E_a[i], R)
            conversion_increment = np.exp(-A/Hr * I_E_T)
            current_fracs[i] *= conversion_increment

        current_conversion = 1 - current_fracs.sum()
        conversions.append(current_conversion)

        new_row = pd.DataFrame({
            'Myr': [Myr_intervals[j]],  # Ensure it's a single-element list or a Series
            'Temperature Interval': [f"{T_intervals[j][0]} - {T_intervals[j][1]}"],
            'Current Conversion': [f"{current_conversion:.4f}"]
        })

        df = pd.concat([df, new_row], ignore_index=True)  # Concatenate new row

        # print(f'Current conversion after temperature increase up to {T_intervals[j][1]}: ', current_conversion, end='\n\n')

    total_conversion = conversions[-1]
    print(f"Total conversion (\u03B1_total): {total_conversion:.4f}")

    return total_conversion, temperatures, conversions, df



#  105 Myr

A = 3.17 * 10**11 * 5171.9 # [s-1] [1e25 1/Myr in 1/s] #??????????????????????????????????????
S_i0  = [
    0.05, 0.09, 0.17, 0.27, 0.21, 0.5, 0.62, 2.11, 6.45,
    11.32, 24.04, 19.16, 18.09, 6.12, 2.92, 2.05, 0.72,
    1.06, 0.53, 0.37, 1.17, 1.98
]
S_i0  = np.array(S_i0) / 100

E_a = np.array([
    47, 48, 49, 50, 51, 52, 53, 54, 55,
    56, 57, 58, 59, 60, 61, 62, 63,
    64, 65, 66, 67, 72
])

T_intervals = [
    (0, 20),(20, 24.62),(24.62, 20.03),(20.03, 46.19),(46.19, 20.07),\
    (20.07, 64.62),(64.62, 20.12),(20.12, 82.98),(82.98, 20.18),\
    (20.18, 100.51),(100.51, 20.23),(20.23, 113.36),(113.36, 20.27),(20.27, 120.17),(120.17, 20.32),\
    (20.32, 128.3),(128.3, 20.38),(20.38, 136.41),(136.41, 20.43),(20.43, 144.47),(144.47, 20.5),(20.5, 152.48)

]  # Temperature intervals in Celsius

Myr_intervals = [105, \
                 100, 95, 90, 85, \
                 80, 75, 70, 65, 60, 55,\
                 50, 45, 40, 35,\
                 30, 25, 20, 15,\
                 10, 5, 0]


time_steps_interval = np.abs(np.diff(Myr_intervals))
time_steps_interval = np.insert(time_steps_interval, 0, 0) # for Hr

Hrs_interval = []
for i in range(len(T_intervals)):
    T_start, T_end = T_intervals[i]
    time_step = time_steps_interval[i]
    Hrs_interval.append((T_end - T_start)/(time_step * 3.15 * 1e13))
    print(T_start, T_end,  time_step)
    print(Hrs_interval)
total_conversion, temperatures, conversions, df = calculate_conversion_non_isothermal(S_i0, A, E_a, T_intervals, np.array(Hrs_interval))

plt.figure(figsize=(8, 6))
plt.plot(Myr_intervals, conversions, marker='o', label="Cumulative Conversion")
plt.xlabel("Time")
plt.ylabel("Cumulative Conversion")
plt.title("Cumulative Conversion vs Temperature for Non-Isothermal Kinetics, 105 Myr")
plt.legend()
plt.grid(True)
plt.show()