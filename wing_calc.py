import pandas as pd
import numpy as np
import sympy as sp

# === CONFIG ===
input_file = '/content/drive/MyDrive/intern/aoa_vs_movMassPos.csv'
output_file = '/content/drive/MyDrive/intern/wing_calc.csv'

# === DATA LOADING ===
df = pd.read_csv(input_file)

rho = 1000  # Water density in kg/m³
V = 0.5     # Glider speed in m/s
q = 0.5 * rho * V**2  # Dynamic pressure
print('q:', q)

K_D_filtered = df.iloc[:, 4]
K_D0_filtered = df.iloc[:, 5]
K_L_filtered = df.iloc[:, 6]
alpha_filtered = df.iloc[:, 2]

def wingCalc():
    K_L, K_D, K_L0, K_D0, alpha_d = sp.symbols('K_L K_D K_L0 K_D0 alpha_d')

    C_L1 = K_L / q
    C_D0 = K_D0 / q
    C_D1 = K_D / q
    LD_max = C_L1 / (2 * sp.sqrt(C_D0 * C_D1))
    FL = 1.0

    C_L1_values = []
    C_D0_values = []
    C_D1_values = []
    S_values = []
    LD_max_values = []
    CL_values = []

    for kd, kd0, kl, a in zip(K_D_filtered, K_D0_filtered, K_L_filtered, alpha_filtered):
        substitution_dict = {
            K_D: kd,
            K_D0: kd0,
            K_L: kl,
            alpha_d: a
        }

        c1 = C_L1.subs(substitution_dict)
        c2 = C_D0.subs(substitution_dict)
        c3 = C_D1.subs(substitution_dict)
        c4 = LD_max.subs(substitution_dict)
        cl_val = c1 * a
        c5 = FL / (q * cl_val) if cl_val != 0 else np.nan  # Avoid division by zero

        C_L1_values.append(float(c1))
        C_D0_values.append(float(c2))
        C_D1_values.append(float(c3))
        LD_max_values.append(float(c4))
        S_values.append(float(c5))
        CL_values.append(float(cl_val))

    return C_L1_values, C_D0_values, C_D1_values, LD_max_values, S_values, CL_values

# === CALCULATIONS ===
CL1_val, CD0_val, CD1_val, LD_maxval, S_val, CL_val = wingCalc()

sVal = np.array(S_val)

# Calculate span (b) and chord (c)
b_values = []
c_values = []
AspectRatio = 6

for s in sVal:
    if np.isnan(s) or s <= 0:
        b = np.nan
        c = np.nan
    else:
        b = float(round(np.sqrt(AspectRatio * s), 3))
        c = float(round(s / b, 3)) if b != 0 else np.nan
    b_values.append(b)
    c_values.append(c)

# === ENSURE EQUAL LENGTHS ===
min_length = min(
    len(CL_val), len(CD0_val), len(CD1_val), len(LD_maxval),
    len(sVal), len(b_values), len(c_values), len(alpha_filtered)
)

CL_val = np.array(CL_val)[:min_length]
CD0_val = np.array(CD0_val)[:min_length]
CD1_val = np.array(CD1_val)[:min_length]
LD_maxval = np.array(LD_maxval)[:min_length]
sVal = np.array(sVal)[:min_length]
b_values = np.array(b_values)[:min_length]
c_values = np.array(c_values)[:min_length]
alpha_filtered = np.array(alpha_filtered)[:min_length]

# Calculate total drag coefficient
CD = CD0_val + CD1_val * np.power(alpha_filtered, 2)

# === BUILD DATAFRAME ===
res_dict = {
    'C_L': CL_val,
    'C_D': CD,
    'alpha_d': alpha_filtered,
    'L/D Max': LD_maxval,
    'S (wing area)': sVal,
    'Span (b)': b_values,
    'Chord (c)': c_values
}
result = pd.DataFrame(res_dict)
print(result)
print('------------------------------------------------------------------------------')
result_filtered = result.dropna(subset=['Span (b)', 'Chord (c)'])
print(result_filtered)
result_filtered.to_csv(output_file, index=False)





