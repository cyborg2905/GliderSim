import numpy as np
import sympy as sp
import pandas as pd

filepath = 'corrected_glider_results.csv'

def alphaFun():
    K_L, K_D, K_L0, K_D0, alpha_d, zeta_d = sp.symbols('K_L K_D K_L0 K_D0 alpha_d zeta_d')

    left_hand_side = alpha_d**2 + (K_L/K_D)*sp.tan(sp.rad(zeta_d))*alpha_d + (K_D0/K_D)

    input = pd.read_csv(filepath)
    K_L_values = input.iloc[:,1].to_numpy()
    K_D0_values = input.iloc[:,2].to_numpy()
    K_D_values = input.iloc[:,3].to_numpy()
    n = len(K_D_values)
    K_L0_values = np.full(n, 0)
    zeta_d_values = np.linspace(-45, 45, n)

    alpha_values = []
    glide_values = []
    k_d_values = []
    k_d0_values = []
    k_l_values = []

    for kd, kd0, kl0, kl, zd in zip(K_D_values, K_D0_values, K_L0_values, K_L_values, zeta_d_values):
        substitution_dict = {
            K_L: kl,
            K_D: kd,
            K_L0: kl0,
            K_D0: kd0,
            zeta_d: zd
        }

        final_lhs = left_hand_side.subs(substitution_dict)
        solutions = sp.solve(final_lhs, alpha_d)

        numeric_solutions = []
        for sol in solutions:
            numeric_sol = sol.evalf()
            if numeric_sol.is_real and abs(float(numeric_sol)) > 1:
                numeric_solutions.append(float(numeric_sol))

        # Only append glide angle and alpha values if we found valid solutions
        if numeric_solutions:  # If the list is not empty
            for alpha_val in numeric_solutions:
                glide_values.append(zd)
                k_d_values.append(kd)
                k_d0_values.append(kd0)
                k_l_values.append(kl)
  # Append glide angle for each valid alpha
            alpha_values.extend(numeric_solutions)  # Add all valid solutions to flat list

    glide_angle = np.array(glide_values)

    print(f"Found {len(alpha_values)} alpha_d values.")
    return alpha_values, glide_angle, k_d_values, k_d0_values, k_l_values


def MovableMassPos(alpha_d_values, glide_angle, N):
    g = 9.81  # m/s^2
    r_p1d, r_p3d, theta_d, m_bar, m_f3, m_f1, v1_d, v3_d, K_M0, K_M, V_d, alpha_d = sp.symbols(
        'r_p1d r_p3d theta_d m_bar m_f3 m_f1 v1_d v3_d K_M0 K_M V_d alpha_d'
    )

    right_hand_side = -r_p3d * sp.tan(theta_d) + (
        (m_f3 - m_f1) * v1_d * v3_d + (K_M0 + K_M * alpha_d) * pow(V_d, 2)
    ) / (m_bar * g * sp.cos(theta_d))

    r_p3d_values = np.full(N, 0.04)
    theta_d_values = glide_angle[:N] - alpha_d_values[:N]
    m_bar_values = np.full(N, 2)
    m_f3_values = np.full(N, 14)
    m_f1_values = np.full(N, 2)
    V_d_values = np.linspace(0.1, 0.37, N)
    K_M0_values = np.full(N, 0)
    K_M_values = np.full(N, -36.5)

    alpha_d_rad = np.radians(alpha_d_values[:N])
    v1_d_values = V_d_values * np.cos(alpha_d_rad)
    v3_d_values = V_d_values * np.sin(alpha_d_rad)

    position_values = []

    for rp3d, td, mbar, mf3, mf1, v1d, v3d, KM0, KM, Vd, ad in zip(
        r_p3d_values, theta_d_values, m_bar_values, m_f3_values, m_f1_values,
        v1_d_values, v3_d_values, K_M0_values, K_M_values, V_d_values, alpha_d_values[:N]
    ):
        substitution_dict = {
            r_p3d: rp3d,
            theta_d: td,
            m_bar: mbar,
            m_f3: mf3,
            m_f1: mf1,
            v1_d: v1d,
            v3_d: v3d,
            K_M0: KM0,
            K_M: KM,
            V_d: Vd,
            alpha_d: ad
        }

        final_rhs = right_hand_side.subs(substitution_dict).evalf()
        position_values.append(float(final_rhs))

    print("Position values computed.")
    return position_values


# === Run ===
alpha_values,glide_angle,K_d,K_d0,K_l = alphaFun()
print(len(glide_angle))
N = len(alpha_values)
positions = MovableMassPos(alpha_values, glide_angle, N)

# DataFrame creation
d = {'Glide Angle':glide_angle,'Angle of Attack': alpha_values, 'Moving Mass Positions': positions,'K_D':K_d,'K_D0':K_d0,'K_L':K_l}
data = pd.DataFrame(d)
data.to_csv('aoa_vs_movMassPos.csv')
print(data)

