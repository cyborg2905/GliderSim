import numpy as np
import sympy as sp
import pandas as pd

def alphaFun():
    K_L, K_D, K_L0, K_D0, alpha_d, zeta_d = sp.symbols('K_L K_D K_L0 K_D0 alpha_d zeta_d')

    left_hand_side = alpha_d**2 + (K_L/K_D)*sp.tan(sp.rad(zeta_d))*alpha_d + (K_D0/K_D)

    K_D_values = np.full(500, 18)
    K_D0_values = np.full(500, 109)
    K_L0_values = np.full(500, 0)
    K_L_values = np.full(500, 306)
    zeta_d_values = np.linspace(-45, 45, 500)

    alpha_values = []

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

        alpha_values.append(numeric_solutions)

    flat_list = [item for sublist in alpha_values for item in sublist]
    print(f"Found {len(flat_list)} alpha_d values.")
    return flat_list


def MovableMassPos(alpha_d_values, N):
    g = 9.81  # m/s^2
    r_p1d, r_p3d, theta_d, m_bar, m_f3, m_f1, v1_d, v3_d, K_M0, K_M, V_d, alpha_d = sp.symbols(
        'r_p1d r_p3d theta_d m_bar m_f3 m_f1 v1_d v3_d K_M0 K_M V_d alpha_d'
    )

    right_hand_side = -r_p3d*sp.tan(theta_d) + (
        (m_f3 - m_f1)*v1_d*v3_d + (K_M0 + K_M*alpha_d)*pow(V_d, 2)
    ) / (m_bar*g*sp.cos(theta_d))

    # Input arrays
    r_p3d_values = np.full(N, 0.04)
    theta_d_values = np.full(N, 0)
    m_bar_values = np.full(N, 2)
    m_f3_values = np.full(N, 14)
    m_f1_values = np.full(N, 2)
    V_d_values = np.linspace(0.1, 0.37, N)
    K_M0_values = np.full(N, 0)
    K_M_values = np.full(N, -36.5)

    # Velocity components
    alpha_d_rad = np.radians(alpha_d_values)
    v1_d_values = V_d_values * np.cos(alpha_d_rad)
    v3_d_values = V_d_values * np.sin(alpha_d_rad)

    position_values = []

    for rp3d, td, mbar, mf3, mf1, v1d, v3d, KM0, KM, Vd, ad in zip(
        r_p3d_values, theta_d_values, m_bar_values, m_f3_values, m_f1_values,
        v1_d_values, v3_d_values, K_M0_values, K_M_values, V_d_values, alpha_d_values
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
alpha_d_values = alphaFun()
N = len(alpha_d_values)
positions = MovableMassPos(alpha_d_values, N)

# DataFrame creation
d = {'Angle of Attack': alpha_d_values, 'Moving Mass Positions': positions}
data = pd.DataFrame(d)
data.to_csv('C:/Users/mukul/Downloads/alpha_vs_movable_mass.csv', index=False)



print(data)
