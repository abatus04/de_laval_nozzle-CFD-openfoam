import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# ─── Load CFD data ────────────────────────────────────────────────────────────
df = pd.read_csv('/home/shoyo/nozzleCFD/centerline_t0020.csv')

x_cfd  = df['Points:0'].values
Ma_cfd = df['Ma'].values
p_cfd  = df['p'].values
T_cfd  = df['T'].values

# ─── Nozzle geometry from exact blockMeshDict spline points ───────────────────
x_spline = np.array([  0,  20,  40,  60,  80, 100,
                      110, 120, 130,
                      150, 170, 190, 210, 230, 250, 280])

r_spline = np.array([ 50, 47.2, 42.8, 37.2, 31.3, 26.0,
                      25.2, 25.0, 25.0,
                      26.8, 30.2, 34.5, 38.8, 43.0, 46.8, 50.0])

radius_interp = interp1d(x_spline, r_spline, kind='cubic')

def nozzle_radius(x):
    x_mm = np.clip(x * 1000, 0, 280)
    return radius_interp(x_mm) / 1000

# Throat area
r_throat = 0.025
A_throat = np.pi * r_throat**2

# Area ratio at each CFD point
A_ratio = np.array([
    (np.pi * nozzle_radius(x)**2) / A_throat
    for x in x_cfd
])

# ─── Isentropic theory ────────────────────────────────────────────────────────
gamma = 1.4
p0    = 500000.0
T0    = 600.0

def area_ratio_from_mach(M, g=1.4):
    return (1.0/M) * ((2/(g+1)) * (1 + (g-1)/2 * M**2)) ** ((g+1)/(2*(g-1)))

def mach_from_area_ratio(AR, supersonic=False):
    if AR < 1.0:
        AR = 1.0
    if supersonic:
        M_lo, M_hi = 1.0, 10.0
    else:
        M_lo, M_hi = 0.001, 0.9999
    for _ in range(200):
        M_mid = (M_lo + M_hi) / 2
        f = area_ratio_from_mach(M_mid) - AR
        if supersonic:
            if f > 0:
                M_hi = M_mid
            else:
                M_lo = M_mid
        else:
            if f < 0:
                M_hi = M_mid
            else:
                M_lo = M_mid
    return (M_lo + M_hi) / 2

x_throat = 0.130

Ma_theory = np.array([
    mach_from_area_ratio(ar, supersonic=(x >= x_throat))
    for x, ar in zip(x_cfd, A_ratio)
])

p_theory = p0 * (1 + (gamma-1)/2 * Ma_theory**2) ** (-gamma/(gamma-1))
T_theory = T0 * (1 + (gamma-1)/2 * Ma_theory**2) ** (-1)

# ─── Plots ────────────────────────────────────────────────────────────────────
x_mm = x_cfd * 1000

fig, axes = plt.subplots(3, 1, figsize=(12, 10))
fig.suptitle('De Laval Nozzle — CFD vs Isentropic Theory', fontsize=14, fontweight='bold')

# Mach
axes[0].plot(x_mm, Ma_cfd,    'b-',  linewidth=2,   label='OpenFOAM CFD')
axes[0].plot(x_mm, Ma_theory, 'r--', linewidth=1.5, label='Isentropic Theory')
axes[0].axvline(x=130, color='gray',  linestyle=':',  label='Throat (x=130mm)')
axes[0].axhline(y=1.0, color='green', linestyle=':', alpha=0.5, label='Mach = 1')
axes[0].set_ylabel('Mach Number')
axes[0].set_title('Mach Number Distribution')
axes[0].legend()
axes[0].grid(True, alpha=0.3)
axes[0].set_ylim(0, 3)

# Pressure
axes[1].plot(x_mm, p_cfd/1000,    'b-',  linewidth=2,   label='OpenFOAM CFD')
axes[1].plot(x_mm, p_theory/1000, 'r--', linewidth=1.5, label='Isentropic Theory')
axes[1].axvline(x=130, color='gray', linestyle=':')
axes[1].set_ylabel('Pressure (kPa)')
axes[1].set_title('Static Pressure Distribution')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

# Temperature
axes[2].plot(x_mm, T_cfd,    'b-',  linewidth=2,   label='OpenFOAM CFD')
axes[2].plot(x_mm, T_theory, 'r--', linewidth=1.5, label='Isentropic Theory')
axes[2].axvline(x=130, color='gray', linestyle=':')
axes[2].set_xlabel('Axial Position (mm)')
axes[2].set_ylabel('Temperature (K)')
axes[2].set_title('Static Temperature Distribution')
axes[2].legend()
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/shoyo/nozzleCFD/results/nozzle_validation.png', dpi=150, bbox_inches='tight')
plt.show()

# ─── Error metrics ────────────────────────────────────────────────────────────
mask = x_cfd < x_throat
mae_mach = np.mean(np.abs(Ma_cfd[mask] - Ma_theory[mask]))
print(f"\n── Validation Metrics ──────────────────────")
print(f"Converging section MAE (Mach):  {mae_mach:.4f}")
print(f"Max CFD Mach:                   {Ma_cfd.max():.3f}")
print(f"Max Theory Mach (exit):         {Ma_theory[-1]:.3f}")
throat_idx = np.argmin(np.abs(x_cfd - x_throat))
print(f"Throat Mach (CFD):              {Ma_cfd[throat_idx]:.3f}")
print(f"Throat Mach (Theory):           {Ma_theory[throat_idx]:.3f}")
print(f"────────────────────────────────────────────")
