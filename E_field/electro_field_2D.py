# Retry: faster solver using sparse linear system (scipy.sparse) for Laplace with Dirichlet electrodes.
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
import os

# --- Parameters (smaller grid for fast solve) ---
Lx = 0.02
Ly = 0.02
nx = 101
ny = 101
V_plate = 0.0
V_line = 200.0
line_count = 4
line_x_margin = 0.002
line_spacing = (Lx - 2*line_x_margin) / (line_count-1)
line_radius = 0.00025
line_y = 0.018

x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
dx = x[1]-x[0]
dy = y[1]-y[0]

phi = np.zeros((ny, nx), dtype=float)
fixed = np.zeros_like(phi, dtype=bool)

# Bottom plate Dirichlet
phi[0, :] = V_plate
fixed[0, :] = True

# Line electrodes
line_centers_x = [line_x_margin + k*line_spacing for k in range(line_count)]
line_half_cells = max(1, int(np.ceil(line_radius / dx)))
i_line = int(round((line_y - y[0]) / dy))
for xc in line_centers_x:
    j_center = int(round((xc - x[0]) / dx))
    j0 = max(0, j_center - line_half_cells)
    j1 = min(nx-1, j_center + line_half_cells)
    phi[i_line, j0:j1+1] = V_line
    fixed[i_line, j0:j1+1] = True

# Build sparse matrix A and RHS b for Laplace: A * phi_vec = b
N = nx * ny
A = lil_matrix((N, N), dtype=float)
b = np.zeros(N, dtype=float)

def idx(i, j):
    return i * nx + j

for i in range(ny):
    for j in range(nx):
        k = idx(i,j)
        if fixed[i,j]:
            A[k,k] = 1.0
            b[k] = phi[i,j]
        else:
            # interior or boundary with Neumann treatment on outer boundaries (copy neighbor)
            # Use 5-point Laplacian: (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1] - 4phi)/dx^2 = 0
            # But if neighbor is outside domain, approximate Neumann by using same value (zero normal) -> mirror
            coef = 0.0
            # up (i+1)
            if i+1 < ny:
                A[k, idx(i+1,j)] = 1.0
            else:
                A[k, idx(i,j)] += 1.0  # mirror
            # down (i-1)
            if i-1 >= 0:
                A[k, idx(i-1,j)] = 1.0
            else:
                A[k, idx(i,j)] += 1.0
            # right (j+1)
            if j+1 < nx:
                A[k, idx(i,j+1)] = 1.0
            else:
                A[k, idx(i,j)] += 1.0
            # left (j-1)
            if j-1 >= 0:
                A[k, idx(i,j-1)] = 1.0
            else:
                A[k, idx(i,j)] += 1.0
            A[k,k] = -4.0 + A[k,k]  # adjust diagonal if mirrored
            b[k] = 0.0

A = csc_matrix(A)
# Solve
phi_vec = spsolve(A, b)
phi = phi_vec.reshape((ny, nx))

# Compute E-field
Ey, Ex = np.gradient(phi, dy, dx)
Ex = -Ex
Ey = -Ey

# interpolation functions
interp_Ex = interpolate.RegularGridInterpolator((y, x), Ex)
interp_Ey = interpolate.RegularGridInterpolator((y, x), Ey)

def electric_field_at(xp, yp):
    pt = np.array([yp, xp])
    Exv = interp_Ex(pt)
    Eyv = interp_Ey(pt)
    return float(Exv), float(Eyv)

def coulomb_force(q_charge, xp, yp):
    Exv, Eyv = electric_field_at(xp, yp)
    return q_charge * Exv, q_charge * Eyv

# plot
plt.figure(figsize=(7,6))
X, Y = np.meshgrid(x, y)
levels = np.linspace(np.min(phi), np.max(phi), 50)
plt.contourf(X, Y, phi, levels=levels)
plt.colorbar(label='Potential [V]')
step = max(1, int(nx/30))
plt.quiver(X[::step, ::step], Y[::step, ::step], Ex[::step, ::step], Ey[::step, ::step], scale=3e4)
plt.scatter(line_centers_x, [line_y]*len(line_centers_x), marker='s', s=40, label='Line electrodes')
plt.hlines(0, 0, Lx, colors='k', linewidth=2, label='Bottom plate (V=0)')
plt.xlabel('x [m]'); plt.ylabel('y [m]')
plt.title('Potential & E-field (sparse Laplace solver)')
plt.legend(loc='upper right', fontsize='small')
plt.tight_layout()
plt.show()

# Show sample field/force values
sample_positions = [(0.005, 0.004), (0.01, 0.01), (0.015, 0.012)]
q_test = -1e-12
for xp, yp in sample_positions:
    Exv, Eyv = electric_field_at(xp, yp)
    Fx, Fy = coulomb_force(q_test, xp, yp)
    print(f"pos ({xp:.4f},{yp:.4f})  E=({Exv:.2e},{Eyv:.2e}) V/m  F=({Fx:.2e},{Fy:.2e}) N for q={q_test} C")

# save for DEM
os.makedirs('data', exist_ok=True)
np.savez('data/electro_field_2D.npz', x=x, y=y, phi=phi, Ex=Ex, Ey=Ey)
print("Saved data/electro_field_2D.npz")
