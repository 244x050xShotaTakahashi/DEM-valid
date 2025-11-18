import numpy as np
data = np.load("./E_field/data/electro_field_2D.npz")
x, y, Ex, Ey = data["x"], data["y"], data["Ex"], data["Ey"]

# 出力: 1行に x, y, Ex, Ey
with open("./E_field/data/electro_field_2D.dat", "w") as f:
    f.write(f"{len(x)} {len(y)}\n")
    for i in range(len(y)):
        for j in range(len(x)):
            f.write(f"{x[j]:.6e} {y[i]:.6e} {Ex[i,j]:.6e} {Ey[i,j]:.6e}\n")
print("Saved electro_field_2D.dat")
