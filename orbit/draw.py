import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

data = np.loadtxt('orbit_diff.txt')

t = data[:, 0]      # Время (часы)
dr = data[:, 1]     # Отклонение координат (км)
dv = data[:, 2]     # Отклонение скоростей (м/с)

# Моделируемая траектория MESSENGER относительно Меркурия (км)
mx, my, mz = data[:, 3], data[:, 4], data[:, 5]

# Референсная траектория MESSENGER из SPICE относительно Меркурия (км)
rx, ry, rz = data[:, 9], data[:, 10], data[:, 11]

mask = t < 60
t = t[mask]
dr = dr[mask]
dv = dv[mask]
mx, my, mz = mx[mask], my[mask], mz[mask]
rx, ry, rz = rx[mask], ry[mask], rz[mask]


fig = plt.figure(figsize=(16, 7))

# левый график
ax1 = fig.add_subplot(121, projection='3d')

ax1.plot(mx, my, mz, color='tab:red', lw=2.5, label='Модель EIH')
ax1.plot(rx, ry, rz, color='tab:blue', linestyle='--', lw=1.5, label='Референс SPICE (DE442)')

ax1.scatter(0, 0, 0, color='darkorange', s=200, marker='o', edgecolors='black', label='Меркурий')

ax1.set_xlabel('X (км)', fontweight='bold')
ax1.set_ylabel('Y (км)', fontweight='bold')
ax1.set_zlabel('Z (км)', fontweight='bold')
ax1.set_title('Траектория MESSENGER относительно Меркурия', fontsize=12, fontweight='bold')
ax1.grid(True)
ax1.legend()

# Делаем равный масштаб осей для корректного отображения 3D-сферы орбиты
max_range = np.array([mx.max()-mx.min(), my.max()-my.min(), mz.max()-mz.min()]).max() / 2.0
mid_x = (mx.max()+mx.min()) * 0.5
mid_y = (my.max()+my.min()) * 0.5
mid_z = (mz.max()+mz.min()) * 0.5
ax1.set_xlim(mid_x - max_range, mid_x + max_range)
ax1.set_ylim(mid_y - max_range, mid_y + max_range)
ax1.set_zlim(mid_z - max_range, mid_z + max_range)


# Правый график
ax2 = fig.add_subplot(122)

ax2.set_xlabel('Время (часы)', fontweight='bold')
ax2.set_ylabel('Отклонение |dr| (км)', color='tab:red', fontweight='bold')
ax2.plot(t, dr, color='tab:red', lw=2, label='|dr|')
ax2.tick_params(axis='y', labelcolor='tab:red')
ax2.grid(True, linestyle=':')

ax3 = ax2.twinx() 
ax3.set_ylabel('Отклонение |dv| (м/с)', color='tab:blue', fontweight='bold')
ax3.plot(t, dv, color='tab:blue', alpha=0.7, lw=1.5, label='|dv|')
ax3.tick_params(axis='y', labelcolor='tab:blue')

ax2.set_title('Накопление расхождения траекторий со временем', fontsize=12, fontweight='bold')

plt.suptitle('Сравнение моделей движения MESSENGER (EIH Релятивизм vs SPICE)', fontsize=14, fontweight='bold', y=0.98)
fig.tight_layout()


plt.savefig('result_plot.png', dpi=150)
plt.show()