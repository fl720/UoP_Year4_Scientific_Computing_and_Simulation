# # Velocity Verlet algorithm
# import matplotlib.pyplot as plt
# import matplotlib.animation
# import numpy as np

# m1 = 1
# (x1, y1) = (1.0, 0.0)
# (vx1, vy1) = (0, 0.2)

# m2 = 2
# (x2, y2) = (-0.5, 0.0)
# (vx2, vy2) = (0, -0.1)

# def grav(x1, y1, x2, y2):   
#     # inverse square law between two particles mass m1, m2
#     dx = x2 - x1
#     dy = y2 - y1
#     r2 = dx * dx + dy * dy
#     r3 = r2 ** 1.5
#     fx = dx / r3
#     fy = dy / r3
#     return (fx/m1, fy/m1, -fx/m2, -fy/m2)

# t=0
# tend=50
# dt = 0.01
# NITER = int((tend - t)/dt + 0.5)

# ys = np.zeros((4, NITER + 1))
# ys [:, 0] = [x1, y1, x2, y2]

# (ax1, ay1, ax2, ay2) = grav(x1, y1, x2, y2)
# for i in range(NITER):
#     vx1 += 0.5 * dt * ax1
#     vy1 += 0.5 * dt * ay1
#     vx2 += 0.5 * dt * ax2
#     vy2 += 0.5 * dt * ay2
#     x1 += dt * vx1
#     y1 += dt * vy1
#     x2 += dt * vx2
#     y2 += dt * vy2
#     (ax1, ay1, ax2, ay2) = grav(x1, y1, x2, y2)
#     vx1 += 0.5 * dt * ax1
#     vy1 += 0.5 * dt * ay1
#     vx2 += 0.5 * dt * ax2
#     vy2 += 0.5 * dt * ay2
#     ys [:, i+1] = [x1, y1, x2, y2]

# fig = plt.figure(figsize=[6,6])
# ax = plt.axes([0., 0., 1., 1.], xlim=(-1.5,1.5), ylim=(-1.5,1.5))
# ax.set_aspect('equal')
# ax.axis('off')
# l = ax.scatter([x1, x2], [y1, y2], c=['blue', 'red'], s=[4, 8])

# OUTPUT_FREQ=10 # display every 10 steps only

# def animate(i):
#     l.set_offsets(np.reshape(ys [0:4, OUTPUT_FREQ * i], (2,2)))

# ani = matplotlib.animation.FuncAnimation(fig, animate, frames=int(NITER/OUTPUT_FREQ))

# # from IPython.display import HTML
# # HTML(ani.to_jshtml())

# plt.show()

# Velocity Verlet algorithm
import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np

m1 = 1
(x1, y1) = (1.0, 0.0)
(vx1, vy1) = (0.0, 0.2)

m2 = 2
(x2, y2) = (-0.5, 0.0)
(vx2, vy2) = (0.0, -0.1)

def grav(x1, y1, x2, y2):   
    # inverse square law between two particles mass m1, m2
    dx = x2 - x1
    dy = y2 - y1
    r2 = dx * dx + dy * dy
    r3 = r2 ** 1.5
    fx = dx / r3
    fy = dy / r3
    return (fx/m1, fy/m1, -fx/m2, -fy/m2)

t=0
tend=50
dt = 0.01
NITER = int((tend - t)/dt + 0.5)

ys = np.zeros((4, NITER + 1))
ys [:, 0] = [x1, y1, x2, y2]

(ax1, ay1, ax2, ay2) = grav(x1, y1, x2, y2)
for i in range(NITER):
    vx1 += 0.5 * dt * ax1
    vy1 += 0.5 * dt * ay1
    vx2 += 0.5 * dt * ax2
    vy2 += 0.5 * dt * ay2
    x1 += dt * vx1
    y1 += dt * vy1
    x2 += dt * vx2
    y2 += dt * vy2
    (ax1, ay1, ax2, ay2) = grav(x1, y1, x2, y2)
    vx1 += 0.5 * dt * ax1
    vy1 += 0.5 * dt * ay1
    vx2 += 0.5 * dt * ax2
    vy2 += 0.5 * dt * ay2
    ys [:, i+1] = [x1, y1, x2, y2]

fig = plt.figure(figsize=[6,6])
# 稍微调整边距，确保轨迹不被切掉
ax = plt.axes([0.1, 0.1, 0.8, 0.8], xlim=(-1.5,1.5), ylim=(-1.5,1.5))
ax.set_aspect('equal')
# ax.axis('off') # 如果需要网格辅助观察，可以注释掉这行

# 1. 创建两条线对象，用于显示历史轨迹
# alpha 设置透明度，lw 设置线宽
trail1, = ax.plot([], [], 'b-', alpha=0.4, lw=1) 
trail2, = ax.plot([], [], 'r-', alpha=0.4, lw=1)

# 保持原有的散点对象（代表当前位置的“头”）
l = ax.scatter([x1, x2], [y1, y2], c=['blue', 'red'], s=[15, 30])

OUTPUT_FREQ=10 

def animate(i):
    # 计算当前步骤的索引
    current_idx = OUTPUT_FREQ * i
    
    # 2. 更新历史轨迹：传入从开始到当前步的所有坐标
    trail1.set_data(ys[0, :current_idx], ys[1, :current_idx])
    trail2.set_data(ys[2, :current_idx], ys[3, :current_idx])
    
    # 更新当前点的位置
    l.set_offsets(np.reshape(ys[0:4, current_idx], (2,2)))
    
    return l, trail1, trail2

# blit=True 可以提高渲染效率
ani = matplotlib.animation.FuncAnimation(fig, animate, frames=int(NITER/OUTPUT_FREQ), blit=True, interval=20)

plt.show()