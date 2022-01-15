import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

nx = 21
nt = 1001

dx = 1/(nx-1)
dt = 0.5/(nt-1)
x = np.array([i * dx for i in range(nx)])

# On choisit le profil obtenu avec la méthode d'Euler implicite et les matrices creuses par exemple.

arr = np.loadtxt("QBonus2_implicite.txt")

fig, ax = plt.subplots()

def animate(i):
    ax.clear()
    ax.plot(x,arr[i])
    ax.set_title(f"Temperature distribution at t = {i * dt}")
    plt.ylim(-0.5,1.75)
    plt.xlabel("x")
    plt.ylabel("Temperature")
    return ax

ani = animation.FuncAnimation(fig,animate,frames = nt,repeat = False)
plt.show()