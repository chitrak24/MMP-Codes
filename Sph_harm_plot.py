import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm

# Grids of polar and azimuthal angles
theta = np.linspace(0, np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)

# Azimithul and Magnetic Quantum numbers
l = int(input("Azimithul quantum number: "))

if l < 0:
    print(" l can't be a negative number")

m = int(input("Magnetic quantum number: "))

if m > l or m < (-l):
    print("quantum no. m belongs to [-l,l]")


# Create a 2-D meshgrid of (theta, phi) angles.
theta, phi = np.meshgrid(theta, phi)

# Calculate the Cartesian coordinates of each point in the mesh.
xyz = np.array([np.sin(theta) * np.sin(phi),
                np.sin(theta) * np.cos(phi),
                np.cos(theta)])

def plot_Y(ax, l, m):
    
    Y = sph_harm(abs(m), l, phi, theta)

    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real
    Yx, Yy, Yz = abs(Y) * xyz

    # Colour the plotted surface a
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('seismic'))
    cmap.set_clim(-0.5, 0.5)

    ax.plot_surface(Yx, Yy, Yz,facecolors=cmap.to_rgba(Y.real),rstride=2, cstride=2)

    # Set the Axes limits and title
    ax.set_title("Spherical Harmonics $\quad$" r'$Y_{{{},{}}}$'.format(l, m))
    ax_lim = 0.5
    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)

    ax.set_xlabel('$X \longrightarrow$', fontsize='12', fontweight='normal')
    ax.set_ylabel('$Y \longrightarrow$', fontsize='12', fontweight='normal')
    ax.set_zlabel('$Z$'                , fontsize='12', fontweight='normal')


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plot_Y(ax, l, m)
plt.savefig("plot_sph_harm.pdf")
plt.show()
