import matplotlib.pyplot as plt
import numpy as np
separation = 4.5 #Atomic lengthscales of separation - needs Angstrom factor
composition = '2Cu'  

x_dist = (np.linspace(1.75,separation, 50))*Ang #Generates a list of positions: #Starts at 0+delta(d), atomic units of separatiom

atoms = Atoms(composition, positions = [(0.0,0.0,0.0), (0.0,0.0, 0.00001)]) #Initially at 0 separation

atoms.set_calculator(calc)
pot_energies = []
forces = []

for sep in x_dist:  #Iterates over regular increments of positions to to compute varying potential energy
    atoms.positions[1,2] = sep
    pot_energies.append(atoms.get_potential_energy())
    forces.append(atoms.get_forces()[0,2]) #Forces plotted in axis parallel to separation axis


ax1 = plt.subplot(111) #Generates a plotting environment to display the data 
ax2 = plt.subplot(111)
ax1.plot(x_dist, pot_energies, label = "Potential energy curve over distance", color = "blue")
ax2.plot(x_dist, forces, label = "Force curve over distance", color = "red")

plt.legend()
plt.show()

#Finding the  larest magntiude of potential energy 

print("From using different values for separation, we we found the repulsive potential energy dominated at separatons \n", end = "")
print("around 1.5 Atomic units. Therefore, extrapolating the lowest value and it's corresponding index gives separation at lowest PE\n", end = "")

position, value = pot_energies.index(min(pot_energies)),min(pot_energies)
print("Minimum potential energy of {0:.3f} eV/A^2, at {1:.3f} atomic units separation".format(value, x_dist[position]))