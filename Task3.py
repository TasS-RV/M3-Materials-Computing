def compute_pressure_pe(strain):
    '''
    This function regenerates the unit cell of copper, measuring PE and 
    pressures for varying Volumetric strain states.
    '''
    #2,2,2 scaling not applied - 'specific' volume used, assuming bulk properties cascade down to unit cell level

    lat_constant = 3.6
    cu = bulk("Cu", "fcc", a=lat_constant, cubic=True)
    cu.cell
    cu.set_calculator(calc) #Set a calculator for atomistic calculations

    sf = 1 + strain #Volumetric strain applied
    cu.set_cell(cu.get_cell()*sf) #Scales cell 
    
    pe_per_atom = cu.get_potential_energy()/cu.get_global_number_of_atoms()
    pressure = (-np.trace(cu.get_stress(voigt = False))/3)/ GPa #Requires conversion to Pa units
    
    volume_per_atom = (lat_constant*sf)**3*Ang
    
    return volume_per_atom, pressure, pe_per_atom #Following values required for each specific volume

#Trying between sensible range of +- 10% axial strain
strains = np.linspace(-0.1, 0.1, 20)

#Appending to individual lists for plotting
vol_atomic = [compute_pressure_pe(strain)[0]/4 for strain in strains]
potential_energies = [compute_pressure_pe(strain)[2] for strain in strains]
pressure = [compute_pressure_pe(strain)[1] for strain in strains]


plt.figure(1)
ax_1 = plt.subplot(2,1,1)
ax_1.plot(vol_atomic, potential_energies)
ax_1.set_title('Potential Energy/ Atom against Specific Volume (A^3)')
plt.show()

ax_2 = plt.subplot(2,1,2)
ax_2.plot(vol_atomic, pressure)
ax_2.set_title('Pressure against Specific Volume (A^3)')
plt.show()

'''
In the formula -VdP/dV, all volumes and changes in volumes are specific (per atom)
therefore can continue using values from arrays in the graphs as /Atom part will cancel.


The calculated value is 91.15 GPa - sufficiently close to experimentally obtained value
'''

vol_lowest_pe = vol_atomic[potential_energies.index(min(potential_energies))]
stable_sep = (vol_lowest_pe*4)**(1/3)
print("Corresponding lattice constant to equilibrium volume: {0:.3f} Angstrom".format(stable_sep))

#Obtain pressure and volume differences wtih infinitesimal linearisation:
lin_strain = 0.01

V0, P0, PE0 = compute_pressure_pe(0)
V1, P1, PE1 = compute_pressure_pe(lin_strain)

#Using given formula for Bulk Modulus:
K = -1*V0*((P1-P0)/(V1-V0))
print("Bulk Modulus calculated: {0:.3f} GPa".format(K))
