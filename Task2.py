###### 1- Unit Test to check accuracy of force function ######


#Gradient assumed to be reference correct force value, use for proportion of error
def test_force_error(atoms_obj, epsilon = 0.05, test = True):
    spacing = 1.8
    grad_force = compute_force(atoms_obj, epsilon, spacing)
    ase_force = atoms_obj.get_forces()[0,2]  
    #Position after .get_forces() represents: (n-1)th atom,  force in (n-1)th axis (in x,y,z) 
    force_error = abs((ase_force-grad_force)/grad_force) 
    
    #If executed as validation, assert, if not return error 
    if test == True:
        assert force_error <= 0.05, "Error higher than 5%, incorrect force computed."
    elif test == False:
        return force_error

'''Ideally, the threshold for error would be tabulated and compared against a reference table 
for varying epsilon values, but for now we can neglect this complication as we intend to investigate
this anyways.'''

#if __name__ == '__main__':
#   test_force_error(atoms, 0.01)

###### 2- Force plots: manual potential gradient against output from ASE method #####

atoms = Atoms(composition, positions = [(0.0,0.0,0.0), (0.0,0.0, 0.0)]) #Initially at 0 separation
atoms.set_calculator(calc)

epsilon = 0.05

def compute_force(atoms_obj, epsilon, spacing):
    #Small change in potential energy
    atoms_obj.positions[1,2] = spacing
    P0 = atoms_obj.get_potential_energy() 

    atoms_obj.positions[1,2] = spacing+epsilon
    P1 = atoms_obj.get_potential_energy() 
    
    #Small perturbation in distance - epsilon
    return (P1-P0)/epsilon

linearised_forces = [compute_force(atoms, epsilon, sep) for sep in x_dist]

plt.figure(1)
ax1 = plt.subplot(2,1,1)

ax1.plot(x_dist, linearised_forces, color = "green")
ax1.plot(x_dist, forces, color = "red")
ax1.set_title("Force against separation curve(atomic units)")
plt.show()


#Checking error for varying epsilon (at an atomic spacing of 1.8 Angstrom)
eps_variation = np.linspace(0.1, 0.6, 20)
error = [100*test_force_error(atoms,eps,False) for eps in eps_variation]

ax2 = plt.subplot(2,1,2)
ax2.plot(eps_variation, error)
ax2.set_title("Epsilon against % error")
plt.show()