######  Task 4 Part 2  - incorrect method ######
'''
from ase.build import bulk
lat_constant = 3.6 #Need to adjust lattice constant to obtain equilbrium state for unit cell (lowest pot. energy)
cu = bulk("Cu", "fcc", a=lat_constant, cubic=True)
cu.cell
cu.set_calculator(calc)

cu_unit = cu.copy() #Preserves a copy of the unit cell

#print(cu.get_stress(voigt = False)) #Checking unstrained stresses

cu.set_cell(cu.get_cell()*(4,4,4))
cu_mat = np.array(cu.get_cell())

xx_strain = 0.06
#Compressive strain applied in xx axis (normal)
cu_mat[0][0] = (1-xx_strain)*cu_mat[0][0]  
cu.set_cell(cu_mat, scale_atoms=True) #ust scale atoms for correct stress values

#print(cu_mat) #Uncomment out to check atomic positions

print(cu.get_stress(voigt = False)) #For a constrained unit cube - would expect applied, AND reaction stresses to equalise
E_Cu = 0.118 #Units GPa - source https://www.engineeringtoolbox.com/young-modulus-d_773.html
E_Cu_convert = E_Cu*GPa

print(f"Stress matrix: (after imposed {xx_strain} strain).")
stress_yy = cu.get_stress(voigt = False)[1][1]
stress_zz = cu.get_stress(voigt = False)[2][2]
stress_xx = cu.get_stress(voigt = False)[0][0]

#-1*strain value used: as compressive strain
#Using the following formula, niu made subject: strain_xx = 1/E(sigma_xx-niu(sigma_yy+sigma_zz))
print("Poisson ratio with lattice constant: {}, applied uniaxial {} compressive strain.".format(lat_constant, xx_strain))
print(((-1*xx_strain)*E_Cu_convert - stress_xx)/(-1*(stress_yy+stress_zz)))
'''



######  Task 4 Part 2  - corrected method with Iteration ######


from ase.build import bulk
lat_constant = stable_sep #Need to adjust lattice constant to obtain equilbrium state for unit cell (lowest pot. energy)
cu = bulk("Cu", "fcc", a=lat_constant, cubic=True)
cu.cell
cu.set_calculator(calc)

cu_unit = cu.copy() #Preserves a copy of the unit cell

print("Unstrained stress state: \n {}".format(cu.get_stress(voigt = False))) #Checking unstrained stresses

#cu.set_cell(cu.get_cell()*(4,4,4))
cu_mat = np.array(cu.get_cell())


xx_strain = 0.2
#Compressive strain applied in xx axis (normal)
cu_mat[0][0] = (1-xx_strain)*cu_mat[0][0] 
cu.set_cell(cu_mat, scale_atoms = True)

'''
-Iteratively apply tensile strain in yy and zz - until 'reaction' strain reduces to 0
-Given Poisson ratio < 1, only need to iterate maximum strain applies in xx.
'''

#Copy of single strained dimensions to apply strains (and restore)
cu_single_strain = cu.copy()
cu_ss_mat = np.array(cu_single_strain.get_cell())

zz_stress = []
yy_stress = []
xx_stress = []

iterations = 1000
strains = np.linspace(0, xx_strain, iterations)

for strain in strains:

    cu_mat[1][1] = (1+strain)*cu_ss_mat[1][1] #Stress appied in yy    
    cu_mat[2][2] = (1+strain)*cu_ss_mat[2][2] #Stress appied in zz

    cu.set_cell(cu_mat, scale_atoms = True) 
    
    #Normal stresses:
    xx_stress.append(cu.get_stress(voigt = False)[0][0])
    yy_stress.append(cu.get_stress(voigt = False)[1][1])
    zz_stress.append(cu.get_stress(voigt = False)[2][2]) 
    
    #print(cu.get_stress(voigt=False)) #Check appropriateness of stress matrix

#Generates plots to visually inspect stress changes
plt.figure(1)

ax1 = plt.subplot(3,1,1)
ax1.plot(strains, np.linspace(0, 0, iterations), color = 'red')
ax1.plot(strains, yy_stress, label = "YY stress against varying strain.")
plt.legend()

ax2 = plt.subplot(3,1,2)
ax2.plot(strains, zz_stress, label = "ZZ stress against varying strain.")
plt.legend()

ax3 = plt.subplot(3,1,3)
ax3.plot(strains, xx_stress, label = "XX stress against varying strain.")
plt.legend()

plt.show()


for n, s in enumerate(yy_stress):
    if s > 0:
        ps_strain = strains[n]
        break; 

print("Poisson ratio where yy and zz strains equate to 0: {}".format(ps_strain/xx_strain))


