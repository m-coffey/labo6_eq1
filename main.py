import os
import numpy as np
import matplotlib.pyplot as plt

naca_profil = '0012'
n_panel_nodes = 160
te_le_density_ratio = .2
p = .4
n_iterations = 2000


alpha_start = 0
alpha_end = 17
alpha_step = 1
alphas = str(alpha_start) + ' ' + str(alpha_end) + ' ' + str(alpha_step)

# Angles d'attaque et pas de discretisation qui serviront a determiner le pas de discretisation a utiliser plus tard
alpha1 = 0
alpha2 = 15
pas_sensibilite_alpha = .05
n_pas_sensibilite_alpha = 30

# Nombre de Reynolds utilise pour les calculs, ainsi que le pas de discretisation
re = 6e6
n_re_sensibilite = 15
pas_sensibilite_re = 50000
re_sensibilite = np.arange(re-n_re_sensibilite*pas_sensibilite_re, re+n_re_sensibilite*pas_sensibilite_re, pas_sensibilite_re)


# Creation des noms pour les fichier d'input et de resultats pour determiner le pas de discretisation
pas_sensibilite_results_alpha1 = 'pas_sensibilite_results_alpha1.txt'
pas_sensibilite_input_alpha1 = 'pas_sensibilite_input_alpha1.txt'
pas_sensibilite_results_alpha2 = 'pas_sensibilite_results_alpha2.txt'
pas_sensibilite_input_alpha2 = 'pas_sensibilite_input_alpha2.txt'


#for i = 1 : n_pas_sensibilite_alpha:

# Si les fichiers de resultats existent deja il faut les effacer
if os.path.exists(pas_sensibilite_results_alpha1):
    os.remove(pas_sensibilite_results_alpha1)
    
if os.path.exists(pas_sensibilite_results_alpha2):
    os.remove(pas_sensibilite_results_alpha2)
    



# Ecriture des fichiers de input
# Pour alpha = 0 deg
f = open(pas_sensibilite_input_alpha1, 'w')
f.write('naca' + naca_profil + '\n')
f.write('ppar \n')
f.write('n' + str(n_panel_nodes) + '\n')
f.write('p' + str(p) + '\n')
f.write('t' + str(te_le_density_ratio) + '\n')
f.write('\n\n')
f.write('oper\n')
f.write('visc\n')
f.write(str(re) + '\n')
f.write('iter' + str(n_iterations) + '\n')
f.write('pacc\n')
f.write(pas_sensibilite_results_alpha1 + '\n')
f.write('\n')
f.write('aseq\n' + str(alpha1 - n_pas_sensibilite_alpha * pas_sensibilite_alpha) + '\n' + str(alpha1 + n_pas_sensibilite_alpha * pas_sensibilite_alpha) + '\n' + str(pas_sensibilite_alpha) + '\n')
f.write('\n')
f.write('quit\n')
f.close()

# Pour alpha = 15 deg
f = open(pas_sensibilite_input_alpha2, 'w')
f.write('naca' + naca_profil + '\n')
f.write('ppar \n')
f.write('n' + str(n_panel_nodes) + '\n')
f.write('p' + str(p) + '\n')
f.write('t' + str(te_le_density_ratio) + '\n')
f.write('\n\n')
f.write('oper\n')
f.write('visc\n')
f.write(str(re) + '\n')
f.write('iter' + str(n_iterations) + '\n')
f.write('pacc\n')
f.write(pas_sensibilite_results_alpha2 + '\n')
f.write('\n')
f.write('aseq\n' + str(alpha2 - n_pas_sensibilite_alpha * pas_sensibilite_alpha) + '\n' + str(alpha2 + n_pas_sensibilite_alpha * pas_sensibilite_alpha) + '\n' + str(pas_sensibilite_alpha) + '\n')
f.write('\n')
f.write('quit\n')
f.close()

# Execution de XFoil
os.system('xfoil.exe < ' + pas_sensibilite_input_alpha1)
os.system('xfoil.exe < ' + pas_sensibilite_input_alpha2)

# Lecture des fichiers de resultats
pas_sensibilite_alpha_data1 = np.loadtxt(pas_sensibilite_results_alpha1, skiprows=12)
pas_sensibilite_alpha_data2 = np.loadtxt(pas_sensibilite_results_alpha2, skiprows=12)

angle_sensibilite_alpha1 = pas_sensibilite_alpha_data1[:,0]
angle_sensibilite_alpha2 = pas_sensibilite_alpha_data2[:,0]

cl_sensibilite_alpha1 = pas_sensibilite_alpha_data1[:,1]
cl_sensibilite_alpha2 = pas_sensibilite_alpha_data2[:,1]

cd_sensibilite_alpha1 = pas_sensibilite_alpha_data1[:,2]
cd_sensibilite_alpha2 = pas_sensibilite_alpha_data2[:,2]

delta_alpha = np.arange(pas_sensibilite_alpha, (n_pas_sensibilite_alpha + 1) * pas_sensibilite_alpha, pas_sensibilite_alpha)
sensibilite_cl_alpha1 = np.zeros((n_pas_sensibilite_alpha))
sensibilite_cl_alpha2 = np.zeros((n_pas_sensibilite_alpha))
sensibilite_cd_alpha1 = np.zeros((n_pas_sensibilite_alpha))
sensibilite_cd_alpha2 = np.zeros((n_pas_sensibilite_alpha))

for i in range(n_pas_sensibilite_alpha):
    sensibilite_cl_alpha1[-(i+1)] = (cl_sensibilite_alpha1[-(i+1)] - cl_sensibilite_alpha1[i]) / (2 * delta_alpha[-(i+1)])
    sensibilite_cl_alpha2[-(i+1)] = (cl_sensibilite_alpha2[-(i+1)] - cl_sensibilite_alpha2[i]) / (2 * delta_alpha[-(i+1)])
    sensibilite_cd_alpha1[-(i+1)] = (cd_sensibilite_alpha1[-(i+1)] - cd_sensibilite_alpha1[i]) / (2 * delta_alpha[-(i+1)])
    sensibilite_cd_alpha2[-(i+1)] = (cd_sensibilite_alpha2[-(i+1)] - cd_sensibilite_alpha2[i]) / (2 * delta_alpha[-(i+1)])
    

fig1= plt.figure()
plt.plot(delta_alpha , sensibilite_cl_alpha1)
plt.xlabel("d_alpha")
plt.ylabel("d_CL/d_alpha")
plt.title('Variation de la sensibilité du CL en fonction du pas de alpha autour de alpha = 0')
plt.show()
    
fig2= plt.figure()
plt.plot(delta_alpha , sensibilite_cl_alpha2)
plt.xlabel("d_alpha")
plt.ylabel("d_CL/d_alpha")
plt.title('Variation de la sensibilité du CL en fonction du pas de alpha autour de alpha = 15')
plt.show()    

fig3= plt.figure()
plt.plot(delta_alpha , sensibilite_cd_alpha1)
plt.xlabel("d_alpha")
plt.ylabel("d_CD/d_alpha")
plt.title('Variation de la sensibilité du CD en fonction du pas de alpha autour de alpha = 0')
plt.show()
    
fig4= plt.figure()
plt.plot(delta_alpha , sensibilite_cd_alpha2)
plt.xlabel("d_alpha")
plt.ylabel("d_CD/d_alpha")
plt.title('Variation de la sensibilité du CD en fonction du pas de alpha autour de alpha = 15')
plt.show()   

# Creation des noms de fichiers textes pour le pas en Reynolds
pas_sensibilite_results_re1 = 'pas_sensibilite_results_re1.txt'
pas_sensibilite_input_re1 = 'pas_sensibilite_input_re1.txt'
pas_sensibilite_results_re2 = 'pas_sensibilite_results_re2.txt'
pas_sensibilite_input_re2 = 'pas_sensibilite_input_re2.txt'

i_re = 0
sensibilite_data_cl_re1 = np.zeros(len(re_sensibilite)) 
sensibilite_data_cd_re1 = np.zeros(len(re_sensibilite)) 
sensibilite_data_cl_re2 = np.zeros(len(re_sensibilite)) 
sensibilite_data_cd_re2 = np.zeros(len(re_sensibilite)) 

n_panel_nodes = 160
te_le_density_ratio = .2
p = 1
n_iterations = 2000

for r in re_sensibilite:
    if os.path.exists(pas_sensibilite_results_re1):
        os.remove(pas_sensibilite_results_re1)
    
    if os.path.exists(pas_sensibilite_results_re2):
        os.remove(pas_sensibilite_results_re2)
        
    # Pour alpha = 0 deg
    f = open(pas_sensibilite_input_re1, 'w')
    f.write('naca' + naca_profil + '\n')
    f.write('ppar \n')
    f.write('n' + str(n_panel_nodes) + '\n')
    f.write('p' + str(p) + '\n')
    f.write('t' + str(te_le_density_ratio) + '\n')
    f.write('\n\n')
    f.write('oper\n')
    f.write('visc\n')
    f.write(str(r) + '\n')
    f.write('iter' + str(n_iterations) + '\n')
    f.write('pacc\n')
    f.write(pas_sensibilite_results_re1 + '\n')
    f.write('\n')
    f.write('alfa' + str(alpha1) + '\n')
    f.write('\n')
    f.write('quit\n')
    f.close()
    
    # Pour alpha = 15 deg
    f = open(pas_sensibilite_input_re2, 'w')
    f.write('naca' + naca_profil + '\n')
    f.write('ppar \n')
    f.write('n' + str(n_panel_nodes) + '\n')
    f.write('p' + str(p) + '\n')
    f.write('t' + str(te_le_density_ratio) + '\n')
    f.write('\n\n')
    f.write('oper\n')
    f.write('visc\n')
    f.write(str(r) + '\n')
    f.write('iter' + str(n_iterations) + '\n')
    f.write('pacc\n')
    f.write(pas_sensibilite_results_re2 + '\n')
    f.write('\n')
    f.write('alfa' + str(alpha2) + '\n')
    f.write('\n')
    f.write('quit\n')
    f.close()
        
        
    os.system('xfoil.exe < ' + pas_sensibilite_input_re1)
    pas_sensibilite_data_re1 = np.loadtxt(pas_sensibilite_results_re1, skiprows=12)
    sensibilite_data_cl_re1[i_re] = pas_sensibilite_data_re1[1]    
    sensibilite_data_cd_re1[i_re] = pas_sensibilite_data_re1[2] 
    
    os.system('xfoil.exe < ' + pas_sensibilite_input_re2)
    pas_sensibilite_data_re2 = np.loadtxt(pas_sensibilite_results_re2, skiprows=12)
    sensibilite_data_cl_re2[i_re] = pas_sensibilite_data_re2[1]    
    sensibilite_data_cd_re2[i_re] = pas_sensibilite_data_re2[2] 
    
    i_re += 1



delta_re = np.arange(pas_sensibilite_re, (n_pas_sensibilite_alpha + 1) * pas_sensibilite_re, pas_sensibilite_re)
    
sensibilite_cl_re1 = np.zeros((n_pas_sensibilite_alpha))
sensibilite_cl_re2 = np.zeros((n_pas_sensibilite_alpha))
sensibilite_cd_re1 = np.zeros((n_pas_sensibilite_alpha))
sensibilite_cd_re2 = np.zeros((n_pas_sensibilite_alpha))
        
for i in range(n_pas_sensibilite_alpha):
    sensibilite_cl_re1[-(i+1)] = (sensibilite_data_cl_re1[-(i+1)] - sensibilite_data_cl_re1[i]) / (2 * delta_re[-(i+1)])
    sensibilite_cl_re2[-(i+1)] = (sensibilite_data_cl_re2[-(i+1)] - sensibilite_data_cl_re2[i]) / (2 * delta_re[-(i+1)])
    sensibilite_cd_re1[-(i+1)] = (sensibilite_data_cd_re1[-(i+1)] - sensibilite_data_cd_re1[i]) / (2 * delta_re[-(i+1)])
    sensibilite_cd_re2[-(i+1)] = (sensibilite_data_cd_re2[-(i+1)] - sensibilite_data_cd_re2[i]) / (2 * delta_re[-(i+1)])
        
fig1= plt.figure()
plt.plot(delta_re , sensibilite_cl_re1)
plt.xlabel("d_alpha")
plt.ylabel("d_CL/d_Re")
plt.title('Variation de la sensibilité du CL en fonction du pas de Re autour de alpha = 0')
plt.show()
    
fig2= plt.figure()
plt.plot(delta_re , sensibilite_cl_re2)
plt.xlabel("d_alpha")
plt.ylabel("d_CL/d_re")
plt.title('Variation de la sensibilité du CL en fonction du pas de Re autour de alpha = 15')
plt.show()    

fig3= plt.figure()
plt.plot(delta_re , sensibilite_cd_re1)
plt.xlabel("d_alpha")
plt.ylabel("d_CD/d_re")
plt.title('Variation de la sensibilité du CD en fonction du pas de Re autour de alpha = 0')
plt.show()
    
fig4= plt.figure()
plt.plot(delta_re , sensibilite_cd_re2)
plt.xlabel("d_alpha")
plt.ylabel("d_CD/d_re")
plt.title('Variation de la sensibilité du CD en fonction du pas de Re autour de alpha = 15')
plt.show()   

d_alpha = 0.2
d_re = 5e5

sensibilite_results_alpha_min_re_min = 'sensibilite_results_alpha_min_re_min.txt'
sensibilite_input_alpha_min_re_min = 'sensibilite_input_alpha_min_re_min.txt'

sensibilite_results_alpha_min_re_max = 'sensibilite_results_alpha_min_re_max.txt'
sensibilite_input_alpha_min_re_max = 'sensibilite_input_alpha_min_re_max.txt'

sensibilite_results_alpha_max_re_min = 'sensibilite_results_alpha_max_re_min.txt'
sensibilite_input_alpha_max_re_min = 'sensibilite_input_alpha_max_re_min.txt'

sensibilite_results_alpha_max_re_max = 'sensibilite_results_alpha_max_re_max.txt'
sensibilite_input_alpha_max_re_max = 'sensibilite_input_alpha_max_re_max.txt'

n_panel_nodes = 180
te_le_density_ratio = .2
p = 1
n_iterations = 2000

if os.path.exists(sensibilite_results_alpha_min_re_min):
    os.remove(sensibilite_results_alpha_min_re_min)
        
if os.path.exists(sensibilite_results_alpha_min_re_max):
    os.remove(sensibilite_results_alpha_min_re_max)
    
if os.path.exists(sensibilite_results_alpha_max_re_min):
    os.remove(sensibilite_results_alpha_max_re_min)
    
if os.path.exists(sensibilite_results_alpha_max_re_max):
    os.remove(sensibilite_results_alpha_max_re_max)
        

# Alpha min, Re min
f = open(sensibilite_input_alpha_min_re_min, 'w')
f.write('naca' + naca_profil + '\n')
f.write('ppar \n')
f.write('n' + str(n_panel_nodes) + '\n')
f.write('p' + str(p) + '\n')
f.write('t' + str(te_le_density_ratio) + '\n')
f.write('\n\n')
f.write('oper\n')
f.write('visc\n')
f.write(str(re-d_re) + '\n')
f.write('iter' + str(n_iterations) + '\n')
f.write('pacc\n')
f.write(sensibilite_results_alpha_min_re_min + '\n')
f.write('\n')
f.write('aseq\n' + str(alpha_start - d_alpha) + '\n' + str(alpha_end - d_alpha) + '\n' + str(alpha_step) + '\n')
f.write('quit\n')
f.close() 
os.system('xfoil.exe < ' + sensibilite_input_alpha_min_re_min)

n_panel_nodes = 180
te_le_density_ratio = .15
p = .9
n_iterations = 2000

# Alpha min, Re max
f = open(sensibilite_input_alpha_min_re_max, 'w')
f.write('naca' + naca_profil + '\n')
f.write('ppar \n')
f.write('n' + str(n_panel_nodes) + '\n')
f.write('p' + str(p) + '\n')
f.write('t' + str(te_le_density_ratio) + '\n')
f.write('\n\n')
f.write('oper\n')
f.write('visc\n')
f.write(str(re+d_re) + '\n')
f.write('iter' + str(n_iterations) + '\n')
f.write('pacc\n')
f.write(sensibilite_results_alpha_min_re_max + '\n')
f.write('\n')
f.write('aseq\n' + str(alpha_start - d_alpha) + '\n' + str(alpha_end - d_alpha) + '\n' + str(alpha_step) + '\n')
f.write('quit\n')
f.close()  
os.system('xfoil.exe < ' + sensibilite_input_alpha_min_re_max)

# Alpha max, Re min
f = open(sensibilite_input_alpha_max_re_min, 'w')
f.write('naca' + naca_profil + '\n')
f.write('ppar \n')
f.write('n' + str(n_panel_nodes) + '\n')
f.write('p' + str(p) + '\n')
f.write('t' + str(te_le_density_ratio) + '\n')
f.write('\n\n')
f.write('oper\n')
f.write('visc\n')
f.write(str(re-d_re) + '\n')
f.write('iter' + str(n_iterations) + '\n')
f.write('pacc\n')
f.write(sensibilite_results_alpha_max_re_min + '\n')
f.write('\n')
f.write('aseq\n' + str(alpha_start + d_alpha) + '\n' + str(alpha_end + d_alpha) + '\n' + str(alpha_step) + '\n')
f.write('quit\n')
f.close()
os.system('xfoil.exe < ' + sensibilite_input_alpha_max_re_min)

# Alpha max, Re max
f = open(sensibilite_input_alpha_max_re_max, 'w')
f.write('naca' + naca_profil + '\n')
f.write('ppar \n')
f.write('n' + str(n_panel_nodes) + '\n')
f.write('p' + str(p) + '\n')
f.write('t' + str(te_le_density_ratio) + '\n')
f.write('\n\n')
f.write('oper\n')
f.write('visc\n')
f.write(str(re+d_re) + '\n')
f.write('iter' + str(n_iterations) + '\n')
f.write('pacc\n')
f.write(sensibilite_results_alpha_max_re_max + '\n')
f.write('\n')
f.write('aseq\n' + str(alpha_start + d_alpha) + '\n' + str(alpha_end + d_alpha) + '\n' + str(alpha_step) + '\n')
f.write('quit\n')
f.close()  
os.system('xfoil.exe < ' + sensibilite_input_alpha_max_re_max)

sensibilite_min_min = np.loadtxt(sensibilite_results_alpha_min_re_min, skiprows=12)
sensibilite_min_max = np.loadtxt(sensibilite_results_alpha_min_re_max, skiprows=12)
sensibilite_max_min = np.loadtxt(sensibilite_results_alpha_max_re_min, skiprows=12)
sensibilite_max_max = np.loadtxt(sensibilite_results_alpha_max_re_max, skiprows=12)

d_cl_d_alpha=(sensibilite_max_min[:,1]-sensibilite_min_min[:,1])/(d_alpha*2)
d_cd_d_alpha=(sensibilite_max_min[:,2]-sensibilite_min_min[:,2])/(d_alpha*2)
d_cl_d_re=(sensibilite_max_max[:,1]-sensibilite_max_min[:,1])/(d_re*2)
d_cd_d_re=(sensibilite_max_max[:,2]-sensibilite_max_min[:,2])/(d_re*2)

alpha=np.arange(0,18,1)


fig1= plt.figure()
plt.plot(alpha, d_cl_d_alpha)
plt.xlabel("alpha")
plt.ylabel("d_CL/d_alpha")
plt.title('Sensibilité du CL/d_alpha en fonction de alpha ')
plt.show()
    
fig2= plt.figure()
plt.plot(alpha, d_cd_d_alpha)
plt.xlabel("alpha")
plt.ylabel("d_CD/d_alpha")
plt.title('Sensibilité du CD/d_alpha en fonction de alpha')
plt.show()    

fig3= plt.figure()
plt.plot(alpha, d_cl_d_re)
plt.xlabel("alpha")
plt.ylabel("d_CL/d_re")
plt.title('Sensibilité du CL/d_Re en fonction de alpha')
plt.show()
    
fig4= plt.figure()
plt.plot(alpha, d_cd_d_re)
plt.xlabel("alpha")
plt.ylabel("d_CD/d_re")
plt.title('Sensibilité du CD/d_Re en fonction de alpha')
plt.show()   