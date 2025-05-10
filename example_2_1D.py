"""
Description:
    This script implements the Finite Element Method (FEM) for solving the 
    Modified Poisson–Nernst–Planck/Navier-Stokes (PNP/NS) equations with 
    non-homogeneous boundary conditions and 4 Dirichlet conditions.
    
    The model is based on the paper:
    "Finite Element Method for the Numerical Simulation of Modified PNP/NS Model"
    (https://arxiv.org/abs/2409.08746)

    The script computes electrostatic potential, ion concentrations, 
    and velocity fields for different parameter settings.

Dependencies:
    - FEniCS
    - NumPy
    - SciPy
    - Matplotlib (optional, for plotting)

Author: Ankur, IIT Roorkee
Date: 06-05-2025
"""

from __future__ import print_function
from fenics import *
import numpy as np

def q_1(y_c, y_a, alpha_c, alpha_a):    # Y_c = u_1, Y_a = u_2
    "Return nonlinear coefficient"
    return -( (1/(alpha_c*y_c))+(1/(1-y_c-y_a))    )
def q_2(y_c, y_a):
    "Return nonlinear coefficient"
    return -( 1/(1-y_c-y_a)   )
def q_2_1(y_c, y_a, alpha_c, alpha_a, z_c, z_a, z_s, psi):
    "Return nonlinear coefficient"
    return psi*(     ((1/alpha_c) -1 )   *  (z_c*y_c+z_a*y_a) - ((z_c/alpha_c) -z_s )    )
def q_3(y_c, y_a):
    "Return nonlinear coefficient"
    return -( 1/(1-y_c-y_a)   )
def q_4(y_c, y_a, alpha_c, alpha_a):
    "Return nonlinear coefficient"
    return -( (1/(alpha_a*y_a))+(1/(1-y_c-y_a))    )
def q_4_1(y_c, y_a, alpha_c, alpha_a, z_c, z_a, z_s, psi):
    "Return nonlinear coefficient"
    return psi*(     ((1/alpha_a) -1 )   *  (z_c*y_c+z_a*y_a) - ((z_a/alpha_a) -z_s )    )
def q_new(new_constant_new, z_c, z_a,y_c, y_a,aa ):
     return new_constant_new* (z_c*y_c+z_a*y_a)*aa

# Create mesh 
mesh = UnitIntervalMesh(1000)

P1 = FiniteElement('P', mesh.ufl_cell(), 1)
R1 =  FiniteElement('R', mesh.ufl_cell(), 0)

element = MixedElement([P1, R1, P1, R1, P1, R1, P1])
V= FunctionSpace(mesh, element)

#Exact solution and Dirichlet conditions

alpha_c = 0.1
alpha_a = 0.1
z_c =1.0
z_a = -1.0
z_s = 0.0
psi = 1.0

Tilde = 1000.0
chaai = 1.0
new_constant = (Tilde/psi)*(1/(1+chaai))

k_cap = 1.0
new_constant_new = (psi/k_cap)


"""print('ue1, ue2 =', ue1,'and', ue2)
print('u1_neu =', g3_Top,'and', g3_Bottom)
print('f1 and f2 =', f_1,'and', f_2)
print('f1 and f2 =', g1_Total_top,'and', g1_Total_bottom, 'and', g2_Total_top,'and', g2_Total_bottom)
print('f1 and f2 =', g1_Total_right,'and', g1_Total_left, 'and', g2_Total_right,'and', g2_Total_left)"""



tol = 1E-14
boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)

class Left(SubDomain):
     def inside(self, x, on_boundary):
            return (on_boundary and near(x[0], 0, tol) )  #Neumann

class Right(SubDomain):
      def inside(self, x, on_boundary):
            return (on_boundary and near(x[0], 1, tol) ) #Neumann 
        
"""class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0)
# Initialize sub-domain instances"""
left = Left()
#top = Top()
right = Right()
#bottom = Bottom()

boundaries.set_all(0)
left.mark(boundaries, 1)
#top.mark(boundaries, 2)
right.mark(boundaries, 3)
#bottom.mark(boundaries, 4)

# Values for u3_bc1 and u3_bc3
u3_bc1_values = [-0.5, -1.0, -2.0]
u3_bc3_values = [0.5, 1.0, 2.0]

# Define u3_bc1 and u3_bc3 as Constants
u3_bc1 = Constant(0.0)  
u3_bc3 = Constant(0.0) 

# Lists to store solutions
solutions1 = []
solutions2 = []
solutions3 = []
solutions4 = []
# Loop through values and solve the problem
for u3_bc1_val, u3_bc3_val in zip(u3_bc1_values, u3_bc3_values):
    # Update Dirichlet boundary conditions
    u3_bc1.assign(Constant(u3_bc1_val))
    u3_bc3.assign(Constant(u3_bc3_val))

    bc31 = DirichletBC(V.sub(4), u3_bc1, boundaries, 1)
    bc32 = DirichletBC(V.sub(4), u3_bc3, boundaries, 3)
    bcs = [bc31, bc32]
    ########################################################################
    domains = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
    domains.set_all(0)
    #g=Constant(0.0)
    # Define variational problem
    dx = Measure('dx', domain=mesh, subdomain_data=domains)
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

    u1_0 = Constant(0.40)
    u2_0 = Constant(0.40)
    u4_0 = Constant(1.0)
    int_u1 = assemble(u1_0*dx(degree=5))  #integration of u (Extra condition)
    int_u2 = assemble(u2_0*dx(degree=5))  #integration of u (Extra condition)
    int_u4 = assemble(u4_0*dx(degree=5))  #integration of u (Extra condition)

    u = TrialFunction(V)
    u_1, c_1, u_2, c_2, u_3, c_4, u_4  = split(u)
    v_1, d_1, v_2, d_2, v_3, d_4, v_4 = TestFunctions(V)
    ##############################################################################################################
    #############################################################################################
    """a_con=Constant(1.001)
    F =  dot(grad(u_1) + a_con*grad(u_2) +  grad(u_3), grad(v_1))*dx(0) - f_1*v_1*dx(0)\
  -(g1_Total_top)*v_1*ds(2) -(g1_Total_bottom)*v_1*ds(4) -(g1_Total_right)*v_1*ds(3) -(g1_Total_left)*v_1*ds(1) \
  +c_1*v_1*dx(0) + u_1*d_1*dx(0) - int_u1*d_1*dx(0) \
  + dot(grad(u_1) + grad(u_2) +  grad(u_3), grad(v_2))*dx(0)  - f_2*v_2*dx(0)\
  -(g2_Total_top)*v_2*ds(2)- (g2_Total_bottom)*v_2*ds(4) -(g2_Total_right)*v_2*ds(3)- (g2_Total_left)*v_2*ds(1)\
   +c_2*v_2*dx(0) + u_2*d_2*dx(0) - int_u2*d_2*dx(0) \
   + dot(grad(u_3), grad(v_3))*dx(0) - ( new_constant* (z_c*u_1+z_a*u_2) )*v_3*dx(0)- f_3*v_3*dx(0)  \
    - g3_Top*v_3*ds(2) - g3_Bottom*v_3*ds(4) \
    - dot(grad(u_4)+ grad(u_3), grad(v_4))*dx(0) - f_4*v_4*dx(0) \
    +(g4_Total_top)*v_4*ds(2) +(g4_Total_bottom)*v_4*ds(4) +(g4_Total_right)*v_4*ds(3) +(g4_Total_left)*v_4*ds(1) \
     +c_4*v_4*dx(0) + u_4*d_4*dx(0) - int_u4*d_4*dx(0)"""
    F =  dot(q_1(u1_0, u2_0, alpha_c, alpha_a)*grad(u_1) + q_2(u1_0, u2_0)*grad(u_2) +  q_2_1(u1_0, u2_0, alpha_c, alpha_a, z_c, z_a, z_s, psi)*grad(u_3), grad(v_1))*dx(0) \
  +c_1*v_1*dx(0) + u_1*d_1*dx(0) - int_u1*d_1*dx(0) \
  + dot(q_3(u1_0, u2_0)*grad(u_1) + q_4(u1_0, u2_0, alpha_c, alpha_a)*grad(u_2) +  q_4_1(u1_0, u2_0, alpha_c, alpha_a, z_c, z_a, z_s, psi)*grad(u_3), grad(v_2))*dx(0)  \
   +c_2*v_2*dx(0) + u_2*d_2*dx(0) - int_u2*d_2*dx(0) \
   + dot(grad(u_3), grad(v_3))*dx(0) - ( new_constant* (z_c*u_1+z_a*u_2)*u4_0 )*v_3*dx(0)  \
    - dot(grad(u_4)+ q_new(new_constant_new, z_c, z_a, u1_0, u2_0, u4_0)*grad(u_3), grad(v_4))*dx(0)  \
     +c_4*v_4*dx(0) + u_4*d_4*dx(0) - int_u4*d_4*dx(0) 
    a, L = lhs(F), rhs(F)
    A, b = assemble_system(a, L)
    for bc in bcs:
      bc.apply(A, b)
    u_k = Function(V)
    solve(A, u_k.vector(), b, 'lu')
    u_k_1, c_k_1, u_k_2, c_k_2, u_k_3, c_k_4, u_k_4= split(u_k)
    

    F =  dot(q_1(u_k_1, u_k_2, alpha_c, alpha_a)*grad(u_1) + q_2(u_k_1, u_k_2)*grad(u_2) +  q_2_1(u_k_1, u_k_2, alpha_c, alpha_a, z_c, z_a, z_s, psi)*grad(u_3), grad(v_1))*dx(0) \
  +c_1*v_1*dx(0) + u_1*d_1*dx(0) - int_u1*d_1*dx(0) \
  + dot(q_3(u_k_1, u_k_2)*grad(u_1) + q_4(u_k_1, u_k_2, alpha_c, alpha_a)*grad(u_2) +  q_4_1(u_k_1, u_k_2, alpha_c, alpha_a, z_c, z_a, z_s, psi)*grad(u_3), grad(v_2))*dx(0)  \
   +c_2*v_2*dx(0) + u_2*d_2*dx(0) - int_u2*d_2*dx(0) \
   + dot(grad(u_3), grad(v_3))*dx(0) - ( new_constant* (z_c*u_1+z_a*u_2)*u_k_4 )*v_3*dx(0) \
    - dot(grad(u_4)+ q_new(new_constant_new, z_c, z_a, u_k_1, u_k_2,u_k_4)*grad(u_3), grad(v_4))*dx(0)  \
     +c_4*v_4*dx(0) + u_4*d_4*dx(0) - int_u4*d_4*dx(0) 
    a, L = lhs(F), rhs(F)


    # Picard iterations
    u = Function(V)   
    eps = 1.0          
    tol = 1.0E-12       
    iter = 0           
    maxiter = 25       
    while eps > tol and iter < maxiter:
     iter += 1
     solve(a == L, u, bcs)
     diff = np.array(u.vector()) - np.array(u_k.vector())
     eps = np.linalg.norm(diff, ord=np.Inf)
     print ('iter=%d: norm=%g' % (iter, eps))
     u_k.assign(u)   # update for next iteration
    # Save solution to file (VTK)
    (u1, c1, u2, c2, u3, c4, u4) = u.split()
    u1.rename("Solution Cation", "")
    u2.rename("Solution Anion", "")
    u3.rename("Potential", "")
    u4.rename("Total number density", "")
    # Append the solution to the list
    solutions1.append(u1.copy())
    solutions2.append(u2.copy())
    solutions3.append(u3.copy())
    solutions4.append(u4.copy())
#################
import matplotlib.pyplot as plt
# Plotting all solutions together
for i, solution in enumerate(solutions1):
    plot(solution, linewidth=3)

plt.xlabel('x-axis', fontsize=14)
plt.ylabel(r'$y_C$', fontsize=14)
plt.ylim(0, 1)  # Set y-axis limits
plt.legend([fr'$\varphi$={u3_bc1_val} on $\Gamma_D^L$, $\varphi$={u3_bc3_val} on $\Gamma_D^R$' for u3_bc1_val, u3_bc3_val in zip(u3_bc1_values, u3_bc3_values)], loc='best', fontsize=14)
# Adjust the scale of the x-axis and y-axis
plt.xticks(fontsize=12)  # You can adjust the font size as needed
plt.yticks(fontsize=12)  # You can adjust the font size as needed
plt.show()
#################
# Plotting all solutions together
for i, solution in enumerate(solutions2):
    plot(solution, linewidth=3)

plt.xlabel('x-axis', fontsize=14)
plt.ylabel(r'$y_A$', fontsize=14)
plt.ylim(0, 1)  # Set y-axis limits
plt.legend([fr'$\varphi$={u3_bc1_val} on $\Gamma_D^L$, $\varphi$={u3_bc3_val} on $\Gamma_D^R$' for u3_bc1_val, u3_bc3_val in zip(u3_bc1_values, u3_bc3_values)], loc='best', fontsize=14)
# Adjust the scale of the x-axis and y-axis
plt.xticks(fontsize=12)  # You can adjust the font size as needed
plt.yticks(fontsize=12)  # You can adjust the font size as needed
plt.show()
#################
# Plotting all solutions together
for i, solution in enumerate(solutions3):
    plot(solution, linewidth=3)

plt.xlabel('x-axis', fontsize=14)
plt.ylabel(r'$\varphi$', fontsize=16)
plt.ylim(-2, 2)  # Set y-axis limits
plt.legend([fr'$\varphi$={u3_bc1_val} on $\Gamma_D^L$, $\varphi$={u3_bc3_val} on $\Gamma_D^R$' for u3_bc1_val, u3_bc3_val in zip(u3_bc1_values, u3_bc3_values)], loc='best', fontsize=14)
# Adjust the scale of the x-axis and y-axis
plt.xticks(fontsize=12)  # You can adjust the font size as needed
plt.yticks(fontsize=12)  # You can adjust the font size as needed
plt.show()
#################
# Plotting all solutions together
for i, solution in enumerate(solutions4):
    plot(solution, linewidth=3)

plt.xlabel('x-axis', fontsize=14)
plt.ylabel('$n$', fontsize=16)
plt.ylim(0.9, 3)  # Set y-axis limits
plt.legend([fr'$\varphi$={u3_bc1_val} on $\Gamma_D^L$, $\varphi$={u3_bc3_val} on $\Gamma_D^R$' for u3_bc1_val, u3_bc3_val in zip(u3_bc1_values, u3_bc3_values)], loc='best', fontsize=14)
plt.xticks(fontsize=12) 
plt.yticks(fontsize=12) 
plt.show()
