

from __future__ import print_function
from fenics import *
def q_1(y_c, y_a):    # Y_c = u_1, Y_a = u_2
    "Return nonlinear coefficient"
    return -( (1/y_c)+(1/(1-y_c-y_a))    )
def q_2(y_c, y_a):
    "Return nonlinear coefficient"
    return -( 1/(1-y_c-y_a)   )
def q_3(y_c, y_a):
    "Return nonlinear coefficient"
    return -( 1/(1-y_c-y_a)   )
def q_4(y_c, y_a):
    "Return nonlinear coefficient"
    return -( (1/y_a)+(1/(1-y_c-y_a))    )

import sympy as sym

# Create mesh 
mesh = UnitSquareMesh(32, 32)
# Define function space for system of solutions
#P1 = FiniteElement('P', triangle, 1)
element = VectorElement('P', triangle, 1, dim=2)
V = FunctionSpace(mesh, element)

x, y = sym.symbols('x[0], x[1]')
#Exact solution and Dirichlet conditions
#ue1 = (2+x + 0.5*y**2)

#ue2 = (3+ x**2 + 2*y**2)
#ue1 = 1/(2+x + 0.5*y**2)

#ue2 = 1/(3+ x**2 + 2*y**2)

ue1 = (sym.exp(-0.1*x-0.2*y) )

ue2 = (sym.exp(-0.2*x-0.1*y))


# Right hand side
f_1 = - sym.diff(q_1(ue1, ue2)*sym.diff(ue1, x), x) - sym.diff(q_1(ue1, ue2)*sym.diff(ue1, y), y) -\
        sym.diff(q_2(ue1, ue2)*sym.diff(ue2, x), x) - sym.diff(q_2(ue1, ue2)*sym.diff(ue2, y), y)  # -Laplace(u)
f_2 = - sym.diff(q_3(ue1, ue2)*sym.diff(ue1, x), x) - sym.diff(q_3(ue1, ue2)*sym.diff(ue1, y), y) -\
        sym.diff(q_4(ue1, ue2)*sym.diff(ue2, x), x) - sym.diff(q_4(ue1, ue2)*sym.diff(ue2, y), y)  # -Laplace(u)
f_1 = sym.simplify(f_1)  
f_2 = sym.simplify(f_2)   
g1_Total_top= q_1(ue1, ue2)*sym.diff(ue1, y).subs(y, 1) + q_2(ue1, ue2)*sym.diff(ue2, y).subs(y, 1)     # For 1st equatioon i.e for (sum of normal derivative )
g1_Total_bottom= -q_1(ue1, ue2)*sym.diff(ue1, y).subs(y, 0)  -q_2(ue1, ue2)*sym.diff(ue2, y).subs(y, 0)     # For 1st equatioon i.e for (sum of normal derivative )
g2_Total_top= q_3(ue1, ue2)*sym.diff(ue1, y).subs(y, 1) + q_4(ue1, ue2)*sym.diff(ue2, y).subs(y, 1)      # For 2nd equatioon i.e for (u_1 +2*u_2 )
g2_Total_bottom= -q_3(ue1, ue2)*sym.diff(ue1, y).subs(y, 0)  -q_4(ue1, ue2)*sym.diff(ue2, y).subs(y, 0)   # For 2nd equatioon i.e for (u_1 +2*u_2 )

g1_Total_right= q_1(ue1, ue2)*sym.diff(ue1, x).subs(x, 1) + q_2(ue1, ue2)*sym.diff(ue2, x).subs(x, 1)     # For 1st equatioon i.e for (sum of normal derivative )
g1_Total_left= -q_1(ue1, ue2)*sym.diff(ue1, x).subs(x, 0)  -q_2(ue1, ue2)*sym.diff(ue2, x).subs(x, 0)     # For 1st equatioon i.e for (sum of normal derivative )
g2_Total_right= q_3(ue1, ue2)*sym.diff(ue1, x).subs(x, 1) + q_4(ue1, ue2)*sym.diff(ue2, x).subs(x, 1)      # For 2nd equatioon i.e for (u_1 +2*u_2 )
g2_Total_left= -q_3(ue1, ue2)*sym.diff(ue1, x).subs(x, 0)  -q_4(ue1, ue2)*sym.diff(ue2, x).subs(x, 0)   # For 2nd equatioon i.e for (u_1 +2*u_2 )
# Define boundary conditions
# Collect variables
variables = [ue1, ue2, f_1, f_2, g1_Total_top, g1_Total_bottom, g2_Total_top, g2_Total_bottom]
variables1 = [ g1_Total_left, g1_Total_right, g2_Total_left, g2_Total_right]

# Turn into C/C++ code strings
variables = [sym.printing.ccode(var) for var in variables]
variables1 = [sym.printing.ccode(var) for var in variables1]


"""print('ue1, ue2 =', ue1,'and', ue2)
print('u1_Diri =', uD1_Left,'and', uD1_Right)
print('u2_Diri =', uD2_Left,'and', uD2_Right)
print('u1_neu =', g1_Top,'and', g1_Bott)
print('u2_neu =', g2_Top,'and', g2_Bott)
print('f1 and f2 =', f1,'and', f2)
print('f1 and f2 =', g1_Total_top,'and', g1_Total_bottom, 'and', g2_Total_top,'and', g2_Total_bottom)"""

# Turn into FEniCS Expressions
variables = [Expression(var, degree=5) for var in variables]
variables1 = [Expression(var, degree=5) for var in variables1]


# Extract variables
ue1, ue2, f_1, f_2, g1_Total_top, g1_Total_bottom, g2_Total_top, g2_Total_bottom = variables
g1_Total_left, g1_Total_right, g2_Total_left, g2_Total_right = variables1





tol = 1E-14
boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)

class Left(SubDomain):
     def inside(self, x, on_boundary):
            return (on_boundary and near(x[0], 0, tol) and x[1]> 0.5-tol and x[1]< 1.0+tol)  #Dirichlet

class BoundaryX00(SubDomain):
        def inside(self, x, on_boundary):
            return (on_boundary and near(x[0], 0, tol) and x[1]< 0.5+tol and x[1]> 0.0-tol)   #Neumannn

class Right(SubDomain):
      def inside(self, x, on_boundary):
            return (on_boundary and near(x[0], 1, tol) and x[1]< 0.5+tol and x[1]> 0.0-tol) #Dirichlet

class BoundaryX11(SubDomain):
        def inside(self, x, on_boundary):
            return (on_boundary and near(x[0], 1, tol) and x[1]> 0.5-tol and x[1]< 1.0+tol) #Neumannn  
        
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0)
# Initialize sub-domain instances
left = Left()
top = Top()
right = Right()
bottom = Bottom()
bx00 = BoundaryX00()
bx11 = BoundaryX11()

boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)

bx00.mark(boundaries, 5)
bx11.mark(boundaries, 6)

bc11 = DirichletBC(V.sub(0), ue1,  boundaries, 1)
bc12 = DirichletBC(V.sub(0), ue1,  boundaries, 2)
bc13 = DirichletBC(V.sub(0), ue1,  boundaries, 3)
bc14 = DirichletBC(V.sub(0), ue1,  boundaries, 4)
bc15 = DirichletBC(V.sub(0), ue1,  boundaries, 5)
bc16 = DirichletBC(V.sub(0), ue1,  boundaries, 6)



bc21 = DirichletBC(V.sub(1), ue2,  boundaries, 1)
bc22 = DirichletBC(V.sub(1), ue2,  boundaries, 2)
bc23 = DirichletBC(V.sub(1), ue2,  boundaries, 3)
bc24 = DirichletBC(V.sub(1), ue2,  boundaries, 4)
bc25 = DirichletBC(V.sub(1), ue2,  boundaries, 5)
bc26 = DirichletBC(V.sub(1), ue2,  boundaries, 6)

bc = [bc11, bc12, bc13, bc14, bc15, bc16, bc21, bc22, bc23, bc24, bc25, bc26 ]
#Neumann bondary condition

#g_1=Expression('x[1]== 1 ? -3 : x[1]==0 ? 0 : 5', degree=2)
#g_2=Expression('x[1]== 1 ? -7 : x[1]==0 ? 0 : 5', degree=2)


########################################################################
domains = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
domains.set_all(0)
#g=Constant(0.0)
# Define variational problem
dx = Measure('dx', domain=mesh, subdomain_data=domains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)



# Define test functions
v_1, v_2 = TestFunctions(V)
u = interpolate( Expression( ('0.1', '0.1'), degree =0 ), V)
# Split system functions to access components
u_1, u_2 = split(u)

print(u_1)


# Define source terms
"""f_1 = Constant(1)
f_2 = Constant(3)"""

# Define expressions used in variational forms



# Define variational problem
F =  dot(q_1(u_1, u_2)*grad(u_1) + q_2(u_1, u_2)*grad(u_2) , grad(v_1))*dx(0) \
  + dot(q_3(u_1, u_2)*grad(u_1) + q_4(u_1, u_2)*grad(u_2) , grad(v_2))*dx(0) \
  - f_1*v_1*dx(0) - f_2*v_2*dx(0) 



solve(F == 0, u, bc, solver_parameters={"newton_solver": {"relative_tolerance": 1e-15}})
####################################################################################################
# Save solution to file (VTK)
# Create VTK files for visualization output
_u_1, _u_2 = u.split()
_u_1.rename("Solution Cation", "")
_u_2.rename("Solution Anion", "")
vtkfile_u_1 = File('Example of Thesis/Linear coupled/u_1.pvd')
vtkfile_u_2 = File('Example of Thesis/Linear coupled/u_2.pvd')
vtkfile_u_1 << (_u_1)
vtkfile_u_2 << (_u_2)
###################################################################################################

#Plotting 
import matplotlib.pyplot as plt
plot(mesh, title='Finite element mesh')
plot(_u_1, title='Solution Cation y_c')
plt.show()
plot(_u_2, title='Solution Cation y_a')
plt.show()
###################################################################################################
# Compute error in L2 norm
error_L2_u1 = errornorm(ue1, _u_1, 'L2') 
error_L2_u2 = errornorm(ue2, _u_2, 'L2') 
###################################################################################################
# Compute maximum error at vertices
vertex_values_u_D1 = ue1.compute_vertex_values(mesh)
vertex_values_u1 = _u_1.compute_vertex_values(mesh)
import numpy as np
error_max_u1 = np.max(np.abs(vertex_values_u_D1 - vertex_values_u1))

# Compute maximum error at vertices
vertex_values_u_D2 = ue2.compute_vertex_values(mesh)
vertex_values_u2 = _u_2.compute_vertex_values(mesh)
import numpy as np
error_max_u2 = np.max(np.abs(vertex_values_u_D2 - vertex_values_u2))
###################################################################################################
# Print errors
print('error_L2 in cation =', error_L2_u1)
print('error_L2 in anion =', error_L2_u2)
print('error_max in cation =', error_max_u1)
print('error_max in anion=', error_max_u2)

#####################################################################################
# Curve plot along x = 0.5 comparing p and w
tol = 0.001  # avoid hitting points outside the domain
y = np.linspace(0 + tol, 1 - tol, 101)
c1=0.8 #Constant value on x axis
points = [(c1, y_) for y_ in y]  # 2D points
u1_line = np.array([_u_1(point) for point in points])
u2_line = np.array([_u_2(point) for point in points])
plt.plot(y, u1_line, 'k', linewidth=2)  # magnify w
plt.plot(y, u2_line, 'b--', linewidth=2)
plt.grid(True)
plt.xlabel(f'y-aixs')
plt.ylabel(f'Value of solutions along x = {c1}')
plt.legend(['Cation', 'Anion'], loc='upper left')
plt.show()
########################################################################################
#####################################################################################
# Curve plot along x = 0.5 comparing p and w
tol = 0.001  # avoid hitting points outside the domain
x = np.linspace(0 + tol, 1 - tol, 101)
c2=0.6 #Constant value on x axis
points = [(x_,c2 ) for x_ in x]  # 2D points
u1_line = np.array([_u_1(point) for point in points])
u2_line = np.array([_u_2(point) for point in points])
plt.plot(y, u1_line, 'k', linewidth=2)  # magnify w
plt.plot(y, u2_line, 'b--', linewidth=2)
plt.grid(True)
plt.xlabel('$x$')
plt.ylabel(f'Value of solutions along y = {c2}')
plt.legend(['Cation', 'Anion'], loc='upper left')
plt.show()
