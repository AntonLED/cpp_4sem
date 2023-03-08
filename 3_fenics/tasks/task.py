from fenics import *

t_end = 10.0
dt = 0.1
k = 50
u_in = 20
u_out = -20

xml_file = "/home/anton/cpp2/3_fenics/tasks/ddiamond.xml"
mesh = Mesh(xml_file)
fd = MeshFunction('size_t', mesh, "/home/anton/cpp2/3_fenics/tasks/ddiamond_facet_region.xml")

V = FunctionSpace(mesh, 'P', 1)

bc1 = DirichletBC(V, Constant(u_in), fd, 5)
# bc2 = DirichletBC(V, Constant(u_out), fd, 6)
bc = [bc1]

u = TrialFunction(V)
v = TestFunction(V)
u_n = Function(V)

F = u*v*dx + dt*k*dot(grad(u), grad(v))*dx - u_n*v*dx
a, L = lhs(F), rhs(F)

u = Function(V)
t = 0
vtkfile = File('output/output.pvd')

num_steps = int(t_end/dt)
for n in range(num_steps):
    t += dt
    solve(a == L, u, bc)
    u_n.assign(u)
    vtkfile << (u, t)