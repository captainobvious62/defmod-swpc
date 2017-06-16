from pylab    import *
from fenics   import *
from scipy.sparse            import spdiags
from pylab                   import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors       import from_levels_and_colors
from dolfin import *
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['legend.fontsize'] = 'medium'

def plot_matrix(M, ax, title, continuous=False, cmap='Greys'):
  """
  plot a matrix <M> with title <title> and a colorbar on subplot (axes object) 
  <ax>.
  """
  M    = array(M)
  m,n  = shape(M)
  M    = M.round(decimals=9)
  cmap = cm.get_cmap(cmap)
  if not continuous:
    unq  = unique(M)
    num  = len(unq)
  im      = ax.imshow(M, cmap=cmap, interpolation='None')
  divider = make_axes_locatable(ax)
  cax     = divider.append_axes("right", size="5%", pad=0.05)
  dim     = r'$%i \times %i$ ' % (m,n)
  ax.set_title(dim + title)
  ax.axis('off')
  cb = colorbar(im, cax=cax)
  if not continuous:
    cb.set_ticks(unq)
    cb.set_ticklabels(unq)

n    = 2
mesh = UnitSquareMesh(n,n)
#mesh = Mesh('meshes/unit_square_mesh.xml')
#mesh = Mesh('meshes/triangle.xml')

# refine mesh :
origin = Point(0.0,0.5,0.0)
#for i in range(1,8):
#  cell_markers = CellFunction("bool", mesh)
#  cell_markers.set_all(False)
#  for cell in cells(mesh):
#    p = cell.midpoint()
#    if p.distance(origin) < 1.0/i:
#      cell_markers[cell] = True
#  mesh = refine(mesh, cell_markers)

## refine mesh :
#for i in range(1,1):
#  cell_markers = CellFunction("bool", mesh)
#  cell_markers.set_all(False)
#  for cell in cells(mesh):
#    p = cell.midpoint()
#    if p.y() <= 1.0/i:
#      cell_markers[cell] = True
#  mesh = refine(mesh, cell_markers)

Q = FunctionSpace(mesh, 'CG', 1)
N = FacetNormal(mesh)
A = FacetArea(mesh)

w = TrialFunction(Q)
v = TestFunction(Q)
u = Function(Q, name='u')

f = Constant(1.0)

a = inner(grad(w), grad(v)) * dx
l = f * v * dx

def left(x, on_boundary):
  return x[0] == 0 and on_boundary

def right(x, on_boundary):
  return x[0] == 1 and on_boundary

def top(x, on_boundary):
  return x[1] == 1 and on_boundary

def bottom(x, on_boundary):
  return x[1] == 0 and on_boundary

bcl = DirichletBC(Q, 0.0, left)
bcr = DirichletBC(Q, 0.0, right)
bct = DirichletBC(Q, 0.0, top)
bcb = DirichletBC(Q, 0.0, bottom)

bc  = [bcl]

solve(a == l, u, bc)

File('output/u.pvd') << u

uv = u.vector().array()

s  = assemble(dot(grad(u), N) * v * ds)
M  = assemble(w * v * dx)
b  = assemble(l)
K  = assemble(a)

t  = (np.dot(K.array(),uv) - b.array())

q = Function(Q, name='q')
q.vector().set_local(t)
q.vector().apply('insert')

File('output/q.pvd') << q

fx = Function(Q, name='fx')
solve(M, fx.vector(), s)
File('output/fx.pvd') << fx

fig = figure(figsize=(10,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

plot_matrix(M.array(), ax1, r'mass matrix $M$',      continuous=True)
plot_matrix(K.array(), ax2, r'stiffness matrix $K$', continuous=False)

tight_layout()
show()

