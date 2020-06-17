# ultraSEM

![GitHub](https://img.shields.io/github/license/danfortunato/ultraSEM?)

ultraSEM is a software system for solving second-order linear partial differential equations (PDEs) on unstructured meshes in two dimension. The software is an implementation of the ultraspherical spectral element method [1], which allows hp-adaptivity to be efficiently performed in the high-p regime (where h is the mesh size and p is the polynomial order).

## Installation

To install ultraSEM, simply add the top-level directory of this repository to your MATLAB path. This can be done via the `pathtool` command or by executing the command `addpath(ultraSEMroot)`. To keep ultraSEM in your path for future MATLAB sessions, execute the `savepath` command.

## Getting started

Solving PDEs with ultraSEM is easy. Suppose we have the PDE `Lu = f`, where `L` is a partial differential operator (PDO) and `f` is a function. Given a domain `dom` and boundary conditions `bc`, the solution `u` can be approximated using degree-`p` polynomials on each element of `dom` by executing

```
S = ultraSEM(dom, L, f, p);
u = S \ bc;
```

The `ultraSEM` object `S` is a direct solver for the PDE `Lu = f`, which may be invoked for the given boundary conditions via the `\` command. The solution is returned as an `ultraSEM.Sol` object, which overloads a host of functions for plotting (e.g., `plot`, `contour`) and evaluation (e.g., `feval`, `norm`).

### Constructing a domain

A domain is represented in the ultraSEM system as an `ultraSEM.Domain` object. Convenient functions for constructing rectangles, quadrilaterals, triangles, and polygons are available via the commands `ultraSEM.rectangle`, `ultraSEM.quad`, `ultraSEM.triangle`, and `ultraSEM.polygon`, respectively. Elements can be combined to form larger domains by merging them with the `&` operator. More general meshes can be constructed using the `refine(dom)` method, which performs uniform h-refinement on a given domain `dom`, or the `refinePoint(dom, [x,y])` method, which performs adaptive h-refinement on `dom` around the point (x, y).

### Specifying a PDE

A PDO is specified by its coefficients for each derivative in the form  `{{dxx, dxy, dyy}, {dx, dy}, b}`, where each term `dxx`, `dxy`, `dyy`, `dx`, `dy`, and `b` can be a scalar (constant coefficient) or function handle (variable coefficient). The following syntax variations are permitted:
  * `{{dxx, dyy}, __, __}`, in which case `dxy = 0`.
  * `{dxx, __, __}`, in which case `dyy = dxx` and `dxy = 0`.
  * `{__, dx, __}`, in which case `dy = dx`.

The righthand side `f` and boundary conditions `bc` may be scalars or function handles.

## Example

Let's solve a simple problem with ultraSEM. First, construct a pentagonal domain:

```
dom = ultraSEM.polygon(5);
plot(dom)
```

<img src="https://www.dropbox.com/s/bvwg40if7w61ayp/Screenshot%202020-06-12%2015.50.17.png?raw=1" width="200">

Then, refine the elements of the domain to make a finer mesh:

```
dom = refine(dom);
plot(dom)
```

<img src="https://www.dropbox.com/s/5y1v1fr9fkgum65/Screenshot%202020-06-12%2015.50.32.png?raw=1" width="200">

Let's solve the constant-coefficient Helmholtz equation `lap(u) + 1000*u = -1` on this domain with zero Dirichlet boundary conditions:
```
L = {1, 0, 1000};
f = -1;
bc = 0;
```
Construct an `ultraSEM` object for this problem using 20th degree polynomials on each element of the mesh:

```
p = 20; 
S = ultraSEM(dom, L, f, p);
```

Solve the problem using the `\` command:

```
u = S \ bc;
plot(u)
```

<img src="https://www.dropbox.com/s/cpw3vpufnapjkj5/Screenshot%202020-06-12%2015.50.38.png?raw=1" width="200">

## References

[1] Daniel Fortunato, Nicholas Hale, and Alex Townsend, *The ultraspherical spectral element method*, https://arxiv.org/abs/2006.08756 (2020).
