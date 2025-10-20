/**
# Incompressible Navier--Stokes solver (centered formulation)

We wish to approximate numerically the incompressible,
variable-density Navier--Stokes equations
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(2\mu\mathbf{D})\right] + 
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u} = 0
$$
with the deformation tensor 
$\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$.

The scheme implemented here is close to that used in Gerris ([Popinet,
2003](/src/references.bib#popinet2003), [Popinet,
2009](/src/references.bib#popinet2009), [Lagrée et al,
2011](/src/references.bib#lagree2011)).

We will use the generic time loop, a CFL-limited timestep, the
Bell-Collela-Glaz advection scheme and the implicit viscosity
solver. If embedded boundaries are used, a different scheme is used
for viscosity. */

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#if EMBED
# include "viscosity-embed.h"
#else
# include "viscosity.h"
#endif

/**
The primary variables are the centered pressure field $p$ and the
centered velocity field $\mathbf{u}$. The centered vector field
$\mathbf{g}$ will contain pressure gradients and acceleration terms.

We will also need an auxilliary face velocity field $\mathbf{u}_f$ and
the associated centered pressure field $p_f$. */

scalar p[];
vector u[], g[];
scalar pf[];
face vector uf[];

/**
In the case of variable density, the user will need to define both the
face and centered specific volume fields ($\alpha$ and $\alpha_c$
respectively) i.e. $1/\rho$. If not specified by the user, these         
fields are set to one i.e. the density is unity.

Viscosity is set by defining the face dynamic viscosity $\mu$; default
is zero.

The face field $\mathbf{a}$ defines the acceleration term; default is
zero.

The statistics for the (multigrid) solution of the pressure Poisson
problems and implicit viscosity are stored in *mgp*, *mgpf*, *mgu*
respectively. 

If *stokes* is set to *true*, the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$ is omitted. This is a
reference to [Stokes flows](http://en.wikipedia.org/wiki/Stokes_flow)
for which inertia is negligible compared to viscosity. */

(const) face vector mu = zerof, a = zerof, alpha = unityf;
(const) scalar rho = unity;
mgstats mgp = {0}, mgpf = {0}, mgu = {0};
bool stokes = false;

/**
## Boundary conditions

For the default symmetric boundary conditions, we need to ensure that
the normal component of the velocity is zero after projection. This
means that, at the boundary, the acceleration $\mathbf{a}$ must be
balanced by the pressure gradient. Taking care of boundary orientation
and staggering of $\mathbf{a}$, this can be written */

#if EMBED
# define neumann_pressure(i) (alpha.n[i] ? a.n[i]*fm.n[i]/alpha.n[i] :	//alpha.n[i]定义为比体积即1/ρ alpha是面比体积场 alphac是中心比体积场
			      a.n[i]*rho[]/(cm[] + SEPS))
#else
# define neumann_pressure(i) (a.n[i]*fm.n[i]/alpha.n[i])  //定义为∂p/∂n=a/alpha即∂p/∂n=ρa，通过设置压力梯度来抵消加速度，使边界法向速度保持零
#endif

p[right] = neumann (neumann_pressure(ghost));   //使用右边界的外侧面，即ghost面，通过设置压力梯度来抵消加速度，使边界法向速度保持零
p[left]  = neumann (- neumann_pressure(0));    //使用左边界的内侧面，即网格内部面，通过设置压力梯度来抵消加速度，使边界法向速度保持零

#if AXI
uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                             // is zero on the axis of symmetry
p[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
p[top]    = neumann (neumann_pressure(ghost));
p[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
p[front]  = neumann (neumann_pressure(ghost));
p[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI

/**
For [embedded boundaries on trees](/src/embed-tree.h), we need to
define the pressure gradient for prolongation of pressure close to
embedded boundaries. */

#if TREE && EMBED
void pressure_embed_gradient (Point point, scalar p, coord * g)
{
  foreach_dimension()
    g->x = rho[]/(cm[] + SEPS)*(a.x[] + a.x[1])/2.;
}
#endif // TREE && EMBED

/**
## Initial conditions */

event defaults (i = 0)
{

  /**
  We reset the multigrid parameters to their default values. */
  
  mgp = (mgstats){0};
  mgpf = (mgstats){0};
  mgu = (mgstats){0};  
  
  CFL = 0.8;

  /**
  The pressures are never dumped. */

  p.nodump = pf.nodump = true;
  
  /**
  The default density field is set to unity (times the metric and the
  solid factors). */

  if (alpha.x.i == unityf.x.i) {     //默认单相流、均匀密度流体
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {   //如果 α（比体积）不是常数（即流体密度随空间变化），则我们需要更新它，使它在边界或截断面上乘以局部面系数 fm.x[]
    face vector alphav = alpha;
    foreach_face()
      alphav.x[] = fm.x[];
  }

  /**
  On trees, refinement of the face-centered velocity field needs to
  preserve the divergence-free condition. */

#if TREE
  uf.x.refine = refine_face_solenoidal;

  /**
  When using [embedded boundaries](/src/embed.h), the restriction and
  prolongation operators need to take the boundary into account. */

#if EMBED
  uf.x.refine = refine_face;
  foreach_dimension()
    uf.x.prolongation = refine_embed_face_x;
  for (scalar s in {p, pf, u, g}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
    s.depends = list_add (s.depends, cs);
  }
  for (scalar s in {p, pf})
    s.embed_gradient = pressure_embed_gradient;
#endif // EMBED
#endif // TREE

  /**
  We set the dimensions of the velocity field. */

  foreach()
    foreach_dimension()
      dimensional (u.x[] == Delta/t);    //设置速度场 u 的物理量纲等价于 长度/时间
}


/**
We had some objects to display by default. */

event default_display (i = 0)
  display ("squares (color = 'u.x', spread = -1);");  //bview 绘图命令字符串, squares: 绘制方格图（每个网格单元显示为一个方块）color = 'u.x':用颜色显示速度分量 u.x spread = -1:平滑或插值参数（-1 表示启用基于父网格的平滑显示）

/**
After user initialisation, we initialise the face velocity and fluid
properties. */

double dtmax;

event init (i = 0)
{
  trash ({uf}); //跟踪uf
  foreach_face()
    uf.x[] = fm.x[]*face_value (u.x, 0); //u 中心速度 uf 面速度 把单元中心速度 u 插值到面上

  /**
  We update fluid properties. */

  event ("properties");   //执行event properties，更新流体特性

  /**
  We set the initial timestep (this is useful only when restoring from
  a previous run). */

  dtmax = DT;
  event ("stability");
}

/**
## Time integration

The timestep for this iteration is controlled by the CFL condition,
applied to the face centered velocity field $\mathbf{u}_f$; and the
timing of upcoming events. */

event set_dtmax (i++,last) dtmax = DT;   //在每个时间步把 dtmax 设置为 预设的最大时间步 DT

event stability (i++,last) {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));  //设置时间步长dt,若为Stokes 方程，直接使用最大时间步。 否则CFL 条件和面速度 uf 计算安全时间步
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. Note that this assumes that tracer fields
are defined at time $t-\Delta t/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */

event vof (i++,last);
event tracer_advection (i++,last);
event tracer_diffusion (i++,last);

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++,last);

/**
### Predicted face velocity field                //预测面速度场

For second-order in time integration of the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$, we need to define the face
velocity field $\mathbf{u}_f$ at time $t+\Delta t/2$. We use a version
of the Bell-Collela-Glaz [advection scheme](/src/bcg.h) and the
pressure gradient and acceleration terms at time $t$ (stored in vector
$\mathbf{g}$). */

void prediction()  
{
  vector du;
  foreach_dimension() {
    scalar s = new scalar;
    du.x = s;       
  }

  if (u.x.gradient)                 //1 定义了自定义梯度
    foreach()
      foreach_dimension() {
#if EMBED         
        if (!fs.x[] || !fs.x[1]) //fs.x[] 检查当前面是否属于流体，如果面上没有流体，梯度设为 0
	  du.x[] = 0.;
	else
#endif
	  du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1])/Delta;
      }              //1 自定义的梯度函数(u.x.gradient)， 函数形式为 gradient(left, center, right)，返回该单元 x 方向的速度梯度
  else             //2 未定义自定义梯度
    foreach()
      foreach_dimension() {
#if EMBED
        if (!fs.x[] || !fs.x[1])
	  du.x[] = 0.;
	else
#endif
	  du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);  //用二阶中心差分：∂u∂x≈ui+1−ui−1/2Δ
    }        //2

  trash ({uf});
  foreach_face() {             //计算面速度场uf
    double un = dt*(u.x[] + u.x[-1])/(2.*Delta), s = sign(un);   //计算当前面 x 方向的无量纲速度，s = sign(un)：上风方向，s = 1 → 正向流，i = -1，取左侧单元的速度。s = -1 → 反向流，i = 0，取有侧单元的速度
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i] + (g.x[] + g.x[-1])*dt/4. + s*(1. - s*un)*du.x[i]*Delta/2.;  //更新x 方向面速度
    #if dimension > 1
    if (fm.y[i,0] && fm.y[i,1]) {
      double fyy = u.y[i] < 0. ? u.x[i,1] - u.x[i] : u.x[i] - u.x[i,-1];
      uf.x[] -= dt*u.y[i]*fyy/(2.*Delta);          //y 方向对 x 方向速度的对流修正
    }
    #endif
    #if dimension > 2
    if (fm.z[i,0,0] && fm.z[i,0,1]) {
      double fzz = u.z[i] < 0. ? u.x[i,0,1] - u.x[i] : u.x[i] - u.x[i,0,-1];
      uf.x[] -= dt*u.z[i]*fzz/(2.*Delta);             //计算 z 方向对 x 方向速度的修正
    }
    #endif
    uf.x[] *= fm.x[];         //按流体比例缩放速度
  }

  delete ((scalar *){du});
}

/**
### Advection term

We predict the face velocity field $\mathbf{u}_f$ at time $t+\Delta
t/2$ then project it to make it divergence-free. We can then use it to
compute the velocity advection term, using the standard
Bell-Collela-Glaz advection scheme for each component of the velocity
field. */

event advection_term (i++,last)      //对流项计算
{
  if (!stokes) {                    //若为斯托克斯流，跳过计算对流项
    prediction();                   //中间速度预测
    mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);    //半步压力投影 解如下方程：∇⋅(1/ρ∇p)=∇⋅uf/Δt，修正 uf 使其无散度，并更新压力修正场 pf
    advection ((scalar *){u}, uf, dt, (scalar *){g});      //执行真正的速度对流项更新，即离散化：∂u/∂t+∇⋅(u⊗u)=g，显式推进对流项un+1=un−Δt∇⋅(uf​u)+Δtg​,根据面速度计算单元速度变化
/**
### Viscous term

We first define a function which adds the pressure gradient and
acceleration terms. */

static void correction (double dt)    //粘性项  static 定义为静态函数，只在当前源文件中可见
{
  foreach()
    foreach_dimension()
      u.x[] += dt*g.x[];   //在当前时间步 dt 内，给速度场 u 添加重力加速度的影响
}

/**
The viscous term is computed implicitly. We first add the pressure
gradient and acceleration terms, as computed at time $t$, then call
the implicit viscosity solver. We then remove the acceleration and
pressure gradient terms as they will be replaced by their values at
time $t+\Delta t$. */

event viscous_term (i++,last)
{
  if (constant(mu.x) != 0.) {         //检查粘度是否为0，若黏度为 0，流动是无黏性的，就跳过
    correction (dt);             //调用函数correction，将重力加速度加到速度场 上
    mgu = viscosity (u, mu, rho, dt, mgu.nrelax);  //调用 Basilisk 内置函数 viscosity()，隐式求解黏性扩散项：∂u/∂t=1/ρ∇⋅(μ∇u)
    correction (-dt);        //将之前加上的重力修正撤回，这样做到了只在黏性求解“临时”考虑了重力，但不改变速度的物理值
  }

  /**
  We reset the acceleration field (if it is not a constant). */

  if (!is_constant(a.x)) {    //检查 a.x（x 方向的加速度分量）是否恒定
    face vector af = a;  //创建一个 面向量 af，引用加速度向量 a
    trash ({af});  //清理 缓存或临时数据。保证 af 的数据干净，不受之前时间步残留值影响
    foreach_face()
      af.x[] = 0.;  //把 x 方向的加速度 af.x[] 置为 0
  }
}

/**
### Acceleration term

The acceleration term $\mathbf{a}$ needs careful treatment as many
equilibrium solutions depend on exact balance between the acceleration
term and the pressure gradient: for example Laplace's balance for
surface tension or hydrostatic pressure in the presence of gravity.

To ensure a consistent discretisation, the acceleration term is
defined on faces as are pressure gradients and the centered combined
acceleration and pressure gradient term $\mathbf{g}$ is obtained by
averaging. 

The (provisionary) face velocity field at time $t+\Delta t$ is
obtained by interpolation from the centered velocity field. The
acceleration term is added. */

event acceleration (i++,last)
{
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*(face_value (u.x, 0) + dt*a.x[]);
}

/**
## Approximate projection

This function constructs the centered pressure gradient and
acceleration field *g* using the face-centered acceleration field *a*
and the cell-centered pressure field *p*. */

void centered_gradient (scalar p, vector g)
{

  /**
  We first compute a face field $\mathbf{g}_f$ combining both
  acceleration and pressure gradient. */

  face vector gf[];
  foreach_face()
    gf.x[] = fm.x[]*a.x[] - alpha.x[]*(p[] - p[-1])/Delta;

  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({g});
  foreach()
    foreach_dimension()
      g.x[] = (gf.x[] + gf.x[1])/(fm.x[] + fm.x[1] + SEPS);
}

/**
To get the pressure field at time $t + \Delta t$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). Then compute the centered gradient field *g*. */

event projection (i++,last)
{
  mgp = project (uf, p, alpha, dt, mgp.nrelax);
  centered_gradient (p, g);

  /**
  We add the gradient field *g* to the centered velocity field. */

  correction (dt);
}

/**
Some derived solvers need to hook themselves at the end of the
timestep. */

event end_timestep (i++, last);

/**
## Adaptivity

After mesh adaptation fluid properties need to be updated. When using
[embedded boundaries](/src/embed.h) the fluid fractions and face
fluxes need to be checked for inconsistencies. */

#if TREE
event adapt (i++,last) {
#if EMBED
  fractions_cleanup (cs, fs);
  foreach_face()
    if (uf.x[] && !fs.x[])
      uf.x[] = 0.;
#endif
  event ("properties");
}
#endif

/**
## See also

* [Double projection](double-projection.h)
* [Performance monitoring](perfs.h)
*/
