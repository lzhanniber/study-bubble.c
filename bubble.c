/**
# Bubble rising in a large tank

We wish to study the behaviour of a single bubble rising "in a large
tank" i.e. far from any boundaries.

We use the centered Navier--Stokes solver and log performance
statistics. */

#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h" 

/**
We have two phases e.g. air and water. For large viscosity and density
ratios, the harmonic mean for the viscosity tends to work better than
the default arithmetic mean. We "overload" the default by defining the
*mu()* macro before including the code for two phases. */

#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))  //使用粘度的调和平均值代替默认的算术平均值
#include "two-phase.h"

/**
We also need surface tension, and in 3D only we will use the
$\lambda_2$ criterion of [Jeong and Hussain,
1995](/src/references.bib#jeong1995) to display the vortices using
Basilisk View. */

#include "tension.h"
#if dimension == 3
# include "lambda2.h"
#endif
#include "view.h"

/**
We can control the maximum runtime. */

#include "maxruntime.h"

/**
The density ratio is 1000 and the dynamic viscosity ratio 100. */

#define RHOR 1000.
#define MUR 100.

/**
We try to replicate the results of [Cano-Lozano et al,
2016](/src/references.bib#cano2016) (obtained with Gerris). Aside from
the ratios above, there are two independent parameters which can be
described by the Galilei number
$$
Ga^2 = \frac{g D^3}{\nu^2}     //Ga数表征重力驱动力 与 黏性阻力 的相对强度，控制气泡上升速度与流动稳定性
$$
with $g$ the acceleration of gravity, $D$ the diameter of the bubble
and $\nu$ the kinematic viscosity of the outer fluid; and the
[Bond/Eötvös](https://en.wikipedia.org/wiki/E%C3%B6tv%C3%B6s_number)
number
$$
Bo = \frac{\rho g D^2}{\sigma}  //Bo数表征重力（浮力）与表面张力之间的相对强弱，Bo数越大，气泡/液滴趋于保持球形，控制气泡的形变
$$
with $\rho$ the density of the outer fluid and $\sigma$ the surface
tension coefficient.

We consider two bubbles studied by Cano-Lozano et al, 2016. */

#if BUBBLE19
// Bubble 19 of Cano-Lozano et al, P.R.Fluids, 2016
const double Ga = 100.8;
const double Bo = 4.;
const double MAXTIME = 82;
#else
// Bubble 26 of Cano-Lozano et al, P.R.Fluids, 2016
const double Ga = 100.25;
const double Bo = 10.;
const double MAXTIME = 110;
#endif

/**
We choose as length unit the diameter of the bubble. The domain is
$120^3$. *Zi* is the initial position of the bubble relative to the
bottom wall. The acceleration of gravity is set to unity, which gives
a characteristic rise velocity also of order unity, which gives a
maximum time for the simulation comparable to the domain size. */

const double WIDTH = 120. [1];
const double Zi = 3.5;
int LEVEL = 12;               //定义初始条件（包括参数的定义、计算域、气泡的初始位置等）

/**
The main function can take two optional parameters: the maximum level
of adaptive refinement (as well as an optional maximum runtime). */

int main (int argc, char * argv[]) {
  maxruntime (&argc, argv);
  if (argc > 1)
    LEVEL = atoi (argv[1]);             //从命令行读取运行参数，设置模拟最大运行时间，并根据输入设置网格分辨率等级（LEVEL）
  
  /**
  We set the domain geometry and initial refinement. */
  
  size (WIDTH);
  origin (-L0/2, 0, -L0/2);
  init_grid (128);

  /**
  We set the physical parameters: densities, viscosities and surface
  tension. */
  
  rho1 = 1. [0];
  rho2 = rho1/RHOR;
  mu1 = 1./Ga;
  mu2 = 1./(MUR*Ga);
  f.sigma = 1./Bo;

  /**
  We reduce the tolerance on the divergence of the flow. This is
  important to minimise mass conservation errors for these simulations
  which are very long. */
  
  TOLERANCE = 1e-4 [*];
  run();
}

/**
For the initial conditions, we first try to restore the simulation
from a previous "restart", if this fails we refine the mesh locally to
the maximum level, in a sphere of diameter 1.5 around the bubble. We
then initialise the volume fraction for a bubble initially at (0,Zi,0)
of diameter unity. */

event init (t = 0) {
  if (!restore (file = "restart")) {
    refine (sq(x) + sq(y - Zi) + sq(z) - sq(0.75) < 0 && level < LEVEL);   //在气泡区域细化网格，使该区域的分辨率更高
    fraction (f, sq(x) + sq(y - Zi) + sq(z) - sq(.5));                     //根据给定的标量场φ(x) 定义体积分数函数 f,初始化两相流中界面的位置:若 φ < 0，则 f = 0（气泡内部）φ > 0，则 f = 1（外部流体）
  }
}

/**
We add the acceleration of gravity (unity) in the downward (-y)
direction. */

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 1.;     //遍历所有 y 方向的网格面，在 y 方向上施加 -1 的加速度（沿 负 y 方向、大小为 1 ）
}

/**
We adapt the mesh by controlling the error on the volume fraction and
velocity field. */

event adapt (i++) {
  double uemax = 1e-2;
  adapt_wavelet ({f,u}, (double[]){0.01,uemax,uemax,uemax}, LEVEL, 5);   //体积分数和速度场的变化来判断是否细化网格
}

/**
## Outputs

Every ten timesteps, we output the time, volume, position, and
velocity of the bubble. */

event logfile (i += 10) {
  double xb = 0., yb = 0., zb = 0., sb = 0.;
  double vbx = 0., vby = 0., vbz = 0.;
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb)
	  reduction(+:vbx) reduction(+:vby) reduction(+:vbz)
	  reduction(+:sb)) {               //计算气泡位置（质心坐标）、体积、速度
    double dv = (1. - f[])*dv();     
    xb += x*dv;
    yb += y*dv;
    zb += z*dv;
    vbx += u.x[]*dv;
    vby += u.y[]*dv;
    vbz += u.z[]*dv;
    sb += dv; 
  }
  fprintf (stderr,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb,
	   xb/sb, yb/sb, zb/sb,
	   vbx/sb, vby/sb, vbz/sb);
  fflush (stderr);
}

/**
Every time unit, we output a full snapshot of the simulation, to be
able to restart and for visualisation. In three dimensions, we compute
the value of the $\lambda_2$ field which will be used for
visualisation of vortices, as well as the streamwise vorticity
$\omega_y = \partial_x u_z - \partial_z u_x$. */

event snapshot (t = 1; t <= MAXTIME; t++)
{
  scalar l2[], omegay[];
#if dimension == 3
  lambda2 (u, l2);     //用于涡旋可视化的 λ2字段的值
  foreach()
    omegay[] = (u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);   //流向涡度 ωy=∂xuz−∂zux 的值
#endif
  
  char name[80];                              //定义字符串缓冲区，用来存放文件名
  sprintf (name, "dump-%03d", (int) t);      //生成格式化文件名，如 "dump-001"
  dump (file = name);                        //将当前模拟的全部数据保存到一个“快照”文件中
}

/**
We make a movie of the trailing wake and bubble shape. */

event movie (t = 21; t <= MAXTIME; t += 0.25)
{

#if BUBBLE19
  
  view (fov = 5.0278, quat = {-0.132839,0.513023,0.0748175,0.844727},
	tx = 0.00149469, ty = -0.355489, width = 300, height = 800);
  
  /**
  For $50 < t < 60$ we do a tracking shot from the initial camera
  position (above) to a new position, following and zooming on the
  bubble. */
  
  travelling (50, 60, fov = 2.07254, tx = 0.00524944, ty = -0.513744);

  /**
  For $60 < t < 82$ we just follow the bubble. */
  
  travelling (60, 82, tx = 0.00898601, ty = -0.703841);

#else // BUBBLE26

  /**
  The viewing parameters are different for bubble 26. */

  view (fov = 5.0278, quat = {-0.0487657,0.654023,0.0506598,0.7532},
	tx = -0.00353642, ty = -0.302285, width = 300, height = 800);

  travelling (50, 60, fov = 2.07254, tx = -0.0022909, ty = -0.402237);

  travelling (60, 110, tx = -0.00644264, ty = -0.741016);
#endif
  
  /**
  We use a different color for the bubble interface and for the
  isosurface. */
  
  clear();
  draw_vof ("f", fc = {0.13,0.47,0.77});
#if dimension == 3
  scalar l2[];
  lambda2 (u, l2);
  isosurface ("l2", -0.0002);
#endif

  save ("bubble.mp4");
}

/**
## Parallel runs

These simulations are expensive (in 3D), mostly because of the
timestep restriction due to surface tension and of the long evolution
required to reach an established quasi-stationary regime. The results
presented below were obtained using 96 cores on [occigen at
CINES](https://www.cines.fr/calcul/materiels/occigen/) in 12 hours for
each bubble. To run the simulation in 3D use

~~~bash
local% qcc -source -grid=octree -D_MPI=1 bubble.c
local% scp _bubble.c occigen.cines.fr:
~~~

and see the [isotropic turbulence
example](isotropic.c#running-with-mpi-on-occigen) for details.

## Results

The evolution of the bubble Reynols number for bubble 19 and 26 can be
compared to Figures 20.c and 20.f respectively of Cano-Lozano et al,
2016. The results are close, although the transition to established
regime seems to be faster in Basilisk. The final Reynolds numbers are
comparable (within a few percent) but the amplitude of oscillation
seems to be somewhat larger in Basilisk. Note that the Basilisk
results presented here use a different resolution for the bubble and
its wake than the Gerris results of Cano-Lozano et al, 2016.

~~~gnuplot Evolution of the Reynolds number
set grid
set xlabel 'Time'
set ylabel 'Reynolds'
set key bottom
plot [0:80]'bubble.19' u 1:($7*100.8) w l t 'bubble 19', \
           'bubble.26' u 1:($7*100.25) w l t 'bubble 26'
~~~

The trajectories of the center of gravity of the bubble reveal the two
regimes.

<table>
<tr>
<td>
~~~gnuplot
unset grid
set term svg size 480,640 font ",10"
set xlabel 'x'
set ylabel 'y'
set zlabel 'z'
set xyplane 0
splot [-0.7:0.7][-0.7:0.7]'bubble.19' u 3:5:4 w l t ''
~~~
</td>
<td>
~~~gnuplot
splot [-0.7:0.7][-0.7:0.7]'bubble.26' u 3:5:4 w l t ''
~~~
</td>
</tr>
<tr>
<td><center>Zig-zag trajectory for bubble 19</center></td>
<td><center>Spiralling trajectory for bubble 26</center></td>
</tr>
</table>

On the movies below, the bubble interface is in dark blue and the
white surface is the $\lambda_2$ isosurface showing vortical
structures.

<center>
<table>
<tr>
<td>
![Zig-zag regime. $Ga=100.8$, $Bo=4$.](bubble/bubble.19.mp4)(autoplay loop)
</td>
<td>
![Spiralling regime. $Ga=100.25$, $Bo=10$.](bubble/bubble.26.mp4)(autoplay loop)
</td>
</tr>
</table>
</center>
*/
