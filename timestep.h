// note: u is weighted by fm
double timestep (const face vector u, double dtmax)
{
  static double previous = 0.;
  if (t == 0.) previous = 0.;
  dtmax /= CFL;
  foreach_face(reduction(min:dtmax))
    if (u.x[] != 0.) {  //跳过静止面（速度为0的地方）
      double dt = Delta/fabs(u.x[]); //计算局部 CFL 步长限制：Δt = Δx /u.x
      assert (fm.x[]);  //检查该面存在有效流体（体积分数不为0）
      dt *= fm.x[];  //按面上流体分数调整步长
      if (dt < dtmax) dtmax = dt;
    }
  dtmax *= CFL; //恢复实际时间步值： Δt=CFL×min⁡(Δx∣u∣)
  if (dtmax > previous)                    //
    dtmax = (previous + 0.1*dtmax)/1.1;    // 这里是一个时间步“缓升”机制：如果当前计算得到的 dtmax 比上一步更大，就不要立刻放大，而是逐渐增加一点
    previous = dtmax;                      //具体公式：dtnew​=(dtold​+0.1dtnew)​/1.1​
    return dtmax;                          //相当于平滑过渡，避免由于速度突变造成时间步突然变大而导致不稳定
}
在显式时间积分中，为了数值稳定性，必须满足 Courant–Friedrichs–Lewy (CFL) 条件：

CFL=uΔt/Δx≤1

即流体不能在一个时间步跨越超过一个网格。

因此时间步长应满足：
Δt≤CFL×Δx/∣u∣
