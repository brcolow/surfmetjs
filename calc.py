import sympy as sym

(xc, yc, vx, vy, r0, rr, x1, y1, x2, y2, t, s) = sym.symbols("xc yc vx vy r0 rr x1 y1 x2 y2 t s")
phi = sym.Symbol("phi", real=True)

eqx1 = xc + vx * t + (r0 + rr * t) * sym.cos(phi)
eqy1 = yc + vy * t + (r0 + rr * t) * sym.sin(phi)

eqx2 = x1 + (x2 - x1) * s
eqy2 = y1 + (y2 - y1) * s

sol = sym.solve([eqx1 - eqx2, eqy1 - eqy2], [t, s])
solS = sol[s]
print('s(phi) =', sym.together(solS))
solT = sol[t]
print('t(phi) =', sym.together(solT))

print(80 * '-')
print('Solutions for s(phi) = 0 ...')
sol_s0 = sym.solve(solS, phi)
print(sol_s0)

print('Solutions for s(phi) = 1 ...')
sol_s1 = sym.solve(solS - 1, phi)
print(sol_s1)

print(80 * '-')

print('The values of phi for the extremes of t(phi) ...')
dt = sym.diff(solT, phi)
sol_t = sym.solve(dt, phi)
print(sol_t)