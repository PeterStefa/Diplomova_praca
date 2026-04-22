import gurobipy as gp
from gurobipy import GRB
import numpy as np

# Parametre
VARIANTS = [5]
CUSTOMERS = 24
ALPHA1, ALPHA2 = 1.0, 1.0
DTL, SL, SR    = 20.0, 1.0, 1.0
GAMMA_P        = 5
K              = [0,1]

for VARIANT in VARIANTS:

    BASE = f"C:/Users/petos/OneDrive/Počítač/Diplomovka1/tspdl_input/TSPDL-MC25/Demand/{VARIANT}/"
    nodes    = np.loadtxt(BASE + f"nodes{VARIANT}_2L_with_demand.csv", delimiter=",")
    tau      = np.loadtxt(BASE + f"tau{VARIANT}_2L.csv",   delimiter=",")
    tauprime = np.loadtxt(BASE + f"tauprime{VARIANT}_2L.csv", delimiter=",")

    # Množiny
    C = list(range(1, CUSTOMERS + 1))
    P = list(range(CUSTOMERS + 1, CUSTOMERS + 3))
    mu = {c: int(nodes[c][3]) for c in C}
    lam   = {c: 1 for c in C}
    kappa = {(c, p): 0 for c in C for p in P}
    gam_p = {p: GAMMA_P for p in P}
    q_c   = {c: nodes[c][4] for c in C}
    l_cp  = {(c, p): tau[c][p] for c in C for p in P}
    Q = max(-(-sum(q_c.values()) // len(K)), max(q_c.values()) + 10)

    # graf
    S, T    = 's', 't'
    V_no_t  = [S] + C + P
    V_no_s  = C + P + [T]
    V_inner = C + P

    def midx(n):    return 0 if n in (S, T) else int(n)
    def t_ij(i, j): return tau[midx(i)][midx(j)]
    def d_ij(i, j): return tauprime[midx(i)][midx(j)]

    A = gp.tuplelist(
        (i, j) for i in V_no_t for j in V_no_s if i != j and not (i == S and j == T)
    )

    # Fikácia premenných (37-39)
    EL = DTL - SR
    forbidden_x      = {(i,j) for (i,j) in A if d_ij(i,j) > EL}
    forbidden_g_tr   = {(i,j) for (i,j) in A if t_ij(i,j) > EL}
    forbidden_g_sort = {(c,i,j) for c in C for (i,j) in A
                        if d_ij(i,c) + d_ij(c,j) > 2*EL - t_ij(i,j)}

    # Model nastavenie
    with gp.Model("VRP_DL") as m:
        m.Params.LazyConstraints = 1
        m.Params.TimeLimit = 3600

        # Premenné
        y     = m.addVars(K, A, vtype=GRB.BINARY, name="y")
        x     = m.addVars(((k,i,j) for k in K for (i,j) in A if (i,j) not in forbidden_x),
                          vtype=GRB.BINARY, name="x")
        g     = m.addVars(((k,c,i,j) for k in K for c in C for (i,j) in A
                           if (i,j) not in forbidden_g_tr and (c,i,j) not in forbidden_g_sort),
                          vtype=GRB.BINARY, name="g")
        theta = m.addVars(K, C, vtype=GRB.BINARY, name="th")
        omega = m.addVars(((k,c,i) for k in K for c in C for i in V_no_t if i != c),
                          vtype=GRB.BINARY, name="om")
        delta = m.addVars(((k,c,j) for k in K for c in C for j in V_no_s if j != c),
                          vtype=GRB.BINARY, name="de")
        phi   = m.addVars(K, C, vtype=GRB.BINARY, name="phi")
        z     = m.addVars(K, C, P, vtype=GRB.BINARY, name="z")
        sigma = m.addVars(K, C, vtype=GRB.CONTINUOUS, lb=0, name="sig")
        sigT  = m.addVars(K, vtype=GRB.CONTINUOUS, lb=0, name="sigT")

        def xv(k,i,j):   return x.get((k,i,j), 0)
        def gv(k,c,i,j): return g.get((k,c,i,j), 0)
        def omv(k,c,i):  return omega.get((k,c,i), 0)
        def dev(k,c,j):  return delta.get((k,c,j), 0)

        # Účelová funkcia
        m.setObjective(
            ALPHA1 * gp.quicksum(
                gp.quicksum(t_ij(i,j)*y[(k,i,j)] for (i,j) in A)
                + gp.quicksum(sigma[(k,c)] for c in C) + sigT[k]
                + gp.quicksum((SL+SR)*theta[(k,c)] for c in C)
                for k in K)
            + ALPHA2 * gp.quicksum(l_cp[(c,p)]*z[(k,c,p)]
                                for k in K for c in C for p in P),
            GRB.MINIMIZE)

        # Ohraničenia
# (2)
        for k in K:
            m.addConstr(gp.quicksum(y[(k,i,j)] for i,j in A.select(S, '*')) == 1)
            m.addConstr(gp.quicksum(y[(k,i,j)] for i,j in A.select('*', T)) == 1)
        # (3)
        m.addConstrs((gp.quicksum(y[(k,i,j)] for i,j in A.select(v, '*')) ==
                      gp.quicksum(y[(k,i,j)] for i,j in A.select('*', v))
                      for k in K for v in V_inner))
        m.addConstrs((gp.quicksum(y[(k,i,j)] for i,j in A.select(v, '*')) <= 1
                      for k in K for v in V_inner))
        # (5)
        m.addConstrs((gp.quicksum(gv(k,c,i,j) for i,j in A.select(S, '*')) == omv(k,c,S)
                      for k in K for c in C))
        # (6)
        m.addConstrs((gp.quicksum(gv(k,c,i,j) for i,j in A.select('*', T)) == dev(k,c,T)
                      for k in K for c in C))
        # (7)
        m.addConstrs((gp.quicksum(gv(k,c,i,j) for i,j in A.select(v, '*')) -
                      gp.quicksum(gv(k,c,i,j) for i,j in A.select('*', v)) ==
                      omv(k,c,v) - dev(k,c,v)
                      for k in K for c in C for v in V_inner))
        # (8)
        m.addConstrs((y[(k,i,j)] + xv(k,i,j) <= 1
                      for k in K for i,j in A.select(S, '*')))
        # (9)
        m.addConstrs((y[(k,i,j)] + xv(k,i,j) <= 1
                      for k in K for i,j in A.select('*', T)))
        # (10)
        m.addConstrs((gp.quicksum(
                          gp.quicksum(y[(k,i,j)] for i,j in A.select(c, '*')) + theta[(k,c)] + phi[(k,c)]
                          for k in K) == 1
                      for c in C))
        # (11)
        m.addConstrs((y[(k,i,j)] + z[(k,j,i)] + xv(k,i,j) + xv(k,j,i) <= 1
                      for k in K for p in P for i,j in A.select(p, '*') if j in C))
        # (12)
        m.addConstrs((y[(k,i,p)] + z[(k,i,p)] + xv(k,i,p) + xv(k,p,i) <= 1
                      for k in K for i in C for p in P if (i,p) in A))
        # (13)
        m.addConstrs((y[(k,i,j)] + xv(k,i,j) + xv(k,j,i) <= 1
                      for k in K for i in C for ii,j in A.select(i, '*') if j in C and i != j))
        # (14)
        m.addConstrs((z[(k,c,p)] <= kappa[(c,p)]
                      for k in K for c in C for p in P))
        # (15)
        m.addConstrs((gp.quicksum(y[(k,i,j)] for i,j in A.select('*', c)) <= lam[c]
                      for k in K for c in C))
        # (16)
        m.addConstrs((theta[(k,c)] <= 1 - mu[c]
                      for k in K for c in C))
        # (17)
        m.addConstrs((phi[(k,c)] == gp.quicksum(z[(k,c,p)] for p in P)
                      for k in K for c in C))
        # (18)
        m.addConstrs((gp.quicksum(z[(k,c,p)] for c in C) <=
                      gam_p[p] * gp.quicksum(y[(k,i,j)] for i,j in A.select(p, '*'))
                      for k in K for p in P))
        # (19)
        m.addConstrs((gp.quicksum(y[(k,i,j)] for i,j in A.select('*', p)) <=
                      gp.quicksum(z[(k,c,p)] for c in C)
                      for k in K for p in P))
        # (20)
        m.addConstrs((gp.quicksum(t_ij(i,j)*gv(k,c,i,j) for (i,j) in A) <= EL*theta[(k,c)]
                      for k in K for c in C))
        # (21)
        m.addConstrs((gp.quicksum(d_ij(i,c)*omv(k,c,i) for i,j in A.select('*', c) if i != c) +
                      gp.quicksum(d_ij(c,j)*dev(k,c,j) for i,j in A.select(c, '*') if j != c) <= EL*theta[(k,c)]
                      for k in K for c in C))
        # (22)
        m.addConstrs((sigma[(k,c)] >=
                      gp.quicksum(d_ij(i,c)*omv(k,c,i) for i,j in A.select('*', c) if i != c) +
                      gp.quicksum(d_ij(c,j)*dev(k,c,j) for i,j in A.select(c, '*') if j != c) -
                      gp.quicksum(t_ij(i,j)*gv(k,c,i,j) for (i,j) in A)
                      for k in K for c in C))
        # (23)
        m.addConstrs((gp.quicksum(g[(k,c,i,j)] for c in C if (k,c,i,j) in g) <= y[(k,i,j)]
                      for k in K for (i,j) in A if any((k,c,i,j) in g for c in C)))
        # (24)
        m.addConstrs((gp.quicksum(omv(k,c,i) for i,j in A.select('*', c) if i != c) == theta[(k,c)]
                      for k in K for c in C))
        m.addConstrs((gp.quicksum(dev(k,c,j) for i,j in A.select(c, '*') if j != c) == theta[(k,c)]
                      for k in K for c in C))
        # (25)
        m.addConstrs((xv(k,i,c) <= theta[(k,c)]
                      for k in K for c in C for i in P + [S] if (i,c) in A))
        # (26)
        m.addConstrs((xv(k,c,j) <= theta[(k,c)]
                      for k in K for c in C for j in P + [T] if (c,j) in A))
        # (27)
        m.addConstrs((xv(k,i,j) <= theta[(k,i)] + theta[(k,j)]
                      for k in K for i in C for ii,j in A.select(i, '*') if j in C and i != j))
        # (28)
        m.addConstrs((xv(k,i,c) <= omv(k,c,i)
                      for k in K for c in C for i in P + [S] if (i,c) in A))
        # (29)
        m.addConstrs((xv(k,h,j) <= dev(k,h,j)
                      for k in K for h in C for j in P + [T] if (h,j) in A))
        # (30)
        m.addConstrs((xv(k,i,j) <= omv(k,j,i) + dev(k,i,j)
                      for k in K for i in C for ii,j in A.select(i, '*') if j in C and i != j))
        # (31)
        m.addConstrs((gp.quicksum(xv(k,i,j) for i,j in A.select(i, '*')) ==
                      gp.quicksum(omega[(k,c,i)] for c in C if (k,c,i) in omega and c != i) + theta[(k,i)]
                      for k in K for i in C))
        m.addConstrs((gp.quicksum(omega[(k,c,i)] for c in C if (k,c,i) in omega and c != i) + theta[(k,i)] <= 1
                      for k in K for i in C))
        # (32)
        m.addConstrs((gp.quicksum(xv(k,p,j) for i,j in A.select(p, '*')) ==
                      gp.quicksum(omega[(k,c,p)] for c in C if (k,c,p) in omega)
                      for k in K for p in P))
        m.addConstrs((gp.quicksum(omega[(k,c,p)] for c in C if (k,c,p) in omega) <= 1
                      for k in K for p in P))
        # (33)
        m.addConstrs((gp.quicksum(xv(k,i,j) for i,j in A.select('*', j)) ==
                      gp.quicksum(delta[(k,c,j)] for c in C if (k,c,j) in delta and c != j) + theta[(k,j)]
                      for k in K for j in C))
        m.addConstrs((gp.quicksum(delta[(k,c,j)] for c in C if (k,c,j) in delta and c != j) + theta[(k,j)] <= 1
                      for k in K for j in C))
        # (34)
        m.addConstrs((gp.quicksum(xv(k,i,p) for i,j in A.select('*', p)) ==
                      gp.quicksum(delta[(k,c,p)] for c in C if (k,c,p) in delta)
                      for k in K for p in P))
        m.addConstrs((gp.quicksum(delta[(k,c,p)] for c in C if (k,c,p) in delta) <= 1
                      for k in K for p in P))
        # (35)
        m.addConstrs((gp.quicksum(z[(k,c,p)] for k in K for c in C) <= gam_p[p]
                      for p in P))
        # (36)
        m.addConstrs((gp.quicksum(q_c[c] * (gp.quicksum(y[(k,i,j)] for i,j in A.select(c, '*'))
                                  + theta[(k,c)] + phi[(k,c)]) for c in C) <= Q
                      for k in K))

        # -- Lazy subtour callback
        def subtour_callback(model, where):
            if where != GRB.Callback.MIPSOL:
                return
            for k in K:
                y_vars = [y[(k,i,j)] for (i,j) in A]
                y_vals = model.cbGetSolution(y_vars)
                y_val  = {(i,j): v for (i,j), v in zip(A, y_vals)}

                succ = {}
                for (i,j), val in y_val.items():
                    if val > 0.5:
                        succ.setdefault(i, []).append(j)
                visited = set()
                stack = [S]
                while stack:
                    node = stack.pop()
                    if node in visited: continue
                    visited.add(node)
                    for nb in succ.get(node, []):
                        stack.append(nb)

                unreachable = [v for v in V_inner if v not in visited]
                if not unreachable:
                    continue

                remaining = set(unreachable)
                while remaining:
                    start = next(iter(remaining))
                    comp = set()
                    bfs = [start]
                    while bfs:
                        node = bfs.pop()
                        if node in comp: continue
                        comp.add(node)
                        for nb in succ.get(node, []):
                            if nb in remaining:
                                bfs.append(nb)
                    remaining -= comp

                    sub = list(comp)
                    internal = [(i,j) for i in sub for j in sub if i != j and (i,j) in A]
                    if not internal:
                        continue
                    custs_in_sub = [c for c in sub if c in C]
                    lkrs_in_sub  = [p for p in sub if p in P]
                    if not custs_in_sub and not lkrs_in_sub:
                        continue
                    q = custs_in_sub[0] if custs_in_sub else lkrs_in_sub[0]
                    lhs = gp.quicksum(y[(k,i,j)] for (i,j) in internal)
                    rhs = (gp.quicksum(1 - theta[(k,c)] - phi[(k,c)]
                                        for c in custs_in_sub if c != q)
                        + len([p for p in lkrs_in_sub if p != q]))
                    model.cbLazy(lhs <= rhs)
        # Výsledky
        m.optimize(subtour_callback)

        for k in K:
            route = [S]
            while route[-1] != T:
                route.append(next(j for i,j in A.select(route[-1], '*') if y[(k,route[-1],j)].X > 0.5))
            route = [0 if n in (S,T) else n for n in route]
            print(f"Vehicle {k}: {route}")

        drones  = {k: [c for c in C if any(z[(k,c,p)].X > 0.5 for p in P)] for k in K}
        lockers = {p: [c for c in C for k in K if z[(k,c,p)].X > 0.5] for p in P}
        print(f"Drones: {[drones[k] for k in K]}")
        print(f"Lockers: { {p: lockers[p] for p in P if lockers[p]} }")
