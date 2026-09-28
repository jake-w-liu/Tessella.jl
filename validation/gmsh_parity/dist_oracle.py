# Oracle values for model_distance differential vs gmsh 4.15.2 OCC kernel.
# Prints "key|d|x1 y1 z1 x2 y2 z2" per query; consumed by eyeball/Julia diff.
import gmsh

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)


def fresh():
    gmsh.model.add("m")
    gmsh.model.occ.removeAllDuplicates()


def q(key, d1, t1, d2, t2):
    r = gmsh.model.occ.getDistance(d1, t1, d2, t2)
    print("{}|{:.17g}|{:.17g} {:.17g} {:.17g}|{:.17g} {:.17g} {:.17g}".format(
        key, r[0], r[1], r[2], r[3], r[4], r[5], r[6]))


# pt-pt: (0,0,0),(3,4,0) -> 5
fresh()
a = gmsh.model.occ.addPoint(0, 0, 0)
b = gmsh.model.occ.addPoint(3, 4, 0)
gmsh.model.occ.synchronize()
q("pt-pt", 0, a, 0, b)

# pt-line: pt(1,4,0) vs x-axis segment (0,0,0)-(2,0,0) -> 4 at (1,0,0)
fresh()
p1 = gmsh.model.occ.addPoint(0, 0, 0)
p2 = gmsh.model.occ.addPoint(2, 0, 0)
ln = gmsh.model.occ.addLine(p1, p2)
pq = gmsh.model.occ.addPoint(1, 4, 0)
gmsh.model.occ.synchronize()
q("pt-line", 0, pq, 1, ln)

# pt-arc: circle arc quarter; query at center-ish
fresh()
c0 = gmsh.model.occ.addPoint(1, 0, 0)
cc = gmsh.model.occ.addPoint(0, 0, 0)
c1 = gmsh.model.occ.addPoint(0, 1, 0)
arc = gmsh.model.occ.addCircleArc(c0, cc, c1)
pq = gmsh.model.occ.addPoint(0, 0, 3)
gmsh.model.occ.synchronize()
q("pt-arc", 0, pq, 1, arc)

# pt-plane(surface): unit square z=0, pt above corner region -> project inside
fresh()
pa = gmsh.model.occ.addPoint(0, 0, 0)
pb = gmsh.model.occ.addPoint(1, 0, 0)
pc = gmsh.model.occ.addPoint(1, 1, 0)
pd = gmsh.model.occ.addPoint(0, 1, 0)
l1 = gmsh.model.occ.addLine(pa, pb)
l2 = gmsh.model.occ.addLine(pb, pc)
l3 = gmsh.model.occ.addLine(pc, pd)
l4 = gmsh.model.occ.addLine(pd, pa)
cl = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
sf = gmsh.model.occ.addPlaneSurface([cl])
pq = gmsh.model.occ.addPoint(0.5, 0.5, 2)
gmsh.model.occ.synchronize()
q("pt-plane", 0, pq, 2, sf)

# pt-box outside: box [0,1]^3, pt (3,0,0) -> 2
fresh()
bx = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
pq = gmsh.model.occ.addPoint(3, 0, 0)
gmsh.model.occ.synchronize()
q("pt-box-out", 0, pq, 3, bx)

# pt-box inside: pt (0.2,0.5,0.5) inside box -> 0
fresh()
bx = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
pq = gmsh.model.occ.addPoint(0.2, 0.5, 0.5)
gmsh.model.occ.synchronize()
q("pt-box-in", 0, pq, 3, bx)

# line-line skew: x-axis seg (0,0,0)-(1,0,0) vs vertical (0,1,3)-(0,1,4)
fresh()
a1 = gmsh.model.occ.addPoint(0, 0, 0)
a2 = gmsh.model.occ.addPoint(1, 0, 0)
b1 = gmsh.model.occ.addPoint(0, 1, 3)
b2 = gmsh.model.occ.addPoint(0, 1, 4)
la = gmsh.model.occ.addLine(a1, a2)
lb = gmsh.model.occ.addLine(b1, b2)
gmsh.model.occ.synchronize()
q("line-line-skew", 1, la, 1, lb)

# line-line crossing -> 0
fresh()
a1 = gmsh.model.occ.addPoint(-1, 0, 0)
a2 = gmsh.model.occ.addPoint(1, 0, 0)
b1 = gmsh.model.occ.addPoint(0, -1, 0)
b2 = gmsh.model.occ.addPoint(0, 1, 0)
la = gmsh.model.occ.addLine(a1, a2)
lb = gmsh.model.occ.addLine(b1, b2)
gmsh.model.occ.synchronize()
q("line-line-x", 1, la, 1, lb)

# line-plane: x-axis seg at z=2 vs square z=0 -> 2
fresh()
pa = gmsh.model.occ.addPoint(0, 0, 0)
pb = gmsh.model.occ.addPoint(1, 0, 0)
pc = gmsh.model.occ.addPoint(1, 1, 0)
pd = gmsh.model.occ.addPoint(0, 1, 0)
l1 = gmsh.model.occ.addLine(pa, pb)
l2 = gmsh.model.occ.addLine(pb, pc)
l3 = gmsh.model.occ.addLine(pc, pd)
l4 = gmsh.model.occ.addLine(pd, pa)
cl = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
sf = gmsh.model.occ.addPlaneSurface([cl])
q1 = gmsh.model.occ.addPoint(0.5, 0.5, 2)
q2 = gmsh.model.occ.addPoint(0.5, 0.6, 2)
lq = gmsh.model.occ.addLine(q1, q2)
gmsh.model.occ.synchronize()
q("line-plane", 1, lq, 2, sf)

# face-face parallel squares 2 apart -> 2
fresh()
def square(z):
    p1 = gmsh.model.occ.addPoint(0, 0, z)
    p2 = gmsh.model.occ.addPoint(1, 0, z)
    p3 = gmsh.model.occ.addPoint(1, 1, z)
    p4 = gmsh.model.occ.addPoint(0, 1, z)
    ls = [gmsh.model.occ.addLine(a, b)
          for a, b in ((p1, p2), (p2, p3), (p3, p4), (p4, p1))]
    return gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop(ls)])
s1 = square(0)
s2 = square(2)
gmsh.model.occ.synchronize()
q("face-face", 2, s1, 2, s2)

# box-box separated along x by 3: [0,1]^3 and [4,5]x[0,1]^2 -> 3
fresh()
b1 = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
b2 = gmsh.model.occ.addBox(4, 0, 0, 1, 1, 1)
gmsh.model.occ.synchronize()
q("box-box", 3, b1, 3, b2)

# box-sphere: box [0,1]^3, sphere c=(4,0.5,0.5) r=1 -> gap 2
fresh()
b1 = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
sp = gmsh.model.occ.addSphere(4, 0.5, 0.5, 1.0)
gmsh.model.occ.synchronize()
q("box-sphere", 3, b1, 3, sp)

# sphere-sphere: r1=r2=1 centers 3 apart -> 1
fresh()
s1 = gmsh.model.occ.addSphere(0, 0, 0, 1.0)
s2 = gmsh.model.occ.addSphere(3, 0, 0, 1.0)
gmsh.model.occ.synchronize()
q("sphere-sphere", 3, s1, 3, s2)

# overlapping sphere-box -> 0
fresh()
b1 = gmsh.model.occ.addBox(0, 0, 0, 2, 2, 2)
s1 = gmsh.model.occ.addSphere(1, 1, 3, 1.0)
gmsh.model.occ.synchronize()
q("box-sphere-ovl", 3, b1, 3, s1)

# cylinder: c=(0,0,0) axez h=3 r=1; box far
fresh()
cy = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, 3, 1.0)
b1 = gmsh.model.occ.addBox(4, 0, 0, 1, 1, 1)
gmsh.model.occ.synchronize()
q("cyl-box", 3, cy, 3, b1)

# torus: z-axis R=3 r=1; point on axis above hole -> nearest rim
fresh()
to = gmsh.model.occ.addTorus(0, 0, 0, 3.0, 1.0)
pq = gmsh.model.occ.addPoint(0, 0, 0)
gmsh.model.occ.synchronize()
q("pt-torus", 0, pq, 3, to)

# curve-face of same box: box edge 1 to box face 6 (should be 0 - shared)
fresh()
b1 = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
gmsh.model.occ.synchronize()
q("edge-face-self", 1, 1, 2, 6)

gmsh.finalize()
print("ORACLE-DONE")
