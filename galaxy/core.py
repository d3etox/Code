import numpy as np
import math

G = 1.0  # gravitational constant (scaled)

class Body:
    def __init__(self, pos, vel, mass, name=None):
        self.pos = np.array(pos, dtype=float)
        self.vel = np.array(vel, dtype=float)
        self.mass = float(mass)
        self.name = name if name else "body"

    def kinetic(self):
        return 0.5 * self.mass * np.dot(self.vel, self.vel)

class Simulation:
    def __init__(self, bodies, dt=0.01, softening=1e-2):
        self.bodies = list(bodies)
        self.dt = dt
        self.softening = softening

    def step_naive(self):
        """Naive O(N^2) force computation and symplectic Euler update."""
        n = len(self.bodies)
        positions = np.array([b.pos for b in self.bodies])
        masses = np.array([b.mass for b in self.bodies])
        forces = np.zeros((n,3))
        for i in range(n):
            for j in range(i+1, n):
                r = positions[j] - positions[i]
                dist2 = np.dot(r,r) + self.softening*self.softening
                invdist3 = 1.0 / (dist2 * math.sqrt(dist2))
                f = G * masses[i] * masses[j] * invdist3 * r
                forces[i] += f
                forces[j] -= f
        for i, b in enumerate(self.bodies):
            acc = forces[i] / b.mass
            b.vel += acc * self.dt
            b.pos += b.vel * self.dt

    def step_quadrant_approx(self):
        """
        Very simple approximate method: split bodies into 4 quadrants in xy-plane,
        compute each quadrant barycenter (mass, position), treat distant bodies
        by their quadrant barycenters. This is a crude approximation inspired by
        hierarchical methods (Barnes-Hut).
        """
        # partition
        qs = {0:[],1:[],2:[],3:[]}  # 4 quadrants
        for b in self.bodies:
            x,y = b.pos[0], b.pos[1]
            if x>=0 and y>=0:
                qs[0].append(b)
            elif x<0 and y>=0:
                qs[1].append(b)
            elif x<0 and y<0:
                qs[2].append(b)
            else:
                qs[3].append(b)
        bary = {}
        for k, lst in qs.items():
            if not lst:
                bary[k] = None
            else:
                total_m = sum(b.mass for b in lst)
                pos = sum((b.mass * b.pos for b in lst)) / total_m
                bary[k] = {'mass': total_m, 'pos': pos}
        # compute forces using barycenters for other quadrants
        n = len(self.bodies)
        forces = [np.zeros(3) for _ in range(n)]
        for i,b in enumerate(self.bodies):
            # intra-quadrant use naive pairwise for small lists
            # first add contributions from other quadrants via barycenters
            x,y = b.pos[0], b.pos[1]
            if x>=0 and y>=0:
                myq = 0
            elif x<0 and y>=0:
                myq = 1
            elif x<0 and y<0:
                myq = 2
            else:
                myq = 3
            # other quadrants as barycenters
            for k, info in bary.items():
                if info is None or k==myq:
                    continue
                r = info['pos'] - b.pos
                dist2 = np.dot(r,r) + self.softening*self.softening
                invdist3 = 1.0 / (dist2 * math.sqrt(dist2))
                f = G * b.mass * info['mass'] * invdist3 * r
                forces[i] += f
        # add intra-quadrant exact pairwise
        for k, lst in qs.items():
            L = len(lst)
            for a in range(L):
                for c in range(a+1,L):
                    bi = lst[a]; bj = lst[c]
                    r = bj.pos - bi.pos
                    dist2 = np.dot(r,r) + self.softening*self.softening
                    invdist3 = 1.0 / (dist2 * math.sqrt(dist2))
                    f = G * bi.mass * bj.mass * invdist3 * r
                    forces[self.bodies.index(bi)] += f
                    forces[self.bodies.index(bj)] -= f
        # update
        for i, b in enumerate(self.bodies):
            acc = forces[i] / b.mass
            b.vel += acc * self.dt
            b.pos += b.vel * self.dt

    def step(self, method='naive'):
        if method == 'naive':
            self.step_naive()
        else:
            self.step_quadrant_approx()

    def energy(self):
        ke = sum(b.kinetic() for b in self.bodies)
        pe = 0.0
        n = len(self.bodies)
        for i in range(n):
            for j in range(i+1,n):
                r = np.linalg.norm(self.bodies[i].pos - self.bodies[j].pos)
                pe -= G * self.bodies[i].mass * self.bodies[j].mass / (r + self.softening)
        return ke + pe
