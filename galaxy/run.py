#!/usr/bin/env python3
import random, csv, os, math
import numpy as np
import matplotlib.pyplot as plt
from galaxy.core import Body, Simulation

def random_galaxy(n, radius=50.0, mass_scale=5.0):
    bodies = []
    for i in range(n):
        r = radius * (random.random()**0.5)
        theta = random.random() * 2 * math.pi
        x = r * math.cos(theta)
        y = r * math.sin(theta)
        z = (random.random()-0.5) * radius * 0.05
        speed = math.sqrt(1.0 * mass_scale / (r+1e-3))
        vx = -speed * math.sin(theta)
        vy = speed * math.cos(theta)
        vz = 0.0
        m = random.random()*mass_scale + 0.1
        bodies.append(Body([x,y,z],[vx,vy,vz],m, name=f"body{i}"))
    return bodies

if __name__ == '__main__':
    n = 200
    bodies = random_galaxy(n)
    sim = Simulation(bodies, dt=0.05, softening=0.5)
    steps = 400
    outdir = '/Users/charlesfabre/Desktop/Code/galaxy/out'
    os.makedirs(outdir, exist_ok=True)
    for s in range(steps):
        # use quadrant approx for speed on larger N
        sim.step(method='quadrant')
        if s%50==0:
            print(f"step {s}, energy {sim.energy():.6f}")
    # save CSV
    with open(os.path.join(outdir,'positions.csv'),'w',newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['name','x','y','z','vx','vy','vz','mass'])
        for b in sim.bodies:
            writer.writerow([b.name]+list(b.pos)+list(b.vel)+[b.mass])
    # save scatter XY
    xs = [b.pos[0] for b in sim.bodies]
    ys = [b.pos[1] for b in sim.bodies]
    plt.figure(figsize=(6,6))
    plt.scatter(xs, ys)
    plt.title('Galaxy simulation (projection xy)')
    plt.xlabel('x'); plt.ylabel('y')
    plt.savefig(os.path.join(outdir,'galaxy_xy.png'), dpi=150)
    print('Saved outputs in', outdir)
