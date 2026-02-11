#!/usr/bin/env python3
import math
import numpy as np
from .core import Sphere, Plane

def complex_scene(grid_size=8, spacing=1.2, wobble=0.2):
    """Create a more complex scene: grid of spheres with varied colors,
    a few large feature spheres, and a checker-like floor plane.
    """
    objects = []
    half = grid_size // 2
    for i in range(-half, half):
        for j in range(-half, half):
            x = i * spacing + (np.random.rand() - 0.5) * wobble
            z = j * spacing + (np.random.rand() - 0.5) * wobble - 6.0
            y = np.random.rand() * 0.4 - 0.2
            r = 0.25 + np.random.rand() * 0.2
            color = (np.random.rand()*0.8+0.2, np.random.rand()*0.8+0.2, np.random.rand()*0.8+0.2)
            spec = np.random.choice([5, 10, 50, 150])
            objects.append(Sphere([x, y, z], r, color=color, spec=spec))

    # Add a few larger spheres for focal points
    objects.append(Sphere([0.0, 0.0, -3.0], 1.0, color=(0.9, 0.3, 0.2), spec=200))
    objects.append(Sphere([-2.5, 0.2, -4.5], 0.8, color=(0.2, 0.6, 0.9), spec=50))
    objects.append(Sphere([2.5, 0.1, -5.0], 0.9, color=(0.8, 0.8, 0.2), spec=30))

    # Checker-like floor (plane)
    floor = Plane([0, -1.2, 0], [0, 1, 0], color=(0.6, 0.6, 0.6), spec=5)
    objects.append(floor)

    return objects
