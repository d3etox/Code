"""Numba-accelerated renderer for the raytracer.

This module converts the scene to Structure-of-Arrays (SoA) and uses Numba
to JIT-compile the hot inner loop. If Numba is not available, the module
raises ImportError so callers can fall back to other renderers.
"""
from typing import List
import numpy as np

try:
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except Exception:
    NUMBA_AVAILABLE = False


def scene_to_soa(scene_objects):
    """Convert scene objects (Sphere, Plane) into contiguous numpy arrays.

    Returns a dict with sphere arrays and plane arrays.
    """
    spheres = [o for o in scene_objects if o.__class__.__name__ == 'Sphere']
    planes = [o for o in scene_objects if o.__class__.__name__ == 'Plane']

    so = {}
    if spheres:
        centers = np.array([s.center for s in spheres], dtype=np.float32)
        radii = np.array([s.radius for s in spheres], dtype=np.float32)
        colors = np.array([s.color for s in spheres], dtype=np.float32)
        specs = np.array([getattr(s, 'spec', 10) for s in spheres], dtype=np.float32)
    else:
        centers = np.zeros((0, 3), dtype=np.float32)
        radii = np.zeros((0,), dtype=np.float32)
        colors = np.zeros((0, 3), dtype=np.float32)
        specs = np.zeros((0,), dtype=np.float32)

    if planes:
        p_points = np.array([p.point for p in planes], dtype=np.float32)
        p_normals = np.array([p.normal for p in planes], dtype=np.float32)
        p_colors = np.array([p.color for p in planes], dtype=np.float32)
        p_specs = np.array([getattr(p, 'spec', 10) for p in planes], dtype=np.float32)
    else:
        p_points = np.zeros((0, 3), dtype=np.float32)
        p_normals = np.zeros((0, 3), dtype=np.float32)
        p_colors = np.zeros((0, 3), dtype=np.float32)
        p_specs = np.zeros((0,), dtype=np.float32)

    so['centers'] = centers
    so['radii'] = radii
    so['colors'] = colors
    so['specs'] = specs
    so['p_points'] = p_points
    so['p_normals'] = p_normals
    so['p_colors'] = p_colors
    so['p_specs'] = p_specs
    return so


if NUMBA_AVAILABLE:
    @njit(parallel=True, fastmath=True)
    def _render_numba(width, height,
                      centers, radii, colors, specs,
                      p_points, p_normals, p_colors, p_specs,
                      cam_pos, cam_right, cam_up, cam_forward, fov):
        out = np.zeros((height, width, 3), dtype=np.uint8)
        aspect = width / height
        scale = np.tan(np.radians(fov * 0.5))

        ns = centers.shape[0]
        np_ = p_points.shape[0]

        for j in prange(height):
            for i in range(width):
                x = (2.0 * (i + 0.5) / width - 1.0) * scale * aspect
                y = (1.0 - 2.0 * (j + 0.5) / height) * scale
                rx = cam_right[0] * x + cam_up[0] * y + cam_forward[0]
                ry = cam_right[1] * x + cam_up[1] * y + cam_forward[1]
                rz = cam_right[2] * x + cam_up[2] * y + cam_forward[2]
                # normalize ray
                rl = (rx*rx + ry*ry + rz*rz)**0.5
                rx /= rl; ry /= rl; rz /= rl

                best_t = 1e30
                col_r = 20; col_g = 20; col_b = 40

                # Sphere intersections
                for s in range(ns):
                    cx = centers[s, 0]; cy = centers[s, 1]; cz = centers[s, 2]
                    ocx = cam_pos[0] - cx
                    ocy = cam_pos[1] - cy
                    ocz = cam_pos[2] - cz
                    b = 2.0 * (rx * ocx + ry * ocy + rz * ocz)
                    c = ocx*ocx + ocy*ocy + ocz*ocz - radii[s]*radii[s]
                    disc = b*b - 4.0*c
                    if disc <= 0:
                        continue
                    sq = disc**0.5
                    t0 = (-b - sq) / 2.0
                    t1 = (-b + sq) / 2.0
                    t = t0 if t0 > 1e-6 else (t1 if t1 > 1e-6 else 1e30)
                    if t < best_t:
                        best_t = t
                        # compute color roughly (diffuse only)
                        # compute hit point and normal
                        hx = cam_pos[0] + rx * t
                        hy = cam_pos[1] + ry * t
                        hz = cam_pos[2] + rz * t
                        nx = hx - cx; ny = hy - cy; nz = hz - cz
                        nlen = (nx*nx + ny*ny + nz*nz)**0.5
                        nx /= nlen; ny /= nlen; nz /= nlen
                        # simple directional light
                        lx = 0.577; ly = -0.577; lz = 0.577
                        llen = (lx*lx + ly*ly + lz*lz)**0.5
                        lx /= llen; ly /= llen; lz /= llen
                        diff = nx*lx + ny*ly + nz*lz
                        if diff < 0.0:
                            diff = 0.0
                        cr = colors[s, 0]; cg = colors[s, 1]; cb = colors[s, 2]
                        ir = 10 + int(255 * (0.1 * cr + diff * cr))
                        ig = 10 + int(255 * (0.1 * cg + diff * cg))
                        ib = 10 + int(255 * (0.1 * cb + diff * cb))
                        col_r = ir if ir < 255 else 255
                        col_g = ig if ig < 255 else 255
                        col_b = ib if ib < 255 else 255

                # Plane intersections (simple shading)
                for p in range(np_):
                    nx = p_normals[p, 0]; ny = p_normals[p, 1]; nz = p_normals[p, 2]
                    denom = nx*rx + ny*ry + nz*rz
                    if abs(denom) < 1e-6:
                        continue
                    t = ((p_points[p,0] - cam_pos[0]) * nx + (p_points[p,1] - cam_pos[1]) * ny + (p_points[p,2] - cam_pos[2]) * nz) / denom
                    if t > 1e-6 and t < best_t:
                        best_t = t
                        # compute simple shading
                        hx = cam_pos[0] + rx * t
                        hy = cam_pos[1] + ry * t
                        hz = cam_pos[2] + rz * t
                        lx = 0.577; ly = -0.577; lz = 0.577
                        diff = nx*lx + ny*ly + nz*lz
                        if diff < 0.0:
                            diff = 0.0
                        cr = p_colors[p, 0]; cg = p_colors[p, 1]; cb = p_colors[p, 2]
                        ir = 10 + int(255 * (0.1 * cr + diff * cr))
                        ig = 10 + int(255 * (0.1 * cg + diff * cg))
                        ib = 10 + int(255 * (0.1 * cb + diff * cb))
                        col_r = ir if ir < 255 else 255
                        col_g = ig if ig < 255 else 255
                        col_b = ib if ib < 255 else 255

                out[j, i, 0] = col_r
                out[j, i, 1] = col_g
                out[j, i, 2] = col_b

        return out


    def render_numba(scene_objects, camera, width=300, height=200):
        """High-level wrapper: convert scene to SoA and call Numba renderer.

        Returns a numpy uint8 array shape (height, width, 3).
        """
        if not NUMBA_AVAILABLE:
            raise ImportError('Numba is not available')

        so = scene_to_soa(scene_objects)
        centers = so['centers']
        radii = so['radii']
        colors = so['colors']
        specs = so['specs']
        p_points = so['p_points']
        p_normals = so['p_normals']
        p_colors = so['p_colors']
        p_specs = so['p_specs']

        cam_pos = np.array(camera.position, dtype=np.float32)
        cam_right = np.array(camera.right, dtype=np.float32)
        cam_up = np.array(camera.up, dtype=np.float32)
        cam_forward = np.array(camera.forward, dtype=np.float32)

        out = _render_numba(width, height,
                            centers, radii, colors, specs,
                            p_points, p_normals, p_colors, p_specs,
                            cam_pos, cam_right, cam_up, cam_forward, camera.fov)
        return out

else:
    def render_numba(*args, **kwargs):
        raise ImportError('Numba is not installed')
