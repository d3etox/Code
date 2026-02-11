#!/usr/bin/env python3
import math
import os
import numpy as np
from raytracer.core import Ray, Sphere, Plane, shade, Camera, RayCapture
from PIL import Image, ImageDraw
import json
import sys
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import numpy as _np
import subprocess
from typing import Tuple
from raytracer import numba_renderer

def get_scene_bounds(scene_objects):
    """Calculate bounding box of all objects in scene"""
    min_bounds = np.array([float('inf')] * 3)
    max_bounds = np.array([float('-inf')] * 3)
    
    for obj in scene_objects:
        if isinstance(obj, Sphere):
            min_bounds = np.minimum(min_bounds, obj.center - obj.radius)
            max_bounds = np.maximum(max_bounds, obj.center + obj.radius)
        elif isinstance(obj, Plane):
            # For planes, just expand bounds reasonably
            min_bounds = np.minimum(min_bounds, obj.point - 10)
            max_bounds = np.maximum(max_bounds, obj.point + 10)
    
    return min_bounds, max_bounds

def render(scene_objects, camera, width=300, height=200, capture_rays=None):
    """Render scene with given camera"""
    aspect = width / height
    img = Image.new('RGB', (width, height))
    pixels = img.load()
    
    scale = math.tan(math.radians(camera.fov * 0.5))
    
    for j in range(height):
        for i in range(width):
            x = (2 * (i + 0.5) / width - 1) * scale * aspect
            y = (1 - 2 * (j + 0.5) / height) * scale
            
            # Ray in camera space
            ray_dir = camera.right * x + camera.up * y + camera.forward
            ray_dir = ray_dir / np.linalg.norm(ray_dir)
            
            ray = Ray(camera.position, ray_dir)
            
            # Capture specific rays if requested
            if capture_rays is not None and (i, j) in capture_rays:
                capture_rays[(i, j)].capture_ray(ray, scene_objects)
            
            nearest = None
            for obj in scene_objects:
                hit = obj.intersect(ray) if hasattr(obj, 'intersect') else None
                if hit:
                    if nearest is None or hit['t'] < nearest['t']:
                        nearest = hit
            
            if nearest:
                col = shade(nearest, scene_objects)
                r, g, b = (int(255 * col[0]), int(255 * col[1]), int(255 * col[2]))
                pixels[i, j] = (r, g, b)
            else:
                pixels[i, j] = (20, 20, 40)
    
    return img

def visualize_ray_paths(img, ray_paths, width, height):
    """Draw captured ray paths onto image"""
    draw = ImageDraw.Draw(img, 'RGBA')
    
    for ray_path in ray_paths:
        segments = ray_path.get_segments()
        for start, end in segments:
            # Simple projection (would need proper camera matrix for accuracy)
            p1 = tuple(int(c) for c in [start[0] * 100 + width/2, start[1] * 100 + height/2][:2])
            p2 = tuple(int(c) for c in [end[0] * 100 + width/2, end[1] * 100 + height/2][:2])
            draw.line([p1, p2], fill=(255, 255, 0, 128), width=2)


def _render_tile(args):
    """Render a rectangular tile (SoA inputs). Runs in worker process."""
    (centers, radii, colors, specs,
     p_points, p_normals, p_colors, p_specs,
     cam_pos, cam_right, cam_up, cam_forward, fov,
     x0, x1, y0, y1, width, height) = args

    tile_h = y1 - y0
    tile_w = x1 - x0
    tile = _np.zeros((tile_h, tile_w, 3), dtype=_np.uint8)

    aspect = width / height
    scale = math.tan(math.radians(fov * 0.5))

    ns = centers.shape[0]
    np_ = p_points.shape[0]

    for j in range(y0, y1):
        for i in range(x0, x1):
            x = (2 * (i + 0.5) / width - 1) * scale * aspect
            y = (1 - 2 * (j + 0.5) / height) * scale
            rx = cam_right[0] * x + cam_up[0] * y + cam_forward[0]
            ry = cam_right[1] * x + cam_up[1] * y + cam_forward[1]
            rz = cam_right[2] * x + cam_up[2] * y + cam_forward[2]
            rlen = (rx*rx + ry*ry + rz*rz) ** 0.5
            rx /= rlen; ry /= rlen; rz /= rlen

            best_t = 1e30
            col = (20,20,40)

            # spheres
            for s in range(ns):
                cx = centers[s,0]; cy = centers[s,1]; cz = centers[s,2]
                ocx = cam_pos[0] - cx; ocy = cam_pos[1] - cy; ocz = cam_pos[2] - cz
                b = 2.0 * (rx * ocx + ry * ocy + rz * ocz)
                c = ocx*ocx + ocy*ocy + ocz*ocz - radii[s]*radii[s]
                disc = b*b - 4.0*c
                if disc <= 0:
                    continue
                sq = math.sqrt(disc)
                t0 = (-b - sq) / 2.0
                t1 = (-b + sq) / 2.0
                t = t0 if t0 > 1e-6 else (t1 if t1 > 1e-6 else 1e30)
                if t < best_t:
                    best_t = t
                    # shading (diffuse-like)
                    cr = colors[s,0]; cg = colors[s,1]; cb = colors[s,2]
                    # compute normal
                    hx = cam_pos[0] + rx * t; hy = cam_pos[1] + ry * t; hz = cam_pos[2] + rz * t
                    nx = (hx - cx); ny = (hy - cy); nz = (hz - cz)
                    nlen = (nx*nx + ny*ny + nz*nz) ** 0.5
                    nx /= nlen; ny /= nlen; nz /= nlen
                    lx = 0.577; ly = -0.577; lz = 0.577
                    diff = nx*lx + ny*ly + nz*lz
                    if diff < 0.0: diff = 0.0
                    ir = int(255 * (0.1 * cr + diff * cr))
                    ig = int(255 * (0.1 * cg + diff * cg))
                    ib = int(255 * (0.1 * cb + diff * cb))
                    col = (ir, ig, ib)

            # planes
            for p in range(np_):
                nx = p_normals[p,0]; ny = p_normals[p,1]; nz = p_normals[p,2]
                denom = nx*rx + ny*ry + nz*rz
                if abs(denom) < 1e-6:
                    continue
                t = ((p_points[p,0] - cam_pos[0]) * nx + (p_points[p,1] - cam_pos[1]) * ny + (p_points[p,2] - cam_pos[2]) * nz) / denom
                if t > 1e-6 and t < best_t:
                    best_t = t
                    cr = p_colors[p,0]; cg = p_colors[p,1]; cb = p_colors[p,2]
                    lx = 0.577; ly = -0.577; lz = 0.577
                    diff = nx*lx + ny*ly + nz*lz
                    if diff < 0.0: diff = 0.0
                    ir = int(255 * (0.1 * cr + diff * cr))
                    ig = int(255 * (0.1 * cg + diff * cg))
                    ib = int(255 * (0.1 * cb + diff * cb))
                    col = (ir, ig, ib)

            tile[j - y0, i - x0, 0] = int(max(0, min(255, col[0])))
            tile[j - y0, i - x0, 1] = int(max(0, min(255, col[1])))
            tile[j - y0, i - x0, 2] = int(max(0, min(255, col[2])))

    return (x0, y0, tile)


def _int255(v):
    return int(max(0, min(255, v * 255)))


def render_parallel(scene_objects, camera, width=300, height=200, tile_size=64, workers=None):
    """Render using multiple processes by splitting the image into tiles.

    Returns a PIL Image.
    """
    # If Numba is available, use the threaded numba renderer (no copy overhead)
    try:
        if getattr(numba_renderer, 'NUMBA_AVAILABLE', False):
            arr = numba_renderer.render_numba(scene_objects, camera, width, height)
            # arr is uint8 numpy (h,w,3)
            return Image.fromarray(arr, 'RGB')
    except Exception:
        # If Numba renderer fails, fall back to process-pool implementation below
        pass

    if workers is None:
        workers = max(1, multiprocessing.cpu_count() - 1)

    # Convert scene to SoA once for the process-pool fallback
    so = numba_renderer.scene_to_soa(scene_objects)
    centers = so['centers']
    radii = so['radii']
    colors = so['colors']
    specs = so['specs']
    p_points = so['p_points']
    p_normals = so['p_normals']
    p_colors = so['p_colors']
    p_specs = so['p_specs']

    cam_pos = _np.array(camera.position, dtype=_np.float32)
    cam_right = _np.array(camera.right, dtype=_np.float32)
    cam_up = _np.array(camera.up, dtype=_np.float32)
    cam_forward = _np.array(camera.forward, dtype=_np.float32)

    args = []
    for y in range(0, height, tile_size):
        for x in range(0, width, tile_size):
            x1 = min(width, x + tile_size)
            y1 = min(height, y + tile_size)
            args.append((centers, radii, colors, specs,
                         p_points, p_normals, p_colors, p_specs,
                         cam_pos, cam_right, cam_up, cam_forward, camera.fov,
                         x, x1, y, y1, width, height))

    final = _np.zeros((height, width, 3), dtype=_np.uint8)

    # Use ProcessPoolExecutor to render tiles in parallel
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for res in ex.map(_render_tile, args):
            x0, y0, tile = res
            h, w, _ = tile.shape
            final[y0:y0+h, x0:x0+w, :] = tile

    # Convert to PIL
    img = Image.fromarray(final, 'RGB')
    return img


def _get_cache_size(level: str) -> int:
    """Try to read cache size in bytes for 'L1', 'L2', or 'L3'.

    On macOS use `sysctl hw.l1dcachesize` etc. Fall back to common defaults.
    """
    level = level.upper()
    try:
        if sys.platform == 'darwin':
            key = {'L1':'hw.l1dcachesize','L2':'hw.l2cachesize','L3':'hw.l3cachesize'}[level]
            out = subprocess.check_output(['sysctl', '-n', key], text=True).strip()
            return int(out)
    except Exception:
        pass
    # Fallback defaults (bytes)
    defaults = {'L1': 32*1024, 'L2': 256*1024, 'L3': 8*1024*1024}
    return defaults.get(level, 256*1024)


def estimate_tile_size(scene_objects, camera, target_cache='L2', bytes_per_ray=32):
    """Estimate a tile size (square) that keeps working set within target cache.

    bytes_per_ray: approximate bytes needed per ray (SoA dirs + origin + temp)
    Returns tile_size in pixels (int).
    """
    so = numba_renderer.scene_to_soa(scene_objects)
    # geometry size
    geom_bytes = so['centers'].nbytes + so['radii'].nbytes + so['colors'].nbytes + so['specs'].nbytes
    geom_bytes += so['p_points'].nbytes + so['p_normals'].nbytes + so['p_colors'].nbytes + so['p_specs'].nbytes

    cache_bytes = _get_cache_size(target_cache)
    # leave margin
    usable = int(cache_bytes * 0.7)
    if usable <= 0:
        usable = cache_bytes

    # rays budget = usable - geometry
    if geom_bytes >= usable:
        # geometry doesn't fit; choose small tile
        return max(8, min(64, int(_np.sqrt(usable / bytes_per_ray))))

    ray_budget = usable - geom_bytes
    rays = int(ray_budget / bytes_per_ray)
    # choose square tile
    tile = int(_np.floor(_np.sqrt(rays)))
    tile = max(8, min(256, tile))
    return tile

def interactive_mode(scene_objects):
    """Interactive camera exploration mode"""
    camera = Camera()
    min_bounds, max_bounds = get_scene_bounds(scene_objects)
    camera.focus_on_bounds(min_bounds, max_bounds)
    
    ray_capture = RayCapture()
    
    print("=== Ray Engine - Interactive Camera Mode ===")
    print("Commands:")
    print("  orbit <theta> <phi>    - Orbit camera (radians)")
    print("  pan <dx> <dy>          - Pan camera")
    print("  dolly <distance>       - Move forward/backward")
    print("  zoom <fov>             - Change field of view")
    print("  reset                  - Reset to default view")
    print("  focus                  - Focus on all objects")
    print("  capture <x> <y>        - Capture ray at pixel")
    print("  render <filename>      - Render and save image")
    print("  info                   - Show camera info")
    print("  quit                   - Exit")
    print()
    
    outdir = '/Users/charlesfabre/Desktop/Code/raytracer/out'
    os.makedirs(outdir, exist_ok=True)
    
    while True:
        try:
            cmd = input("raytracer> ").strip().split()
            if not cmd:
                continue
            
            action = cmd[0].lower()
            
            if action == 'quit':
                print("Exiting...")
                break
            
            elif action == 'orbit' and len(cmd) == 3:
                theta = float(cmd[1])
                phi = float(cmd[2])
                camera.orbit(theta, phi)
                print(f"Orbited by theta={theta}, phi={phi}")
            
            elif action == 'pan' and len(cmd) == 3:
                dx = float(cmd[1])
                dy = float(cmd[2])
                camera.pan(dx, dy)
                print(f"Panned by dx={dx}, dy={dy}")
            
            elif action == 'dolly' and len(cmd) == 2:
                dist = float(cmd[1])
                camera.dolly(dist)
                print(f"Dollied by {dist}")
            
            elif action == 'zoom' and len(cmd) == 2:
                fov = float(cmd[1])
                camera.fov = fov
                print(f"FOV set to {fov}°")
            
            elif action == 'reset':
                camera = Camera()
                camera.focus_on_bounds(min_bounds, max_bounds)
                print("Camera reset")
            
            elif action == 'focus':
                camera.focus_on_bounds(min_bounds, max_bounds)
                print("Focused on all objects")
            
            elif action == 'capture' and len(cmd) == 3:
                x, y = int(cmd[1]), int(cmd[2])
                capture_dict = {(x, y): ray_capture}
                img = render(scene_objects, camera, 800, 600, capture_dict)
                ray_path = ray_capture.captured_rays[-1] if ray_capture.captured_rays else None
                if ray_path:
                    print(f"Captured ray at ({x}, {y})")
                    print(f"  Hits: {len(ray_path.hits)}")
                    for i, (point, obj) in enumerate(ray_path.hits):
                        print(f"    {i+1}. {type(obj).__name__} at {point}")
            
            elif action == 'render' and len(cmd) >= 2:
                filename = cmd[1]
                filepath = os.path.join(outdir, filename)
                img = render(scene_objects, camera, 800, 600)
                img.save(filepath)
                print(f"Rendered to {filepath}")
            
            elif action == 'info':
                print(f"Position: {camera.position}")
                print(f"Look at:  {camera.look_at}")
                print(f"Distance: {camera.distance:.2f}")
                print(f"FOV:      {camera.fov}°")
                print(f"Scene bounds: {min_bounds} to {max_bounds}")
            
            else:
                print("Unknown command or incorrect arguments")
        
        except ValueError as e:
            print(f"Error: {e}")
        except KeyboardInterrupt:
            print("\nExiting...")
            break

if __name__ == '__main__':
    outdir = '/Users/charlesfabre/Desktop/Code/raytracer/out'
    os.makedirs(outdir, exist_ok=True)
    
    # Build a more complex scene (grid of spheres + floor). Fall back to a small sample if import fails.
    try:
        from raytracer.scenes import complex_scene
        scene_objs = complex_scene(grid_size=10, spacing=1.25)
    except Exception:
        s1 = Sphere([0.0, -0.5, -3.0], 0.7, color=(0.8, 0.2, 0.2), spec=100)
        s2 = Sphere([1.0, 0.0, -4.0], 0.9, color=(0.2, 0.8, 0.2), spec=10)
        floor = Plane([0, -1, 0], [0, 1, 0], color=(0.7, 0.7, 0.7), spec=5)
        scene_objs = [s1, s2, floor]

    if len(sys.argv) > 1 and sys.argv[1] == '--batch':
        # Batch render from multiple angles
        camera = Camera()
        min_bounds, max_bounds = get_scene_bounds(scene_objs)
        camera.focus_on_bounds(min_bounds, max_bounds)

        # Render from multiple angles
        for i, (theta, phi) in enumerate([(0, 0), (np.pi/4, 0), (np.pi/2, np.pi/6)]):
            camera.orbit(theta, phi)
            img = render(scene_objs, camera, 600, 400)
            filename = f'raytracer_angle_{i}.png'
            img.save(os.path.join(outdir, filename))
            print(f'Saved {filename}')
    else:
        # Try to launch GUI if available, otherwise fall back to interactive CLI
        try:
            from raytracer import gui as gui_module
            print('Launching GUI...')
            gui_module.main()
        except Exception as e:
            print('GUI unavailable or failed to start, falling back to interactive mode.')
            print('Reason:', e)
            interactive_mode(scene_objs)
