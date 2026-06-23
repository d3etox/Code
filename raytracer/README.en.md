[🇫🇷 Version française](README.md)

# 3D rendering engine — Ray tracing

> A 3D rendering engine built from scratch using ray tracing, with Numba acceleration and a graphical interface.

**Context**: personal project (extension of my preparatory-class 3D graphics-engine project).

## What it does
Photorealistic rendering of 3D scenes:
- Ray/object intersection (**spheres**, **planes**), diffuse + specular lighting (Phong), shadows
- Steerable **camera** (orbit, automatic framing on the scene), multi-angle rendering
- **Numba acceleration** (JIT compilation) for the per-pixel intensive computation
- **Graphical interface** plus interactive / batch modes
- Generation of complex scenes (sphere grids, floor)

## Stack
Python (NumPy, **Numba**, Pillow), GUI. Package architecture.

## Files
- `core.py` — primitives, rays, shading, camera
- `numba_renderer.py` — accelerated rendering
- `scenes.py` — predefined scenes · `gui.py` — interface · `run.py` — entry point
- `out/` — rendered images · `requirements.txt`

## Run
```bash
pip install -r requirements.txt
python -m raytracer.run            # GUI / interactive
python -m raytracer.run --batch    # multi-angle rendering
```
