[🇬🇧 English version](README.en.md)

# Moteur de rendu 3D — Lancer de rayons

> Moteur de rendu 3D par lancer de rayons (ray tracing) écrit de zéro, avec accélération Numba et interface graphique.

**Contexte** : projet personnel (prolongement de mon projet de moteur graphique 3D en prépa).

## Ce que fait le projet
Rendu photoréaliste de scènes 3D :
- Intersection rayon/objets (**sphères**, **plans**), éclairage diffus + spéculaire (Phong), ombres
- **Caméra** orientable (orbite, cadrage automatique sur la scène), rendu multi-angles
- **Accélération Numba** (compilation JIT) pour le calcul intensif par pixel
- **Interface graphique** et mode interactif / mode batch
- Génération de scènes complexes (grilles de sphères, sol)

## Stack
Python (NumPy, **Numba**, Pillow), GUI. Architecture en package.

## Fichiers
- `core.py` — primitives, rayons, shading, caméra
- `numba_renderer.py` — rendu accéléré
- `scenes.py` — scènes prédéfinies · `gui.py` — interface · `run.py` — point d'entrée
- `out/` — images rendues · `requirements.txt`

## Lancer
```bash
pip install -r requirements.txt
python -m raytracer.run            # GUI / interactif
python -m raytracer.run --batch    # rendu multi-angles
```
