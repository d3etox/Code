import numpy as np
import math
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple

EPS = 1e-6

def normalize(v):
    v = np.array(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-12:
        return v
    return v / n

class Ray:
    def __init__(self, orig, dir):
        self.orig = np.array(orig, dtype=float)
        self.dir = normalize(dir)

class Sphere:
    def __init__(self, center, radius, color=(1,1,1), spec=50):
        self.center = np.array(center, dtype=float)
        self.radius = float(radius)
        self.color = np.array(color, dtype=float)
        self.spec = spec

    def intersect(self, ray):
        oc = ray.orig - self.center
        b = 2.0 * np.dot(ray.dir, oc)
        c = np.dot(oc, oc) - self.radius*self.radius
        disc = b*b - 4*c
        if disc < 0:
            return None
        sq = math.sqrt(disc)
        t0 = (-b - sq)/2.0
        t1 = (-b + sq)/2.0
        t = t0 if t0>EPS else (t1 if t1>EPS else None)
        if t is None:
            return None
        hit = ray.orig + ray.dir * t
        normal = normalize(hit - self.center)
        return {'t':t, 'point':hit, 'normal':normal, 'object':self}

class Plane:
    def __init__(self, point, normal, color=(1,1,1), spec=10):
        self.point = np.array(point, dtype=float)
        self.normal = normalize(normal)
        self.color = np.array(color, dtype=float)
        self.spec = spec

    def intersect(self, ray):
        denom = np.dot(self.normal, ray.dir)
        if abs(denom) < EPS:
            return None
        t = np.dot(self.point - ray.orig, self.normal) / denom
        if t < EPS:
            return None
        hit = ray.orig + ray.dir * t
        return {'t':t, 'point':hit, 'normal':self.normal, 'object':self}

def ray_triangle_intersect(ray, v0, v1, v2):
    # Möller–Trumbore
    v0, v1, v2 = map(np.array, (v0,v1,v2))
    edge1 = v1 - v0
    edge2 = v2 - v0
    h = np.cross(ray.dir, edge2)
    a = np.dot(edge1, h)
    if abs(a) < EPS:
        return None
    f = 1.0 / a
    s = ray.orig - v0
    u = f * np.dot(s, h)
    if u < 0.0 or u > 1.0:
        return None
    q = np.cross(s, edge1)
    v = f * np.dot(ray.dir, q)
    if v < 0.0 or u+v > 1.0:
        return None
    t = f * np.dot(edge2, q)
    if t > EPS:
        point = ray.orig + ray.dir * t
        normal = normalize(np.cross(edge1, edge2))
        return {'t':t, 'point':point, 'normal':normal}
    return None

def get_scene_bounds(scene_objects):
    """Calculate bounding box of all objects in scene"""
    min_bounds = np.array([float('inf')] * 3)
    max_bounds = np.array([float('-inf')] * 3)

    for obj in scene_objects:
        if isinstance(obj, Sphere):
            min_bounds = np.minimum(min_bounds, obj.center - obj.radius)
            max_bounds = np.maximum(max_bounds, obj.center + obj.radius)
        elif isinstance(obj, Plane):
            # For planes, expand bounds around the plane point
            min_bounds = np.minimum(min_bounds, obj.point - 10)
            max_bounds = np.maximum(max_bounds, obj.point + 10)
    
    return min_bounds, max_bounds

def shade(hit, scene, light_dir=np.array([1,-1,1])):
    obj = hit['object']
    n = hit['normal']
    l = normalize(-light_dir)
    diff = max(0.0, np.dot(n, l))
    view = normalize(-hit['point'])  # camera at origin convention
    r = 2*np.dot(n,l)*n - l
    spec = max(0.0, np.dot(r, view))**(obj.spec if hasattr(obj,'spec') else 10)
    color = 0.1*obj.color + diff*obj.color + 0.5*spec
    return np.clip(color,0,1)

class Camera:
    """Interactive camera with orbit and tracking capabilities"""
    def __init__(self, position=None, look_at=None, fov=50):
        self.position = np.array(position if position is not None else [0, 0, 5], dtype=float)
        self.look_at = np.array(look_at if look_at is not None else [0, 0, 0], dtype=float)
        self.fov = fov
        self.distance = np.linalg.norm(self.position - self.look_at)
        self.update_vectors()

    def update_vectors(self):
        """Compute forward, right, and up vectors"""
        self.forward = normalize(self.look_at - self.position)
        world_up = np.array([0, 1, 0], dtype=float)
        self.right = normalize(np.cross(self.forward, world_up))
        self.up = normalize(np.cross(self.right, self.forward))

    def pan(self, dx, dy):
        """Move camera left/right and up/down"""
        self.position += self.right * dx
        self.position += self.up * dy
        self.look_at += self.right * dx
        self.look_at += self.up * dy
        self.update_vectors()

    def dolly(self, amount):
        """Move camera forward/backward"""
        direction = normalize(self.look_at - self.position)
        self.position += direction * amount
        self.distance = np.linalg.norm(self.position - self.look_at)
        self.update_vectors()

    def orbit(self, theta, phi):
        """Orbit around look_at point (theta: horizontal, phi: vertical)"""
        rel_pos = self.position - self.look_at
        
        # Convert to spherical coords
        r = np.linalg.norm(rel_pos)
        old_theta = np.arctan2(rel_pos[0], rel_pos[2])
        old_phi = np.arccos(np.clip(rel_pos[1] / r, -1, 1))
        
        # Apply rotation
        new_theta = old_theta + theta
        new_phi = np.clip(old_phi + phi, 0.1, np.pi - 0.1)
        
        # Convert back to cartesian
        self.position = self.look_at + np.array([
            r * np.sin(new_phi) * np.sin(new_theta),
            r * np.cos(new_phi),
            r * np.sin(new_phi) * np.cos(new_theta)
        ])
        self.update_vectors()

    def focus_on_bounds(self, min_bounds, max_bounds):
        """Position camera to view all objects within bounds"""
        center = (min_bounds + max_bounds) / 2
        size = np.linalg.norm(max_bounds - min_bounds)
        distance = size / (2 * np.tan(np.radians(self.fov / 2)))
        
        self.look_at = center
        self.position = center + np.array([1, 1, 1]) * (distance / np.sqrt(3))
        self.distance = np.linalg.norm(self.position - self.look_at)
        self.update_vectors()

@dataclass
class RayPath:
    """Represents a captured ray path for visualization"""
    origin: np.ndarray
    direction: np.ndarray
    hits: List[Tuple[np.ndarray, object]]  # list of (point, object) tuples
    
    def get_segments(self) -> List[Tuple[np.ndarray, np.ndarray]]:
        """Get line segments for drawing the ray path"""
        segments = []
        current_point = self.origin
        
        for hit_point, obj in self.hits:
            segments.append((current_point.copy(), hit_point.copy()))
            current_point = hit_point
        
        # Add ray continuation if no final hit
        if not self.hits:
            end_point = self.origin + self.direction * 100
            segments.append((self.origin.copy(), end_point.copy()))
        
        return segments

class RayCapture:
    """Captures and traces rays through the scene"""
    def __init__(self):
        self.captured_rays: List[RayPath] = []

    def capture_ray(self, ray: Ray, scene_objects, max_bounces: int = 3) -> RayPath:
        """Capture a ray and trace its path through the scene"""
        hits = []
        current_ray = ray
        
        for _ in range(max_bounces):
            nearest = None
            for obj in scene_objects:
                if hasattr(obj, 'intersect'):
                    hit = obj.intersect(current_ray)
                    if hit:
                        if nearest is None or hit['t'] < nearest['t']:
                            nearest = hit
            
            if nearest is None:
                break
            
            hit_point = nearest['point']
            hits.append((hit_point, nearest['object']))
            
            # Reflect ray (simple specular reflection)
            normal = nearest['normal']
            current_ray = Ray(hit_point + normal * EPS, 
                            current_ray.dir - 2 * np.dot(current_ray.dir, normal) * normal)
        
        ray_path = RayPath(ray.orig.copy(), ray.dir.copy(), hits)
        self.captured_rays.append(ray_path)
        return ray_path

    def clear(self):
        """Clear captured rays"""
        self.captured_rays = []
