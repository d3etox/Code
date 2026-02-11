#!/usr/bin/env python3
import sys
import math
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QLabel, QSlider, QPushButton, 
                             QSpinBox, QDoubleSpinBox, QGroupBox, QCheckBox)
from PyQt5.QtGui import QImage, QPixmap, QFont
from PyQt5.QtCore import Qt, QTimer, pyqtSignal, QThread
from PyQt5.QtWidgets import QScrollArea
from PIL import Image
import numpy as np
from raytracer.core import Ray, Sphere, Plane, shade, Camera, RayCapture, get_scene_bounds

class RenderThread(QThread):
    """Render in background thread to keep UI responsive"""
    image_ready = pyqtSignal(QImage)
    
    def __init__(self, scene_objects, camera, width, height):
        super().__init__()
        self.scene_objects = scene_objects
        self.camera = camera
        self.width = width
        self.height = height
        self.running = True
    
    def run(self):
        # Use the process-pool tiled renderer for faster frames
        if not self.running:
            return

        # Try Numba renderer first for max speed, then parallel renderer, then fallback
        try:
            try:
                from raytracer import numba_renderer
                if hasattr(numba_renderer, 'render_numba'):
                    # render_numba returns a uint8 numpy array (h,w,3)
                    arr = numba_renderer.render_numba(self.scene_objects, self.camera, self.width, self.height)
                    h, w, ch = arr.shape
                    bytes_per_line = 3 * w
                    qt_img = QImage(arr.data, w, h, bytes_per_line, QImage.Format_RGB888)
                    self.image_ready.emit(qt_img)
                    return
            except Exception:
                # Numba not available or failed; continue to parallel renderer
                pass

            from raytracer.run import render_parallel, estimate_tile_size
            try:
                tile_size = estimate_tile_size(self.scene_objects, self.camera, target_cache='L2')
            except Exception:
                tile_size = 64
            pil_img = render_parallel(self.scene_objects, self.camera, self.width, self.height, tile_size=tile_size)
            pil_img = pil_img.convert('RGB')
            arr = np.array(pil_img)
            h, w, ch = arr.shape
            bytes_per_line = 3 * w
            qt_img = QImage(arr.data, w, h, bytes_per_line, QImage.Format_RGB888)
            self.image_ready.emit(qt_img)
        except Exception as e:
            # Final fallback to synchronous Python renderer if all else fails
            print('Render failed (numba/parallel), falling back to local render:', e)
            aspect = self.width / self.height
            img_array = np.zeros((self.height, self.width, 3), dtype=np.uint8)
            scale = math.tan(math.radians(self.camera.fov * 0.5))
            for j in range(self.height):
                for i in range(self.width):
                    x = (2 * (i + 0.5) / self.width - 1) * scale * aspect
                    y = (1 - 2 * (j + 0.5) / self.height) * scale
                    ray_dir = self.camera.right * x + self.camera.up * y + self.camera.forward
                    ray_dir = ray_dir / np.linalg.norm(ray_dir)
                    ray = Ray(self.camera.position, ray_dir)
                    nearest = None
                    for obj in self.scene_objects:
                        hit = obj.intersect(ray) if hasattr(obj, 'intersect') else None
                        if hit:
                            if nearest is None or hit['t'] < nearest['t']:
                                nearest = hit
                    if nearest:
                        col = shade(nearest, self.scene_objects)
                        img_array[j, i] = [int(255 * col[0]), int(255 * col[1]), int(255 * col[2])]
                    else:
                        img_array[j, i] = [20, 20, 40]
            h, w, ch = img_array.shape
            bytes_per_line = 3 * w
            qt_img = QImage(img_array.data, w, h, bytes_per_line, QImage.Format_RGB888)
            self.image_ready.emit(qt_img)
    
    def stop(self):
        self.running = False

class RaytracerGUI(QMainWindow):
    def __init__(self, scene_objects):
        super().__init__()
        self.scene_objects = scene_objects
        self.camera = Camera()
        min_bounds, max_bounds = get_scene_bounds(scene_objects)
        self.min_bounds = min_bounds
        self.max_bounds = max_bounds
        self.camera.focus_on_bounds(min_bounds, max_bounds)
        
        self.ray_capture = RayCapture()
        self.render_thread = None
        self.is_rendering = False
        
        self.init_ui()
        self.render_scene()
    
    def init_ui(self):
        self.setWindowTitle("Raytracer - Interactive Viewer")
        self.setGeometry(100, 100, 1400, 800)
        
        # Main layout
        main_widget = QWidget()
        main_layout = QHBoxLayout()
        
        # Canvas (left side)
        self.canvas = QLabel()
        self.canvas.setMinimumSize(800, 600)
        self.canvas.setStyleSheet("background-color: #1a1a2e; border: 1px solid #0f3460;")
        main_layout.addWidget(self.canvas)
        
        # Control panel (right side)
        control_panel = QWidget()
        control_layout = QVBoxLayout()
        
        # Camera position controls
        camera_group = QGroupBox("Camera Controls")
        camera_layout = QVBoxLayout()
        
        # Orbit controls
        orbit_box = QGroupBox("Orbit")
        orbit_layout = QVBoxLayout()
        
        # Theta slider
        orbit_layout.addWidget(QLabel("Horizontal (θ):"))
        self.theta_slider = QSlider(Qt.Horizontal)
        self.theta_slider.setRange(-180, 180)
        self.theta_slider.setValue(0)
        self.theta_slider.sliderMoved.connect(self.on_theta_changed)
        self.theta_slider.setTickPosition(QSlider.TicksBelow)
        self.theta_slider.setTickInterval(45)
        orbit_layout.addWidget(self.theta_slider)
        
        # Phi slider
        orbit_layout.addWidget(QLabel("Vertical (φ):"))
        self.phi_slider = QSlider(Qt.Horizontal)
        self.phi_slider.setRange(10, 170)
        self.phi_slider.setValue(45)
        self.phi_slider.sliderMoved.connect(self.on_phi_changed)
        self.phi_slider.setTickPosition(QSlider.TicksBelow)
        self.phi_slider.setTickInterval(30)
        orbit_layout.addWidget(self.phi_slider)
        
        orbit_box.setLayout(orbit_layout)
        camera_layout.addWidget(orbit_box)
        
        # Distance/Zoom
        zoom_box = QGroupBox("Distance & Zoom")
        zoom_layout = QVBoxLayout()
        
        zoom_layout.addWidget(QLabel("Distance from object:"))
        self.distance_spinbox = QDoubleSpinBox()
        self.distance_spinbox.setRange(0.1, 100)
        self.distance_spinbox.setValue(self.camera.distance)
        self.distance_spinbox.setSingleStep(0.5)
        self.distance_spinbox.valueChanged.connect(self.on_distance_changed)
        zoom_layout.addWidget(self.distance_spinbox)
        
        zoom_layout.addWidget(QLabel("Field of View:"))
        self.fov_spinbox = QDoubleSpinBox()
        self.fov_spinbox.setRange(10, 120)
        self.fov_spinbox.setValue(self.camera.fov)
        self.fov_spinbox.setSingleStep(1)
        self.fov_spinbox.valueChanged.connect(self.on_fov_changed)
        zoom_layout.addWidget(self.fov_spinbox)
        
        zoom_box.setLayout(zoom_layout)
        camera_layout.addWidget(zoom_box)
        
        camera_group.setLayout(camera_layout)
        control_layout.addWidget(camera_group)
        
        # Pan controls
        pan_box = QGroupBox("Pan (Left/Right & Up/Down)")
        pan_layout = QVBoxLayout()
        
        pan_layout.addWidget(QLabel("Horizontal:"))
        self.pan_x_slider = QSlider(Qt.Horizontal)
        self.pan_x_slider.setRange(-50, 50)
        self.pan_x_slider.setValue(0)
        self.pan_x_slider.sliderMoved.connect(self.on_pan_x_changed)
        pan_layout.addWidget(self.pan_x_slider)
        
        pan_layout.addWidget(QLabel("Vertical:"))
        self.pan_y_slider = QSlider(Qt.Horizontal)
        self.pan_y_slider.setRange(-50, 50)
        self.pan_y_slider.setValue(0)
        self.pan_y_slider.sliderMoved.connect(self.on_pan_y_changed)
        pan_layout.addWidget(self.pan_y_slider)
        
        pan_box.setLayout(pan_layout)
        control_layout.addWidget(pan_box)
        
        # Buttons
        button_box = QGroupBox("Actions")
        button_layout = QVBoxLayout()
        
        focus_btn = QPushButton("Focus on Objects")
        focus_btn.clicked.connect(self.on_focus)
        button_layout.addWidget(focus_btn)
        
        reset_btn = QPushButton("Reset View")
        reset_btn.clicked.connect(self.on_reset)
        button_layout.addWidget(reset_btn)
        
        capture_btn = QPushButton("Capture Ray (Center)")
        capture_btn.clicked.connect(self.on_capture)
        button_layout.addWidget(capture_btn)
        
        export_btn = QPushButton("Export Image")
        export_btn.clicked.connect(self.on_export)
        button_layout.addWidget(export_btn)
        
        button_box.setLayout(button_layout)
        control_layout.addWidget(button_box)
        
        # Camera info
        info_box = QGroupBox("Camera Info")
        info_layout = QVBoxLayout()
        
        self.info_label = QLabel()
        self.info_label.setFont(QFont("Courier", 9))
        info_layout.addWidget(self.info_label)
        
        info_box.setLayout(info_layout)
        control_layout.addWidget(info_box)
        
        # Stretch to fill space
        control_layout.addStretch()
        
        control_panel.setLayout(control_layout)
        control_panel.setMaximumWidth(350)
        
        main_layout.addWidget(control_panel)
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)
        
        self.update_info_label()
    
    def render_scene(self):
        """Render the scene asynchronously"""
        if self.is_rendering:
            return
        
        self.is_rendering = True
        
        # Stop previous thread
        if self.render_thread is not None:
            self.render_thread.stop()
            self.render_thread.wait()
        
        # Start new render thread
        self.render_thread = RenderThread(self.scene_objects, self.camera, 800, 600)
        self.render_thread.image_ready.connect(self.on_image_ready)
        self.render_thread.start()
    
    def on_image_ready(self, qt_img):
        """Display rendered image"""
        pixmap = QPixmap.fromImage(qt_img)
        self.canvas.setPixmap(pixmap)
        self.is_rendering = False
        self.update_info_label()
    
    def on_theta_changed(self, value):
        """Handle horizontal orbit slider"""
        if hasattr(self, '_updating'):
            return
        # Convert slider value to radians
        theta_rad = math.radians(value)
        # Recalculate position
        rel_pos = self.camera.position - self.camera.look_at
        r = np.linalg.norm(rel_pos)
        old_phi = np.arccos(np.clip(rel_pos[1] / r, -1, 1))
        
        self.camera.position = self.camera.look_at + np.array([
            r * np.sin(old_phi) * np.sin(theta_rad),
            r * np.cos(old_phi),
            r * np.sin(old_phi) * np.cos(theta_rad)
        ])
        self.camera.update_vectors()
        self.render_scene()
    
    def on_phi_changed(self, value):
        """Handle vertical orbit slider"""
        if hasattr(self, '_updating'):
            return
        # Convert slider value to radians
        phi_rad = math.radians(value)
        # Recalculate position
        rel_pos = self.camera.position - self.camera.look_at
        r = np.linalg.norm(rel_pos)
        theta_rad = np.arctan2(rel_pos[0], rel_pos[2])
        
        self.camera.position = self.camera.look_at + np.array([
            r * np.sin(phi_rad) * np.sin(theta_rad),
            r * np.cos(phi_rad),
            r * np.sin(phi_rad) * np.cos(theta_rad)
        ])
        self.camera.update_vectors()
        self.render_scene()
    
    def on_distance_changed(self, value):
        """Handle distance spinbox"""
        rel_pos = self.camera.position - self.camera.look_at
        direction = rel_pos / np.linalg.norm(rel_pos)
        self.camera.position = self.camera.look_at + direction * value
        self.camera.distance = value
        self.camera.update_vectors()
        self.render_scene()
    
    def on_fov_changed(self, value):
        """Handle FOV spinbox"""
        self.camera.fov = value
        self.render_scene()
    
    def on_pan_x_changed(self, value):
        """Handle horizontal pan"""
        self.camera.pan(value * 0.01, 0)
        self.render_scene()
    
    def on_pan_y_changed(self, value):
        """Handle vertical pan"""
        self.camera.pan(0, value * 0.01)
        self.render_scene()
    
    def on_focus(self):
        """Focus camera on all objects"""
        self.camera.focus_on_bounds(self.min_bounds, self.max_bounds)
        self.update_sliders()
        self.render_scene()
    
    def on_reset(self):
        """Reset to default view"""
        self.camera = Camera()
        self.camera.focus_on_bounds(self.min_bounds, self.max_bounds)
        self.update_sliders()
        self.render_scene()
    
    def on_capture(self):
        """Capture ray from center"""
        # Create ray from center of viewport
        x, y = 0, 0
        ray_dir = self.camera.right * x + self.camera.up * y + self.camera.forward
        ray_dir = ray_dir / np.linalg.norm(ray_dir)
        ray = Ray(self.camera.position, ray_dir)
        
        self.ray_capture.capture_ray(ray, self.scene_objects, max_bounces=5)
        ray_path = self.ray_capture.captured_rays[-1]
        
        print(f"\n=== Captured Ray ===")
        print(f"Origin: {ray_path.origin}")
        print(f"Direction: {ray_path.direction}")
        print(f"Hits: {len(ray_path.hits)}")
        for i, (point, obj) in enumerate(ray_path.hits):
            print(f"  {i+1}. {type(obj).__name__} at {point}")
    
    def on_export(self):
        """Export current view as image"""
        import os
        outdir = '/Users/charlesfabre/Desktop/Code/raytracer/out'
        os.makedirs(outdir, exist_ok=True)
        
        # Save the current pixmap
        if self.canvas.pixmap():
            count = len([f for f in os.listdir(outdir) if f.startswith('export_')])
            filename = os.path.join(outdir, f'export_{count}.png')
            self.canvas.pixmap().save(filename)
            print(f"Exported to {filename}")
    
    def update_sliders(self):
        """Update slider positions to match camera state"""
        self._updating = True
        rel_pos = self.camera.position - self.camera.look_at
        r = np.linalg.norm(rel_pos)
        theta = np.arctan2(rel_pos[0], rel_pos[2])
        phi = np.arccos(np.clip(rel_pos[1] / r, -1, 1))
        
        self.theta_slider.setValue(int(math.degrees(theta)))
        self.phi_slider.setValue(int(math.degrees(phi)))
        self.distance_spinbox.setValue(r)
        self.fov_spinbox.setValue(self.camera.fov)
        
        delattr(self, '_updating')
    
    def update_info_label(self):
        """Update camera info display"""
        pos = self.camera.position
        look = self.camera.look_at
        info_text = f"""Position: ({pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f})
Look At: ({look[0]:.2f}, {look[1]:.2f}, {look[2]:.2f})
Distance: {self.camera.distance:.2f}
FOV: {self.camera.fov:.1f}°
"""
        self.info_label.setText(info_text)
    
    def closeEvent(self, event):
        """Clean up on close"""
        if self.render_thread is not None:
            self.render_thread.stop()
            self.render_thread.wait()
        event.accept()

def main():
    # Build a more complex scene if available
    try:
        from raytracer.scenes import complex_scene
        scene_objs = complex_scene(grid_size=10, spacing=1.25)
    except Exception:
        s1 = Sphere([0.0, -0.5, -3.0], 0.7, color=(0.8, 0.2, 0.2), spec=100)
        s2 = Sphere([1.0, 0.0, -4.0], 0.9, color=(0.2, 0.8, 0.2), spec=10)
        floor = Plane([0, -1, 0], [0, 1, 0], color=(0.7, 0.7, 0.7), spec=5)
        scene_objs = [s1, s2, floor]

    app = QApplication(sys.argv)
    gui = RaytracerGUI(scene_objs)
    gui.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
