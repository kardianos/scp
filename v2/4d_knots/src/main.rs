use bytemuck::{Pod, Zeroable};
use egui::Slider;
use egui_wgpu::ScreenDescriptor;
use egui_winit::State as EguiState;
use nalgebra::{Matrix4, Point3, UnitQuaternion, Vector3};
use std::sync::Arc;
use wgpu::util::DeviceExt;
use winit::event::{ElementState, KeyEvent, MouseButton, MouseScrollDelta};
use winit::event_loop::{ControlFlow, EventLoop};
use winit::keyboard::{KeyCode, PhysicalKey};
use winit::window::Window;
use winit::{dpi, event};

const KNOT_RESOLUTION: usize = 100;
const TUBE_SEGMENTS: usize = 8;
const FIELD_RESOLUTION: usize = 16;

#[derive(Clone, Copy, Pod, Zeroable)]
#[repr(C)]
struct Vertex {
    position: [f32; 3],
    normal: [f32; 3],
    color: [f32; 4],
}

#[derive(Clone, Copy, Pod, Zeroable)]
#[repr(C)]
struct Uniforms {
    mvp: [[f32; 4]; 4],
    light_dir: [f32; 3],
    slice_w: f32,
}

#[derive(Clone)]
struct KnotPoint {
    pos: Vector3<f32>,
    w: f32,
}

#[derive(Clone)]
struct Knot {
    points: Vec<KnotPoint>,
    w_velocity: f32,
    center: Vector3<f32>,      // Center position in 3D space
    velocity: Vector3<f32>,    // 3D velocity
    color: [f32; 3],
}

impl Knot {
    fn torus_knot(p: i32, q: i32, r_major: f32, r_minor: f32, center: Vector3<f32>, color: [f32; 3], initial_w: f32) -> Self {
        let points: Vec<KnotPoint> = (0..KNOT_RESOLUTION)
            .map(|i| {
                let t = i as f32 / KNOT_RESOLUTION as f32 * std::f32::consts::TAU;
                let phi = p as f32 * t;
                let theta = q as f32 * t;
                
                let x = (r_major + r_minor * phi.cos()) * theta.cos();
                let y = (r_major + r_minor * phi.cos()) * theta.sin();
                let z = r_minor * phi.sin();
                
                KnotPoint {
                    pos: Vector3::new(x, y, z),  // Local coordinates (relative to center)
                    w: initial_w,
                }
            })
            .collect();
        
        Self { points, w_velocity: 0.0, center, velocity: Vector3::zeros(), color }
    }
    
    fn trefoil_link_pair() -> (Self, Self) {
        let k1 = Knot::torus_knot(2, 3, 0.4, 0.12, Vector3::new(-0.3, 0.0, 0.0), [0.9, 0.2, 0.2], 0.0);
        let mut k2 = Knot::torus_knot(2, 3, 0.4, 0.12, Vector3::new(0.3, 0.0, 0.0), [0.2, 0.9, 0.2], 0.6);
        
        // Give k2 initial velocities
        k2.velocity = Vector3::new(-0.15, 0.05, 0.02);  // Moving toward k1
        k2.w_velocity = -0.2;  // Moving toward k1 in W
        
        (k1, k2)
    }
    
    fn get_4d_point(&self, i: usize) -> (Vector3<f32>, f32) {
        let p = &self.points[i % self.points.len()];
        (self.center + p.pos, p.w)  // World position = center + local
    }
    
    fn average_w(&self) -> f32 {
        self.points.iter().map(|p| p.w).sum::<f32>() / self.points.len() as f32
    }
    
    fn bounding_radius(&self) -> f32 {
        self.points.iter().map(|p| p.pos.magnitude()).fold(0.0, f32::max)
    }
}

struct EnergyField {
    data: Vec<f32>,
    size: usize,
}

impl EnergyField {
    fn new() -> Self {
        Self {
            data: vec![0.0; FIELD_RESOLUTION * FIELD_RESOLUTION * FIELD_RESOLUTION],
            size: FIELD_RESOLUTION,
        }
    }
    
    fn compute(&mut self, knots: &[Knot], w_slice: f32, barrier: f32, binding: f32) {
        let half_extent = 1.5;
        
        for xi in 0..self.size {
            for yi in 0..self.size {
                for zi in 0..self.size {
                    let x = (xi as f32 / (self.size - 1) as f32 - 0.5) * 2.0 * half_extent;
                    let y = (yi as f32 / (self.size - 1) as f32 - 0.5) * 2.0 * half_extent;
                    let z = (zi as f32 / (self.size - 1) as f32 - 0.5) * 2.0 * half_extent;
                    
                    let pos = Vector3::new(x, y, z);
                    let mut energy = 0.0;
                    
                    for knot in knots {
                        for point in &knot.points {
                            let world_pos = knot.center + point.pos;
                            let dist_3d = (pos - world_pos).magnitude();
                            let dw = (point.w - w_slice).abs();
                            
                            if dist_3d < 1.5 {
                                let spatial = (-dist_3d * dist_3d / 0.1).exp();
                                let w_factor = (-dw * dw / 0.3).exp();
                                
                                energy += barrier * spatial * w_factor * 0.1;
                                
                                let mid_w = (-dw * dw / 0.8).exp();
                                energy -= binding * spatial * mid_w * 0.05;
                            }
                        }
                    }
                    
                    let idx = xi * self.size * self.size + yi * self.size + zi;
                    self.data[idx] = energy.max(0.0);
                }
            }
        }
    }
}

struct Simulation {
    knots: Vec<Knot>,
    energy_field: EnergyField,
    energy_barrier: f32,
    binding_strength: f32,
    repulsion_3d: f32,         // 3D repulsion strength
    w_damping: f32,
    spring_constant: f32,
    w_range: f32,
    bounds: f32,               // Spatial bounds
    paused: bool,
}

impl Simulation {
    fn new_linked() -> Self {
        let (k1, k2) = Knot::trefoil_link_pair();
        let knots = vec![k1, k2];
        
        Self {
            knots,
            energy_field: EnergyField::new(),
            energy_barrier: 0.8,
            binding_strength: 0.4,
            repulsion_3d: 0.5,
            w_damping: 0.02,
            spring_constant: 2.0,
            w_range: 3.0,
            bounds: 1.5,
            paused: false,
        }
    }
    
    fn compute_3d_force(&self, knot_idx: usize) -> Vector3<f32> {
        let knot = &self.knots[knot_idx];
        let mut force = Vector3::zeros();
        let avg_w = knot.average_w();
        
        // Repulsion from other knots (modulated by W distance)
        for (j, other) in self.knots.iter().enumerate() {
            if j == knot_idx { continue; }
            
            let delta = knot.center - other.center;
            let dist = delta.magnitude();
            let other_avg_w = other.average_w();
            let dw = (avg_w - other_avg_w).abs();
            
            // Repulsion is strong when W is similar, weak when W differs
            let w_factor = (-dw * dw / 0.3).exp();
            
            if dist < 1.2 && dist > 0.01 {
                // Repulsive force (inverse square-ish)
                let repulsion = self.repulsion_3d * w_factor / (dist * dist + 0.1);
                force += delta.normalize() * repulsion;
            }
        }
        
        // Boundary force (soft walls)
        for axis in 0..3 {
            let pos = knot.center[axis];
            if pos > self.bounds {
                force[axis] -= (pos - self.bounds) * 2.0;
            } else if pos < -self.bounds {
                force[axis] -= (pos + self.bounds) * 2.0;
            }
        }
        
        force
    }
    
    fn compute_local_w_force(&self, knot_idx: usize, point_idx: usize) -> f32 {
        let knot = &self.knots[knot_idx];
        let (pos1, w1) = knot.get_4d_point(point_idx);
        
        let mut force = 0.0;
        
        for (j, other) in self.knots.iter().enumerate() {
            if j == knot_idx { continue; }
            
            for k in 0..other.points.len() {
                let (pos2, w2) = other.get_4d_point(k);
                let dist_3d = (pos1 - pos2).magnitude();
                
                if dist_3d < 0.6 && dist_3d > 0.01 {
                    let dw = w1 - w2;
                    let w_prox = (-dw * dw / 0.2).exp();
                    let s_prox = (-dist_3d * dist_3d / 0.08).exp();
                    
                    // Barrier in W when 3D-overlapping
                    let barrier = -self.energy_barrier * dw * w_prox * s_prox * 1.5;
                    
                    // Binding at medium W
                    let mid_w = (-dw * dw / 0.6).exp() * (1.0 - w_prox);
                    let binding = self.binding_strength * (-dw) * mid_w * s_prox;
                    
                    force += barrier + binding;
                }
            }
        }
        
        // W cohesion spring
        let n = knot.points.len();
        let prev_w = knot.points[(point_idx + n - 1) % n].w;
        let next_w = knot.points[(point_idx + 1) % n].w;
        force -= self.spring_constant * (w1 - (prev_w + next_w) / 2.0);
        
        force
    }
    
    fn step(&mut self, dt: f32, w_slice: f32) {
        if self.paused { return; }
        
        // Compute W forces for each point
        let w_forces: Vec<Vec<f32>> = self.knots.iter().enumerate()
            .map(|(i, _)| {
                (0..self.knots[i].points.len())
                    .map(|j| self.compute_local_w_force(i, j))
                    .collect()
            })
            .collect();
        
        // Compute 3D forces for each knot
        let forces_3d: Vec<Vector3<f32>> = self.knots.iter().enumerate()
            .map(|(i, _)| self.compute_3d_force(i))
            .collect();
        
        // Update knot positions
        for (i, knot) in self.knots.iter_mut().enumerate() {
            // Update W for each point
            for (j, point) in knot.points.iter_mut().enumerate() {
                point.w += w_forces[i][j] * dt * 0.25;
                point.w = point.w.clamp(-self.w_range, self.w_range);
            }
            
            // Global W drift
            let avg_w = knot.average_w();
            knot.w_velocity += -0.05 * avg_w * dt;
            knot.w_velocity *= 1.0 - self.w_damping;
            for point in &mut knot.points {
                point.w += knot.w_velocity * dt;
                point.w = point.w.clamp(-self.w_range, self.w_range);
            }
            
            // Update 3D position
            knot.velocity += forces_3d[i] * dt;
            knot.velocity *= 0.98; // Damping
            knot.center += knot.velocity * dt;
            
            // Clamp to bounds
            for axis in 0..3 {
                knot.center[axis] = knot.center[axis].clamp(-self.bounds * 1.2, self.bounds * 1.2);
            }
        }
        
        self.energy_field.compute(&self.knots, w_slice, self.energy_barrier, self.binding_strength);
    }
    
    fn compute_total_energy(&self) -> f32 {
        let mut energy = 0.0;
        
        for (i, knot) in self.knots.iter().enumerate() {
            for (pi, point) in knot.points.iter().enumerate() {
                for (j, other) in self.knots.iter().enumerate() {
                    if j <= i { continue; }
                    for other_point in &other.points {
                        let world_pos = knot.center + point.pos;
                        let other_world = other.center + other_point.pos;
                        let dist_3d = (world_pos - other_world).magnitude();
                        if dist_3d < 0.6 && dist_3d > 0.01 {
                            let dw = (point.w - other_point.w).abs();
                            let w_prox = (-dw * dw / 0.2).exp();
                            let s_prox = (-dist_3d * dist_3d / 0.08).exp();
                            energy += self.energy_barrier * w_prox * s_prox * 0.3;
                        }
                    }
                }
                
                let next_w = knot.points[(pi + 1) % knot.points.len()].w;
                energy += 0.5 * self.spring_constant * (point.w - next_w).powi(2) * 0.1;
            }
        }
        
        energy
    }
}

struct FieldMesh {
    vertices: Vec<Vertex>,
    indices: Vec<u16>,
}

impl FieldMesh {
    fn from_field(field: &EnergyField, _w_slice: f32, threshold: f32, _slice_width: f32) -> Self {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();
        
        let half_extent = 1.5;
        let cell_size = 2.0 * half_extent / (field.size - 1) as f32;
        let max_cells = 3000; // Limit to prevent buffer overflow
        
        let mut cell_count = 0;
        
        // Generate field visualization as small cubes at high energy locations
        for xi in 0..field.size {
            for yi in 0..field.size {
                for zi in 0..field.size {
                    if cell_count >= max_cells { break; }
                    
                    let idx = xi * field.size * field.size + yi * field.size + zi;
                    let energy = field.data[idx];
                    
                    if energy > threshold {
                        cell_count += 1;
                        let x = (xi as f32 / (field.size - 1) as f32 - 0.5) * 2.0 * half_extent;
                        let y = (yi as f32 / (field.size - 1) as f32 - 0.5) * 2.0 * half_extent;
                        let z = (zi as f32 / (field.size - 1) as f32 - 0.5) * 2.0 * half_extent;
                        
                        let alpha = ((energy - threshold) / (1.0 - threshold)).min(1.0) * 0.4;
                        
                        let offset = vertices.len() as u16;
                        let s = cell_size * 0.35;
                        
                        let corners = [
                            [x-s, y-s, z-s], [x+s, y-s, z-s], [x+s, y+s, z-s], [x-s, y+s, z-s],
                            [x-s, y-s, z+s], [x+s, y-s, z+s], [x+s, y+s, z+s], [x-s, y+s, z+s],
                        ];
                        
                        for c in &corners {
                            vertices.push(Vertex { position: *c, normal: [0.0, 1.0, 0.0], color: [0.3, 0.5, 1.0, alpha] });
                        }
                        
                        let faces: [[u16; 4]; 6] = [
                            [0,1,2,3], [4,5,6,7], [0,1,5,4], [2,3,7,6], [0,3,7,4], [1,2,6,5],
                        ];
                        
                        for f in &faces {
                            indices.extend_from_slice(&[offset + f[0], offset + f[1], offset + f[2]]);
                            indices.extend_from_slice(&[offset + f[0], offset + f[2], offset + f[3]]);
                        }
                    }
                }
            }
        }
        
        Self { vertices, indices }
    }
}

struct AppState {
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    surface: wgpu::Surface<'static>,
    render_pipeline: wgpu::RenderPipeline,
    vertex_buffer: wgpu::Buffer,
    index_buffer: wgpu::Buffer,
    num_indices: u32,
    uniform_buffer: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    depth_texture: wgpu::Texture,
    depth_view: wgpu::TextureView,
    egui_ctx: egui::Context,
    egui_state: EguiState,
    egui_renderer: egui_wgpu::Renderer,
    window: Arc<Window>,
    
    camera_rotation: UnitQuaternion<f32>,
    camera_distance: f32,
    w_slice: f32,
    slice_width: f32,
    show_field: bool,
    field_threshold: f32,
    
    sim: Simulation,
    
    mouse_pos: Option<(f64, f64)>,
    dragging: bool,
}

impl AppState {
    fn new(window: Window) -> Self {
        let window = Arc::new(window);
        let size = window.inner_size();
        
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::PRIMARY,
            ..Default::default()
        });
        
        let surface = instance.create_surface(window.clone()).unwrap();
        
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: Some(&surface),
            force_fallback_adapter: false,
        })).unwrap();
        
        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: None,
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::default(),
            },
            None,
        )).unwrap();
        
        let caps = surface.get_capabilities(&adapter);
        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: caps.formats.iter().copied().find(|f| f.is_srgb()).unwrap_or(caps.formats[0]),
            width: size.width,
            height: size.height,
            present_mode: wgpu::PresentMode::AutoVsync,
            alpha_mode: wgpu::CompositeAlphaMode::Auto,
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);
        
        let shader_source = include_str!("shader.wgsl");
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Knot Shader"),
            source: wgpu::ShaderSource::Wgsl(shader_source.into()),
        });
        
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("Uniform Layout"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });
        
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Pipeline Layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });
        
        let render_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Knot Pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: "vs_main",
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<Vertex>() as wgpu::BufferAddress,
                    step_mode: wgpu::VertexStepMode::Vertex,
                    attributes: &[
                        wgpu::VertexAttribute { offset: 0, shader_location: 0, format: wgpu::VertexFormat::Float32x3 },
                        wgpu::VertexAttribute { offset: 12, shader_location: 1, format: wgpu::VertexFormat::Float32x3 },
                        wgpu::VertexAttribute { offset: 24, shader_location: 2, format: wgpu::VertexFormat::Float32x4 },
                    ],
                }],
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: "fs_main",
                targets: &[Some(wgpu::ColorTargetState {
                    format: config.format,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None,
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: wgpu::TextureFormat::Depth32Float,
                depth_write_enabled: true,
                depth_compare: wgpu::CompareFunction::Less,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState::default(),
            multiview: None,
        });
        
        let sim = Simulation::new_linked();
        let (vertices, indices) = Self::generate_geometry(&sim, 0.0, 0.5, false, 0.1);
        
        // Pre-allocate large buffers for dynamic geometry
        let max_vertices = 50000;
        let max_indices = 200000;
        
        let vertex_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Vertex Buffer"),
            size: (max_vertices * std::mem::size_of::<Vertex>()) as u64,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let index_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Index Buffer"),
            size: (max_indices * std::mem::size_of::<u16>()) as u64,
            usage: wgpu::BufferUsages::INDEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        
        let uniforms = Uniforms {
            mvp: Self::identity_matrix(),
            light_dir: [0.5, 0.8, 0.3],
            slice_w: 0.0,
        };
        let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Uniform Buffer"),
            contents: bytemuck::bytes_of(&uniforms),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Uniform Bind Group"),
            layout: &bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            }],
        });
        
        let depth_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Depth Texture"),
            size: wgpu::Extent3d { width: config.width, height: config.height, depth_or_array_layers: 1 },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });
        let depth_view = depth_texture.create_view(&wgpu::TextureViewDescriptor::default());
        
        let egui_ctx = egui::Context::default();
        let egui_state = EguiState::new(egui_ctx.clone(), egui::ViewportId::from_hash_of("main"), &window, None, None);
        let egui_renderer = egui_wgpu::Renderer::new(&device, config.format, None, 1);
        
        Self {
            device,
            queue,
            config,
            surface,
            render_pipeline,
            vertex_buffer,
            index_buffer,
            num_indices: indices.len() as u32,
            uniform_buffer,
            bind_group,
            depth_texture,
            depth_view,
            egui_ctx,
            egui_state,
            egui_renderer,
            window,
            camera_rotation: UnitQuaternion::from_axis_angle(&Vector3::y_axis(), 0.3),
            camera_distance: 5.0,
            w_slice: 0.0,
            slice_width: 0.5,
            show_field: true,
            field_threshold: 0.3,
            sim,
            mouse_pos: None,
            dragging: false,
        }
    }
    
    fn identity_matrix() -> [[f32; 4]; 4] {
        [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    }
    
    fn matrix_to_array(m: &Matrix4<f32>) -> [[f32; 4]; 4] {
        let cols = m.column_iter().collect::<Vec<_>>();
        [
            [cols[0].x, cols[0].y, cols[0].z, cols[0].w],
            [cols[1].x, cols[1].y, cols[1].z, cols[1].w],
            [cols[2].x, cols[2].y, cols[2].z, cols[2].w],
            [cols[3].x, cols[3].y, cols[3].z, cols[3].w],
        ]
    }
    
    fn generate_geometry(sim: &Simulation, w_slice: f32, slice_width: f32, show_field: bool, field_threshold: f32) -> (Vec<Vertex>, Vec<u16>) {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();
        
        // Knots
        for knot in &sim.knots {
            Self::generate_knot(knot, w_slice, slice_width, &mut vertices, &mut indices);
        }
        
        // Energy field
        if show_field {
            let field_mesh = FieldMesh::from_field(&sim.energy_field, w_slice, field_threshold, slice_width);
            let offset = vertices.len() as u16;
            vertices.extend(field_mesh.vertices);
            for i in field_mesh.indices {
                indices.push(offset + i);
            }
        }
        
        (vertices, indices)
    }
    
    fn generate_knot(knot: &Knot, w_slice: f32, slice_width: f32, vertices: &mut Vec<Vertex>, indices: &mut Vec<u16>) {
        let radius = 0.03;
        let n = knot.points.len();
        
        for i in 0..n {
            let (p1_pos, p1_w) = knot.get_4d_point(i);
            let (p2_pos, p2_w) = knot.get_4d_point(i + 1);
            
            let w_dist = ((p1_w + p2_w) / 2.0 - w_slice).abs();
            let alpha = (-w_dist * w_dist / (slice_width * slice_width)).exp().clamp(0.05, 1.0);
            
            if alpha < 0.05 { continue; }
            
            let color = [knot.color[0], knot.color[1], knot.color[2], alpha];
            
            let dir = p2_pos - p1_pos;
            let len = dir.magnitude();
            if len < 1e-6 { continue; }
            let dir = dir.normalize();
            
            let up = if dir.y.abs() < 0.9 { Vector3::y() } else { Vector3::x() };
            let right = dir.cross(&up).normalize();
            let up = right.cross(&dir);
            
            let start_idx = vertices.len() as u16;
            
            for j in 0..TUBE_SEGMENTS {
                let angle = (j as f32 / TUBE_SEGMENTS as f32) * std::f32::consts::TAU;
                let offset = right * (angle.cos() * radius) + up * (angle.sin() * radius);
                
                vertices.push(Vertex {
                    position: [p1_pos.x + offset.x, p1_pos.y + offset.y, p1_pos.z + offset.z],
                    normal: [offset.x, offset.y, offset.z],
                    color,
                });
            }
            
            for j in 0..TUBE_SEGMENTS {
                let angle = (j as f32 / TUBE_SEGMENTS as f32) * std::f32::consts::TAU;
                let offset = right * (angle.cos() * radius) + up * (angle.sin() * radius);
                
                vertices.push(Vertex {
                    position: [p2_pos.x + offset.x, p2_pos.y + offset.y, p2_pos.z + offset.z],
                    normal: [offset.x, offset.y, offset.z],
                    color,
                });
            }
            
            for j in 0..TUBE_SEGMENTS {
                let i0 = start_idx + j as u16;
                let i1 = start_idx + ((j + 1) % TUBE_SEGMENTS) as u16;
                let i2 = i0 + TUBE_SEGMENTS as u16;
                let i3 = i1 + TUBE_SEGMENTS as u16;
                indices.extend_from_slice(&[i0, i1, i2, i1, i3, i2]);
            }
        }
    }
    
    fn resize(&mut self, new_size: winit::dpi::PhysicalSize<u32>) {
        if new_size.width == 0 || new_size.height == 0 { return; }
        self.config.width = new_size.width;
        self.config.height = new_size.height;
        self.surface.configure(&self.device, &self.config);
        
        self.depth_texture = self.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Depth Texture"),
            size: wgpu::Extent3d { width: new_size.width, height: new_size.height, depth_or_array_layers: 1 },
            mip_level_count: 1, sample_count: 1, dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float, usage: wgpu::TextureUsages::RENDER_ATTACHMENT, view_formats: &[],
        });
        self.depth_view = self.depth_texture.create_view(&wgpu::TextureViewDescriptor::default());
    }
    
    fn update(&mut self, dt: f32) {
        self.sim.step(dt, self.w_slice);
    }
    
    fn handle_event(&mut self, event: &event::WindowEvent) {
        let egui_response = self.egui_state.on_window_event(&self.window, event);
        if egui_response.consumed { return; }
        
        match event {
            event::WindowEvent::KeyboardInput { event: KeyEvent { physical_key: PhysicalKey::Code(code), state: ElementState::Pressed, .. }, .. } => {
                match code {
                    KeyCode::Space => { self.sim.paused = !self.sim.paused; }
                    KeyCode::KeyR => { self.sim = Simulation::new_linked(); }
                    KeyCode::KeyF => { self.show_field = !self.show_field; }
                    _ => {}
                }
            }
            event::WindowEvent::MouseWheel { delta, .. } => {
                let scroll = match delta {
                    MouseScrollDelta::LineDelta(_, y) => *y,
                    MouseScrollDelta::PixelDelta(p) => p.y as f32 * 0.01,
                };
                self.camera_distance = (self.camera_distance - scroll * 0.3).clamp(2.0, 12.0);
            }
            event::WindowEvent::MouseInput { state, button, .. } => {
                if *button == MouseButton::Left { self.dragging = *state == ElementState::Pressed; }
            }
            event::WindowEvent::CursorMoved { position, .. } => {
                if self.dragging {
                    if let Some((prev_x, prev_y)) = self.mouse_pos {
                        let dx = position.x - prev_x;
                        let dy = position.y - prev_y;
                        let rot_x = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), (dx * 0.005) as f32);
                        let rot_y = UnitQuaternion::from_axis_angle(&Vector3::x_axis(), (dy * 0.005) as f32);
                        self.camera_rotation = rot_x * self.camera_rotation * rot_y;
                    }
                }
                self.mouse_pos = Some((position.x, position.y));
            }
            _ => {}
        }
    }
    
    fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        let aspect = self.config.width as f32 / self.config.height as f32;
        let proj = Matrix4::new_perspective(aspect, std::f32::consts::FRAC_PI_4, 0.1, 100.0);
        
        let base_pos = Vector3::new(0.0, 0.5, self.camera_distance);
        let cam_pos = self.camera_rotation.transform_vector(&base_pos);
        let view = Matrix4::look_at_rh(&Point3::from(cam_pos), &Point3::new(0.0, 0.0, 0.0), &Vector3::y());
        let mvp = proj * view;
        
        let uniforms = Uniforms { mvp: Self::matrix_to_array(&mvp), light_dir: [0.5, 0.8, 0.3], slice_w: self.w_slice };
        self.queue.write_buffer(&self.uniform_buffer, 0, bytemuck::bytes_of(&uniforms));
        
        let (vertices, indices) = Self::generate_geometry(&self.sim, self.w_slice, self.slice_width, self.show_field, self.field_threshold);
        self.queue.write_buffer(&self.vertex_buffer, 0, bytemuck::cast_slice(&vertices));
        self.queue.write_buffer(&self.index_buffer, 0, bytemuck::cast_slice(&indices));
        self.num_indices = indices.len() as u32;
        
        let output = self.surface.get_current_texture()?;
        let view = output.texture.create_view(&wgpu::TextureViewDescriptor::default());
        
        let mut encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("Render Encoder") });
        
        {
            let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view, resolve_target: None,
                    ops: wgpu::Operations { load: wgpu::LoadOp::Clear(wgpu::Color { r: 0.02, g: 0.02, b: 0.05, a: 1.0 }), store: wgpu::StoreOp::Store },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.depth_view,
                    depth_ops: Some(wgpu::Operations { load: wgpu::LoadOp::Clear(1.0), store: wgpu::StoreOp::Store }),
                    stencil_ops: None,
                }),
                timestamp_writes: None, occlusion_query_set: None,
            });
            
            render_pass.set_pipeline(&self.render_pipeline);
            render_pass.set_bind_group(0, &self.bind_group, &[]);
            render_pass.set_vertex_buffer(0, self.vertex_buffer.slice(..));
            render_pass.set_index_buffer(self.index_buffer.slice(..), wgpu::IndexFormat::Uint16);
            render_pass.draw_indexed(0..self.num_indices, 0, 0..1);
        }
        
        let energy = self.sim.compute_total_energy();
        let egui_input = self.egui_state.take_egui_input(&self.window);
        let egui::FullOutput { platform_output, shapes, textures_delta, .. } = 
            self.egui_ctx.run(egui_input, |ctx| {
                egui::Window::new("4D Knots").default_pos([10.0, 10.0]).show(ctx, |ui| {
                    ui.heading("Linked Knots in 4D");
                    ui.label("Two linked trefoils slip through in W");
                    ui.separator();
                    
                    ui.label("W-Slice (4D cross-section):");
                    ui.add(Slider::new(&mut self.w_slice, -2.0..=2.0).text(""));
                    ui.add(Slider::new(&mut self.slice_width, 0.1..=1.5).text("width"));
                    
                    ui.checkbox(&mut self.show_field, "Show Energy Field");
                    if self.show_field {
                        ui.add(Slider::new(&mut self.field_threshold, 0.0..=1.0).text("threshold"));
                    }
                    
                    ui.separator();
                    ui.label("Physics:");
                    ui.add(Slider::new(&mut self.sim.energy_barrier, 0.0..=3.0).text("W-Barrier"));
                    ui.add(Slider::new(&mut self.sim.binding_strength, 0.0..=2.0).text("W-Binding"));
                    ui.add(Slider::new(&mut self.sim.repulsion_3d, 0.0..=2.0).text("3D-Repulsion"));
                    
                    ui.separator();
                    if ui.button(if self.sim.paused { "Resume" } else { "Pause" }).clicked() { self.sim.paused = !self.sim.paused; }
                    if ui.button("Reset").clicked() { self.sim = Simulation::new_linked(); }
                    
                    ui.separator();
                    ui.label(format!("Energy: {:.2}", energy));
                    for (i, k) in self.sim.knots.iter().enumerate() {
                        ui.label(format!("Knot {}: pos=({:.1},{:.1},{:.1}) w={:.2} v=({:.2},{:.2},{:.2})", 
                            i+1, k.center.x, k.center.y, k.center.z, k.average_w(),
                            k.velocity.x, k.velocity.y, k.velocity.z));
                    }
                    
                    ui.separator();
                    ui.label("Space:Pause R:Reset F:Field");
                });
            });
        
        self.egui_state.handle_platform_output(&self.window, platform_output);
        
        let screen_descriptor = ScreenDescriptor {
            size_in_pixels: [self.config.width, self.config.height],
            pixels_per_point: self.window.scale_factor() as f32,
        };
        
        for (id, delta) in textures_delta.set.iter() { self.egui_renderer.update_texture(&self.device, &self.queue, *id, delta); }
        for id in textures_delta.free.iter() { self.egui_renderer.free_texture(id); }
        
        let clipped = self.egui_ctx.tessellate(shapes, screen_descriptor.pixels_per_point);
        self.egui_renderer.update_buffers(&self.device, &self.queue, &mut encoder, &clipped, &screen_descriptor);
        
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Egui"), color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view, resolve_target: None,
                    ops: wgpu::Operations { load: wgpu::LoadOp::Load, store: wgpu::StoreOp::Store },
                })],
                depth_stencil_attachment: None, timestamp_writes: None, occlusion_query_set: None,
            });
            self.egui_renderer.render(&mut pass, &clipped, &screen_descriptor);
        }
        
        self.queue.submit(std::iter::once(encoder.finish()));
        output.present();
        Ok(())
    }
}

fn main() {
    env_logger::init();
    let event_loop = EventLoop::new().unwrap();
    let window = Window::new(&event_loop).unwrap();
    window.set_title("4D Knot Visualizer");
    let _ = window.request_inner_size(dpi::LogicalSize::new(1280, 720));
    
    let mut state = AppState::new(window);
    event_loop.set_control_flow(ControlFlow::Poll);
    
    let mut last_time = instant::Instant::now();
    event_loop.run(move |event, window_target| {
        match event {
            event::Event::WindowEvent { event, window_id } if window_id == state.window.id() => {
                match event {
                    event::WindowEvent::CloseRequested => { window_target.exit(); }
                    event::WindowEvent::Resized(size) => { state.resize(size); }
                    event::WindowEvent::RedrawRequested => {
                        let now = instant::Instant::now();
                        let dt = (now - last_time).as_secs_f32().min(0.1);
                        last_time = now;
                        state.update(dt);
                        match state.render() {
                            Ok(_) => {}
                            Err(wgpu::SurfaceError::Lost) => state.resize(state.window.inner_size()),
                            Err(wgpu::SurfaceError::OutOfMemory) => window_target.exit(),
                            Err(e) => eprintln!("{:?}", e),
                        }
                    }
                    _ => { state.handle_event(&event); }
                }
            }
            event::Event::AboutToWait => { state.window.request_redraw(); }
            _ => {}
        }
    }).unwrap();
}
