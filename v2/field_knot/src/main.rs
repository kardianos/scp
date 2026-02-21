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

const GRID_SIZE: usize = 48;
const FIELD_VOXEL_SCALE: usize = 2;

#[derive(Clone, Copy, Pod, Zeroable)]
#[repr(C)]
struct Vertex {
    position: [f32; 3],
    color: [f32; 4],
}

#[derive(Clone, Copy, Pod, Zeroable)]
#[repr(C)]
struct Uniforms {
    mvp: [[f32; 4]; 4],
    time: f32,
    _pad: [f32; 7],
}

struct ComplexField {
    real: Vec<f32>,
    imag: Vec<f32>,
    n: usize,
    extent: f32,
}

impl ComplexField {
    fn new(n: usize, extent: f32) -> Self {
        let size = n * n * n;
        Self {
            real: vec![0.0; size],
            imag: vec![0.0; size],
            n,
            extent,
        }
    }

    fn initialize_random(&mut self, amplitude: f32) {
        use std::f64::consts::TAU;
        let n = self.n;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let idx = i * n * n + j * n + k;
                    let phase = rand_phase();
                    self.real[idx] = amplitude * phase.cos();
                    self.imag[idx] = amplitude * phase.sin();
                }
            }
        }
    }

    fn initialize_dense(&mut self, center_amplitude: f32) {
        let n = self.n;
        let half = self.extent;

        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let x = (i as f32 / (n - 1) as f32 - 0.5) * 2.0 * half;
                    let y = (j as f32 / (n - 1) as f32 - 0.5) * 2.0 * half;
                    let z = (k as f32 / (n - 1) as f32 - 0.5) * 2.0 * half;

                    let r = (x * x + y * y + z * z).sqrt();
                    // Use a larger radius for more spread
                    let envelope = (-r * r / 4.0).exp();

                    let phase = rand_phase();
                    let amp = center_amplitude * envelope;

                    let idx = i * n * n + j * n + k;
                    self.real[idx] = amp * phase.cos();
                    self.imag[idx] = amp * phase.sin();
                }
            }
        }
    }

    fn density(&self) -> Vec<f32> {
        self.real
            .iter()
            .zip(self.imag.iter())
            .map(|(r, i)| r * r + i * i)
            .collect()
    }

    fn phase(&self) -> Vec<f32> {
        self.real
            .iter()
            .zip(self.imag.iter())
            .map(|(r, i)| i.atan2(*r))
            .collect()
    }

    fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        let n = self.n;
        ((i + n) % n) * n * n + ((j + n) % n) * n + ((k + n) % n)
    }

    fn laplacian(&self, i: usize, j: usize, k: usize) -> (f32, f32) {
        let n = self.n;
        let idx_center = self.idx(i, j, k);

        let (r_c, i_c) = (self.real[idx_center], self.imag[idx_center]);

        let neighbors = [
            self.idx(i + 1, j, k),
            self.idx(i - 1, j, k),
            self.idx(i, j + 1, k),
            self.idx(i, j - 1, k),
            self.idx(i, j, k + 1),
            self.idx(i, j, k - 1),
        ];

        let mut r_sum = 0.0f32;
        let mut i_sum = 0.0f32;

        for &idx in &neighbors {
            r_sum += self.real[idx];
            i_sum += self.imag[idx];
        }

        (r_sum - 6.0 * r_c, i_sum - 6.0 * i_c)
    }

    fn amplitude_at(&self, i: usize, j: usize, k: usize) -> f32 {
        let idx = self.idx(i, j, k);
        (self.real[idx] * self.real[idx] + self.imag[idx] * self.imag[idx]).sqrt()
    }

    fn phase_at(&self, i: usize, j: usize, k: usize) -> f32 {
        let idx = self.idx(i, j, k);
        self.imag[idx].atan2(self.real[idx])
    }
}

fn rand_phase() -> f32 {
    use std::f64::consts::TAU;
    (fastrand::f64() * TAU) as f32
}

struct FieldSimulation {
    field: ComplexField,
    temperature: f32,
    temperature_initial: f32,
    temperature_final: f32,
    cooling_steps: usize,
    step_count: usize,
    dt: f32,
    diffusion: f32,
    noise_scale: f32,
    paused: bool,

    density: Vec<f32>,
    phase: Vec<f32>,
    winding: Vec<f32>,
    defects: Vec<(f32, f32, f32)>,

    history_defects: Vec<usize>,
    history_temperature: Vec<f32>,
}

impl FieldSimulation {
    fn new() -> Self {
        let mut field = ComplexField::new(GRID_SIZE, 2.0);
        field.initialize_dense(0.8);

        let density = field.density();
        let phase = field.phase();
        let winding = vec![0.0; density.len()];

        Self {
            field,
            temperature: 0.5,
            temperature_initial: 0.5,
            temperature_final: 0.01,
            cooling_steps: 400,
            step_count: 0,
            dt: 0.02,
            diffusion: 0.15,
            noise_scale: 0.01,
            paused: false,
            density,
            phase,
            winding,
            defects: Vec::new(),
            history_defects: Vec::new(),
            history_temperature: Vec::new(),
        }
    }

    fn reset(&mut self) {
        self.field.initialize_dense(0.8);
        self.temperature = self.temperature_initial;
        self.step_count = 0;
        self.defects.clear();
        self.history_defects.clear();
        self.history_temperature.clear();
        self.update_derived();
    }

    fn update_derived(&mut self) {
        self.density = self.field.density();
        self.phase = self.field.phase();
        self.compute_winding();
        self.detect_defects();
    }

    fn compute_winding(&mut self) {
        let n = self.field.n;

        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let idx = self.field.idx(i, j, k);
                    let rho = self.density[idx];

                    if rho < 0.1 {
                        self.winding[idx] = 0.0;
                        continue;
                    }

                    let ip = self.field.idx(i + 1, j, k);
                    let jp = self.field.idx(i, j + 1, k);
                    let kp = self.field.idx(i, j, k + 1);

                    let mut dpx = self.phase[ip] - self.phase[idx];
                    let mut dpy = self.phase[jp] - self.phase[idx];
                    let mut dpz = self.phase[kp] - self.phase[idx];

                    if dpx > std::f32::consts::PI {
                        dpx -= std::f32::consts::TAU;
                    }
                    if dpx < -std::f32::consts::PI {
                        dpx += std::f32::consts::TAU;
                    }
                    if dpy > std::f32::consts::PI {
                        dpy -= std::f32::consts::TAU;
                    }
                    if dpy < -std::f32::consts::PI {
                        dpy += std::f32::consts::TAU;
                    }
                    if dpz > std::f32::consts::PI {
                        dpz -= std::f32::consts::TAU;
                    }
                    if dpz < -std::f32::consts::PI {
                        dpz += std::f32::consts::TAU;
                    }

                    self.winding[idx] = (dpx.abs() + dpy.abs() + dpz.abs()) * rho.sqrt();
                }
            }
        }
    }

    fn detect_defects(&mut self) {
        self.defects.clear();
        let n = self.field.n;
        let half = self.field.extent;

        let threshold = {
            let mut w_sorted: Vec<f32> = self.winding.iter().copied().collect();
            w_sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());
            w_sorted[w_sorted.len() / 20]
        };

        for i in 1..n - 1 {
            for j in 1..n - 1 {
                for k in 1..n - 1 {
                    let idx = self.field.idx(i, j, k);
                    let w = self.winding[idx];

                    if w > threshold
                        && w >= self.winding[self.field.idx(i - 1, j, k)]
                        && w >= self.winding[self.field.idx(i + 1, j, k)]
                        && w >= self.winding[self.field.idx(i, j - 1, k)]
                        && w >= self.winding[self.field.idx(i, j + 1, k)]
                        && w >= self.winding[self.field.idx(i, j, k - 1)]
                        && w >= self.winding[self.field.idx(i, j, k + 1)]
                    {
                        let x = (i as f32 / (n - 1) as f32 - 0.5) * 2.0 * half;
                        let y = (j as f32 / (n - 1) as f32 - 0.5) * 2.0 * half;
                        let z = (k as f32 / (n - 1) as f32 - 0.5) * 2.0 * half;
                        self.defects.push((x, y, z));
                    }
                }
            }
        }
    }

    fn step(&mut self) {
        if self.paused {
            return;
        }

        self.step_count += 1;

        if self.step_count < self.cooling_steps {
            let t = self.step_count as f32 / self.cooling_steps as f32;
            self.temperature = self.temperature_initial * (1.0 - t) + self.temperature_final * t;
        } else {
            self.temperature = self.temperature_final;
        }

        let n = self.field.n;
        let dt = self.dt;
        let diff = self.diffusion;
        let T = self.temperature;
        let noise = self.noise_scale * T.sqrt();

        let mut new_real = vec![0.0f32; n * n * n];
        let mut new_imag = vec![0.0f32; n * n * n];

        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let idx = self.field.idx(i, j, k);

                    let rho = self.density[idx];
                    let (lap_r, lap_i) = self.field.laplacian(i, j, k);

                    let nonlinear_r = (1.0 - rho) * self.field.real[idx];
                    let nonlinear_i = (1.0 - rho) * self.field.imag[idx];

                    let thermal_r = -T * self.field.real[idx];
                    let thermal_i = -T * self.field.imag[idx];

                    let noise_r = noise * (fastrand::f64() * 2.0 - 1.0) as f32;
                    let noise_i = noise * (fastrand::f64() * 2.0 - 1.0) as f32;

                    new_real[idx] = self.field.real[idx]
                        + dt * (diff * lap_r + nonlinear_r + thermal_r + noise_r);
                    new_imag[idx] = self.field.imag[idx]
                        + dt * (diff * lap_i + nonlinear_i + thermal_i + noise_i);
                }
            }
        }

        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let idx = self.field.idx(i, j, k);
                    let amp = (new_real[idx] * new_real[idx] + new_imag[idx] * new_imag[idx])
                        .sqrt()
                        .min(3.0);
                    let phase = new_imag[idx].atan2(new_real[idx]);
                    self.field.real[idx] = amp * phase.cos();
                    self.field.imag[idx] = amp * phase.sin();
                }
            }
        }

        self.update_derived();

        if self.step_count % 10 == 0 {
            self.history_defects.push(self.defects.len());
            self.history_temperature.push(self.temperature);
        }
    }

    fn mean_density(&self) -> f32 {
        self.density.iter().sum::<f32>() / self.density.len() as f32
    }

    fn max_density(&self) -> f32 {
        self.density.iter().cloned().fold(0.0, f32::max)
    }
}

struct FieldMesh {
    vertices: Vec<Vertex>,
    indices: Vec<u32>,
}

impl FieldMesh {
    fn add_quad(
        vertices: &mut Vec<Vertex>,
        indices: &mut Vec<u32>,
        x: f32,
        y: f32,
        z: f32,
        size: f32,
        color: [f32; 4],
    ) {
        let base = vertices.len() as u32;

        // Create a small axis-aligned box instead of a quad
        let s = size;
        let corners = [
            [x - s, y - s, z - s],
            [x + s, y - s, z - s],
            [x + s, y + s, z - s],
            [x - s, y + s, z - s],
            [x - s, y - s, z + s],
            [x + s, y - s, z + s],
            [x + s, y + s, z + s],
            [x - s, y + s, z + s],
        ];

        for c in &corners {
            vertices.push(Vertex {
                position: *c,
                color,
            });
        }

        // 12 triangles for a cube
        let faces: [[u32; 4]; 6] = [
            [0, 1, 2, 3], // front
            [5, 4, 7, 6], // back
            [4, 0, 3, 7], // left
            [1, 5, 6, 2], // right
            [3, 2, 6, 7], // top
            [4, 5, 1, 0], // bottom
        ];

        for f in &faces {
            indices.extend_from_slice(&[base + f[0], base + f[1], base + f[2]]);
            indices.extend_from_slice(&[base + f[0], base + f[2], base + f[3]]);
        }
    }

    fn from_defects(defects: &[(f32, f32, f32)], _point_size: f32) -> Self {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        for &(x, y, z) in defects {
            Self::add_quad(
                &mut vertices,
                &mut indices,
                x,
                y,
                z,
                0.03,
                [1.0, 0.3, 0.1, 0.9],
            );
        }

        Self { vertices, indices }
    }

    fn from_density_field(
        density: &[f32],
        phase: &[f32],
        n: usize,
        extent: f32,
        threshold: f32,
        max_points: usize,
    ) -> Self {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        let half = extent;
        let mut count = 0;

        // Sample more uniformly through 3D space
        let step = if n > 30 { 2 } else { 1 };

        for i in (0..n).step_by(step) {
            for j in (0..n).step_by(step) {
                for k in (0..n).step_by(step) {
                    if count >= max_points {
                        break;
                    }

                    let idx = i * n * n + j * n + k;
                    let rho = density[idx];

                    if rho > threshold {
                        count += 1;
                        let x = (i as f32 / (n - 1) as f32 - 0.5) * 2.0 * half;
                        let y = (j as f32 / (n - 1) as f32 - 0.5) * 2.0 * half;
                        let z = (k as f32 / (n - 1) as f32 - 0.5) * 2.0 * half;

                        let intensity = ((rho - threshold) / (1.5 - threshold)).min(1.0);
                        let phi = phase[idx];

                        // Color varies with phase AND z-position for 3D visibility
                        let z_norm = (z / extent + 1.0) / 2.0; // 0 to 1
                        let r = 0.3 + 0.5 * z_norm + 0.2 * (phi / std::f32::consts::PI + 1.0) / 2.0;
                        let g = 0.2 + 0.4 * intensity;
                        let b = 0.8 - 0.3 * z_norm - 0.4 * (phi / std::f32::consts::PI + 1.0) / 2.0;
                        let alpha = 0.4 + 0.4 * intensity;

                        let size = 0.02 + 0.025 * intensity;
                        Self::add_quad(
                            &mut vertices,
                            &mut indices,
                            x,
                            y,
                            z,
                            size,
                            [r, g, b, alpha],
                        );
                    }
                }
            }
        }

        Self { vertices, indices }
    }

    fn from_winding_field(
        winding: &[f32],
        n: usize,
        extent: f32,
        threshold_fraction: f32,
        max_points: usize,
    ) -> Self {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        let mut w_sorted: Vec<f32> = winding.iter().copied().collect();
        w_sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());
        let threshold = w_sorted[(w_sorted.len() as f32 * threshold_fraction) as usize];

        let half = extent;
        let mut count = 0;

        let step = if n > 30 { 2 } else { 1 };

        for i in (0..n).step_by(step) {
            for j in (0..n).step_by(step) {
                for k in (0..n).step_by(step) {
                    if count >= max_points {
                        break;
                    }

                    let idx = i * n * n + j * n + k;
                    let w = winding[idx];

                    if w > threshold {
                        count += 1;
                        let x = (i as f32 / (n - 1) as f32 - 0.5) * 2.0 * half;
                        let y = (j as f32 / (n - 1) as f32 - 0.5) * 2.0 * half;
                        let z = (k as f32 / (n - 1) as f32 - 0.5) * 2.0 * half;

                        let intensity = (w / (threshold * 3.0)).min(1.0);
                        let z_norm = (z / half + 1.0) / 2.0;
                        let size = 0.02 + 0.03 * intensity;
                        // Orange-yellow color varying with z
                        Self::add_quad(
                            &mut vertices,
                            &mut indices,
                            x,
                            y,
                            z,
                            size,
                            [1.0, 0.3 + 0.4 * z_norm, 0.0, 0.8],
                        );
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
    render_pipeline_points: wgpu::RenderPipeline,
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
    point_size: f32,
    show_density: bool,
    show_defects: bool,
    show_winding: bool,
    density_threshold: f32,

    sim: FieldSimulation,

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
        }))
        .unwrap();

        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: None,
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::default(),
            },
            None,
        ))
        .unwrap();

        let caps = surface.get_capabilities(&adapter);
        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: caps
                .formats
                .iter()
                .copied()
                .find(|f| f.is_srgb())
                .unwrap_or(caps.formats[0]),
            width: size.width,
            height: size.height,
            present_mode: wgpu::PresentMode::AutoVsync,
            alpha_mode: wgpu::CompositeAlphaMode::Auto,
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        let shader_source = include_str!("shader_field.wgsl");
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Field Shader"),
            source: wgpu::ShaderSource::Wgsl(shader_source.into()),
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("Uniform Layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Pipeline Layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let render_pipeline_points =
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("Field Pipeline"),
                layout: Some(&pipeline_layout),
                vertex: wgpu::VertexState {
                    module: &shader,
                    entry_point: "vs_main",
                    buffers: &[wgpu::VertexBufferLayout {
                        array_stride: std::mem::size_of::<Vertex>() as wgpu::BufferAddress,
                        step_mode: wgpu::VertexStepMode::Vertex,
                        attributes: &[
                            wgpu::VertexAttribute {
                                offset: 0,
                                shader_location: 0,
                                format: wgpu::VertexFormat::Float32x3,
                            },
                            wgpu::VertexAttribute {
                                offset: 12,
                                shader_location: 1,
                                format: wgpu::VertexFormat::Float32x4,
                            },
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

        let sim = FieldSimulation::new();

        let max_vertices = 500000;
        let vertex_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Vertex Buffer"),
            size: (max_vertices * std::mem::size_of::<Vertex>()) as u64,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let max_indices = 2000000;
        let index_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Index Buffer"),
            size: (max_indices * std::mem::size_of::<u32>()) as u64,
            usage: wgpu::BufferUsages::INDEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let uniforms = Uniforms {
            mvp: Self::identity_matrix(),
            time: 0.0,
            _pad: [0.0; 7],
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
            size: wgpu::Extent3d {
                width: config.width,
                height: config.height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });
        let depth_view = depth_texture.create_view(&wgpu::TextureViewDescriptor::default());

        let egui_ctx = egui::Context::default();
        let egui_state = EguiState::new(
            egui_ctx.clone(),
            egui::ViewportId::from_hash_of("main"),
            &window,
            None,
            None,
        );
        let egui_renderer = egui_wgpu::Renderer::new(&device, config.format, None, 1);

        Self {
            device,
            queue,
            config,
            surface,
            render_pipeline_points,
            vertex_buffer,
            index_buffer,
            num_indices: 0,
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
            point_size: 6.0,
            show_density: true,
            show_defects: true,
            show_winding: false,
            density_threshold: 0.005,
            sim,
            mouse_pos: None,
            dragging: false,
        }
    }

    fn identity_matrix() -> [[f32; 4]; 4] {
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
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

    fn generate_geometry(&mut self) -> (Vec<Vertex>, Vec<u32>) {
        let mut vertices = Vec::new();
        let mut indices = Vec::new();

        if self.show_defects {
            let defect_mesh = FieldMesh::from_defects(&self.sim.defects, 0.05);
            let offset = vertices.len() as u32;
            vertices.extend(defect_mesh.vertices);
            for i in defect_mesh.indices {
                indices.push(offset + i);
            }
        }

        if self.show_winding {
            let winding_mesh = FieldMesh::from_winding_field(
                &self.sim.winding,
                self.sim.field.n,
                self.sim.field.extent,
                0.97,
                15000,
            );
            let offset = vertices.len() as u32;
            vertices.extend(winding_mesh.vertices);
            for i in winding_mesh.indices {
                indices.push(offset + i);
            }
        }

        if self.show_density {
            let density_mesh = FieldMesh::from_density_field(
                &self.sim.density,
                &self.sim.phase,
                self.sim.field.n,
                self.sim.field.extent,
                self.density_threshold,
                10000,
            );
            let offset = vertices.len() as u32;
            vertices.extend(density_mesh.vertices);
            for i in density_mesh.indices {
                indices.push(offset + i);
            }
        }

        (vertices, indices)
    }

    fn resize(&mut self, new_size: winit::dpi::PhysicalSize<u32>) {
        if new_size.width == 0 || new_size.height == 0 {
            return;
        }
        self.config.width = new_size.width;
        self.config.height = new_size.height;
        self.surface.configure(&self.device, &self.config);

        self.depth_texture = self.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Depth Texture"),
            size: wgpu::Extent3d {
                width: new_size.width,
                height: new_size.height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });
        self.depth_view = self
            .depth_texture
            .create_view(&wgpu::TextureViewDescriptor::default());
    }

    fn update(&mut self, _dt: f32) {
        self.sim.step();
    }

    fn handle_event(&mut self, event: &event::WindowEvent) {
        // Check if egui wants this input BEFORE processing
        let egui_wants_pointer = self.egui_ctx.wants_pointer_input();
        let egui_wants_keyboard = self.egui_ctx.wants_keyboard_input();

        let egui_response = self.egui_state.on_window_event(&self.window, event);

        match event {
            event::WindowEvent::KeyboardInput {
                event:
                    KeyEvent {
                        physical_key: PhysicalKey::Code(code),
                        state: ElementState::Pressed,
                        ..
                    },
                ..
            } => {
                // Only handle if egui doesn't want keyboard
                if !egui_wants_keyboard && !egui_response.consumed {
                    match code {
                        KeyCode::Space => {
                            self.sim.paused = !self.sim.paused;
                        }
                        KeyCode::KeyR => {
                            self.sim.reset();
                        }
                        KeyCode::KeyD => {
                            self.show_density = !self.show_density;
                        }
                        KeyCode::KeyF => {
                            self.show_defects = !self.show_defects;
                        }
                        KeyCode::KeyW => {
                            self.show_winding = !self.show_winding;
                        }
                        _ => {}
                    }
                }
            }
            event::WindowEvent::MouseWheel { delta, .. } => {
                // Only zoom if egui doesn't want pointer
                if !egui_wants_pointer {
                    let scroll = match delta {
                        MouseScrollDelta::LineDelta(_, y) => *y,
                        MouseScrollDelta::PixelDelta(p) => p.y as f32 * 0.01,
                    };
                    self.camera_distance = (self.camera_distance - scroll * 0.3).clamp(2.0, 12.0);
                }
            }
            event::WindowEvent::MouseInput { state, button, .. } => {
                // Only start drag if egui doesn't want pointer
                if *button == MouseButton::Left && !egui_wants_pointer {
                    self.dragging = *state == ElementState::Pressed;
                }
            }
            event::WindowEvent::CursorMoved { position, .. } => {
                // Only rotate if dragging AND egui doesn't want pointer
                if self.dragging && !egui_wants_pointer {
                    if let Some((prev_x, prev_y)) = self.mouse_pos {
                        let dx = position.x - prev_x;
                        let dy = position.y - prev_y;
                        let rot_x = UnitQuaternion::from_axis_angle(
                            &Vector3::y_axis(),
                            (dx * 0.005) as f32,
                        );
                        let rot_y = UnitQuaternion::from_axis_angle(
                            &Vector3::x_axis(),
                            (dy * 0.005) as f32,
                        );
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

        let base_pos = Vector3::new(0.0, 0.3, self.camera_distance);
        let cam_pos = self.camera_rotation.transform_vector(&base_pos);
        let view = Matrix4::look_at_rh(
            &Point3::from(cam_pos),
            &Point3::new(0.0, 0.0, 0.0),
            &Vector3::y(),
        );
        let mvp = proj * view;

        let uniforms = Uniforms {
            mvp: Self::matrix_to_array(&mvp),
            time: self.sim.step_count as f32 * self.sim.dt,
            _pad: [0.0; 7],
        };
        self.queue
            .write_buffer(&self.uniform_buffer, 0, bytemuck::bytes_of(&uniforms));

        let (vertices, indices) = self.generate_geometry();
        self.queue
            .write_buffer(&self.vertex_buffer, 0, bytemuck::cast_slice(&vertices));
        self.queue
            .write_buffer(&self.index_buffer, 0, bytemuck::cast_slice(&indices));
        self.num_indices = indices.len() as u32;

        let output = self.surface.get_current_texture()?;
        let view = output
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Render Encoder"),
            });

        {
            let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.02,
                            g: 0.02,
                            b: 0.05,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });

            render_pass.set_pipeline(&self.render_pipeline_points);
            render_pass.set_bind_group(0, &self.bind_group, &[]);
            render_pass.set_vertex_buffer(0, self.vertex_buffer.slice(..));
            render_pass.set_index_buffer(self.index_buffer.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..self.num_indices, 0, 0..1);
        }

        let egui_input = self.egui_state.take_egui_input(&self.window);
        let egui::FullOutput {
            platform_output,
            shapes,
            textures_delta,
            ..
        } = self.egui_ctx.run(egui_input, |ctx| {
            egui::Window::new("Emergent Field Knots")
                .default_pos([10.0, 10.0])
                .show(ctx, |ui| {
                    ui.heading("Phase Transition Dynamics");
                    ui.label(format!(
                        "Step: {}, T={:.3}",
                        self.sim.step_count, self.sim.temperature
                    ));
                    ui.separator();

                    ui.label("Visualization:");
                    ui.checkbox(&mut self.show_defects, "Defects (Knots)");
                    ui.checkbox(&mut self.show_density, "Field Density");
                    ui.checkbox(&mut self.show_winding, "Phase Winding");

                    if self.show_density {
                        ui.add(
                            Slider::new(&mut self.density_threshold, 0.001..=0.5)
                                .text("density threshold"),
                        );
                    }

                    ui.separator();
                    ui.label("Simulation Parameters:");
                    ui.add(Slider::new(&mut self.sim.diffusion, 0.01..=0.5).text("diffusion"));
                    ui.add(Slider::new(&mut self.sim.noise_scale, 0.0..=0.1).text("noise"));
                    ui.add(
                        Slider::new(&mut self.sim.temperature_initial, 0.1..=1.0).text("T initial"),
                    );

                    ui.separator();
                    if ui
                        .button(if self.sim.paused { "Resume" } else { "Pause" })
                        .clicked()
                    {
                        self.sim.paused = !self.sim.paused;
                    }
                    if ui.button("Reset").clicked() {
                        self.sim.reset();
                    }

                    ui.separator();
                    ui.label(format!("Defects: {}", self.sim.defects.len()));
                    ui.label(format!("ρ_mean: {:.3}", self.sim.mean_density()));
                    ui.label(format!("ρ_max: {:.3}", self.sim.max_density()));

                    ui.separator();
                    ui.label("Space:Pause R:Reset");
                    ui.label("D:Density F:Defects W:Winding");
                });
        });

        self.egui_state
            .handle_platform_output(&self.window, platform_output);

        let screen_descriptor = ScreenDescriptor {
            size_in_pixels: [self.config.width, self.config.height],
            pixels_per_point: self.window.scale_factor() as f32,
        };

        for (id, delta) in textures_delta.set.iter() {
            self.egui_renderer
                .update_texture(&self.device, &self.queue, *id, delta);
        }
        for id in textures_delta.free.iter() {
            self.egui_renderer.free_texture(id);
        }

        let clipped = self
            .egui_ctx
            .tessellate(shapes, screen_descriptor.pixels_per_point);
        self.egui_renderer.update_buffers(
            &self.device,
            &self.queue,
            &mut encoder,
            &clipped,
            &screen_descriptor,
        );

        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Egui"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            self.egui_renderer
                .render(&mut pass, &clipped, &screen_descriptor);
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
    window.set_title("Emergent Field Knots - Phase Transition Simulation");
    let _ = window.request_inner_size(dpi::LogicalSize::new(1280, 720));

    let mut state = AppState::new(window);
    event_loop.set_control_flow(ControlFlow::Poll);

    let mut last_time = instant::Instant::now();
    event_loop
        .run(move |event, window_target| match event {
            event::Event::WindowEvent { event, window_id } if window_id == state.window.id() => {
                match event {
                    event::WindowEvent::CloseRequested => {
                        window_target.exit();
                    }
                    event::WindowEvent::Resized(size) => {
                        state.resize(size);
                    }
                    event::WindowEvent::RedrawRequested => {
                        let now = instant::Instant::now();
                        let dt = (now - last_time).as_secs_f32().min(0.05);
                        last_time = now;
                        state.update(dt);
                        match state.render() {
                            Ok(_) => {}
                            Err(wgpu::SurfaceError::Lost) => {
                                state.resize(state.window.inner_size())
                            }
                            Err(wgpu::SurfaceError::OutOfMemory) => window_target.exit(),
                            Err(e) => eprintln!("{:?}", e),
                        }
                    }
                    _ => {
                        state.handle_event(&event);
                    }
                }
            }
            event::Event::AboutToWait => {
                state.window.request_redraw();
            }
            _ => {}
        })
        .unwrap();
}
