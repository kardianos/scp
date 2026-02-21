struct Uniforms {
    mvp: mat4x4<f32>,
    light_dir: vec3<f32>,
    time: f32,
}

@group(0) @binding(0) var<uniform> uniforms: Uniforms;

struct VertexInput {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) color: vec4<f32>,
}

struct VertexOutput {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) world_pos: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) color: vec4<f32>,
}

@vertex
fn vs_main(input: VertexInput) -> VertexOutput {
    var output: VertexOutput;
    let world_pos = vec4<f32>(input.position, 1.0);
    output.clip_pos = uniforms.mvp * world_pos;
    output.world_pos = input.position;
    output.normal = input.normal;
    output.color = input.color;
    return output;
}

@fragment
fn fs_main(input: VertexOutput) -> @location(0) vec4<f32> {
    let normal = normalize(input.normal);
    let light_dir = normalize(uniforms.light_dir);
    
    let ambient = 0.3;
    let diffuse = max(dot(normal, light_dir), 0.0) * 0.5;
    
    let lighting = ambient + diffuse;
    let final_color = input.color.rgb * lighting;
    
    return vec4<f32>(final_color, input.color.a);
}
