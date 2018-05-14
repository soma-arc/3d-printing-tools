#version 150 core

in vec4 uv;
in vec4 position;
out vec2 coord;
out vec3 vPosition;

void main() {
    const vec2 one = vec2(1., 1.);
    coord = (uv.xy * 2 - one);
    vPosition = position.xyz;
    gl_Position = vec4(coord, uv.zw);
} 
