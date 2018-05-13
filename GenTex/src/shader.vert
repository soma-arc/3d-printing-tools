#version 150 core

in vec4 uv;
in vec4 position;
out vec2 coord;
out vec3 vPosition;

void main() {
    coord = uv.xy;
    vPosition = position.xyz;
    gl_Position = uv;
} 
