#version 150 core
in vec4 position;
out vec2 coord;

void main() {
    coord = position.xy;
    gl_Position = position;
} 
