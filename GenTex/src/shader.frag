#version 150 core
in vec2 coord;
out vec4 fragment;

void main() {
    fragment = vec4(0.5 + coord.x, 0.5 + coord.y, 0.0, 1.0);
} 
