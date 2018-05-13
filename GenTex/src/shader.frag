#version 150 core
in vec2 coord;
in vec3 vPosition;
out vec4 fragment;

void main() {
    fragment = vec4(vPosition, 1.0);
} 
