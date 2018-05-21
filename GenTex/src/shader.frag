#version 150 core
in vec2 coord;
in vec3 vPosition;
out vec4 fragment;

struct Sphere {
    vec3 center;
    vec2 r;
};

struct Plane {
    vec3 origin;
    vec3 normal;
};

uniform Sphere u_prismSpheres[4];
uniform int u_numPrismSpheres;
uniform Plane u_prismPlanes[4];
uniform int u_numPrismPlanes;

uniform Sphere u_boundingSphere;
uniform Sphere u_finiteSpheres[8];
uniform int u_numFiniteSpheres;

uniform vec3 u_bboxMin;
uniform bool u_isFinite;

vec3 Hsv2rgb(float h, float s, float v){
    vec3 c = vec3(h, s, v);
    const vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

const float DISPLAY_GAMMA_COEFF = 1. / 2.2;
vec4 GammaCorrect(vec4 rgba) {
    return vec4((min(pow(rgba.r, DISPLAY_GAMMA_COEFF), 1.)),
                (min(pow(rgba.g, DISPLAY_GAMMA_COEFF), 1.)),
                (min(pow(rgba.b, DISPLAY_GAMMA_COEFF), 1.)),
                rgba.a);
}

vec3 SphereInvert(inout vec3 pos, vec3 center, vec2 r) {
    vec3 diff = pos - center;
    float lenSq = dot(diff, diff);
    float k = r.y / lenSq;
    pos = (diff * k) + center;
    return pos;
}

const int MAX_ITERATIONS = 100;
float IIS(vec3 pos) {
    int invNum = 0;
    float d;
    for(int i = 0; i < MAX_ITERATIONS; i++) {
        bool inFund = true;
		for(int n = 0; n < u_numPrismSpheres; n++) {
            if(distance(pos, u_prismSpheres[n].center) < u_prismSpheres[n].r.x) {
                invNum++;
                SphereInvert(pos,
                             u_prismSpheres[n].center,
                             u_prismSpheres[n].r);
                inFund = false;
            }
		}

        for (int n = 0; n < u_numPrismPlanes; n++){
            pos -= u_prismPlanes[n].origin;
            d = dot(pos, u_prismPlanes[n].normal);
            if(d > 0.) {
                invNum++;
                pos -= 2. * d * u_prismPlanes[n].normal;
                inFund = false;
            }
            pos += u_prismPlanes[n].origin;
		}

        if(inFund) break;
    }
    return float(invNum);
}

void main() {
    float n = IIS(vPosition + u_bboxMin);
    fragment = vec4(Hsv2rgb(n * 0.01, 1., 1.), 1.0);
}
