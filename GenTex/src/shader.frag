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
uniform Plane u_boundingPlanes[8];
uniform int u_numBoundingPlanes;
uniform Sphere u_inversionSphere;
uniform Sphere u_spheirahedraSpheres[8];
uniform int u_numDividePlanes;
uniform Plane u_dividePlanes[1];

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

void SphereInvert(inout vec3 pos, vec3 center, vec2 r) {
    vec3 diff = pos - center;
    float lenSq = dot(diff, diff);
    float k = r.y / lenSq;
    pos = (diff * k) + center;
}

const int MAX_ITERATIONS = 30;
float IIS(vec3 pos) {
    float invNum = 0.;

    float d;
    for(int i = 0; i < MAX_ITERATIONS; i++) {
        bool inFund = true;
		for(int n = 0; n < u_numPrismSpheres; n++) {
            if(distance(pos, u_prismSpheres[n].center) < u_prismSpheres[n].r.x) {
                invNum++;
                SphereInvert(pos,
                             u_prismSpheres[n].center,
                             u_prismSpheres[n].r);
                continue;
            }
		}

        for (int n = 0; n < u_numPrismPlanes; n++){
            pos -= u_prismPlanes[n].origin;
            d = dot(pos, u_prismPlanes[n].normal);
            if(d > 0.) {
                invNum++;
                pos -= 2. * d * u_prismPlanes[n].normal;
                pos += u_prismPlanes[n].origin;
                continue;
            }
            pos += u_prismPlanes[n].origin;
		}

        if(inFund) break;
    }
    return  invNum;
}

void main() {
    float n = float(IIS(vPosition));
    fragment = vec4(Hsv2rgb(-0.13 + n * 0.01, 1., 1.), 1.0);
}
