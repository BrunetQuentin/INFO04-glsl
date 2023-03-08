// screen resolution
#define R iResolution.xy

// minimum distance to objects
#define DIST_MIN 1.

// maximum distance to objects
#define DIST_MAX 50.0

// max number of steps for the ray-marching
#define RAY_MARCH_STEPS 100

// consider hit if we reach this distance
#define RAY_MARCH_PRECI 0.0001

// for ray direction computation
#define PI 3.14159265359

// ray structure
struct Ray {
    vec3 o; // origin
    vec3 d; // direction
};

struct Surface {
    float t; // surface distance
    vec3 c; // surface color
};

float sdSphere( vec3 p, float s ) {
  return length(p)-s;
}

float sdBox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

Surface add(in Surface s1, in Surface s2){
    if(s1.t<s2.t)
        return s1;
    return s2;
}


Surface opSmoothUnion( in Surface s1, in Surface s2, float k ) {
    float h = clamp( 0.5 + 0.5*(s2.t-s1.t)/k, 0.0, 1.0 );
    return Surface(mix( s2.t, s1.t, h ) - k*h*(1.0-h), mix(s1.c, s2.c, 0.5))  - vec3(k*h*(1.0-h)); // TODO mix color only when collapsing
}

Surface scene(in vec3 p) {
    Surface final;

    vec3 t = (p+iTime)*6.;
    float d = (cos(t.x)*cos(t.y)*cos(t.z))/5.0;

    Surface cubeSurface = Surface(sdBox(p,vec3(1.0)), vec3(0.5,0.5,0.5));
    Surface torusSurface = Surface(sdTorus(p, vec2(1.0,0.6)), vec3(1.0, 0.0, 0.0));

    final = opSmoothUnion(cubeSurface,torusSurface, 0.2);

    return final;
}

Surface march(in Ray r) {
    float t = DIST_MIN;

    for(int i=0;i<RAY_MARCH_STEPS,t<=DIST_MAX;++i) {
        Surface s = scene(r.o+t*r.d);

        if(s.t<RAY_MARCH_PRECI) {
            return Surface(t+s.t,s.c);
        }

        t = t+s.t;
    }

    return Surface(DIST_MAX,vec3(0));
}

vec3 normalAt(in Surface s,in Ray r) {
    const float e = 0.01;
    vec3 p = r.o+s.t*r.d;
    float nx = scene(vec3(p.x+e,p.y,p.z)).t-scene(vec3(p.x-e,p.y,p.z)).t;
    float ny = scene(vec3(p.x,p.y+e,p.z)).t-scene(vec3(p.x,p.y-e,p.z)).t;
    float nz = scene(vec3(p.x,p.y,p.z+e)).t-scene(vec3(p.x,p.y,p.z-e)).t;

    return normalize(vec3(nx,ny,nz));
}

vec3 shade(in Surface s,in Ray r) {
    vec3 n = normalAt(s,r);
    vec3 l = normalize(vec3(1.,1.,-1.));
    
    float d = dot(n,l);
    return d*s.c;
}

Ray camRay(in vec2 p) {
    // p is the current pixel coord, in [-1,1]

    // normalized mouse position
    vec2 m = iMouse.xy/R.y;
    
    // camera position
    float DP = 10.;
    float d = DP/2.;
    vec3 ro = vec3(d*cos(6.0*m.x),DP/5.0,d*sin(6.0*m.x) );

    //vec3 ro = vec3(0.,0.,-7.);

    // target point
    vec3 ta = vec3(0.0,0.0,0.0);

    // camera view vector
    vec3 cw = normalize(ta-ro);

    // camera up vector
    vec3 cp = vec3(0.0,1.0,0.0);

    // camera right vector
    vec3 cu = normalize(cross(cw,cp));

    // camera (normalized) up vector
    vec3 cv = normalize(cross(cu,cw));
    
    float fovDeg = 45.;
    float fovRad = (fovDeg/360.)*2.*PI;
    float zf = 1./tan(fovRad/2.);
    
    // view vector, including perspective (the more you multiply cw, the less fovy)
    vec3 rd = normalize(p.x*cu + p.y*cv*(R.y/R.x) + zf*cw);

    return Ray(ro,rd);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (fragCoord / R.xy)*2.-1.;
    

    Ray r = camRay(uv);
    Surface s = march(r);
    vec3 c = vec3(0.5);
    
    if(s.t<DIST_MAX) {
        c = shade(s,r);
    }
    
    fragColor = vec4(c,1.0);
}