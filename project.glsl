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


struct Ray {
    vec3 o; // origin
    vec3 d; // direction
};

struct Specular {
    float intensity;
    float size;
    float sharpness;
};

struct Surface {
    float t; // surface distance
    vec3 c; // surface color
    Specular spec; // specular
    float diff; // diffuse
};

float sdSphere( vec3 p, float s ) {
  return length(p)-s;
}

float sdPlane( vec3 p, vec3 n, float h ){
  // n must be normalized
  return dot(p,n) + h;
}

Surface add(in Surface s1, in Surface s2){
    if(s1.t<s2.t)
        return s1;
    return s2;
}


float hash13(vec3 p) {
    return fract(sin(dot(p,vec3(12.9898,78.233,45.5432)))*43758.5453123);
}

float vnoise(in vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);

    return mix(
                mix(mix(hash13(p+vec3(0.,0.,0.)),hash13(p+vec3(1.,0.,0.)),f.x),
                      mix(hash13(p+vec3(0.,1.,0.)),hash13(p+vec3(1.,1.,0.)),f.x),f.y),
                mix(mix(hash13(p+vec3(0.,0.,1.)),hash13(p+vec3(1.,0.,1.)),f.x),
                      mix(hash13(p+vec3(0.,1.,1.)),hash13(p+vec3(1.,1.,1.)),f.x),f.y),f.z);
}

float fnoise(in vec3 p,in float amplitude,in float frequency,in float persistence, in int nboctaves) {
    float a = amplitude;
    float f = frequency;
    float n = 0.0;

    for(int i=0;i<nboctaves;++i) {
        n = n+a*vnoise(p*f);
        f = f*2.;
        a = a*persistence;
    }
    
    return n;
}

float noiseTexture(in vec3 p) {

    vec3 t = (p+iTime)*4.;

    return fnoise(t,0.5,1.0,0.5,4);
}

Surface scene(in vec3 p) {
    Surface final;

    vec3 t = (p+iTime)*6.;
    float d = (cos(t.x)*cos(t.y)*cos(t.z))/5.0;

    float plane = sdPlane(p, vec3(0.0,1.0,0.0), 1.0);
    plane = plane + 1.0;
    Surface planeSurface = Surface(sdPlane(p, vec3(0.0,1.0,0.0), 1.0), vec3(1.0, 0.0, 0.0), Specular(0.0,0.0,0.0), 0.0);

    Surface rondSurface = Surface(sdSphere(p, 0.5), vec3(noiseTexture(p)), Specular(1.0,0.0,0.0), 1.0);

    final = add(rondSurface, planeSurface);

    return final;
}

Surface march(in Ray r) {
    float t = DIST_MIN;

    for(int i=0;i<RAY_MARCH_STEPS,t<=DIST_MAX;++i) {
        Surface s = scene(r.o+t*r.d);

        if(s.t<RAY_MARCH_PRECI) {
            s.t = s.t + t;
            return s;
        }

        t = t+s.t;
    }

    return Surface(DIST_MAX,vec3(0.0, 0.0, 0.0), Specular(0.0,0.0,0.0), 0.0);
}

vec3 normalAt(in Surface s,in Ray r) {
    const float e = 0.01;
    vec3 p = r.o+s.t*r.d;
    float nx = scene(vec3(p.x+e,p.y,p.z)).t-scene(vec3(p.x-e,p.y,p.z)).t;
    float ny = scene(vec3(p.x,p.y+e,p.z)).t-scene(vec3(p.x,p.y-e,p.z)).t;
    float nz = scene(vec3(p.x,p.y,p.z+e)).t-scene(vec3(p.x,p.y,p.z-e)).t;

    return normalize(vec3(nx,ny,nz));
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

vec3 shade(in Surface s,in Ray r) {
    vec3 n = normalAt(s,r);
    vec3 l = normalize(vec3(1.,1.,-1.));
    vec3 v = -r.d;
    vec3 e = reflect(-l,n);
    
    vec3 Ks = vec3(1.);
    vec3 Ka = vec3(0.);
    float sh = 50.;
    
    float diff = max(dot(n,l),0.) ;
    float spec = pow(max(dot(e,v),0.),sh) * s.spec.intensity;

    vec3 color = (Ka + diff);
    
    return color * s.c + Ks*spec;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (fragCoord / R.xy)*2.-1.;
    

    Ray r = camRay(uv);
    Surface s = march(r);
    s = march(r);
    vec3 c = vec3(0.5);
    
    if(s.t<DIST_MAX) {
        c = shade(s,r);
    }
    
    fragColor = vec4(c,1.0);
}