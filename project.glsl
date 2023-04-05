// screen resolution
#define R iResolution.xy

// minimum distance to objects
#define DIST_MIN 1.

// maximum distance to objects
#define DIST_MAX 2000.0

// max number of steps for the ray-marching
#define RAY_MARCH_STEPS 100

// consider hit if we reach this distance
#define RAY_MARCH_PRECI 0.001

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

// function random
vec2 random (vec2 st) {
    st = vec2(dot(st,vec2(127.1,311.7)),
              dot(st,vec2(269.5,183.3)));
    return -1.0 + 2.0*fract(sin(st)*43758.5453123);
}

// function of a noise function
float noiseFunction(vec2 p) {
    vec2 i = floor(p);
    vec2 f = fract(p);
    vec2 u = f*f*(3.0-2.0*f);
    return mix(mix(dot(random(i + vec2(0.0, 0.0)), f - vec2(0.0, 0.0)),
                   dot(random(i + vec2(1.0, 0.0)), f - vec2(1.0, 0.0)), u.x),
               mix(dot(random(i + vec2(0.0, 1.0)), f - vec2(0.0, 1.0)),
                   dot(random(i + vec2(1.0, 1.0)), f - vec2(1.0, 1.0)), u.x), u.y);
}

// function that generate perlin noise  in 2D
float perlinNoise(vec2 p, float frequency, float amplitude, int octaves, float lacunarity, float persistence) {
    float total = 0.0;
    for (int i = 0; i < octaves; i++) {
        total += noiseFunction(p * frequency) * amplitude;
        frequency *= lacunarity;
        amplitude *= persistence;
    }
    return total;
}

float rand(vec2 n) { 
	return fract(sin(dot(n, vec2(12.9898, 4.1414))) * 43758.5453);
}
float noise(vec2 p){
	vec2 ip = floor(p);
	vec2 u = fract(p);
	u = u*u*(3.0-2.0*u);
	
	float res = mix(
		mix(rand(ip),rand(ip+vec2(1.0,0.0)),u.x),
		mix(rand(ip+vec2(0.0,1.0)),rand(ip+vec2(1.0,1.0)),u.x),u.y);
	return res*res;
}
float fbm2(vec2 p) {
    float r = 0.0;
    float amp = 1.0;
    float freq = 1.0;
    for(int i = 0; i < 5; i++) {
        r += amp * noise(freq*p);
        amp *= 0.5; 
        freq *= 2.0;
    }
    return r;
}

//https://www.shadertoy.com/view/MsBBWK
vec3 marble(float u, float v) {
    vec2 p = vec2(u, v);
    vec2 q = vec2(0);
    vec2 r = vec2(0);
    q.x = fbm2(p + vec2(1, 1));
    q.y = fbm2(p + vec2(2, 2));
    r.x = fbm2(p + 2.0*q + vec2(1, 1));
    r.y = fbm2(p + 2.0*q + vec2(2, 2));
    return fbm2(p + 4.0*r) * vec3(0, q.x + r.x, q.y + r.y);
}

Surface sdMountain(in vec3 p, in bool calculateTexture) { 

    float plane = sdSphere(p+vec3(0.0,2005.0,0.0), 2000.0);

    if(plane < DIST_MAX){
        float mountain = perlinNoise(vec2(p.x, p.z - (iTime * 0.5)), 0.2, 4.0, 5, 2.0, 0.5);
        // multiply by a factor to make the mountain at the same level as the plane when approaching on the x axis
        mountain = mountain * smoothstep(0.5, 3.0, abs(p.x)) * - 1.0;
        plane = plane + mountain;
    }
    
    return Surface(plane, vec3(0.6078, 0.3294, 0.1059), Specular(0.3,0.01,1.0), 0.2);
}

Surface sdPlanet(in vec3 p, in bool calculateTexture) {
vec3 marbleNoise;
    if(calculateTexture){
        marbleNoise = marble((p.x+(iTime))/50.0, p.y/50.0);
    }else {
        marbleNoise = vec3(0.0);
    }

    Surface planet = Surface(sdSphere(p+vec3(0.0, 30.0, 1000.0), 100.0), marbleNoise, Specular(marbleNoise.b,0.1,0.5), 5.0); // normalise marblenoise.b
    return planet;
}

Surface scene(in vec3 p, in bool calculateTexture) {
    Surface final;

    Surface planeSurface = sdMountain(p, calculateTexture);

    Surface planet = sdPlanet(p, calculateTexture);
 
    final = add(planet, planeSurface);
    // final = planeSurface;
    // final = planet;

    return final;
}

Surface march(in Ray r) {
    float t = DIST_MIN;

    for(int i=0;i<RAY_MARCH_STEPS,t<=DIST_MAX;++i) {
        Surface s = scene(r.o+t*r.d, true);

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
    float nx = scene(vec3(p.x+e,p.y,p.z), false).t-scene(vec3(p.x-e,p.y,p.z), false).t;
    float ny = scene(vec3(p.x,p.y+e,p.z), false).t-scene(vec3(p.x,p.y-e,p.z), false).t;
    float nz = scene(vec3(p.x,p.y,p.z+e), false).t-scene(vec3(p.x,p.y,p.z-e), false).t;

    return normalize(vec3(nx,ny,nz));
}

Ray camRay(in vec2 p) {
    // p is the current pixel coord, in [-1,1]

    // normalized mouse position
    vec2 m = iMouse.xy/R.y;
    
    // camera position
    float DP = 10.;
    float d = DP/2.;
    vec3 ro = vec3(d*cos(6.0*m.x),DP/5.0 - 1.6,d*sin(6.0*m.x) );

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
    vec3 l = normalize(vec3(1.,1.,-1));
    vec3 v = -r.d;
    vec3 e = reflect(-l,n);
    
    vec3 Ka = vec3(0.0);
    
    float diff = min(max(dot(n,l + s.diff),0.), 1.0);
    float spec = pow(max(dot(e,v),0.),s.spec.size * 50.) * s.spec.intensity;

    vec3 color = (Ka + diff);
    
    return color * s.c + s.spec.sharpness * spec;
}

//https://www.shadertoy.com/view/sddSzX
float n11(float p) {
	return fract(sin(p*154.101)*313.019);
}

float n21(vec2 p) {
	return sin(dot(p, vec2(7., 157.)));
}


float star(vec3 p) {
	float z = 1.;
	vec2 gv = fract(p.xy*z) - 0.5;
	vec2 id = floor(p.xy*z);
	gv.x += sin(n21(id)*354.23) * 0.3;
	gv.y += sin(n11(n21(id))*914.19) * 0.3;
	float r = n11(n21(id));
	return 0.1*n11(r)*abs(sin(p.z+r*133.12))*0.4/length(gv)*0.1;
}

float stars(in vec3 p) {
	float z = 1., m = 0.;
	for(int i=1; i<=6;i++){
		vec3 t = vec3(0., 0., p.z + iTime * 0.3);
		z *= 2.;
		m += star(vec3(p.xy*z, 1.)+t);
	}
	return m;
}
//https://www.shadertoy.com/view/sddSzX

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (fragCoord / R.xy)*2.-1.;
    

    Ray r = camRay(uv);
    Surface s = march(r);

    vec3 c;
    
    if(s.t<DIST_MAX) {
        c = shade(s,r);
    }else {
        vec3 sphere =  r.d;
        c = vec3(stars(sphere));
    }
    
    fragColor = vec4(c,1.0);
}