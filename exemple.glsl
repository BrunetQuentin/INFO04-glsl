#define PI 3.14159265
#define TAU (2*PI)
#define PHI (1.618033988749895)

float smin( float a, float b, float k )
{
    float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k ;
}


float hash(float x) {
 	return fract(sin(x * 4341.1) * 3123.8);   
}

float hash11(float p)
{
	p = fract(p * .1031);
	p *= p + 19.19;
	p *= p + p;
	return fract(p);
}

// method by fizzer
vec3 hashHs(vec3 n, float seed)
{
	float u = hash11(78.233+seed);
	float v = hash11(10.873+seed);
	float a = 6.2831853 * v;
	u = 2.0*u - 1.0;
    return normalize(n+vec3(sqrt(1.0-u*u) * vec2(cos(a), sin(a)), u));    
}

void pR(inout vec2 p, float a) {
	p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

float sdBox( vec3 p, vec3 b )
{
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float sphere(vec3 p, float r) {
 	return length(p) - r;   
}

float crate(vec3 p, float d) {
    float sdf = max(-sdBox(p, vec3(11., d, d)), sdBox(p, vec3(10.)));
    sdf = max(-sdBox(p.zyx, vec3(11., d, d)), sdf);

    return max(-sdBox(p.yxz, vec3(11., d, d)), sdf);;   
}


float c = 0.;
float vols[4] = float[](0., 0.,0., 0.);

float fSphere(vec3 p, float r) {
	return length(p) - r;
}

float blob(vec3 bp, float r) {
    vec3 p = bp;
    
    bp.xzy *= 0.4;
    bp.x -= sin(vols[2] * 10.) * 5.;
    pR(bp.xz, iTime + vols[3] * 2.);
    pR(bp.yx, iTime * .8 + vols[3] * 12.0);
    
    bp.xyz += sin(1.5*bp.yzx + vec2(0., iTime + vols[0]).xxy * 10.);
    bp.xyz -= sin(1.5*bp.yzx + vec2(0., iTime + vols[0]).yxx * 10.) * 1.5;
    
    
    bp.xzy *= .3;
    return fSphere(bp, r);
    
}

float map(vec3 p) {
    vec3 op = p, op2=p;
    
    
    
    float sdf = max(0., -sdBox(p, vec3(45.)));

    float a = sdBox(op2 + vec3(0., 31., 0.), vec3(20., 20., 20.));

    sdf = mix(fSphere(op, 7.), blob(op, 1.4), vols[2]);

    sdf = min(sdf, a);
    
    pR(op.xz, -iTime * .25);
    
    sdf = smin(sdf, crate(op * (1. + vols[0] * .1), 8.6 + vols[0]), 1.9);
    sdf = min(sdf, -sdBox(p, vec3(45.)));
    
    return sdf;
}

float trace(vec3 ro, vec3 rd, int steps) {
    float d = 1e32;
    float s = 0.;
    
    for (int i = 0; i < steps; i++) {
    	d = map(ro + rd * s);
        
        if (d < 0.001 || d > 30.) break;
        s += d;
    }
    
    return s;    
}

float calcSoftshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax)
{
	float res = 1.0;
    float t = mint;
    float ph = 1e10; // big, such that y = 0 on the first iteration
    
    for( int i=0; i<32; i++ )
    {
		float h = map( ro + rd*t );
        float y = h*h/(2.0*ph);
        float d = sqrt(h*h-y*y);
        res = min( res, 10.0*d/max(0.0,t-y) );
        ph = h;        
        t += h;        
        if( res<0.0001 || t>tmax ) break;
        
    }
    return clamp( res, 0.8, 1.0 );
}

vec3 randomSphereDir(vec2 rnd)
{
	float s = rnd.x*PI*2.;
	float t = rnd.y*2.-1.;
	return vec3(sin(s), cos(s), t) / sqrt(1.0 + t * t);
}

vec3 randomHemisphereDir(vec3 dir, float i)
{
	vec3 v = randomSphereDir( vec2(hash(i+1.), hash(i+2.)) );
	return v * sign(dot(v, dir));
}

float ambientOcclusion( in vec3 p, in vec3 n, in float maxDist, in float falloff )
{
	const int nbIte = 32;
    const float nbIteInv = 1./float(nbIte);
    const float rad = 1.-1.*nbIteInv; //Hemispherical factor (self occlusion correction)
    
	float ao = 0.;
    
    for( int i=0; i<nbIte; i++ )
    {
        float l = hash(float(i)) * maxDist;
        vec3 rd = normalize(n+randomHemisphereDir(n, l )*rad)*l;
        
        ao += (l - max(map( p + rd ),0.)) / maxDist * falloff;
    }
	
    return clamp( 1.-ao*nbIteInv, 0., 1.);
}

vec3 normals(vec3 p) {
    vec3 d = vec3(0.001, 0, 0);
 	return normalize(
        vec3(
        	map(p + d.xyy) - map(p - d.xyy),
            map(p + d.yxy) - map(p - d.yxy),
            map(p + d.yyx) - map(p - d.yyx)
        )
    );
}

void image(out vec4 color, in vec2 uv)
{
	
    vols = float[](
        texture(iChannel1, vec2(0.1, .2)).r,
        texture(iChannel1, vec2(0.6, 0.2)).r,
        mix(0., texture(iChannel1, vec2(0.2, 1.)).r * 1., min(67., iTime) / 67.), // mix sphere / blob
        texture(iChannel1, vec2(0.6, 0.2)).r
    );
    
    
    vec3 
        vuv = normalize(vec3(0., 1., 0.)),
    	eyePosition = vec3(sin(iMouse.x / 300.) * 40., 1.,  cos(iMouse.x / 300.) * 40.);

    pR(eyePosition.xz, +iTime* .1);
    pR(eyePosition.yx, - blob(eyePosition + vols[2], 1.) * .01);
    vec3
        vrp = vec3(0.),
    	vpn = normalize(vrp - eyePosition),
    	u = normalize(cross(vuv, vpn)),
    	v = cross(vpn, u),
    	vcv = (eyePosition + vpn),
    	scrCoord = (vcv + uv.x * u + uv.y * v),
    	viewDirection = normalize(scrCoord - eyePosition);   
    
    vec3 ro = eyePosition;
    vec3 rd = viewDirection;
    vec4 c, oc = vec4(0.);
    
    vec3 LIGHT = vec3(sin(-iTime) * 15., -10., cos(-iTime) * 15.);
    vec3 nLight = normalize(LIGHT);

    int initialTracingSteps = 99;
    int tracingSteps = initialTracingSteps;
    int bounces = 4;
    
    for (int bounce = 1; bounce < bounces; bounce++) {
        
        float march = trace(ro, rd, tracingSteps);
		tracingSteps -= initialTracingSteps / bounces;
        
        ro += rd * march;
        vec3 n = normals(ro);
		
        c.rgb = vec3(1.);    

        if (ro.y > 10.) {
            vec3 ro2 = ro;
            ro2.x += sin(iTime + (vols[1] * ro2.x * .1 + ro2.z * .2) + vols[2] ) * 5.;
            float stripes = (n.y < .0 && ro.y > 0.) ? (1. * ceil(fract((ro2.z-ro2.x) * .1 + iTime + vols[2]) - .5) * 3.) : 0.;
            c += stripes * vols[2];            
        }
        
                
        vec3 op = ro;
        
        float bl = mix(fSphere(op, 7.), blob(op, 1.4), vols[2]); // to blob
        
        pR(op.xz, iTime * .25);
        
        float cr = crate(op * (1. + vols[0] * .1), 8.6 + vols[0]); // to box
        
        c.b += mix(1., 0., clamp(bl + cr, 0., 1.)) * 5.;
        
        
        c += pow(max(0., dot(n, nLight)), 5.) * 2.5;
        c = mix(vec4(.5), c, (1.-(march/ 125.)) * .5);

        //c += (1. / distance(LIGHT, ro)) * 1.;
        
        c *= ambientOcclusion(ro, n, 4.2, 1.2);
		
        c.b += mix(.1, 0., clamp((bl + cr) * .5, 0., 1.)) * 1.;
                
        oc = mix(oc, c, .2 / float(bounce));
        
        rd = reflect(rd, n);
        //rd = refract(rd, normalize(ro), .99);
        //ro += mix(ro, hashHs(rd / 1., float(99 + bounce)), .1) * .0001;
        
        rd.x += hash(ro.x)* .01;
        rd.y += hash(ro.y)* .01;
        rd.z += hash(ro.z)* .05;
        
        rd = normalize(rd);
        ro += rd * .01;
    }
    
    uv += 3.;
    
    //oc *= 1.+ max(0., pow(abs(dot(viewDirection, nLight)) * 1.001, 124.));
    //oc = pow(abs(oc), 1. / vec4(1.75)) * .45;
    
    color = 1. * max(
        vec4(0.), 
        vec4(
            vec3((oc * 3. - length(uv/ 12.)) + hash(rd.x / uv.y + mod(iTime, 1.)) * .001), 1.));
    
    //color = oc;
    color = pow(color * 1., vec4(1.7)) * 3.;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (2.*fragCoord.xy-iResolution.xy) / iResolution.x;
    
    image(fragColor, uv);
    
}