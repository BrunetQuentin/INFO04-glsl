#define PI 3.14159265359

// fonctions de distance signées
float sdCircle(vec2 p,float r) {
    return length(p)-r;
}

float sdBox( in vec2 p, in vec2 b ) {
    vec2 d = abs(p)-b;
    return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}


float sdSegment( in vec2 p, in vec2 a, in vec2 b ) {
    vec2 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}

// operateurs sur les fonstions de distance signées
float substract(float sd1,float sd2) {
    return max(sd1,-sd2);
}

float add(float sd1,float sd2) {
    return min(sd1,sd2);
}

float intersect(float sd1,float sd2) {
    return max(sd1,sd2);
}

// transformations (du domaine)
vec2 translate(in vec2 p,in vec2 t) {
    return p-t;
}

vec2 scale(in vec2 p, vec2 f) {
    return p*f;
}

vec2 rotate(in vec2 p, float a) {
    float ca = cos(a);
    float sa = sin(a);
    return vec2(ca*p.x-sa*p.y, sa*p.x+ca*p.y);
}

// des objets et/ou scenes
float boite(in vec2 uv) {
    float L = 1.5;
    float c1 = sdBox(uv,vec2(L,0.5));
    vec2 uvCircle = translate(uv,vec2(cos(iTime)*(L-0.5),0.0));
    float c2 = sdCircle(uvCircle,0.4);
    float c = substract(c1,c2);
    return c;
}

float roue(in vec2 p) {
    float c = sdCircle(p,0.7);
    c = substract(c,sdCircle(p,0.68));
    
    vec2 a = rotate(vec2(-0.7,0.0),iTime*3.0);
    vec2 b = rotate(vec2( 0.7,0.0),iTime*3.0);
    int nbSeg = 20;
    float angle = ((180.0/float(nbSeg))/360.0)*2.0*PI;
    for(int i=0;i<nbSeg;++i) {
        vec2 ar = rotate(a,float(i)*angle);
        vec2 br = rotate(b,float(i)*angle);
        c = add(c,sdSegment(p,ar,br));
    }
    return c;
}

float other(in vec2 uv) {
    uv *= 2.0;
    vec2 a = vec2(-1.0,0.0);
    vec2 b = vec2( 1.0,0.0);
    float s = sdSegment(uv,a,b);
    
    vec2 tmp1 = translate(uv,vec2(cos(iTime),0.0));
    vec2 tmp2 = translate(uv,vec2(cos(iTime)*-1.0,0.0));
    
    s = add(s,roue(tmp1));
    s = add(s,roue(tmp2));
    
    return s;
}

#define ROUES
//#define BORDURE
//#define BOITE

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    //vec2 uv = fragCoord/iResolution.y;

    vec2 uv = (2.*fragCoord - iResolution.xy) / iResolution.y;  
    
    uv = rotate(uv,iTime*0.5);

#ifdef ROUES
    float c = other(uv);
    c = smoothstep(0.0,0.02,c);
#endif

#ifdef BORDURE
    // affchage de la bordure d'un cercle
    float c1 = sdCircle(uv,0.5);
    float c2 = sdCircle(uv,0.4);
    float c = substract(c1,c2);
    
    c = smoothstep(0.0,0.02,c);
#endif

#ifdef BOITE
    // une boite trouée qui tourne et se repete...
    uv = mod(uv*1.0,vec2(1.0))*2.0-1.0;
    uv = rotate(uv,iTime*0.5);
    uv = scale(uv,vec2(2.0,2.0));
    
    float c = boite(uv);
    c = smoothstep(0.0,0.02,c);
#endif
    
    
    vec3 col = vec3(c);
       
    // Output to screen
    fragColor = vec4(col,1.0);
}