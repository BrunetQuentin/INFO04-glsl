#define R vec2(4,4)

float sdCircle( vec2 p, float r )
{
    return length(p) - r;
}

float sdPentagon(in vec2 p, in float r){
    const vec3 k = vec3(0.809016994,0.587785252,0.726542528);
    p.x = abs(p.x);
    p -= 2.0*min(dot(p,vec2(-k.x,k.y)),0.0)*vec2(-k.x,k.y);
    p -= 2.0*min(dot(p,vec2(k.y,k.z)),0.0)*vec2(k.y,k.z);
    p -= vec2(clamp(p.x, -r*k.z, r*k.z),r);
    return length(p)*sign(p.y);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord){
    vec3 cFond = vec3(1.0);
    vec3 cCercle = vec3(0.0,0.0,1.0);
    vec2 uv = (fragCoord*2.0 -R.xy)/R.y;

    float d = sdCircle(fragCoord, 0.2);

    float test = smoothstep(-0.01,0.001,d);
    vec3 col = mix(cCercle, cFond, test);

    fragColor = vec4(vec3(d), 1.0);
}