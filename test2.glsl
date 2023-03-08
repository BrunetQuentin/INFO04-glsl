const float PI = 3.14;
const int resolution = 16;
vec2 period = vec2( 10.0, 10.0 );
float power = 1.0;
float size = 20.0;

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

// make a random function that ruturn a number between 0 and 1 given a time value
float random21(float t) {
    return fract(sin(t) * 43758.5453123);
}

float turbulence( vec2 pos, in float size ) {

    float value = 0.0, initialSize = size;
    float x = pos.x;
    float y = pos.y;
    
    for ( int i = 0; i<resolution; i++ ) {
    	value += noiseFunction(vec2(x / size, y / size)) * size;
        size /= 2.0;
    }
    
    return( 128.0 * value / initialSize );
}

float marble( in vec2 p ) {

  	float x = p.x;
    float y = p.y;

    float xy = x / iResolution.y * period.x;

    xy += y * period.y / iResolution.x;
    xy += power * turbulence( p, size ) / 256.0;

    return sin( 256.0 * xy * PI );
    
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/iResolution.xy * 20.0;

    // number for all perlinNoise parameter progressively with cosinus
    float persistence = cos(iTime * 0.1) * 0.5;
    float lacunarity = cos(iTime * 0.1) * 0.5;
    

    float noiseR = perlinNoise(uv, 2.0, 10.0, 10, lacunarity, persistence);
    float noiseG = perlinNoise(uv, 2.0, 10.0, 10, lacunarity / 2.0, persistence / 2.0);
    float noiseB = perlinNoise(uv, 2.0, 10.0, 10, lacunarity / 4.0, persistence / 4.0);

    // Output to screen
    fragColor = vec4(vec3(marble(uv), marble(uv), marble(uv)),1.0);

}