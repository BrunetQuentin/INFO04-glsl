
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

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/iResolution.xy * 20.0;

    // number for all perlinNoise parameter progressively with cosinus
    float persistence = cos(iTime * 0.1) * 0.5 + 0.5;
    float lacunarity = cos(iTime * 0.1) * 0.5 + 0.5;
    

    float noise = perlinNoise(uv, 1.0, 1.0, 10, lacunarity, persistence);

    // Output to screen
    fragColor = vec4(vec3(noise),1.0);
}