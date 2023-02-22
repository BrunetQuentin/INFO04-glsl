//make sdSphere
float sdSphere( vec3 p, float s )
{
    return length(p)-s;
}

// function that translate the pixel in 2D
vec2 translate(vec2 uv, float x, float y)
{
    return vec2(uv.x + x, uv.y + y);
}

// function that round the pixel to a certain level
float roundPixel( float d, float w )
{
    if(d < w)
    {
        return d;
    }
    else
    {
        return 1.0;
    }
}

// function to anti-alias the pixel
float antiAlias(float d, float w)
{
    float a = smoothstep(-w, w, d);
    return a;
}

// function that make atmospheric effect
vec4 makeAtmospheric(vec2 uv)
{
    float atmospheric = sdSphere(vec3(uv,0.0),0.3);

    float roundAtmospheric = roundPixel(atmospheric, 0.15);

    float roundAtmosphericAA = antiAlias(atmospheric, 0.02);

    return vec4(vec3(roundAtmosphericAA, 0.0, 0.0), 0.95);
}

vec4 makePlanet(vec2 uv)
{
    float planet = sdSphere(vec3(uv,0.0),0.2);

    float roundPlanet = roundPixel(planet, 0.1);

    float roundPlanetAA = antiAlias(planet, 0.002);

    return vec4(vec3(roundPlanetAA), 0.90);
}

// combine twopixel acording to the alpha value
vec4 combinePixel(vec4 pixel1, vec4 pixel2)
{
    vec4 finalColor = vec4(0.0);

    if(pixel1.w > pixel2.w)
    {
        finalColor = pixel1;
    }
    else
    {
        finalColor = pixel2;
    }

    return vec4(finalColor.xyz, 1.0);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = fragCoord.xy / iResolution.xy;
    vec2 uv2 = translate(uv, -0.5, -0.5);
    vec4 finalColor = vec4(0.0);

    vec4 sky = vec4(0.3, 0.0, 0.0, 1.0);
    
    vec4 planet = makePlanet(uv2);

    vec4 atmospheric = makeAtmospheric(uv2);

    vec4 planetPixels = combinePixel(planet, atmospheric);
    
    // finalColor = combinePixel(sky, planetPixels);

    fragColor = atmospheric;

}