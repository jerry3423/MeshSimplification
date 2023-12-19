#version 330 core

out vec4 FragColor;

in vec2 texCoord;
in vec3 Normal;

uniform sampler2D modelTexture;
uniform bool enableTexture;

void main()
{
    // Set light parameters
    vec3 lightColor = vec3(1, 1, 1);
    vec3 lightDir = vec3(0, 0, 1);

    // Calculate diffuse color
    float diff = max(dot(Normal, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    vec3 result = diffuse;
    // If use texture
	if (enableTexture) {
        result = texture(modelTexture, texCoord).rgb;
    } 
    FragColor = vec4(result, 1.0);
}