#version 330 core

layout (location = 0) in vec3 aPos; // vertex position
layout (location = 1) in vec2 aTexCoords; // vertex texture coordinates
layout (location = 2) in vec3 aNormal; // vertex normal

out vec2 texCoord;
out vec3 Normal;

uniform mat4 model; // Model matrix
uniform mat4 view; // View matrix
uniform mat4 projection; // Projection matrix

void main()
{
    gl_Position = projection * view * model * vec4(aPos, 1.0);

    // Pass texture coordinates to the fragment shader
    texCoord = aTexCoords;

    // Transform vertex normal to model space coordinate system
    Normal =  transpose(inverse(mat3(model))) * aNormal;
}