#version 430

// Variables coming from the VBO
layout( location = 0 ) in vec3 inputObjectSpacePosition;
layout( location = 1 ) in vec3 inputObjectSpaceNormal;
layout( location = 2 ) in vec2 inputTextureCoords;

// Values to be passed on the pipeline
out vec3 worldPosition;
out vec3 worldNormal;
out vec2 textureCoords;
out vec4 lightspace_pos;

// External parameters of the shader
uniform mat4 world;
uniform mat4 worldInvTransp;
uniform mat4 viewProj;
uniform mat4 shadowVP;

void main()
{
	gl_Position = viewProj * world * vec4( inputObjectSpacePosition, 1 );
	worldPosition  = (world   * vec4(inputObjectSpacePosition,  1)).xyz;
	worldNormal = (worldInvTransp * vec4(inputObjectSpaceNormal, 0)).xyz;
	textureCoords = inputTextureCoords;
	// Where is the vertex from the light's point of view
	lightspace_pos = shadowVP * world * vec4( inputObjectSpacePosition, 1 );
}