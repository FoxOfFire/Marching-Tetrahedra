#version 430

// Per fragment variables coming from the pipeline
in vec3 worldPosition;
in vec3 worldNormal;
in vec2 textureCoords;
in vec4 lightspace_pos;

// Outgoing values - fragment color
out vec4 outputColor;

// External parameters of the shader
uniform sampler2D textureImage;
uniform sampler2D textureShadow;

uniform vec3 cameraPosition;

uniform vec4 lightPosition;
uniform vec3 La = vec3(0.0, 0.0, 0.0 );
uniform vec3 Ld = vec3(1.0, 1.0, 1.0 );
uniform vec3 Ls = vec3(1.0, 1.0, 1.0 );

uniform vec4 lightPosition2;
uniform vec3 La2 = vec3(0.0, 0.0, 0.0 );
uniform vec3 Ld2 = vec3(1.0, 1.0, 1.0 );
uniform vec3 Ls2 = vec3(1.0, 1.0, 1.0 );


uniform float lightConstantAttenuation    = 1.0;
uniform float lightLinearAttenuation      = 0.0;
uniform float lightQuadraticAttenuation   = 0.0;

// anyag tulajdonsagok

uniform vec3 Ka = vec3( 1.0 );
uniform vec3 Kd = vec3( 1.0 );
uniform vec3 Ks = vec3( 1.0 );

uniform float Shininess = 1.0;
uniform int ignoreNormals = 0;

/* help:
	    - normalize: http://www.opengl.org/sdk/docs/manglsl/xhtml/normalize.xml
	    - dot: http://www.opengl.org/sdk/docs/manglsl/xhtml/dot.xml
	    - clamp: http://www.opengl.org/sdk/docs/manglsl/xhtml/clamp.xml
		- reflect: http://www.opengl.org/sdk/docs/manglsl/xhtml/reflect.xml
				reflect(incident_vector, normal_vector);
		- pow: http://www.opengl.org/sdk/docs/manglsl/xhtml/pow.xml
				pow(base, exponent);
*/

struct LightProperties
{
	vec4 pos;
	vec3 La;
	vec3 Ld;
	vec3 Ls;
	float constantAttenuation;
	float linearAttenuation;
	float quadraticAttenuation;
};

struct MaterialProperties
{
	vec3 Ka;
	vec3 Kd;
	vec3 Ks;
	float Shininess;
};

vec3 lighting(LightProperties light, vec3 position, vec3 normal, MaterialProperties material)
{
	
	vec3 ToLight; // A fényforrásBA mutató vektor
	float LightDistance = 0.0; // A fényforrástól vett távolság
	
	if ( light.pos.w == 0.0 ) // irány fényforrás (directional light)
	{
		// Irányfényforrás esetén minden pont ugyan abbóla az irányból van megvilágítva
		ToLight	= light.pos.xyz;
		// A távolságot 0-n hagyjuk, hogy az attenuáció ne változtassa a fényt
	}
	else				  // pont fényforrás (point light)
	{
		// Pontfényforrás esetén kiszámoljuk a fragment pontból a fényforrásba mutató vektort, ...
		ToLight	= light.pos.xyz - position;
		// ... és a távolságot a fényforrástól
		LightDistance = length(ToLight);
	}
	//  Normalizáljuk a fényforrásba mutató vektort
	ToLight = normalize(ToLight);

	if(ignoreNormals == 1){
		normal = ToLight;
	}
	
	// Attenuáció (fényelhalás) kiszámítása
	float Attenuation = 1.0 / ( light.constantAttenuation + light.linearAttenuation * LightDistance + light.quadraticAttenuation * LightDistance * LightDistance);
	
	// Ambiens komponens
	// Ambiens fény mindenhol ugyanakkora
	vec3 Ambient = light.La * material.Ka;

	// Diffúz komponens
	// A diffúz fényforrásból érkező fény mennyisége arányos a fényforrásba mutató vektor és a normálvektor skaláris szorzatával
	// és az attenuációval
	float DiffuseFactor = max(dot(ToLight,normal), 0.0) * Attenuation;
	vec3 Diffuse = DiffuseFactor * light.Ld * material.Kd;
	
	// Spekuláris komponens
	vec3 viewDir = normalize( cameraPosition - position ); // A fragmentből a kamerába mutató vektor
	vec3 reflectDir = reflect( -ToLight, normal ); // Tökéletes visszaverődés vektora
	
	// A spekuláris komponens a tökéletes visszaverődés iránya és a kamera irányától függ.
	// A koncentráltsága cos()^s alakban számoljuk, ahol s a fényességet meghatározó paraméter.
	// Szintén függ az attenuációtól.
	float SpecularFactor = pow(max( dot( viewDir, reflectDir) ,0.0), material.Shininess) * Attenuation;
	vec3 Specular = SpecularFactor * light.Ls * material.Ks;

	return Ambient + Diffuse + Specular;
}


void main()
{
	// A fragment normálvektora
	// MINDIG normalizáljuk!
	vec3 normal = normalize( worldNormal );

	LightProperties light;

	light.constantAttenuation = lightConstantAttenuation;
	light.linearAttenuation = lightLinearAttenuation;
	light.quadraticAttenuation = lightQuadraticAttenuation;

	MaterialProperties material;
	material.Ka = Ka;
	material.Kd = Kd;
	material.Ks = Ks;
	material.Shininess = Shininess;

	light.La = La;
	light.Ld = Ld;
	light.Ls = Ls;
	light.pos = lightPosition;
	vec3 shadedColor = lighting(light, worldPosition, normal, material);

	light.La = La2;
	light.Ld = Ld2;
	light.Ls = Ls2;
	light.pos = lightPosition2;
	vec3 shadedColor2 = lighting(light, worldPosition, normal, material);


	//outputColor = vec4(shadedColor, 1) * texture(textureImage, textureCoords);
	outputColor = (vec4(shadedColor, 1) +vec4(shadedColor2, 1)) * (vec4(textureCoords.xy,0,0));
}
