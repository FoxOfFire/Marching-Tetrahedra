#pragma once

#include "Experimental.hpp"
// GLM
#include <array>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>

// GLEW
#include <GL/glew.h>

// SDL
#include <SDL3/SDL.h>
#include <SDL3/SDL_opengl.h>

// Utils
#include "Camera.h"
#include "CameraManipulator.h"
#include "GLUtils.hpp"

struct SUpdateInfo
{
	float ElapsedTimeInSec = 0.0f;	// Elapsed time since start of the program
	float DeltaTimeInSec = 0.0f;	// Elapsed time since last update
};

class CMyApp
{
public:
	CMyApp();
	~CMyApp();

	bool Init();
	void Clean();

	void Update(const SUpdateInfo&);
	void Render();
	void RenderGUI();

	void KeyboardDown(const SDL_KeyboardEvent&);
	void KeyboardUp(const SDL_KeyboardEvent&);
	void MouseMove(const SDL_MouseMotionEvent&);
	void MouseDown(const SDL_MouseButtonEvent&);
	void MouseUp(const SDL_MouseButtonEvent&);
	void MouseWheel(const SDL_MouseWheelEvent&);
	void Resize(int, int);

	MeshObject<Vertex> MakeVoxel(const std::array<float, 8>&,const glm::vec3&, bool,bool,float);
	MeshObject<Vertex> MakeTetrahedra(const std::array<glm::vec4,4>&);
	MeshObject<Vertex> MakeTetrahedraFrame(const std::array<glm::vec4,4>&);
	MeshObject<Vertex> MergeGeometry(const MeshObject<Vertex>&, const MeshObject<Vertex>&);
	MeshObject<Vertex> MakeVoxelChunk(const glm::ivec3&,bool,bool);
	void MakeVoxelGrid(const glm::vec3&, bool,const std::initializer_list<VertexAttributeDescriptor>);
	
	void PrintVector(const glm::vec4&);

	void GeneratePerlinNoise();
	void ResetChunkNoise();
	void UpdateChunkNoise(const glm::ivec3&);
	signed char SamplePerlinNoise(const glm::vec3&);

	float GetPixelDepthAt(const glm::vec2&);

	typedef std::array<std::array<std::array<float, 4>, 4>, 4> ChunkData;
	ChunkData ExtractChunkData(const glm::ivec3&,const glm::vec3&);
	
	void FixMeshNormals(MeshObject<Vertex>&);
	

	void OtherEvent(const SDL_Event&);
protected:
	void SetupDebugCallback();
	void RenderGeometry(const glm::mat4&, GLuint);

	//
	// Variables
	//
	float m_ElapsedTimeInSec = 0.0f;
	int m_width = 640, m_height = 480;


	// Camera
	Camera m_camera;
	CameraManipulator m_cameraManipulator;

	//
	// OpenGL
	//

	void DrawAxes();

	// Shader variables
	GLuint m_programID = 0;				// Shader of the objects
	GLuint m_programAxesID = 0;			// Program showing X,Y,Z directions
	GLuint m_programPostprocessID = 0;	// Postprocess program

	// Shader initialization and termination
	void InitShaders();
	void CleanShaders();
	void InitAxesShader();
	void CleanAxesShader();

	// Geometry variables
	OGLObject m_voxelFrame = {};
	glm::mat4 testChunk = glm::mat4(0);
	std::array<std::array<std::array<glm::vec3,4> ,4>,4> m_testPerlinMap;

	

	MeshObject<Vertex> m_frame = {};
	OGLObject m_testVoxel = {};
	std::array<float,8> m_testFrame= {1,1,1,1,-1,-1,-1,-1};

	constexpr static float m_chunkPerlinRatio =  1.5f;
	const static int m_chunkCount = 7;
	const static int m_perlinSize = m_chunkCount/m_chunkPerlinRatio;
	glm::vec3 m_start = -glm::vec3((float)m_chunkCount*3.f/2.f);

	std::array<std::array<std::array<glm::vec3,m_perlinSize> ,m_perlinSize>,m_perlinSize> m_perlinMap;
	std::array<std::array<std::array<glm::mat4,m_chunkCount> ,m_chunkCount>,m_chunkCount> m_chunkMap;

	std::array<std::array<std::array<OGLObject,m_chunkCount> ,m_chunkCount>,m_chunkCount> m_chunkMeshes;


	float m_epsilon = .5;

	bool m_makeVoxelGrid = true;
	bool m_smoothShade = false;
	bool m_ortho = false;
	bool m_5tetra = true;
	bool m_chunkBorders = false;
	// Geometry initialization and termination
	void InitGeometry();
	void CleanGeometry();

	// Texture variables
	GLuint m_SamplerID = 0;
	GLuint m_metalTextureID = 0;
	
	// Texture initialization and termination
	void InitTextures();
	void CleanTextures();

	// Light source
	glm::vec4 m_lightPosition = (glm::vec4(-m_chunkCount*3, m_chunkCount*3, -m_chunkCount*3,1.0f));
	glm::vec3 m_La = glm::vec3(.05);
	glm::vec3 m_Ld = glm::vec3(.5);
	glm::vec3 m_Ls = glm::vec3(.5);

	glm::vec4 m_lightPosition2 = (glm::vec4(m_chunkCount*3, m_chunkCount*3, m_chunkCount*3,1.0f));
	glm::vec3 m_La2 = glm::vec3(.05);
	glm::vec3 m_Ld2 = glm::vec3(.5);
	glm::vec3 m_Ls2 = glm::vec3(.5);


};
