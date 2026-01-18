#include  "MyApp.h"
#include "GLUtils.hpp"
#include "SDL_GLDebugMessageCallback.h"
#include "ObjParser.h"
#include "ProgramBuilder.h"

#include <array>
#include <glm/common.hpp>
#include <glm/ext/quaternion_geometric.hpp>
#include <glm/ext/vector_float4.hpp>
#include <glm/fwd.hpp>
#include <glm/geometric.hpp>
#include <imgui.h>
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>

CMyApp::CMyApp()
{
}

CMyApp::~CMyApp()
{
}

void CMyApp::PrintVector(const glm::vec4& vec){
	std::cout << "x:" << vec.r<<"\ty:"<<vec.g<<"\tz:"<<vec.b<<"\tw:"<<vec.a<<std::endl;
}

void CMyApp::GeneratePerlinNoise(){
	std::cout<< "Generating perlin noise"<<std::endl;

	for(int i = 0; i < m_perlinSize; i++){
		if(!m_makeVoxelGrid && i>3)
			continue;
		for(int j = 0; j < m_perlinSize; j++){
			if(!m_makeVoxelGrid && j>3)
				continue;
			for(int k = 0; k < m_perlinSize; k++){
				if(!m_makeVoxelGrid){
					if(k>3)
						continue;
					m_perlinMap[i][j][k] =glm::vec3(0,1,0)* m_testFrame[(i/2)*4+(j/2)*2+(k/2)];
					continue;
				}
				if( i==0 || j==0 || k==0 ||
					i==m_perlinSize-1 ||
					j==m_perlinSize-1 ||
					k==m_perlinSize-1){
					m_perlinMap[i][j][k] = glm::vec3(0,-1,0);
					continue;
				}
				glm::vec3 rdir = glm::vec3(1);
				while(glm::length(rdir) > 1){
					float x = (float)(rand()%2001-1000)/1000;
					float y = (float)(rand()%2001-1000)/1000;
					float z = (float)(rand()%2001-1000)/1000;
					rdir = glm::vec3(x,y,z);
				}
				rdir = normalize(rdir);
				m_perlinMap[i][j][k] = rdir;
			}
		}
	}
	ResetChunkNoise();
}

void CMyApp::ResetChunkNoise(){
	std::cout<< "Reseting chunk noise values"<<std::endl;

	if(m_makeVoxelGrid){
		for(int i = 0; i < m_chunkCount; i++)
			for(int j = 0; j < m_chunkCount; j++)
				for(int k = 0; k < m_chunkCount; k++){
					UpdateChunkNoise(glm::ivec3(i,j,k));
				}
	}
	else{
		UpdateChunkNoise(glm::ivec3(0,0,0));
	}
		
}
void CMyApp::UpdateChunkNoise(const glm::ivec3& coords){
	std::cout<< "Updating chunk noise values"<<std::endl;

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++){
			signed char chs[] = {'\0','\0','\0','\0'};
			for(int k = 0; k < 4; k++){
				chs[k] = SamplePerlinNoise(glm::vec3(coords)*3.f+ glm::vec3(i,j,k) );
			}
			m_chunkMap[coords.x][coords.y][coords.z][i][j] = *((glm::float32*)chs);
			
		}
}

signed char CMyApp::SamplePerlinNoise(const glm::vec3& vec){
	float scalar;
	glm::vec3 scaledVec;
	if(m_makeVoxelGrid){
		scalar = (float)(m_perlinSize-1)/(float)(m_chunkCount*3+1);
		scaledVec = vec*scalar;
	}
	else{
		scalar = (4.f-m_epsilon)/4.f;
		scaledVec = vec*scalar + glm::vec3(1)*(m_epsilon + 0.01f);
	}

	int x = std::floor(scaledVec.x);
	int y = std::floor(scaledVec.y);
	int z = std::floor(scaledVec.z);
	float xScalar = scaledVec.x - (float)x;
	float yScalar = scaledVec.y - (float)y;
	float zScalar = scaledVec.z - (float)z;

	glm::vec3 a,b,c,d;
	a = mix(m_perlinMap[x][y][z],m_perlinMap[x+1][y][z],xScalar);
	b = mix(m_perlinMap[x][y][z+1],m_perlinMap[x+1][y][z+1],xScalar);
	c = mix(m_perlinMap[x][y+1][z],m_perlinMap[x+1][y+1][z],xScalar);
	d = mix(m_perlinMap[x][y+1][z+1],m_perlinMap[x+1][y+1][z+1],xScalar);

	a = mix(a,c,yScalar);
	b = mix(b,d,yScalar);
	a = mix(a,b,zScalar);

	float res = glm::dot(a,glm::vec3(0,1,0));
	
	//[-1,1]->[0,1]
	res = (res+1)/2;
	//bonus
	res = 6*powf(res,5.f)-15*powf(res,4.f)+10*powf(res,3.f);
	//[0,1]->[-1,1]
	res = res*2-1;
	res *= 128;
	res = glm::clamp(res,-128.f,127.f);
	signed char ch = (signed char)((int)res);

	return ch;
}

CMyApp::ChunkData CMyApp::ExtractChunkData(const glm::ivec3& coords,const glm::vec3& start){
	CMyApp::ChunkData data = {};
	glm::mat4 chunk = m_chunkMap[coords.x][coords.y][coords.z];
	assert(sizeof(chunk[0][0]) == sizeof(glm::float32));


	for(int i = 0 ; i < 4; i++){
		for(int j = 0; j < 4; j++){
			glm::float32 val = chunk[i][j];
			signed char* vals = (signed char*)(&val);
			for(int k = 0; k < 4; k++){
				signed char ch = vals[k];
				glm::vec3 pos = glm::vec3(i,j,k)+start;
				data[i][j][k] = (float)((int)ch)/128.f;
			}
		}
	}
	return data;
}


void CMyApp::MakeVoxelGrid(const glm::vec3& start,bool is_wireframe,const std::initializer_list<VertexAttributeDescriptor> vertexAttribList){
	std::cout << "Making Voxel Grid chunk objects" << std::endl;

	bool odd = 0;
	int size = m_chunkCount;
	bool evenSize = size%2==0;
	
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size;j++){
			for(int k = 0; k<size;k++){
				MeshObject<Vertex> chunk = MakeVoxelChunk(glm::ivec3(i,j,k),is_wireframe,odd);
				if(m_smoothShade){
					FixMeshNormals(chunk);
				}
				m_chunkMeshes[i][j][k] = CreateGLObjectFromMesh( chunk,vertexAttribList);
				if(!m_makeVoxelGrid)
					return;
				odd = !odd;
			}
			odd = !odd && evenSize || odd && !evenSize;
		}
		odd = !odd && evenSize || odd && !evenSize;
		
	}
}

MeshObject<Vertex> CMyApp::MakeVoxelChunk(const glm::ivec3& chunkPos,bool is_wireframe,bool odd){
	MeshObject<Vertex> chunk = {{},{}};
	CMyApp::ChunkData data= ExtractChunkData(chunkPos, glm::vec3(chunkPos)*3.f);
	glm::vec3 start;
	
	if(m_makeVoxelGrid){
		start = m_start + glm::vec3(chunkPos)*3.f;
	}
	else{
		start = glm::vec3(0);
	}
	
	for(int i = 0; i < 3; i++){
		
		for(int j = 0; j < 3; j++){

			for(int k = 0; k < 3; k++){
				odd = !odd;
				if((i!=1||j!=1||k!=1)&&!m_makeVoxelGrid)
					continue;


				std::array<float, 8> corners = {
					data[i  ][j  ][k  ],
					data[i+1][j  ][k  ],
					data[i  ][j  ][k+1],
					data[i+1][j  ][k+1],
					data[i  ][j+1][k  ],
					data[i+1][j+1][k  ],
					data[i  ][j+1][k+1],
					data[i+1][j+1][k+1]
				};
				chunk = MergeGeometry(chunk,MakeVoxel(corners, start+glm::vec3(i,j,k), false, odd,1.f));
				if(is_wireframe){
					m_frame = MergeGeometry(m_frame,MakeVoxel(corners, start+glm::vec3(i,j,k), true, odd,1.f));
				}
			}
		}
	}
	if(m_chunkBorders){
		std::array<float, 8> corners = {
			data[0][0][0],
			data[3][0][0],
			data[0][0][3],
			data[3][0][3],
			data[0][3][0],
			data[3][3][0],
			data[0][3][3],
			data[3][3][3]
		};
		m_frame = MergeGeometry(m_frame,MakeVoxel(corners, start, true, odd,3.f));
	}

	return chunk;
}

MeshObject<Vertex> CMyApp::MakeVoxel(const std::array<float,8>& weights,const glm::vec3& start,bool is_wireframe,bool rotations,float scale ){
	MeshObject<Vertex> voxelMesh = {{},{}};
	std::vector<std::vector<int>> indexes;
	if(!m_5tetra){
		indexes = {{
			{0,1,7,5},
			{0,3,7,1},
			{0,2,7,3},
			{0,6,7,2},
			{0,4,7,6},
			{0,5,7,4}
		}};
	}
	else {
		indexes = {{
			{0,3,5,1},
			{0,3,2,6},
			{0,5,6,4},
			{0,3,5,6},
			{3,5,6,7}
		}};
		if(rotations){
			std::vector<int> order = {2,3,0,1,6,7,4,5};
			std::vector<std::vector<int>> newindexes= {};
			for(int i : {0,1,2,3,4}){
				std::vector<int> vec = {};
				for(int j : {0,1,2,3}){
					vec.push_back( order[indexes[i][j]]);
				}
				newindexes.push_back(vec);
			}
			indexes = newindexes;

		}
		


	}

	for(auto ind : indexes){
		std::array<glm::vec4,4> tetrahedra = {};
		for(int i = 0; i < 4; i++){
			glm::vec3 pos = start + glm::vec3(ind[i]%2,(ind[i]/4)%2,(ind[i]/2)%2)*scale;
			tetrahedra[i] = glm::vec4(pos,weights[ind[i]]);
		}
		if (is_wireframe){
			voxelMesh = MergeGeometry(voxelMesh, MakeTetrahedraFrame(tetrahedra));
		}
		else{
			voxelMesh = MergeGeometry(voxelMesh, MakeTetrahedra(tetrahedra));
		}
	}
	if(!is_wireframe){
		FixMeshNormals(voxelMesh);
	}
	return voxelMesh;
}

MeshObject<Vertex> CMyApp::MakeTetrahedraFrame(const std::array<glm::vec4, 4>& points){
	MeshObject<Vertex> mesh = {{},{}};
	for(auto p :points){
		Vertex v = {glm::vec3(p),glm::vec3(1,0,0),glm::vec2((1.f-p.a)/2,(p.a+1.f)/2.f)};
		mesh.vertexArray.push_back(v);
	}
	mesh.indexArray = {
		0,1,2,
		0,2,3,
		0,3,1,
		1,2,3
	};
	return mesh;
}

void CMyApp::FixMeshNormals(MeshObject<Vertex>& mesh){
	//fixing normals
	for(int i = 0; i < mesh.vertexArray.size();i++){
		glm::vec3 norm = glm::vec3(0);
		glm::vec3 pos = mesh.vertexArray[i].position;
		if(m_makeVoxelGrid)
			pos -= m_start;
		
		glm::vec3 hx = glm::vec3(1,0,0)*m_epsilon;
		float evax = (SamplePerlinNoise(pos+hx)-SamplePerlinNoise(pos-hx))/(2*m_epsilon);

		glm::vec3 hy = glm::vec3(0,1,0)*m_epsilon;
		float evay = (SamplePerlinNoise(pos+hy)-SamplePerlinNoise(pos-hy))/(2*m_epsilon);

		glm::vec3 hz = glm::vec3(0,0,1)*m_epsilon;
		float evaz = (SamplePerlinNoise(pos+hz)-SamplePerlinNoise(pos-hz))/(2*m_epsilon);

		norm = hx*evax+hy*evay+hz*evaz;

		mesh.vertexArray[i].normal = glm::normalize(-norm);


	}

	//fixing index order
	for(int i = 0; i < mesh.indexArray.size();i+= 3){
		int a,b,c;
		a = mesh.indexArray[i];
		b = mesh.indexArray[i+1];
		c = mesh.indexArray[i+2];
		glm::vec3 x,y,z;
		x = mesh.vertexArray[a].position;
		y = mesh.vertexArray[b].position;
		z = mesh.vertexArray[c].position;
		
		glm::vec3 cross = glm::cross(x-y,z-y);
		if (glm::dot(cross,mesh.vertexArray[b].normal) > 0){
			mesh.indexArray[i+1] = c;
			mesh.indexArray[i+2] = b;
		}
	}
	
}
MeshObject<Vertex> CMyApp::MakeTetrahedra(const std::array<glm::vec4, 4>& points){
	MeshObject<Vertex> mesh = {{},{}};

	std::vector<int> inPoints = {};
	std::vector<int> outPoints = {};

	for(int i = 0; i < 4; i++){
		glm::vec4 p = points[i];
		if(p.a > 0.0){
			inPoints.push_back(i);
		}
		else{
			outPoints.push_back(i);
		}
	}
	
	if(inPoints.empty() || outPoints.empty())
		return mesh;

	if((int)inPoints.size() - (int)outPoints.size() == 0){
		mesh.indexArray = {0,2,1,1,2,3};
	}
	else {
		mesh.indexArray = {0,1,2};
	}
	
	for(int i : inPoints)
		for(int o: outPoints){
			glm::vec4 in = points[i];
			glm::vec4 out = points[o];

			float mixer = (-out.a)/(in.a-out.a);
			glm::vec2 texCoords = glm::vec2(mixer,mixer);

			texCoords = glm::vec2(1,1);

			glm::vec3 newVertPos = glm::vec3(glm::mix(out,in,mixer));
			Vertex v = {newVertPos,glm::vec3(0),texCoords};
			mesh.vertexArray.push_back(v);
		}
	return mesh;
}

MeshObject<Vertex> CMyApp::MergeGeometry(const MeshObject<Vertex>& meshA, const MeshObject<Vertex>& meshB){
	std::vector<GLuint> indexmap = {};
	MeshObject<Vertex> meshFinal = meshA;
	if(meshA.vertexArray.size() == 0)
		return meshB;
	if(meshB.vertexArray.size() == 0)
		return meshA;

	for (int i = 0;i < meshB.vertexArray.size();i++){
		Vertex nexVert = meshB.vertexArray[i];
		bool vertFound = false;
		GLuint size = meshFinal.vertexArray.size();

		for(GLuint j= 0 ;j < size && m_smoothShade;j++){
			if(glm::length( meshFinal.vertexArray[j].position- nexVert.position)<0.001){
				vertFound = true;
				size = j;
			}
		}

		indexmap.push_back(size);
		
		if (vertFound){
			meshFinal.vertexArray[size].normal = nexVert.normal+meshFinal.vertexArray[size].normal;
			meshFinal.vertexArray[size].texcoord = (nexVert.texcoord+meshFinal.vertexArray[size].texcoord)/2.f;
		}
		else{
			meshFinal.vertexArray.push_back(nexVert);
			
		}
	}
	for(GLuint i : meshB.indexArray){
		meshFinal.indexArray.push_back(indexmap[i]);
	}
	return meshFinal;
}

float CMyApp::GetPixelDepthAt(const glm::vec2& mouse){
	
	return 0;
}


void CMyApp::SetupDebugCallback()
{
	// Enable and set the debug callback function if we are in debug context
	GLint context_flags;
	glGetIntegerv(GL_CONTEXT_FLAGS, &context_flags);
	if (context_flags & GL_CONTEXT_FLAG_DEBUG_BIT) {
		glEnable(GL_DEBUG_OUTPUT);
		glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
		glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DEBUG_SEVERITY_NOTIFICATION, 0, nullptr, GL_FALSE);
		glDebugMessageControl(GL_DONT_CARE, GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR, GL_DONT_CARE, 0, nullptr, GL_FALSE);
		glDebugMessageCallback(SDL_GLDebugMessageCallback, nullptr);
	}
}

void CMyApp::InitShaders()
{
	m_programID = glCreateProgram();
	ProgramBuilder{ m_programID }
		.ShaderStage(GL_VERTEX_SHADER, "Shaders/Vert_PosNormTexShadow.vert")
		.ShaderStage(GL_FRAGMENT_SHADER, "Shaders/Frag_LightingShadow.frag")
.Link();

	// We don't need to output color, we only need the depth values
	// from vertex shader stage
	m_programPostprocessID = glCreateProgram();
	ProgramBuilder{ m_programPostprocessID }
		.ShaderStage(GL_VERTEX_SHADER, "Shaders/Vert_Shadow.vert")
		.Link();

	InitAxesShader();
}

void CMyApp::CleanShaders()
{
	glDeleteProgram(m_programID);
	glDeleteProgram(m_programPostprocessID);
	CleanAxesShader();
}

void CMyApp::InitAxesShader()
{
	m_programAxesID = glCreateProgram();
	ProgramBuilder{ m_programAxesID }
		.ShaderStage(GL_VERTEX_SHADER, "Shaders/Vert_axes.vert")
		.ShaderStage(GL_FRAGMENT_SHADER, "Shaders/Frag_PosCol.frag")
		.Link();
}

void CMyApp::CleanAxesShader()
{
	glDeleteProgram(m_programAxesID);
}


void CMyApp::InitGeometry()
{
	const std::initializer_list<VertexAttributeDescriptor> vertexAttribList =
	{
		{ 0, offsetof(Vertex, position	), 3, GL_FLOAT },
		{ 1, offsetof(Vertex, normal	), 3, GL_FLOAT },
		{ 2, offsetof(Vertex, texcoord	), 2, GL_FLOAT },
	};

	glm::vec3 start = glm::vec3(-m_chunkCount,-m_chunkCount,-m_chunkCount)*2.f+glm::vec3(.5);
	if(m_makeVoxelGrid){
		MakeVoxelGrid(start,false,vertexAttribList);
		m_voxelFrame = CreateGLObjectFromMesh(m_frame,vertexAttribList);
	}
	else{
		m_frame = {};
		MeshObject<Vertex> voxel = MakeVoxelChunk(glm::ivec3(0,0,0),false, false);

		m_voxelFrame = CreateGLObjectFromMesh(m_frame,vertexAttribList);;
		m_testVoxel = CreateGLObjectFromMesh(voxel,vertexAttribList);;
	}
}

void CMyApp::CleanGeometry()
{
	int size = m_chunkCount;	
	for(int i = 0; i < size; i++)
		for(int j = 0; j < size;j++)
			for(int k = 0; k<size;k++){
				CleanOGLObject(m_chunkMeshes[i][j][k]);
			}
	CleanOGLObject(m_voxelFrame);
	CleanOGLObject(m_testVoxel);
}

void CMyApp::InitTextures()
{
	glCreateSamplers( 1, &m_SamplerID );
	glSamplerParameteri( m_SamplerID, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
	glSamplerParameteri( m_SamplerID, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
	glSamplerParameteri( m_SamplerID, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
	glSamplerParameteri( m_SamplerID, GL_TEXTURE_MAG_FILTER, GL_LINEAR );


	glCreateTextures( GL_TEXTURE_2D, 1, &m_metalTextureID );

	glGenerateTextureMipmap( m_metalTextureID );
}

void CMyApp::CleanTextures()
{
	glDeleteTextures(1, &m_metalTextureID);
    glDeleteSamplers( 1, &m_SamplerID );
}




bool CMyApp::Init()
{
	
	SetupDebugCallback();
	std::cout << std::endl;

	srand(time(0));
	GeneratePerlinNoise();


	// Set a bluish clear color
	// glClear() will use this for clearing the color buffer.
	glClearColor(0.125f, 0.25f, 0.5f, 1.0f);

	InitShaders();
	InitGeometry();
	InitTextures();

	//
	// Other
	//

	glEnable(GL_CULL_FACE);	 // Enable discarding the back-facing faces.
	glCullFace(GL_BACK);     // GL_BACK: facets facing away from camera, GL_FRONT: facets facing towards the camera
	glEnable(GL_DEPTH_TEST); // Enable depth testing. (for overlapping geometry)

	// Camera
	m_camera.SetView(
		glm::vec3(0, 20, 20),	// From where we look at the scene - eye
		glm::vec3(0, 0, 0),		// Which point of the scene we are looking at - at
		glm::vec3(0, 1, 0)		// Upwards direction - up
	);
	m_camera.SetAngle(90);
	m_cameraManipulator.SetCamera(&m_camera);

	return true;
}

void CMyApp::Clean()
{
	CleanShaders();
	CleanGeometry();
	CleanTextures();
}

void CMyApp::Update(const SUpdateInfo& updateInfo)
{
	m_cameraManipulator.Update(updateInfo.DeltaTimeInSec);

	m_ElapsedTimeInSec = updateInfo.ElapsedTimeInSec;

    // World transform of the suzannes
}

void CMyApp::DrawAxes()
{
	glUseProgram(m_programAxesID);
	glm::mat4 axisWorld = glm::translate(m_camera.GetAt());
	glUniformMatrix4fv(ul("viewProj"), 1, GL_FALSE, glm::value_ptr(m_camera.GetViewProj()));
	glUniformMatrix4fv(ul("world"), 1, GL_FALSE, glm::value_ptr(axisWorld));

	// We always want to see it, regardless of whether there is an object in front of it
	glDisable(GL_DEPTH_TEST);
	glDrawArrays(GL_LINES, 0, 6);
	glUseProgram(0);
	glEnable(GL_DEPTH_TEST);
}

void CMyApp::RenderGeometry(const glm::mat4& viewProj, GLuint program)
{
	glUseProgram(program);

	glUniformMatrix4fv(ul("world"), 1, GL_FALSE, glm::value_ptr(glm::mat4(1)));
	glUniformMatrix4fv(ul("viewProj"), 1, GL_FALSE, glm::value_ptr(viewProj));
	glUniformMatrix4fv(ul( "worldInvTransp"), 1, GL_FALSE, glm::value_ptr(glm::mat4(1)));

	glBindTextureUnit( 0, m_metalTextureID );
	glBindSampler( 0, m_SamplerID );
	glUniform1i(ul("textureImage"), 0);

	glUniform4fv( ul( "lightPosition" ), 1, glm::value_ptr( m_lightPosition ) );
	glUniform3fv( ul( "La" ), 1, glm::value_ptr( m_La ) );
	glUniform3fv( ul( "Ld" ), 1, glm::value_ptr( m_Ld ) );
	glUniform3fv( ul( "Ls" ), 1, glm::value_ptr( m_Ls ) );

	glUniform4fv( ul( "lightPosition2" ), 1, glm::value_ptr( m_lightPosition2 ) );
	glUniform3fv( ul( "La2" ), 1, glm::value_ptr( m_La2 ) );
	glUniform3fv( ul( "Ld2" ), 1, glm::value_ptr( m_Ld2 ) );
	glUniform3fv( ul( "Ls2" ), 1, glm::value_ptr( m_Ls2 ) );

	glUniform3fv(ul("cameraPosition"), 1, glm::value_ptr(m_camera.GetEye()));

	
	if(m_makeVoxelGrid){
		int size = m_chunkCount;
		for(int i = 0; i < size; i++)
			for(int j = 0; j < size;j++)
				for(int k = 0; k<size;k++){
					glBindVertexArray(m_chunkMeshes[i][j][k].vaoID );
					glDrawElements(GL_TRIANGLES,m_chunkMeshes[i][j][k].count,GL_UNSIGNED_INT,0);
				}

	}
	else{
		glBindVertexArray(m_testVoxel.vaoID );
		glDrawElements(GL_TRIANGLES,m_testVoxel.count,GL_UNSIGNED_INT,0);

	}

	if(!m_makeVoxelGrid || m_chunkBorders){
		glDisable(GL_CULL_FACE);	 
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
		glUniform1i(ul("ignoreNormals"), 1);

		glBindVertexArray(m_voxelFrame.vaoID);
		glDrawElements(GL_TRIANGLES,m_voxelFrame.count,GL_UNSIGNED_INT,0);

		glEnable(GL_CULL_FACE);	 // Enable discarding the back-facing faces.
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
		glUniform1i(ul("ignoreNormals"), 0);
	}


	// We can unbind them
	//glBindTextureUnit( 0, 0);
    glBindSampler( 0, 0 );
	//glBindTextureUnit(1, 0);
	glBindSampler( 1, 0 );
	glUseProgram(0);
	glBindVertexArray(0);
}

void CMyApp::Render()
{

	// 2.
	// Draw mesh to screen

	glBindFramebuffer(GL_FRAMEBUFFER, 0);				// default framebuffer (the backbuffer)
	glViewport(0, 0, m_width, m_height);				// We need to set the render area back
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// clearing the default fbo
	RenderGeometry(m_camera.GetViewProj(), m_programID);

	DrawAxes();
}

void CMyApp::RenderGUI()
{
	// ImGui::ShowDemoWindow();

	ImGui::Begin("Lighting window");
	{
		ImGui::Text("Light 1");
		ImGui::SliderFloat3("light_pos_1", &m_lightPosition.x, -m_chunkCount*3, m_chunkCount*3);
		ImGui::SliderFloat3("La_1", &m_La.x, 0.f, 1.f);
		ImGui::SliderFloat3("Ld_1", &m_Ld.x, 0.f, 1.f);
		ImGui::SliderFloat3("Ls_1", &m_Ls.x, 0.f, 1.f);

		ImGui::Text("Light 2");
		ImGui::SliderFloat3("light_pos_2", &m_lightPosition2.x, -m_chunkCount*3, m_chunkCount*3);
		ImGui::SliderFloat3("La_2", &m_La2.x, 0.f, 1.f);
		ImGui::SliderFloat3("Ld_2", &m_Ld2.x, 0.f, 1.f);
		ImGui::SliderFloat3("Ls_2", &m_Ls2.x, 0.f, 1.f);

		ImGui::Separator();
		if( ImGui::Checkbox("orthogonal", &m_ortho)){
			m_camera.SetProj(
				m_camera.GetAngle(), 
				m_camera.GetAspect(), 
				m_camera.GetZNear(),
				m_camera.GetZFar(),
				!m_ortho
			);
		}
		if(
			ImGui::Checkbox("Smooth shading" , &m_smoothShade) ||
			ImGui::Checkbox("5 tethra voxels" , &m_5tetra)  ||
			ImGui::Checkbox("Chunk borders" , &m_chunkBorders) ||
			ImGui::SliderFloat("Epsilon",&m_epsilon ,0.001f,3.f)
		){
				CleanGeometry();
				InitGeometry();
		}
 
		ImGui::Separator();
		if(
			ImGui::Checkbox("Large Grid" , &m_makeVoxelGrid) ||
			ImGui::SliderFloat4("Fisrt 4 point",&m_testFrame[0],-1.f,1.f) ||
			ImGui::SliderFloat4("Second 4 point",&m_testFrame[4],-1.f,1.f)
		){
				GeneratePerlinNoise();
				CleanGeometry();
				InitGeometry();
		}
		ImGui::Separator();
		if(ImGui::Button("RefreshNoise")){
			GeneratePerlinNoise();
			CleanGeometry();
			InitGeometry();
		}

	}
	ImGui::End();
}

// https://wiki.libsdl.org/SDL3/SDL_KeyboardEvent
// https://wiki.libsdl.org/SDL3/SDL_Keysym
// https://wiki.libsdl.org/SDL3/SDL_Keycode
// https://wiki.libsdl.org/SDL3/SDL_Keymod

void CMyApp::KeyboardDown(const SDL_KeyboardEvent& key)
{
	if (!key.repeat) // Triggers only once when held
	{
		if (key.key == SDLK_F5 && key.mod & SDL_KMOD_CTRL) // CTRL + F5
		{
			CleanShaders();
			InitShaders();
		}
		if (key.key == SDLK_F1) // F1
		{
			GLint polygonModeFrontAndBack[2] = {};
			// https://registry.khronos.org/OpenGL-Refpages/gl4/html/glGet.xhtml
			glGetIntegerv(GL_POLYGON_MODE, polygonModeFrontAndBack); // Query the current polygon mode. It gives the front and back modes separately.
			GLenum polygonMode = (polygonModeFrontAndBack[0] != GL_FILL ? GL_FILL : GL_LINE); // Switch between FILL and LINE
			// https://registry.khronos.org/OpenGL-Refpages/gl4/html/glPolygonMode.xhtml
			glPolygonMode(GL_FRONT_AND_BACK, polygonMode); // Set the new polygon mode
		}
	}
	m_cameraManipulator.KeyboardDown(key);
}

void CMyApp::KeyboardUp(const SDL_KeyboardEvent& key)
{
	m_cameraManipulator.KeyboardUp(key);
}

// https://wiki.libsdl.org/SDL3/SDL_MouseMotionEvent

void CMyApp::MouseMove(const SDL_MouseMotionEvent& mouse)
{
	m_cameraManipulator.MouseMove(mouse);
}

// https://wiki.libsdl.org/SDL3/SDL_MouseButtonEvent

void CMyApp::MouseDown(const SDL_MouseButtonEvent& mouse)
{
}

void CMyApp::MouseUp(const SDL_MouseButtonEvent& mouse)
{
}

// https://wiki.libsdl.org/SDL3/SDL_MouseWheelEvent

void CMyApp::MouseWheel(const SDL_MouseWheelEvent& wheel)
{
	m_cameraManipulator.MouseWheel(wheel);
}

// New window size
void CMyApp::Resize(int _w, int _h)
{
	glViewport(0, 0, _w, _h);
	m_camera.SetAspect(static_cast<float>(_w) / _h);
	m_width = _w;
	m_height = _h;
}

// Other SDL events
// https://wiki.libsdl.org/SDL3/SDL_Event

void CMyApp::OtherEvent(const SDL_Event& ev)
{
}
