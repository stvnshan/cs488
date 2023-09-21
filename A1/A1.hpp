// Termm--Fall 2023

#pragma once

#include <glm/glm.hpp>

#include "cs488-framework/CS488Window.hpp"
#include "cs488-framework/OpenGLImport.hpp"
#include "cs488-framework/ShaderProgram.hpp"

#include "maze.hpp"

class A1 : public CS488Window {
public:
	A1();
	virtual ~A1();

protected:
	virtual void init() override;
	virtual void appLogic() override;
	virtual void guiLogic() override;
	virtual void draw() override;
	virtual void cleanup() override;

	virtual bool cursorEnterWindowEvent(int entered) override;
	virtual bool mouseMoveEvent(double xPos, double yPos) override;
	virtual bool mouseButtonInputEvent(int button, int actions, int mods) override;
	virtual bool mouseScrollEvent(double xOffSet, double yOffSet) override;
	virtual bool windowResizeEvent(int width, int height) override;
	virtual bool keyInputEvent(int key, int action, int mods) override;
	

private:
	void initGrid();
	bool ifWall();
	void redrawCubes();
	void digg();
	void reset();

	// Fields related to the shader and uniforms.
	ShaderProgram m_shader;
	GLint P_uni; // Uniform location for Projection matrix.
	GLint V_uni; // Uniform location for View matrix.
	GLint M_uni; // Uniform location for Model matrix.
	GLint col_uni;   // Uniform location for cube colour.

	// Fields related to grid geometry.
	GLuint m_grid_vao; // Vertex Array Object
	GLuint m_grid_vbo; // Vertex Buffer Object

	GLuint cube_vao;
	GLuint cube_vbo;

	GLuint sphere_vao;
	GLuint sphere_vbo;

	GLuint s_indices_vao;
	GLuint s_indices_vbo;

	// Matrices controlling the camera and projection.
	glm::mat4 proj;
	glm::mat4 view;

	float colour[3];
	float mazeC[3];
	float floorC[3];
	float sphereC[3];
	int current_col;
	


	Maze m = Maze(16);
	int blockNum;


	int sx = 0;
	int sz = 0;

	int wallS = 1;

	GLfloat fov = 1;
	GLfloat radius = 0.5f;  
	GLfloat camX = 0;
	GLfloat camZ = 1; 
	GLfloat camY = 0;   
	GLfloat lastX = 0; 
	GLfloat lastY = 0;
	GLfloat xoffset;    
	GLfloat camRotDistance = 0.0f;    
	GLboolean isFirstMouse = GL_TRUE;    
	GLboolean isLeftMousePress = GL_FALSE;


	float resistanceD = 0;
	float deltaTime = 0;
	float lastFrame = 0;
	float curFrame = 0;
	float speed = 0;


	bool resis = false;
	bool firsttime = true;

	bool digged = false;

	
};
