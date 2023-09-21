// Termm--Fall 2023

#include "A1.hpp"
#include "cs488-framework/GlErrorCheck.hpp"

#include <iostream>

#include <sys/types.h>
#include <unistd.h>

#include <imgui/imgui.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace glm;
using namespace std;

static const size_t DIM = 16;

//----------------------------------------------------------------------------------------
// Constructor
A1::A1()
	: current_col( 0 )
{
	colour[0] = 0.0f;
	colour[1] = 0.0f;
	colour[2] = 0.0f;
	mazeC[0] = 1;
	mazeC[1] = 0;
	mazeC[2] = 0;
	sphereC[0] = 0;
	sphereC[1] = 1;
	sphereC[2] = 0;
	floorC[0] = 0;
	floorC[1] = 0;
	floorC[2] = 0;

}

//----------------------------------------------------------------------------------------
// Destructor
A1::~A1()
{}

//----------------------------------------------------------------------------------------
/*
 * Called once, at program start.
 */
void A1::init()
{
	// Initialize random number generator
	int rseed=getpid();
	srandom(rseed);
	// Print random number seed in case we want to rerun with
	// same random numbers
	cout << "Random number seed = " << rseed << endl;
	

	// DELETE FROM HERE...
	// Maze m(DIM);
	// m.digMaze();
	// m.printMaze();
	// ...TO HERE
	// m = Maze(DIM);
	blockNum = 0;
	
	// Set the background colour.
	glClearColor( 0.3, 0.5, 0.7, 1.0 );

	// Build the shader
	m_shader.generateProgramObject();
	m_shader.attachVertexShader(
		getAssetFilePath( "VertexShader.vs" ).c_str() );
	m_shader.attachFragmentShader(
		getAssetFilePath( "FragmentShader.fs" ).c_str() );
	m_shader.link();

	// Set up the uniforms
	P_uni = m_shader.getUniformLocation( "P" );
	V_uni = m_shader.getUniformLocation( "V" );
	M_uni = m_shader.getUniformLocation( "M" );
	col_uni = m_shader.getUniformLocation( "colour" );

	initGrid();

	// Set up initial view and projection matrices (need to do this here,
	// since it depends on the GLFW window being set up correctly).
	view = glm::lookAt( 
		glm::vec3( 0.0f, 2.*float(DIM)*2.0*M_SQRT1_2, float(DIM)*2.0*M_SQRT1_2 ),
		glm::vec3( 0.0f, 0.0f, 0.0f ),
		glm::vec3( 0.0f, 1.0f, 0.0f ) );

	proj = glm::perspective( 
		glm::radians( 30.0f ),
		float( m_framebufferWidth ) / float( m_framebufferHeight ),
		1.0f, 1000.0f );

	
}

void A1::initGrid()
{
	// size_t sz = 3 * 2 * 2 * (DIM+3);

	// float *verts = new float[ sz ];
	// size_t ct = 0;
	// for( int idx = 0; idx < DIM+3; ++idx ) {
	// 	verts[ ct ] = -1;
	// 	verts[ ct+1 ] = 0;
	// 	verts[ ct+2 ] = idx-1;
	// 	verts[ ct+3 ] = DIM+1;
	// 	verts[ ct+4 ] = 0;
	// 	verts[ ct+5 ] = idx-1;

	// 	ct += 6;


	// 	verts[ ct ] = idx-1;
	// 	verts[ ct+1 ] = 0;
	// 	verts[ ct+2 ] = -1;
	// 	verts[ ct+3 ] = idx-1;
	// 	verts[ ct+4 ] = 0;
	// 	verts[ ct+5 ] = DIM+1;

	// 	ct += 6;
	// }

	float *verts = new float[ 3 * 6 ];
	verts[0] = -1;
	verts[1] = 0;
	verts[2] = -1;

	verts[3] = DIM+1;
	verts[4] = 0;
	verts[5] = -1;

	verts[6] = -1;
	verts[7] = 0;
	verts[8] = DIM + 1;


	verts[9] = DIM + 1;
	verts[10] = 0;
	verts[11] = DIM + 1;

	verts[12] = DIM+1;
	verts[13] = 0;
	verts[14] = -1;

	verts[15] = -1;
	verts[16] = 0;
	verts[17] = DIM + 1;








	

	//sphere vertices

	// std::vector<float>().swap(vertices);
	// std::vector<float>().swap(normals);
	// std::vector<float>().swap(texCoords);
	int sectorCount = 50;
	int stackCount = 50;
	float PI = 3.14;
	float *sphere = new float[3 * (stackCount + 1) * (sectorCount + 1)];
	float x, y, z, xy;                              // vertex position
	float nx, ny, nz, lengthInv = 1.0f / radius;    // vertex normal
	float s, t;                                     // vertex texCoord

	float sectorStep = 2 * PI / sectorCount;
	float stackStep = PI / stackCount;
	float sectorAngle, stackAngle;
	
	int counter = 0;
	for(int i = 0; i <= stackCount; ++i)
	{
		stackAngle = PI / 2 - i * stackStep;        // starting from pi/2 to -pi/2
		xy = radius * cosf(stackAngle);             // r * cos(u)
		z = radius * sinf(stackAngle) + 0.5;              // r * sin(u)

		// add (sectorCount+1) vertices per stack
		// first and last vertices have same position and normal, but different tex coords
		for(int j = 0; j <= sectorCount; ++j)
		{
			sectorAngle = j * sectorStep;           // starting from 0 to 2pi

			// vertex position (x, y, z)
			x = xy * cosf(sectorAngle) + 0.5;             // r * cos(u) * cos(v)
			y = xy * sinf(sectorAngle);             // r * cos(u) * sin(v)
			sphere[counter] = x;
			sphere[counter + 1] = y;
			sphere[counter + 2] = z;
			counter +=3;

			// // normalized vertex normal (nx, ny, nz)
			// nx = x * lengthInv;
			// ny = y * lengthInv;
			// nz = z * lengthInv;
			// normals.push_back(nx);
			// normals.push_back(ny);
			// normals.push_back(nz);

			// // vertex tex coord (s, t) range between [0, 1]
			// s = (float)j / sectorCount;
			// t = (float)i / stackCount;
			// texCoords.push_back(s);
			// texCoords.push_back(t);
		}
	}

	// std::vector<int> indices;
	// std::vector<int> lineIndices;
	GLuint * indices = new GLuint[6  *(stackCount -1) * sectorCount];
	int k1, k2;
	counter = 0;
	for(int i = 0; i < stackCount; ++i)
	{
		k1 = i * (sectorCount + 1);     // beginning of current stack
		k2 = k1 + sectorCount + 1;      // beginning of next stack

		for(int j = 0; j < sectorCount; ++j, ++k1, ++k2)
		{
			// 2 triangles per sector excluding first and last stacks
			// k1 => k2 => k1+1
			if(i != 0)
			{
				// indices.push_back(k1);
				// indices.push_back(k2);
				// indices.push_back(k1 + 1);
				indices[counter] = k1;
				indices[counter + 1] = k2;
				indices[counter + 2] = k1 + 1;
				counter += 3;
			}

			// k1+1 => k2 => k2+1
			if(i != (stackCount-1))
			{
				// indices.push_back(k1 + 1);
				// indices.push_back(k2);
				// indices.push_back(k2 + 1);
				indices[counter] = k1+1;
				indices[counter + 1] = k2;
				indices[counter + 2] = k2 + 1;
				counter += 3;

			}

			// // store indices for lines
			// // vertical lines for all stacks, k1 => k2
			// lineIndices.push_back(k1);
			// lineIndices.push_back(k2);
			// if(i != 0)  // horizontal lines except 1st stack, k1 => k+1
			// {
			// 	lineIndices.push_back(k1);
			// 	lineIndices.push_back(k1 + 1);
			// }
		}
	}


	



	// Create the vertex array to record buffer assignments.
	glGenVertexArrays( 1, &m_grid_vao );
	glBindVertexArray( m_grid_vao );

	// Create the grid vertex buffer
	glGenBuffers( 1, &m_grid_vbo );
	glBindBuffer( GL_ARRAY_BUFFER, m_grid_vbo );
	glBufferData( GL_ARRAY_BUFFER, 3*6*sizeof(float),
		verts, GL_STATIC_DRAW );


	// Specify the means of extracting the position values properly.
	GLint posAttrib = m_shader.getAttribLocation( "position" );
	glEnableVertexAttribArray( posAttrib );
	glVertexAttribPointer( posAttrib, 3, GL_FLOAT, GL_FALSE, 0, nullptr );







//sphere
	glGenVertexArrays( 1, &sphere_vao);
	glBindVertexArray( sphere_vao);
	glGenBuffers( 1, &sphere_vbo );
	glBindBuffer( GL_ARRAY_BUFFER, sphere_vbo);
	glBufferData( GL_ARRAY_BUFFER, 3*51*51 *sizeof(float),
		sphere, GL_STATIC_DRAW );

//s_indices
	glGenBuffers( 1, &s_indices_vbo );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, s_indices_vbo);
	glBufferData( GL_ELEMENT_ARRAY_BUFFER, 6*49*50 *sizeof(GLuint),
		indices, GL_STATIC_DRAW );


	posAttrib = m_shader.getAttribLocation( "position" );
	glEnableVertexAttribArray( posAttrib );
	glVertexAttribPointer( posAttrib, 3, GL_FLOAT, GL_FALSE, 0, nullptr );

	// Reset state to prevent rogue code from messing with *my* 
	// stuff!
	glBindVertexArray( 0 );
	glBindBuffer( GL_ARRAY_BUFFER, 0 );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );

	// OpenGL has the buffer now, there's no need for us to keep a copy.
	delete [] verts;
	//delete[] cube;
	delete[] sphere;

	CHECK_GL_ERRORS;
}

//----------------------------------------------------------------------------------------
/*
 * Called once per frame, before guiLogic().
 */
void A1::appLogic()
{	

	if(! ImGui::IsMouseDragging() && !firsttime){
		
		// curFrame = glfwGetTime();
		// deltaTime = curFrame - lastFrame;
		// lastFrame = curFrame; 
		if(!resis){
			speed = resistanceD/deltaTime;
		}
		// camX = sin(3.14/2 + resistanceD + speed* (glfwGetTime() - curFrame) ) * radius  ;
		// camZ = cos(3.14/2 + resistanceD + speed* (glfwGetTime() - curFrame)) * radius  ;
		
		camX = cos( 3.14/2 +  camRotDistance + speed /500 ) * radius;
		camZ = sin( 3.14/2 + camRotDistance + speed /500) * radius;


		camRotDistance += speed/500;
		speed = speed * 0.99;

		

		view = glm::lookAt(
                glm::vec3(float(DIM)*M_SQRT1_2 * camX *4 * fov, 2.*float(DIM)*2.0*M_SQRT1_2*fov,  float(DIM)*M_SQRT1_2 * camZ * fov * 4),  
                glm::vec3(0, 0, 0),    
                glm::vec3(0.0f, 1.0f, 0.0f) 
           );
		resis = true;
	}
	

	// Place per frame, application logic here ...
}

//----------------------------------------------------------------------------------------
/*
 * Called once per frame, after appLogic(), but before the draw() method.
 */
void A1::guiLogic()
{
	// We already know there's only going to be one window, so for 
	// simplicity we'll store button states in static local variables.
	// If there was ever a possibility of having multiple instances of
	// A1 running simultaneously, this would break; you'd want to make
	// this into instance fields of A1.
	static bool showTestWindow(false);
	static bool showDebugWindow(true);

	ImGuiWindowFlags windowFlags(ImGuiWindowFlags_AlwaysAutoResize);
	float opacity(0.5f);

	ImGui::Begin("Debug Window", &showDebugWindow, ImVec2(100,100), opacity, windowFlags);
		if( ImGui::Button( "Quit Application" ) ) {
			glfwSetWindowShouldClose(m_window, GL_TRUE);
		}

		// Eventually you'll create multiple colour widgets with
		// radio buttons.  If you use PushID/PopID to give them all
		// unique IDs, then ImGui will be able to keep them separate.
		// This is unnecessary with a single colour selector and
		// radio button, but I'm leaving it in as an example.

		// Prefixing a widget name with "##" keeps it from being
		// displayed.

		ImGui::PushID( 0 );
		
		if( ImGui::RadioButton( "Floor", &current_col, 0 ) ) {
			floorC[0] = colour[0];
			floorC[1] = colour[1];
			floorC[2] = colour[2];
			// Select this colour.
		}
		ImGui::SameLine();
		if( ImGui::RadioButton( "Maze", &current_col, 0 ) ) {
			mazeC[0] = colour[0];
			mazeC[1] = colour[1];
			mazeC[2] = colour[2];
			// Select this colour.
		}
		ImGui::SameLine();
		if( ImGui::RadioButton( "Sphere", &current_col, 0 ) ) {
			sphereC[0] = colour[0];
			sphereC[1] = colour[1];
			sphereC[2] = colour[2];
			// Select this colour.
		}
		ImGui::ColorEdit3( "##Colour", colour );

		ImGui::PopID();

		
		if(ImGui::Button("Dig")){
			
			digg();
		}
		if(ImGui::Button("Reset")){
			reset();
		}

/*
		// For convenience, you can uncomment this to show ImGui's massive
		// demonstration window right in your application.  Very handy for
		// browsing around to get the widget you want.  Then look in 
		// shared/imgui/imgui_demo.cpp to see how it's done.
		if( ImGui::Button( "Test Window" ) ) {
			showTestWindow = !showTestWindow;
		}
*/

		ImGui::Text( "Framerate: %.1f FPS", ImGui::GetIO().Framerate );

	ImGui::End();

	if( showTestWindow ) {
		ImGui::ShowTestWindow( &showTestWindow );
	}
}



//----------------------------------------------------------------------------------------
/*
 * Called once per frame, after guiLogic().
 */
void A1::draw()
{
	// Create a global transformation for the model (centre it).
	mat4 W;
	mat4 W1;
	mat4 W2;
	// W = glm::translate( W, vec3( -float(DIM)/2.0f, 0, -float(DIM)/2.0f ) );
	W = glm::translate( W, vec3( -float(DIM)/2.0f, 0, -float(DIM)/2.0f ) );


	m_shader.enable();
		glEnable( GL_DEPTH_TEST );

		
		// int xpos, ypos;
		// glfwGetGPos(&xpos, &ypos);
		// glfwSetMousePos(1024/2, 768/2);
		// const float radius = 10.0f;
		// float camX = sin(xpos) * radius;
		// float camZ = cos(ypos) * radius;
		// glm::mat4 view;
		// view = glm::lookAt(glm::vec3(camX, 0.0, camZ), glm::vec3(0.0, 0.0, 0.0), glm::vec3(0.0, 1.0, 0.0));
		
		
		
		

		//glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);




		glUniformMatrix4fv( P_uni, 1, GL_FALSE, value_ptr( proj ) );
		glUniformMatrix4fv( V_uni, 1, GL_FALSE, value_ptr( view ) );
		glUniformMatrix4fv( M_uni, 1, GL_FALSE, value_ptr( W ) );

		

		// Just draw the grid for now.

		glBindVertexArray( m_grid_vao );
		glUniform3f( col_uni, floorC[0], floorC[1], floorC[2] );
		glDrawArrays( GL_TRIANGLES, 0, 6 );



		//Draw the cubes

		// Set the uniform's value.

		if(digged){
			W2 = glm::scale( W, vec3(1,wallS,1));
			glUniformMatrix4fv( M_uni, 1, GL_FALSE, value_ptr( W2 ) );
			glBindVertexArray( cube_vao );
			glUniform3f( col_uni, mazeC[0], mazeC[1], mazeC[2] );
			glDrawArrays( GL_TRIANGLES, 0, 3*6*2*blockNum );
		}
		
		

		//Draw sphere
		
		
		
		W1 = glm::translate( W, vec3( sx, 0, sz ) );
		glUniformMatrix4fv( M_uni, 1, GL_FALSE, value_ptr( W1 ) );

		glBindVertexArray( sphere_vao );
		glUniform3f( col_uni, sphereC[0], sphereC[1], sphereC[2]);
		//glDrawArrays( GL_TRIANGLES, 0, 51*51  );


		glDrawElements( GL_TRIANGLES,  6*49*50 , GL_UNSIGNED_INT, 0);




		// Highlight the active square.



		
	m_shader.disable();

	// Restore defaults
	glBindVertexArray( 0 );

	CHECK_GL_ERRORS;
}

//----------------------------------------------------------------------------------------
/*
 * Called once, after program is signaled to terminate.
 */
void A1::cleanup()
{}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles cursor entering the window area events.
 */
bool A1::cursorEnterWindowEvent (
		int entered
) {
	bool eventHandled(false);

	// Fill in with event handling code...

	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles mouse cursor movement events.
 */
bool A1::mouseMoveEvent(double xPos, double yPos) 
{
	bool eventHandled(false);

	if (!ImGui::IsMouseHoveringAnyWindow()) {
		// Put some code here to handle rotations.  Probably need to
		// check whether we're *dragging*, not just moving the mouse.
		// Probably need some instance variables to track the current
		// rotation amount, and maybe the previous X position (so 
		// that you can rotate relative to the *change* in X.
		
		// if(isLeftMousePress){
		// 	if (isFirstMouse){
		// 		lastX = xPos;
		// 		isFirstMouse = GL_FALSE;
		// 	}
		// 	xoffset = xPos - lastX; 
		// 	lastX = xPos;
		// 	camRotDistance += xoffset / 10;    
		// 	camX = sin(camRotDistance) * radius;
		// 	camZ = cos(camRotDistance) * radius;
		// }else{
		// 	lastX = xPos;
		// 	lastY = yPos;
		// }
		if(ImGui::IsMouseDragging()){
			resis = false;
			firsttime = false;
			if (isFirstMouse){
				lastX = xPos;
				isFirstMouse = GL_FALSE;
			}
			xoffset = xPos - lastX; 
			lastX = xPos;
			camRotDistance += xoffset / 70;    
			camX = cos(3.14/2 + camRotDistance) * radius;
			camZ = sin(3.14/2 + camRotDistance) * radius;
			resistanceD = xoffset / 70;
			curFrame = glfwGetTime();
			deltaTime = curFrame - lastFrame;
			lastFrame = curFrame; 

			
			view = glm::lookAt(
                glm::vec3(float(DIM)*M_SQRT1_2 * camX *4 * fov, 2.*float(DIM)*2.0*M_SQRT1_2 * fov,  float(DIM)*M_SQRT1_2 * camZ * fov * 4),  
                glm::vec3(0, 0, 0),    
                glm::vec3(0.0f, 1.0f, 0.0f) 
           );
		}else{
			lastX = xPos;
			lastY = yPos;
			
		}

	}

	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles mouse button events.
 */
bool A1::mouseButtonInputEvent(int button, int actions, int mods) {
	bool eventHandled(false);

	if (!ImGui::IsMouseHoveringAnyWindow()) {
		// The user clicked in the window.  If it's the left
		// mouse button, initiate a rotation.
		
		
		
	}

	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles mouse scroll wheel events.
 */
bool A1::mouseScrollEvent(double xOffSet, double yOffSet) {
	bool eventHandled(false);

	// Zoom in or out.
	fov -= yOffSet;
    if (fov < 0.2f)
        fov = 0.2f;
    if (fov > 5.0f)
        fov = 5.0f;

	view = glm::lookAt(
                glm::vec3(float(DIM)*M_SQRT1_2 * camX *4 * fov, 2.*float(DIM)*2.0*M_SQRT1_2 * fov,  float(DIM)*M_SQRT1_2 * camZ * 4 * fov),  
                glm::vec3(0, 0, 0),    
                glm::vec3(0.0f, 1.0f, 0.0f) 
           );

	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles window resize events.
 */
bool A1::windowResizeEvent(int width, int height) {
	bool eventHandled(false);

	// Fill in with event handling code...

	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles key input events.
 */
bool A1::keyInputEvent(int key, int action, int mods) {
	bool eventHandled(false);

	// Fill in with event handling code...
	if( action == GLFW_PRESS ) {
		// Respond to some key events.
		if(key == GLFW_KEY_RIGHT){
			if(sx+1 > DIM + 1){
				
			}else if(m.getValue(sz,sx+1) != 1){
				sx++;
			}else if(glfwGetKey(m_window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS){
				sx++;
				m.setValue(sz,sx,0);
				redrawCubes();
			}
		
		}else if(key == GLFW_KEY_LEFT){
			if(sx-1 < -1){
				
			}else if(m.getValue(sz,sx-1) != 1){
				sx--;
			}else if(glfwGetKey(m_window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS){
				sx--;
				m.setValue(sz,sx,0);
				redrawCubes();
			}
			
		}else if(key == GLFW_KEY_UP){
			if(sz-1 < -1){
				
			}
			else if(m.getValue(sz-1,sx) != 1){
				sz--;
			}else if(glfwGetKey(m_window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS){
				sz--;
				m.setValue(sz,sx,0);
				redrawCubes();
			}
			
		}else if(key == GLFW_KEY_DOWN){
			if(sz + 1 > DIM + 1){}
			else if(m.getValue(sz + 1, sx) != 1 ){
				sz++;
			}else if(glfwGetKey(m_window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS){
				sz++;
				m.setValue(sz,sx,0);
				redrawCubes();
			}
			
		}else if(key == GLFW_KEY_SPACE){
			wallS++;
		}else if(key == GLFW_KEY_BACKSPACE){
			if(wallS >= 1){
				wallS--;
			}
		}else if(key == GLFW_KEY_Q){
			glfwSetWindowShouldClose(m_window, GL_TRUE);
		}else if(key == GLFW_KEY_D){
			digg();
		}else if(key == GLFW_KEY_R){
			reset();
		}
	}

	return eventHandled;
}



void A1::redrawCubes(){
	blockNum--;
	float *cube = new float[3*3*12*blockNum];
	int counter = 0;
	for(int i = 0; i < DIM; i++){
		for(int j = 0; j < DIM; j++){
			int curV = m.getValue(j,i);
			if(curV == 1){
				int tmp[] = {
					i,0,j,   i,1,j,   i,1,j+1,
					i,0,j,   i,0,j+1,   i,1,j+1,
					i+1,0,j,   i+1,1,j,   i+1,1,j+1,
					i+1,0,j,   i+1,0,j+1,   i+1,1,j+1,

					i+1,0,j,  i,0,j, i+1,1,j,
					i,1,j,  i,0,j, i+1,1,j,
					i+1,0,j+1,  i,0,j+1, i+1,1,j+1,
					i,1,j+1,  i,0,j+1, i+1,1,j+1,
					
					i,0,j,  i,0,j+1, i+1,0,j,
					i,0,j+1, i+1,0,j+1,i+1,0,j,
					i,1,j,  i,1,j+1, i+1,1,j,
					i,1,j+1, i+1,1,j+1,i+1,1,j,

				};
				for(int i = 0; i < 3 * 3 * 12;i++){
					cube[counter + i] = tmp[i];
				}
				counter += 3*3*12;
			}
		}
	}



	glGenVertexArrays( 1, &cube_vao);
	glBindVertexArray( cube_vao);
	glGenBuffers( 1, &cube_vbo );
	glBindBuffer( GL_ARRAY_BUFFER, cube_vbo);
	glBufferData( GL_ARRAY_BUFFER, 3*3*12*blockNum *sizeof(float),
		cube, GL_STATIC_DRAW );



	GLint posAttrib = m_shader.getAttribLocation( "position" );
	glEnableVertexAttribArray( posAttrib );
	glVertexAttribPointer( posAttrib, 3, GL_FLOAT, GL_FALSE, 0, nullptr );



	glBindVertexArray( 0 );
	glBindBuffer( GL_ARRAY_BUFFER, 0 );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );

}



void A1::digg(){
	digged = true;
	m.digMaze();
	blockNum = 0;
	for(int i = 0; i < DIM; i++){
		for(int j = 0; j < DIM; j++){
		
			if(m.getValue(j,i) == 1){
				blockNum++;
			}else if(j == 0){
				sx = i;
			}
	
		}
	}
	float *cube = new float[3*3*12*blockNum];
	int counter = 0;
	for(int i = 0; i < DIM; i++){
		for(int j = 0; j < DIM; j++){
			int curV = m.getValue(j,i);
			if(curV == 1){
				int tmp[] = {
					i,0,j,   i,1,j,   i,1,j+1,
					i,0,j,   i,0,j+1,   i,1,j+1,
					i+1,0,j,   i+1,1,j,   i+1,1,j+1,
					i+1,0,j,   i+1,0,j+1,   i+1,1,j+1,

					i+1,0,j,  i,0,j, i+1,1,j,
					i,1,j,  i,0,j, i+1,1,j,
					i+1,0,j+1,  i,0,j+1, i+1,1,j+1,
					i,1,j+1,  i,0,j+1, i+1,1,j+1,
					
					i,0,j,  i,0,j+1, i+1,0,j,
					i,0,j+1, i+1,0,j+1,i+1,0,j,
					i,1,j,  i,1,j+1, i+1,1,j,
					i,1,j+1, i+1,1,j+1,i+1,1,j,

				};
				for(int i = 0; i < 3 * 3 * 12;i++){
					cube[counter + i] = tmp[i];
				}
				counter += 3*3*12;
			}
		}
	}


	//cube
	glGenVertexArrays( 1, &cube_vao);
	glBindVertexArray( cube_vao);
	glGenBuffers( 1, &cube_vbo );
	glBindBuffer( GL_ARRAY_BUFFER, cube_vbo);
	glBufferData( GL_ARRAY_BUFFER, 3*3*12*blockNum *sizeof(float),
		cube, GL_STATIC_DRAW );



	GLint posAttrib = m_shader.getAttribLocation( "position" );
	glEnableVertexAttribArray( posAttrib );
	glVertexAttribPointer( posAttrib, 3, GL_FLOAT, GL_FALSE, 0, nullptr );

	glBindVertexArray( 0 );
	glBindBuffer( GL_ARRAY_BUFFER, 0 );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );

	//delete cube[];



}

void A1::reset(){
	colour[0] = 0.0f;
	colour[1] = 0.0f;
	colour[2] = 0.0f;
	mazeC[0] = 1;
	mazeC[1] = 0;
	mazeC[2] = 0;
	sphereC[0] = 0;
	sphereC[1] = 1;
	sphereC[2] = 0;
	floorC[0] = 0;
	floorC[1] = 0;
	floorC[2] = 0;

	
	
	Maze m = Maze(16);
	int blockNum = 0;


	int sx = 0;
	int sz = 0;

	int wallS = 1;

	fov = 1;
	radius = 0.5f;  
	camX = 0; 
	camZ = 1; 
	camY = 1;   
	lastX = 0; 
	lastY = 0;
	//xoffset;    
	camRotDistance = 0.0f;    
	isFirstMouse = GL_TRUE;    
	isLeftMousePress = GL_FALSE;


	resistanceD = 0;
	deltaTime = 0;
	lastFrame = 0;
	curFrame = 0;
	speed = 0;


	resis = false;
	firsttime = true;

	digged = false;


	// Set up the uniforms
	P_uni = m_shader.getUniformLocation( "P" );
	V_uni = m_shader.getUniformLocation( "V" );
	M_uni = m_shader.getUniformLocation( "M" );
	col_uni = m_shader.getUniformLocation( "colour" );

	initGrid();

	// Set up initial view and projection matrices (need to do this here,
	// since it depends on the GLFW window being set up correctly).
	view = glm::lookAt( 
		glm::vec3( 0.0f, 2.*float(DIM)*2.0*M_SQRT1_2, float(DIM)*2.0*M_SQRT1_2 ),
		glm::vec3( 0.0f, 0.0f, 0.0f ),
		glm::vec3( 0.0f, 1.0f, 0.0f ) );

	proj = glm::perspective( 
		glm::radians( 30.0f ),
		float( m_framebufferWidth ) / float( m_framebufferHeight ),
		1.0f, 1000.0f );
	
}

