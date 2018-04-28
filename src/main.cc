#include <GL/glew.h>
#include <dirent.h>

// #include "bone_geometry.h"
// #include "procedure_geometry.h"
#include "render_pass.h"
#include "config.h"
#include "gui.h"
#include "tictoc.h"

#include "cloth.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>



#include <glm/gtx/component_wise.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/io.hpp>
#include <debuggl.h>

int window_width = 800, window_height = 600;
const std::string window_title = "Skinning";

const char* vertex_shader =
#include "shaders/default.vert"
;

const char* geometry_shader =
#include "shaders/default.geom"
;

const char* fragment_shader =
#include "shaders/default.frag"
;


const char* cloth_vertex_shader =
#include "shaders/cloth.vert"
;

const char* cloth_fragment_shader =
#include "shaders/cloth.frag"
;

const char* spring_vertex_shader =
#include "shaders/spring.vert"
;

const char* spring_fragment_shader =
#include "shaders/spring.frag"
;

const char* bend_spring_vertex_shader =
#include "shaders/bend_spring.vert"
;

const char* bend_spring_fragment_shader =
#include "shaders/bend_spring.frag"
;
// FIXME: Add more shaders here.

void ErrorCallback(int error, const char* description) {
	std::cerr << "GLFW Error: " << description << "\n";
}

GLFWwindow* init_glefw()
{
	if (!glfwInit())
		exit(EXIT_FAILURE);
	glfwSetErrorCallback(ErrorCallback);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_SAMPLES, 4);
	auto ret = glfwCreateWindow(window_width, window_height, window_title.data(), nullptr, nullptr);
	CHECK_SUCCESS(ret != nullptr);
	glfwMakeContextCurrent(ret);
	glewExperimental = GL_TRUE;
	CHECK_SUCCESS(glewInit() == GLEW_OK);
	glGetError();  // clear GLEW's error for it
	glfwSwapInterval(1);
	const GLubyte* renderer = glGetString(GL_RENDERER);  // get renderer string
	const GLubyte* version = glGetString(GL_VERSION);    // version as a string
	std::cout << "Renderer: " << renderer << "\n";
	std::cout << "OpenGL version supported:" << version << "\n";

	return ret;
}

int main(int argc, char* argv[])
{

	GLFWwindow *window = init_glefw();
	GUI gui(window);


	

	int cloth_x_size = 11;
	int cloth_z_size = 11;

	Cloth cloth(cloth_x_size, cloth_z_size);
	TicTocTimer *timer = new TicTocTimer;
	*timer = tic();


	glm::vec4 light_position = glm::vec4(0.0f, 100.0f, 0.0f, 1.0f);
	MatrixPointers mats; // Define MatrixPointers here for lambda to capture
	/*
	 * In the following we are going to define several lambda functions to bind Uniforms.
	 *
	 * Introduction about lambda functions:
	 *      http://en.cppreference.com/w/cpp/language/lambda
	 *      http://www.stroustrup.com/C++11FAQ.html#lambda
	 */
	/*
	 * The following lambda functions are defined to bind uniforms
	 */
	auto matrix_binder = [](int loc, const void* data) {
		glUniformMatrix4fv(loc, 1, GL_FALSE, (const GLfloat*)data);
	};
	auto vector_binder = [](int loc, const void* data) {
		glUniform4fv(loc, 1, (const GLfloat*)data);
	};
	auto vector3_binder = [](int loc, const void* data) {
		glUniform3fv(loc, 1, (const GLfloat*)data);
	};
	auto float_binder = [](int loc, const void* data) {
		glUniform1fv(loc, 1, (const GLfloat*)data);
	};
	auto int_binder = [](int loc, const void* data) {
		glUniform1iv(loc, 1, (const GLint*)data);
	};
	

	/*
	 * The lambda functions below are used to retrieve data
	 */
	auto std_model_data = [&mats]() -> const void* {
		return mats.model;
	};
	auto std_view_data = [&mats]() -> const void* {
		return mats.view;
	};
	auto std_camera_data  = [&gui]() -> const void* {
		return &gui.getCamera()[0];
	};
	auto std_proj_data = [&mats]() -> const void* {
		return mats.projection;
	};
	auto std_light_data = [&light_position]() -> const void* {
		return &light_position[0];
	};
	auto alpha_data  = [&gui]() -> const void* {
		static const float transparet = 0.5; // Alpha constant goes here
		static const float non_transparet = 1.0;
		if (gui.isTransparent())
			return &transparet;
		else
			return &non_transparet;
	};
	glm::mat4 identity_model_mat(1.0);
	auto identity_model_data = [&identity_model_mat]() -> const void* {
		return &identity_model_mat[0][0];
	};

	// FIXME: add more lambdas for data_source if you want to use RenderPass.
	//        Otherwise, do whatever you like here


	
	ShaderUniform std_model = { "model", matrix_binder, std_model_data };
	ShaderUniform std_view = { "view", matrix_binder, std_view_data };
	ShaderUniform std_camera = { "camera_position", vector3_binder, std_camera_data };
	ShaderUniform std_proj = { "projection", matrix_binder, std_proj_data };
	ShaderUniform std_light = { "light_position", vector_binder, std_light_data };
	ShaderUniform object_alpha = { "alpha", float_binder, alpha_data };
	ShaderUniform identity_model = {"model", matrix_binder, identity_model_data };

	// FIXME: define more ShaderUniforms for RenderPass if you want to use it.
	//        Otherwise, do whatever you like here


	
	// Cloth render pass
	RenderDataInput cloth_pass_input;
	cloth_pass_input.assign(0, "vertex_position", cloth.vertices.data(), cloth.vertices.size(), 3, GL_FLOAT);
	cloth_pass_input.assign(1, "uv", cloth.cloth_uv_coords.data(), cloth.cloth_uv_coords.size(), 2, GL_FLOAT);
	
	RenderPass cloth_pass(-1,
			cloth_pass_input,
			{ cloth_vertex_shader, nullptr, cloth_fragment_shader },
			{ std_model, std_view, std_proj, std_light },
			{ "fragment_color" }
			);

	// structural springs render pass 
	RenderDataInput struct_spring_pass_input;
	struct_spring_pass_input.assign(0, "vertex_position", cloth.struct_spring_vertices.data(), cloth.struct_spring_vertices.size(), 3, GL_FLOAT);
	
	RenderPass struct_spring_pass(-1,
			struct_spring_pass_input,
			{ spring_vertex_shader, nullptr, spring_fragment_shader },
			{ std_model, std_view, std_proj, std_light },
			{ "fragment_color" }
			);

	// bending springs render pass
	RenderDataInput bend_spring_pass_input;
	bend_spring_pass_input.assign(0, "vertex_position", cloth.bend_spring_vertices.data(), cloth.bend_spring_vertices.size(), 3, GL_FLOAT);
	
	RenderPass bend_spring_pass(-1,
			bend_spring_pass_input,
			{ bend_spring_vertex_shader, nullptr, bend_spring_fragment_shader },
			{ std_model, std_view, std_proj, std_light },
			{ "fragment_color" }
			);

	
	bool draw_cloth = true;
	bool draw_struct_spring = true;
	bool draw_bend_spring = true;

	toc(timer);
	while (!glfwWindowShouldClose(window)) {
		// Setup some basic window stuff.
		glfwGetFramebufferSize(window, &window_width, &window_height);
		glViewport(0, 0, window_width, window_height);
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_MULTISAMPLE);
		glEnable(GL_BLEND);
		glEnable(GL_CULL_FACE);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDepthFunc(GL_LESS);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glCullFace(GL_BACK);

		gui.updateMatrices();
		mats = gui.getMatrixPointers();

		if (gui.isPoseDirty()) {
			// update animation
			gui.clearPose();
		}

		if (gui.toResetSystem()) {
			// 
			gui.clearResetFlag();
		}

		if (gui.toRandomDisturb()) {
			//
			gui.clearDisturbFlag();
		}

		if (gui.toTear()) {
			cloth.randomTear();
			gui.clearTearFlag();
		}

		
		float delta_t = (float) toc(timer) * gui.getTimeSpeed();
		delta_t *= 1.0;
		// std::cout << "delta_t = " << delta_t << std::endl;
		cloth.animate(delta_t);
		


		if (draw_cloth) {
			glDisable(GL_CULL_FACE);
			cloth_pass.updateVBO(0, cloth.vertices.data(), cloth.vertices.size());
			cloth_pass.updateVBO(1, cloth.cloth_uv_coords.data(), cloth.cloth_uv_coords.size());
			cloth_pass.setup();

			CHECK_GL_ERROR(glDrawArrays(GL_TRIANGLES,
										0,
		                              	cloth.vertices.size()));
		}

		if (draw_struct_spring) {
			struct_spring_pass.updateVBO(0, cloth.struct_spring_vertices.data(), cloth.struct_spring_vertices.size());
			struct_spring_pass.setup();

			CHECK_GL_ERROR(glDrawArrays(GL_LINES,
										0,
		                              	cloth.struct_spring_vertices.size()));
		}

		if (draw_bend_spring) {
			bend_spring_pass.updateVBO(0, cloth.bend_spring_vertices.data(), cloth.bend_spring_vertices.size());
			bend_spring_pass.setup();

			CHECK_GL_ERROR(glDrawArrays(GL_LINES,
										0,
		                              	cloth.bend_spring_vertices.size()));
		}



		// Poll and swap.
		glfwPollEvents();
		glfwSwapBuffers(window);
	}
	glfwDestroyWindow(window);
	glfwTerminate();
	exit(EXIT_SUCCESS);
}
