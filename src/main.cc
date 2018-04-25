#include <GL/glew.h>
#include <dirent.h>

// #include "bone_geometry.h"
#include "procedure_geometry.h"
#include "render_pass.h"
#include "config.h"
#include "gui.h"
#include "mass_spring.h"
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

const char* floor_fragment_shader =
#include "shaders/floor.frag"
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

	std::vector<glm::vec4> floor_vertices;
	std::vector<glm::uvec3> floor_faces;
	create_floor(floor_vertices, floor_faces);

	int cloth_x_size = 3;
	int cloth_z_size = 3;
	MassSpringSystem ms_system(cloth_x_size, cloth_z_size);
	TicTocTimer *timer = new TicTocTimer;
	*timer = tic();


	Cloth cloth(cloth_x_size, cloth_z_size);


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
	}; // This returns point to model matrix
	glm::mat4 floor_model_matrix = glm::mat4(1.0f);
	auto floor_model_data = [&floor_model_matrix]() -> const void* {
		return &floor_model_matrix[0][0];
	}; // This return model matrix for the floor.
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
	ShaderUniform floor_model = { "model", matrix_binder, floor_model_data };
	ShaderUniform std_view = { "view", matrix_binder, std_view_data };
	ShaderUniform std_camera = { "camera_position", vector3_binder, std_camera_data };
	ShaderUniform std_proj = { "projection", matrix_binder, std_proj_data };
	ShaderUniform std_light = { "light_position", vector_binder, std_light_data };
	ShaderUniform object_alpha = { "alpha", float_binder, alpha_data };
	ShaderUniform identity_model = {"model", matrix_binder, identity_model_data };

	// Floor render pass
	RenderDataInput floor_pass_input;
	floor_pass_input.assign(0, "vertex_position", floor_vertices.data(), floor_vertices.size(), 4, GL_FLOAT);
	floor_pass_input.assignIndex(floor_faces.data(), floor_faces.size(), 3);
	RenderPass floor_pass(-1,
			floor_pass_input,
			{ vertex_shader, geometry_shader, floor_fragment_shader},
			{ floor_model, std_view, std_proj, std_light },
			{ "fragment_color" }
			);
	// FIXME: define more ShaderUniforms for RenderPass if you want to use it.
	//        Otherwise, do whatever you like here


	// for(int i = 0; i < ms_system.node_positions.size(); i++) {
	// 	std::cout << glm::to_string(ms_system.node_positions[i]) << ", ";
	// }
	// std::cout << std::endl;

	// for(int i = 0; i < ms_system.line_indices.size(); i++) {
	// 	std::cout << glm::to_string(ms_system.line_indices[i]) << ", ";
	// }

	// Cloth render pass
	RenderDataInput cloth_pass_input;
	cloth_pass_input.assign(0, "vertex_position", ms_system.node_positions.data(), ms_system.node_positions.size(), 3, GL_FLOAT);
	cloth_pass_input.assignIndex(ms_system.line_indices.data(), ms_system.line_indices.size(), 2);
	
	RenderPass cloth_pass(-1,
			cloth_pass_input,
			{ cloth_vertex_shader, nullptr, cloth_fragment_shader },
			{ std_model, std_view, std_proj, std_light },
			{ "fragment_color" }
			);


	RenderDataInput tri_cloth_pass_input;
	tri_cloth_pass_input.assign(0, "vertex_position", cloth.vertices.data(), cloth.vertices.size(), 3, GL_FLOAT);
	tri_cloth_pass_input.assign(1, "uv", cloth.cloth_uv_coords.data(), cloth.cloth_uv_coords.size(), 2, GL_FLOAT);
	
	RenderPass tri_cloth_pass(-1,
			tri_cloth_pass_input,
			{ cloth_vertex_shader, nullptr, cloth_fragment_shader },
			{ std_model, std_view, std_proj, std_light },
			{ "fragment_color" }
			);


	RenderDataInput spring_pass_input;
	spring_pass_input.assign(0, "vertex_position", cloth.spring_vertices.data(), cloth.spring_vertices.size(), 3, GL_FLOAT);
	
	RenderPass spring_pass(-1,
			spring_pass_input,
			{ spring_vertex_shader, nullptr, spring_fragment_shader },
			{ std_model, std_view, std_proj, std_light },
			{ "fragment_color" }
			);

	bool draw_floor = false;
	bool draw_cloth = false;
	bool draw_tri_cloth = true;
	bool draw_spring = true;
	
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
			ms_system.resetSystem();
			gui.clearResetFlag();
		}

		if (gui.toRandomDisturb()) {
			ms_system.randomDisturb();
			gui.clearDisturbFlag();
		}

		
		float delta_t = (float) toc(timer) * gui.getTimeSpeed();

		cloth.animate(delta_t);
		
		// ms_system.animate(delta_t);

		// std::cout << "delta t: " << delta_t << std::endl;
		// std::cout << "force: " << glm::to_string(ms_system.nodes_[38].force) << std::endl;
		// std::cout << "velocity: " << glm::to_string(ms_system.nodes_[38].velocity) << std::endl;
		// std::cout << "position: " << glm::to_string(ms_system.nodes_[38].position) << std::endl;
		// std::cout << std::endl;



		// Then draw floor.
		if (draw_floor) {
			floor_pass.setup();
			// Draw our triangles.
			CHECK_GL_ERROR(glDrawElements(GL_TRIANGLES,
			                              floor_faces.size() * 3,
			                              GL_UNSIGNED_INT, 0));
		}

		if (draw_cloth) {
			cloth_pass.updateVBO(0, ms_system.node_positions.data(), ms_system.node_positions.size());
			cloth_pass_input.assignIndex(ms_system.line_indices.data(), ms_system.line_indices.size(), 2);
			cloth_pass.updateIndexBuffer(ms_system.line_indices.data(), ms_system.line_indices.size());
			cloth_pass.setup();
			// Draw our triangles.
			CHECK_GL_ERROR(glDrawElements(GL_LINES,
			                              ms_system.line_indices.size() * 2,
			                              GL_UNSIGNED_INT, 0));
		}

		if (draw_tri_cloth) {
			glDisable(GL_CULL_FACE);
			tri_cloth_pass.updateVBO(0, cloth.vertices.data(), cloth.vertices.size());
			tri_cloth_pass.updateVBO(1, cloth.cloth_uv_coords.data(), cloth.cloth_uv_coords.size());
			tri_cloth_pass.setup();

			CHECK_GL_ERROR(glDrawArrays(GL_TRIANGLES,
										0,
		                              	cloth.vertices.size()));
		}

		if (draw_spring) {
			spring_pass.updateVBO(0, cloth.spring_vertices.data(), cloth.spring_vertices.size());
			spring_pass.setup();

			CHECK_GL_ERROR(glDrawArrays(GL_LINES,
										0,
		                              	cloth.spring_vertices.size()));
		}



		// Poll and swap.
		glfwPollEvents();
		glfwSwapBuffers(window);
	}
	glfwDestroyWindow(window);
	glfwTerminate();
	exit(EXIT_SUCCESS);
}
