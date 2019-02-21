#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];
    std::fill(state.image_color, state.image_color + (width * height), make_pixel(0,0,0));
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    // std::cout<<"TODO: implement rendering."<<std::endl;
    for (int i = 0; i < state.num_vertices; i = i + 3) {
		const data_geometry* pass[3];
		for (int j = 0; j < 3; j++) {
			int ind = ((i + j) * state.floats_per_vertex);
			data_geometry* dg = new data_geometry();
      			dg->data = state.vertex_data + ind;
			data_vertex dv;
			dv.data = dg->data;
			state.vertex_shader(dv, *dg, state.uniform_data);
			pass[j] = dg;
		}
		rasterize_triangle(state, pass);
	}
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    // std::cout<<"TODO: implement rasterization"<<std::endl;
    vec2 A, B, C;
    A[0] = (in[0]->gl_Position[0] / in[0]->gl_Position[3]);
    for(int i = 0; i < 3; i++) {
  		for (int j = 0; j < 3; j++) {
			inn[i][j] = in[i]->gl_Position[j] / in[i]->gl_Position[3];
  			//in[i]->gl_Position[j] = in[i]->gl_Position[j] / in[i]->gl_Position[3];
  		}
  		inn[i][0] = (inn[i][0] + 1) * (state.image_width / 2);
  		inn[i][1] = (inn[i][1] + 1) * (state.image_height / 2);	
  		int image_index = inn[i][1] * state.image_width + inn[i][0];
  		state.image_color[image_index] = make_pixel(255, 255, 255);
		/*for (int j = 0; j < state.floats_per_vertex; j++) {
			std::cout << in[i]->data[j] << std::endl;
		}
		std::cout << "----------------------" << std::endl;  */	
  	}
	float area = 0.5 * ((inn[0][0] * (inn[1][1] - inn[2][1])) + 
			(inn[1][0] * (inn[2][1] - inn[0][1])) + 
			(inn[2][0] * (inn[0][1] - inn[1][1])));                        
  int startX = std::min(std::min(inn[0][0], inn[1][0]), inn[2][0]);
  int endX = std::max(std::max(inn[0][0], inn[1][0]), inn[2][0]);
  int startY = std::min(std::min(inn[0][1], inn[1][1]), inn[2][1]);
  int endY = std::max(std::max(inn[0][1], inn[1][1]), inn[2][1]);

  for (int i = startX; i < endX; i++) {
      for (int j = startY; j < endY; j++) {
          	/*float Area = 0.5 * ((i * (in[1]->gl_Position[1] - in[2]->gl_Position[1])) + 
          						          (in[1]->gl_Position[0] * (in[2]->gl_Position[1] - j)) + 
          						          (in[2]->gl_Position[0] * (j - in[1]->gl_Position[1])));*/
          	float Brea = 0.5 * ((inn[0][0] * (j - inn[2][1])) + 
          						          (i * (inn[2][1] - inn[0][1])) + 
          						          (inn[2][0] * (inn[0][1] - j)));
          	float Grea = 0.5 * ((inn[0][0] * (inn[1][1] - j)) + 
          						          (inn[1][0] * (j - inn[0][1])) + 
          						          (i * (inn[0][1] - inn[1][1])));
           //double alpha = Area / area;
           double beta = Brea / area;
           double gamma = Grea / area;
           double alpha = (1 - beta) - gamma;
           if (alpha >= 0 && beta >= 0 && gamma >= 0 && (alpha + beta + gamma) <= 1.001) {
               int image_index = (j * state.image_width) + i;
		state.image_color[image_index] = make_pixel(in[1]->data[state.floats_per_vertex - 3], in[1]->data[state.floats_per_vertex - 2], in[1]->data[state.floats_per_vertex - 1]);
		/*float r = (alpha * in[0]->data[4]) + (beta * in[1]->data[4]) + (gamma * in[2]->data[4]);
		float g = (alpha * in[0]->data[5]) + (beta * in[1]->data[5]) + (gamma * in[2]->data[5]);
		float b = (alpha * in[0]->data[6]) + (beta * in[1]->data[6]) + (gamma * in[2]->data[6]);
		state.image_color[image_index] = make_pixel(r, g, b);*/
           }
      }
  }
}

