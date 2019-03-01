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
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];
    std::fill(state.image_color, state.image_color + (width * height), make_pixel(0,0,0));
    std::fill(state.image_depth, state.image_depth + (width * height), 1);
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
    if (type == render_type::triangle || type == render_type::indexed) {
        int numVal = 0;
        if (type == render_type::triangle) {
  			    numVal = state.num_vertices;
        } else if (type == render_type::indexed) {
            numVal = 3 * state.num_triangles;
        }
        for (int i = 0; i < numVal; i = i + 3) {
      		const data_geometry* pass[3];
      		for (int j = 0; j < 3; j++) {
            int ind = 0;
            if (type == render_type::triangle) {
      			    ind = ((i + j) * state.floats_per_vertex);
            } else if (type == render_type::indexed) {
                ind = (state.index_data[(i + j)] * state.floats_per_vertex);
            }
      			data_geometry* dg = new data_geometry();
      			dg->data = state.vertex_data + ind;
      			data_vertex dv;
      			dv.data = dg->data;
      			state.vertex_shader(dv, *dg, state.uniform_data);
      			pass[j] = dg;
      		}
      		clip_triangle(state, pass, 0);
          for (int j = 0; j < 3; j++) 
              delete pass[j]; 
      	}
    } else if (type == render_type::fan) {
        for (int i = 0; i < state.num_vertices; i++) {
  		    const data_geometry* pass[3];
      		for (int j = 0; j < 3; j++) {
            int ind = ((i + j) * state.floats_per_vertex);
            if (j == 0) {
                ind = 0;
            }
      			data_geometry* dg = new data_geometry();
      			dg->data = state.vertex_data + ind;
      			data_vertex dv;
      			dv.data = dg->data;
      			state.vertex_shader(dv, *dg, state.uniform_data);
      			pass[j] = dg;
      		}
      		clip_triangle(state, pass, 0);
          for (int j = 0; j < 3; j++) 
              delete pass[j]; 
      	}
    } else if (type == render_type::strip) {
        for (int i = 0; i < (state.num_vertices - 2); i++) {
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
      		clip_triangle(state, pass, 0);
          for (int j = 0; j < 3; j++) 
              delete pass[j]; 
      	}
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
    // std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    
    //Depending on the face, set which axis and what value to check against
    float cnt = 0, ind = 0, val = 1;
    ivec3 inVec = {0, 0, 0};
    if (face == 0) { 		// x = 1
		ind = 0;
		val = 1;
	} else if (face == 1) { // x = -1
		ind = 0;
		val = -1;
	} else if (face == 2) { // y = 1
		ind = 1;
		val = 1;
	} else if (face == 3) { // y = -1
		ind = 1;
		val = -1;
	} else if (face == 4) { // z = 1
		ind = 2;
		val = 1;
	} else if (face == 5) { // z = -1
		ind = 2;
		val = -1;
	}
 
	//If vertex is inside, set it's inVec value to 1 and increment a counter (inVec stores which vertex is inside, and the counter tells how many vertices are inside)
	for (int i = 0; i < 3; i++) {
		if (val < 0 && in[i]->gl_Position[ind] >= (val * in[i]->gl_Position[3])) {
			cnt++;
			inVec[i] = 1;
		} else if (val > 0 && in[i]->gl_Position[ind] <= (val * in[i]->gl_Position[3])) {
			cnt++;
			inVec[i] = 1;
		}
	}
 
	//Handle simple cases of all in or all out
	if (cnt == 3) {
		clip_triangle(state,in,face+1); //If all vertices are inside, pass triangle unchanged to next face
    return;
	} else if (cnt == 0) {
		return; //If all vertices are outside, don't rasterize it
	}
 
	//If one vertex is inside, set A to that vertex and if two vertices are inside, set A to the vertex that is outside
	int chk = 0;
	if (cnt == 1) {
		chk = 1;
	} else if (cnt == 2) {
		chk = 0;
	}
	const data_geometry* A = 0;
	const data_geometry* B = 0;
	const data_geometry* C = 0;


	for (int i = 0; i < 3; i++) {
		if (inVec[i] == chk) {
			A = in[i];
      B = in[(i + 1) % 3];
      C = in[(i + 2) % 3];
      break;
		}
	}

	//Calculate barycentric values to calculate ab and ac
	float wA = A->gl_Position[3], wB = B->gl_Position[3], wC = C->gl_Position[3];
  float aV = A->gl_Position[ind], bV = B->gl_Position[ind], cV = C->gl_Position[ind];
	float alphaB = ((val * wB) - bV) / ((aV - (val * wA)) + ((val * wB) - bV));
	float alphaC = ((val * wC) - cV) / ((aV - (val * wA)) + ((val * wC) - cV));
 
  //Create and populate new vertices
	data_geometry* AB = new data_geometry();
	data_geometry* AC = new data_geometry();
	AB->data = new float[state.floats_per_vertex];
  AC->data = new float[state.floats_per_vertex];
  for (int i = 0; i < state.floats_per_vertex; i++) {
      AB->data[i] = (alphaB * A->data[i]) + ((1 - alphaB) * B->data[i]);
      AC->data[i] = (alphaC * A->data[i]) + ((1 - alphaC) * C->data[i]);
  }
  for (int i = 0; i < 4; i++) {
      AB->gl_Position[i] = (alphaB * A->gl_Position[i]) + ((1 - alphaB) * B->gl_Position[i]);
      AC->gl_Position[i] = (alphaC * A->gl_Position[i]) + ((1 - alphaC) * C->gl_Position[i]);
  }

  //Call clip_triangle
  const data_geometry* pass[3];
	if (cnt == 1) {	
   			pass[0] = A;
        pass[1] = AB;
        pass[2] = AC;
      clip_triangle(state,pass,face+1);
      delete AB;
      delete AC;
      return;
	} else if (cnt == 2) {
    			pass[0] = AB;
          pass[1] = B;
          pass[2] = C;
          clip_triangle(state,pass,face+1);
          pass[0] = C;
          pass[1] = AC;
          pass[2] = AB;
          clip_triangle(state,pass,face+1);
      delete AB;
      delete AC;
      return;
	}
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    // std::cout<<"TODO: implement rasterization"<<std::endl;
    vec3 A, B, C;
    A[0] = (((in[0]->gl_Position[0] / in[0]->gl_Position[3]) + 1) * (state.image_width / 2));
    A[1] = (((in[0]->gl_Position[1] / in[0]->gl_Position[3]) + 1) * (state.image_height / 2));
    A[2] =    in[0]->gl_Position[2] / in[0]->gl_Position[3];
    B[0] = (((in[1]->gl_Position[0] / in[1]->gl_Position[3]) + 1) * (state.image_width / 2));
    B[1] = (((in[1]->gl_Position[1] / in[1]->gl_Position[3]) + 1) * (state.image_height / 2));
    B[2] =    in[1]->gl_Position[2] / in[1]->gl_Position[3];
    C[0] = (((in[2]->gl_Position[0] / in[2]->gl_Position[3]) + 1) * (state.image_width / 2));
    C[1] = (((in[2]->gl_Position[1] / in[2]->gl_Position[3]) + 1) * (state.image_height / 2));
    C[2] =    in[2]->gl_Position[2] / in[2]->gl_Position[3];
	float area = 0.5 * ((A[0] * (B[1] - C[1])) + (B[0] * (C[1] - A[1])) + (C[0] * (A[1] - B[1])));                        
  int startX = std::min(std::min(A[0], B[0]), C[0]);
  int endX = std::max(std::max(A[0], B[0]), C[0]);
  int startY = std::min(std::min(A[1], B[1]), C[1]);
  int endY = std::max(std::max(A[1], B[1]), C[1]);
  for (int i = startX; i <= endX; i++) {
      for (int j = startY; j <= endY; j++) {
            //float Area = 0.5 * ((i * (B[1] - C[1])) + (B[0] * (C[1] - j)) + (C[0] * (j - B[1])));
          	float Brea = 0.5 * ((A[0] * (j - C[1])) + (i * (C[1] - A[1])) + (C[0] * (A[1] - j)));
          	float Grea = 0.5 * ((A[0] * (B[1] - j)) + (B[0] * (j - A[1])) + (i * (A[1] - B[1])));
           //double alpha = Area / area;
           double beta = Brea / area;
           double gamma = Grea / area;
           double alpha = (1 - beta) - gamma;
           if (alpha >= -0.001 && beta >= -0.001 && gamma >= -0.001 && (alpha + beta + gamma) <= 1.0001) {
               int image_index = (j * state.image_width) + i;
              float iPol[state.floats_per_vertex];
              for (int k = 0; k < state.floats_per_vertex; k++) {
                 if (state.interp_rules[k] == interp_type::flat) {
                     iPol[k] = in[0]->data[k];
                 } else if (state.interp_rules[k] == interp_type::smooth) {
                     float aP = 1 / in[0]->gl_Position[3], bP = 1 / in[1]->gl_Position[3], gP = 1 / in[2]->gl_Position[3];
                     aP = alpha * aP;
                     bP = beta * bP;
                     gP = gamma * gP;
                     iPol[k] = ((aP * in[0]->data[k]) + (bP * in[1]->data[k]) + (gP * in[2]->data[k])) / (aP + bP + gP);
                 } else if (state.interp_rules[k] == interp_type::noperspective) {
                     iPol[k] = (alpha * in[0]->data[k]) + (beta * in[1]->data[k]) + (gamma * in[2]->data[k]);
                 }
             }
             float zVal = (alpha * A[2]) + (beta * B[2]) + (gamma * C[2]);
             if (zVal <= state.image_depth[image_index]) {
                 data_fragment df; 
                 df.data = iPol;
                 data_output dO;
                 state.fragment_shader(df, dO, state.uniform_data);
                 state.image_color[image_index] = make_pixel(255 * dO.output_color[0], 255 * dO.output_color[1], 255 * dO.output_color[2]);
                 state.image_depth[image_index] = zVal;
             }
           }
      }
  }
}

