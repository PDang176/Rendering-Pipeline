#include "driver_state.h"
#include <cstring>
#include <algorithm>
#include <cfloat>

// Helper Functions
float area(vec2 a, vec2 b, vec2 c);
int get_index(int i, int j, int width);

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
    // Set image width and height
    state.image_width=width;
    state.image_height=height;

    // Allocate space for all pixels in image_color and image_depth
    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];
    // Initialize all the pixel's color to black and depth to max depth
    for(int i = 0; i < width * height; i++){
        state.image_color[i] = make_pixel(0, 0, 0);
        state.image_depth[i] = FLT_MAX;
    }

    // std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
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
    // Allocate an array of data_geometry objects (one for each vertex)
    data_geometry * dg_arr = new data_geometry[state.num_vertices];

    // Allocate each data_geometry object's data array
    for(int i = 0; i < state.num_vertices; i++){
        dg_arr[i].data = new float[MAX_FLOATS_PER_VERTEX];
    }

    // Call vertex_shader for each vertex
    data_vertex dv;
    for(int i = 0; i < state.num_vertices; i++){
        // Create a data_vertex passing in the location of the first float in vertex_data for that data_vertex
        dv.data = &state.vertex_data[i * state.floats_per_vertex];
        state.vertex_shader(dv, dg_arr[i], state.uniform_data);
    }

    switch(type){
        case render_type::triangle:
            for(int i = 0; i < state.num_vertices; i++){
                dg_arr[i].data = &state.vertex_data[i * state.floats_per_vertex];
            }
            for(int i = 0; i < state.num_vertices / 3; i++){
                rasterize_triangle(state, dg_arr[3 * i], dg_arr[3 * i + 1], dg_arr[3 * i + 2]);
            }
            break;

        case render_type::indexed:
            break;

        case render_type::fan:
            break;

        case render_type::strip:
            break;

        default:
            std::cerr << "Error: Invalid render_type" << std::endl;
    }

    // Deallocate data_geometry array
    delete[] dg_arr;

    // std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state, v0, v1, v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    int width = state.image_width;
    int height = state.image_height;

    // Converting NDC to pixel space coordinates
    float x[3];
    float y[3];
    float z[3];

    x[0] = ((width / 2.0) * v0.gl_Position[0] / v0.gl_Position[3]) + (width / 2.0) - 0.5;
    y[0] = ((height / 2.0) * v0.gl_Position[1] / v0.gl_Position[3]) + (height / 2.0) - 0.5;
    z[0] = v0.gl_Position[2] / v0.gl_Position[3];

    x[1] = ((width / 2.0) * v1.gl_Position[0] / v1.gl_Position[3]) + (width / 2.0) - 0.5;
    y[1] = ((height / 2.0) * v1.gl_Position[1] / v1.gl_Position[3]) + (height / 2.0) - 0.5;
    z[1] = v1.gl_Position[2] / v1.gl_Position[3];

    x[2] = ((width / 2.0) * v2.gl_Position[0] / v2.gl_Position[3]) + (width / 2.0) - 0.5;
    y[2] = ((height / 2.0) * v2.gl_Position[1] / v2.gl_Position[3]) + (height / 2.0) - 0.5;
    z[2] = v2.gl_Position[2] / v2.gl_Position[3];

    // Calculate bounding box of the triangle
    int min_x = std::min(std::min(x[0], x[1]), x[2]);
    int max_x = std::max(std::max(x[0], x[1]), x[2]);
    int min_y = std::min(std::min(y[0], y[1]), y[2]);
    int max_y = std::max(std::max(y[0], y[1]), y[2]);

    // Calculate area of triangle
    vec2 p0{x[0], y[0]};
    vec2 p1{x[1], y[1]};
    vec2 p2{x[2], y[2]};
    float ABC = area(p0, p1, p2);

    for(int i = min_x; i<= max_x; i++){
        for(int j = min_y; j <= max_y; j++){
            float alpha, beta, gamma;

            // Current point
            vec2 p(i, j);
            
            // Calculating barycentric coordinates for the current point
            alpha = area(p, p1, p2) / ABC;
            beta = area(p, p2, p0) / ABC;
            gamma = area(p, p0, p1) / ABC;

            // Inside triangle
            if(alpha >= 0 && beta >= 0 && gamma >= 0){
                // Calculating the depth for the current point
                float depth = (alpha * z[0]) + (beta * z[1]) + (gamma * z[2]);

                // Check if we found the closest point
                int index = get_index(i, j, width);
                if(depth < state.image_depth[index]){
                    // Set new depth
                    state.image_depth[index] = depth;

                    // Fragment Shading
                    float* data = new float[MAX_FLOATS_PER_VERTEX];
                    data_fragment df;
                    data_output output;

                    // Loop through all floats in the current point
                    for(int k = 0; k < state.floats_per_vertex; k++){
                        switch(state.interp_rules[k]){
                            case interp_type::flat:
                                data[k] = v0.data[k];
                                break;
                            case interp_type::smooth:
                                break;
                            case interp_type::noperspective:
                                data[k] = (alpha * v0.data[k]) + (beta * v1.data[k]) + (gamma * v2.data[k]);
                                break;
                            default:
                                std::cerr << "Error: Invalid interp_rule" << std::endl;
                        }
                    }
                    df.data = data;
                    state.fragment_shader(df, output, state.uniform_data);

                    // Set pixel to final color
                    state.image_color[index] = make_pixel(output.output_color[0] * 255, output.output_color[1] * 255, output.output_color[2] * 255);

                    // Deallocate data array
                    delete[] data;
                }
            }
        }
    }

    // std::cout<<"TODO: implement rasterization"<<std::endl;
}

float area(vec2 a, vec2 b, vec2 c){
    return 0.5 * (((b[0] * c[1]) - (c[0] * b[1])) + ((c[0] * a[1]) - (a[0] * c[1])) + ((a[0] * b[1]) - (b[0] * a[1])));
}

int get_index(int i, int j, int width){
    return (j * width) + i;
}