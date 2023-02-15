#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
    for (int i = 0; i < sqrt(sample_rate); i++) {
          for (int j = 0; j < sqrt(sample_rate); j++) {
              int newy = y * sqrt(sample_rate) + i;
              int newx = x * sqrt(sample_rate) + j;
              int newW =  width * sqrt(sample_rate);
              sample_buffer[(newy * newW + newx)] = c;
          }
      }
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    // convert to the nearest whole pixel
    x0 = int(x0) * sqrt(sample_rate);
    y0 = int(y0) * sqrt(sample_rate);
    x1 = int(x1) * sqrt(sample_rate);
    y1 = int(y1) * sqrt(sample_rate);
    x2 = int(x2) * sqrt(sample_rate);
    y2 = int(y2) * sqrt(sample_rate);

    // coordinates for bounding box that we will loop over
    int xmax = ceil(max(max(x0, x1), x2));
    int xmin = floor(min(min(x0, x1), x2));
    int ymax = ceil(max(max(y0, y1), y2));
    int ymin = floor(min(min(y0, y1), y2));

    for (int x = xmin; x < xmax; x++) {
        for (int y = ymin; y < ymax; y++) {
            int dx0 = x1 - x0;
            int dy0 = y1 - y0;
            int dx1 = x2 - x1;
            int dy1 = y2 - y1;
            int dx2 = x0 - x2;
            int dy2 = y0 - y2;

            int lx0 = -((x + 0.5) - x0) * dy0 + ((y + 0.5) - y0) * dx0;
            int lx1 = -((x + 0.5) - x1) * dy1 + ((y + 0.5) - y1) * dx1;
            int lx2 = -((x + 0.5) - x2) * dy2 + ((y + 0.5) - y2) * dx2;
            if ((lx0 >= 0 and lx1 >= 0 and lx2 >= 0) or (lx0 <= 0 and lx1 <= 0 and lx2 <= 0)) {
                sample_buffer[(y * width * sqrt(sample_rate) + x)] = color;
            }
        }
    }
  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
      // convert to the nearest whole pixel
      x0 = int(x0) * sqrt(sample_rate);
      y0 = int(y0) * sqrt(sample_rate);
      x1 = int(x1) * sqrt(sample_rate);
      y1 = int(y1) * sqrt(sample_rate);
      x2 = int(x2) * sqrt(sample_rate);
      y2 = int(y2) * sqrt(sample_rate);

      // coordinates for bounding box that we will loop over
      int xmax = ceil(max(max(x0, x1), x2));
      int xmin = floor(min(min(x0, x1), x2));
      int ymax = ceil(max(max(y0, y1), y2));
      int ymin = floor(min(min(y0, y1), y2));

      for (int x = xmin; x < xmax; x++) {
          for (int y = ymin; y < ymax; y++) {
              int dx0 = x1 - x0;
              int dy0 = y1 - y0;
              int dx1 = x2 - x1;
              int dy1 = y2 - y1;
              int dx2 = x0 - x2;
              int dy2 = y0 - y2;

              int lx0 = -((x + 0.5) - x0)*dy0 + ((y + 0.5) - y0)*dx0;
              int lx1 = -((x + 0.5) - x1)*dy1 + ((y + 0.5) - y1)*dx1;
              int lx2 = -((x + 0.5) - x2)*dy2 + ((y + 0.5) - y2)*dx2;
              if ((lx0 >= 0 and lx1 >= 0 and lx2 >= 0) or (lx0 <= 0 and lx1 <= 0 and lx2 <= 0)) {
                  float alpha = (-(x-x1)*(y2-y1)+(y-y1)*(x2-x1)) / (-(x0-x1)*(y2-y1) + (y0-y1)*(x2-x1));
                  float beta = (-(x-x2)*(y0-y2)+(y-y2)*(x0-x2)) / (-(x1-x2)*(y0-y2) + (y1-y2)*(x0-x2));
                  float gamma = 1 - alpha - beta;
                  Color color = (alpha * c0) + (beta * c1) + (gamma * c2);
                  sample_buffer[(y * width * sqrt(sample_rate) + x)] = color;
              }
          }
      }
}


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // convert to the nearest whole pixel
    x0 = int(x0) * sqrt(sample_rate);
    y0 = int(y0) * sqrt(sample_rate);
    x1 = int(x1) * sqrt(sample_rate);
    y1 = int(y1) * sqrt(sample_rate);
    x2 = int(x2) * sqrt(sample_rate);
    y2 = int(y2) * sqrt(sample_rate);

    // coordinates for bounding box that we will loop over
    int xmax = ceil(max(max(x0, x1), x2));
    int xmin = floor(min(min(x0, x1), x2));
    int ymax = ceil(max(max(y0, y1), y2));
    int ymin = floor(min(min(y0, y1), y2));
    SampleParams sp;
    sp.lsm = this->lsm;
    sp.psm = this->psm;

    for (int x = xmin; x < xmax; x++) {
        for (int y = ymin; y < ymax; y++) {
            int dx0 = x1 - x0;
            int dy0 = y1 - y0;
            int dx1 = x2 - x1;
            int dy1 = y2 - y1;
            int dx2 = x0 - x2;
            int dy2 = y0 - y2;

            int lx0 = -((x + 0.5) - x0)*dy0 + ((y + 0.5) - y0)*dx0;
            int lx1 = -((x + 0.5) - x1)*dy1 + ((y + 0.5) - y1)*dx1;
            int lx2 = -((x + 0.5) - x2)*dy2 + ((y + 0.5) - y2)*dx2;
            if ((lx0 >= 0 and lx1 >= 0 and lx2 >= 0) or (lx0 <= 0 and lx1 <= 0 and lx2 <= 0)) {
                float alpha = (-(x-x1)*(y2-y1)+(y-y1)*(x2-x1)) / (-(x0-x1)*(y2-y1) + (y0-y1)*(x2-x1));
                float beta = (-(x-x2)*(y0-y2)+(y-y2)*(x0-x2)) / (-(x1-x2)*(y0-y2) + (y1-y2)*(x0-x2));
                float gamma = 1 - alpha - beta;
                Vector2D uv = {(alpha * u0) + (beta * u1) + (gamma * u2), (alpha * v0) + (beta * v1) + (gamma * v2)};
                sp.p_uv = uv;

                if (x + 1 < width) {
                    alpha = (-((x + 1) - x1) * (y2 - y1) + (y - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                    beta = (-((x + 1) - x2) * (y0 - y2) + (y - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                    gamma = 1 - alpha - beta;
                    Vector2D dx_uv = {(alpha * u0) + (beta * u1) + (gamma * u2), (alpha * v0) + (beta * v1) + (gamma * v2)};
                    sp.p_dx_uv = dx_uv;
                }

                if (y + 1 < height) {
                    alpha = (-(x - x1) * (y2 - y1) + ((y + 1) - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                    beta = (-(x - x2) * (y0 - y2) + ((y + 1) - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                    gamma = 1 - alpha - beta;
                    Vector2D dy_uv = {(alpha * u0) + (beta * u1) + (gamma * u2), (alpha * v0) + (beta * v1) + (gamma * v2)};
                    sp.p_dy_uv = dy_uv;
                }

                sample_buffer[(y * width * sqrt(sample_rate) + x)] = tex.sample(sp);

            }
        }
    }
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle




  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support


      for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Color col;
          for (int i = 0; i < sqrt(sample_rate); i++) {
              for (int j = 0; j < sqrt(sample_rate); j++) {
                  int newy = y* sqrt(sample_rate) + i;
                  int newx = x * sqrt(sample_rate) + j;
                  int newW =  width * sqrt(sample_rate);
                  col += sample_buffer[newy * newW + newx];
              }
          }
        col *= float(1.0/sample_rate);
        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }


      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
