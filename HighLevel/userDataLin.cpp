#include <opencvmex.hpp>

using namespace cv;

struct Time {
    float dt;
    int firstFrame;
    int lastFrame;
};

struct World {
    Vec<float, 3> points[];
    Vec<float, 3> segments[];
};

struct Robot {
    int id;
    const char* name;
    const char* type;
    const char* motion;
    Vec<float, 3> position;
    Vec<float, 3> orientationDegrees;
    Vec<float, 3> positionStd;
    Vec<float, 3> orientationStd;
    Vec<float, 3> dx;
    Vec<float, 3> daDegrees;
    Vec<float, 3> dxStd;
    Vec<float, 3> daStd;
};

struct ImGrid {
    Vec<int, 2> numCells;
    boolean skipOuter;
};

struct Sensor {
  int id;
  const char* name;
  const char* type;
  int robot;
  Vec<float, 3> position;
  Vec<float, 3> orientationDegrees;
  Vec<float, 3> positionStd;
  Vec<float, 3> orientationStd;
  Vec<int, 2> imageSize;
  float pixErrorStd;
  Vec<int, 4> intrinsic;
  Vec<float, 3> distortion[];
  bool frameInMap;
  ImGrid imGrid; 
};

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    
}