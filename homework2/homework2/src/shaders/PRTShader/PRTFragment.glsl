#ifdef GL_ES
precision mediump float;
#endif

uniform mat3 uPrecomputeLr;
uniform mat3 uPrecomputeLg;
uniform mat3 uPrecomputeLb;

varying mat3 vPrecomputeLT;
// varying highp vec3 vColor;

float dotMat(mat3 m1, mat3 m2)
{
  float sum = 0.0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      sum += max(m1[i][j] * m2[i][j], 0.0);
    }
  }
  return sum;
}

void main(void) {  
  vec3 PRTColor = vec3(dotMat(uPrecomputeLr, vPrecomputeLT),
  dotMat(uPrecomputeLg, vPrecomputeLT),
  dotMat(uPrecomputeLb, vPrecomputeLT));
  gl_FragColor = vec4(PRTColor, 1.0);
}