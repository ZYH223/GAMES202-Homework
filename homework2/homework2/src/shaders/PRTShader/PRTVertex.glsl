attribute vec3 aVertexPosition;
attribute vec3 aNormalPosition;
attribute vec2 aTextureCoord;

attribute mat3 aPrecomputeLT;

uniform mat4 uModelMatrix;
uniform mat4 uViewMatrix;
uniform mat4 uProjectionMatrix;

varying highp vec2 vTextureCoord;
varying highp vec3 vFragPos;
varying highp vec3 vNormal;

varying highp mat3 vPrecomputeLT;
// varying highp vec3 vColor;

void main(void) {

  vFragPos = (uModelMatrix * vec4(aVertexPosition, 1.0)).xyz;
  vNormal = (uModelMatrix * vec4(aNormalPosition, 0.0)).xyz;
  vTextureCoord = aTextureCoord;

  gl_Position = uProjectionMatrix * uViewMatrix * uModelMatrix *
                vec4(aVertexPosition, 1.0);
  vPrecomputeLT = aPrecomputeLT;
  // vColor = vec3(dotMat(uPrecomputeLr, aPrecomputeLT),
  // dotMat(uPrecomputeLg, aPrecomputeLT),
  // dotMat(uPrecomputeLb, aPrecomputeLT));
}