const LightCubeVertexShader = `
attribute vec3 aVertexPosition;

uniform mat4 uModelViewMatrix;
uniform mat4 uProjectionMatrix;


void main(void) {

  gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(aVertexPosition, 1.0);

}
`;

const LightCubeFragmentShader = `
#ifdef GL_ES
precision mediump float;
#endif

uniform float uLigIntensity;
uniform vec3 uLightColor;

void main(void) {
    
  //gl_FragColor = vec4(1,1,1, 1.0);
  gl_FragColor = vec4(uLightColor, 1.0);
}
`;
const VertexShader = `
attribute vec3 aVertexPosition;
attribute vec3 aNormalPosition;
attribute vec2 aTextureCoord;

uniform mat4 uModelViewMatrix;
uniform mat4 uProjectionMatrix;

varying highp vec3 vFragPos;
varying highp vec3 vNormal;
varying highp vec2 vTextureCoord;

void main(void) {

  vFragPos = aVertexPosition;
  vNormal = aNormalPosition;

  gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(aVertexPosition, 1.0);

  vTextureCoord = aTextureCoord;

}
`;

const FragmentShader = `
#ifdef GL_ES
precision mediump float;
#endif

uniform int uTextureSample;
uniform vec3 uKd;
uniform sampler2D uSampler;
uniform vec3 uLightPos;
uniform vec3 uCameraPos;

varying highp vec3 vFragPos;
varying highp vec3 vNormal;
varying highp vec2 vTextureCoord;

void main(void) {
  
  if (uTextureSample == 1) {
    gl_FragColor = texture2D(uSampler, vTextureCoord);
  } else {
    gl_FragColor = vec4(uKd,1);
  }
  gl_FragColor = vec4(vNormal, 1);
}
`;

const PhongVertexShader = `
attribute vec3 aVertexPosition;
attribute vec3 aNormalPosition;
attribute vec2 aTextureCoord;

uniform mat4 uModelViewMatrix;
uniform mat4 uProjectionMatrix;

varying highp vec3 vFragPos;
varying highp vec3 vNormal;
varying highp vec2 vTextureCoord;

void main(void) {
  vec4 position = vec4(aVertexPosition, 1.0);
  vec4 normal = vec4(aNormalPosition, 1.0);
  gl_Position = uProjectionMatrix * uModelViewMatrix * position;
  // normal = uProjectionMatrix * uModelViewMatrix * normal;
  vFragPos = (uModelViewMatrix * position).xyz;
  vNormal= normal.xyz;
  vTextureCoord = aTextureCoord;
}
`;

const PhongFragmentShader = `
#ifdef GL_ES
precision mediump float;
#endif
uniform sampler2D uSampler;
uniform vec3 uKd;
uniform vec3 uKs;
uniform vec3 uLightPos;
uniform vec3 uCameraPos;
uniform float uLightIntensity;
uniform int uTextureSample;

varying highp vec3 vFragPos;
varying highp vec3 vNormal;
varying highp vec2 vTextureCoord;

void main(void) {
  vec3 baseColor;
  if (uTextureSample == 1) {
    baseColor = texture2D(uSampler, vTextureCoord).rgb;
    // Gamma Correction
    baseColor = pow(baseColor, vec3(2.2));
  } else {
    baseColor = vec3(1.0);
  }

  vec3 ambient = 0.05 * baseColor;
  vec3 lightDir = normalize(uLightPos - vFragPos);
  vec3 normal = normalize(vNormal);
  vec3 viewDir = normalize(uCameraPos - vFragPos);
  float attenuation = uLightIntensity / length(uLightPos - vFragPos);
  float shinness = 32.0;

  vec3 diffuse = attenuation * baseColor * max(dot(lightDir, normal), 0.0);
  vec3 reflectDir = reflect(-lightDir , normal);
  vec3 specular = uKs * attenuation * baseColor * pow(max(dot(viewDir, reflectDir), 0.0), shinness);
  
  // gl_FragColor = vec4(vTextureCoord.x, vTextureCoord.y, 0.0, 1.0);
  // gl_FragColor = vec4(specular, 1.0);
  // gl_FragColor = vec4(texture2D(uSampler, vec2(0.1, 0.1)).rgb, 1.0);
  gl_FragColor = vec4(pow((ambient + diffuse + specular), vec3(1.0/2.2)), 1.0);
}
`;