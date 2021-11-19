class DirectionalLight {

    constructor(lightIntensity, lightColor, lightPos, focalPoint, lightUp, hasShadowMap, gl) {
        this.mesh = Mesh.cube(setTransform(0, 0, 0, 0.2, 0.2, 0.2, 0));
        this.mat = new EmissiveMaterial(lightIntensity, lightColor);
        this.lightPos = lightPos;
        this.focalPoint = focalPoint;
        this.lightUp = lightUp

        this.hasShadowMap = hasShadowMap;
        this.fbo = new FBO(gl);
        if (!this.fbo) {
            console.log("无法设置帧缓冲区对象");
            return;
        }
    }

    CalcLightMVP(translate, scale) {
        let lightMVP = mat4.create();
        let modelMatrix = mat4.create();
        let viewMatrix = mat4.create();
        let projectionMatrix = mat4.create();
        const lightPos = this.lightPos;
        const focalPoint = this.focalPoint;
        const lightUp = this.lightUp
        // Model transform
        mat4.translate(modelMatrix, modelMatrix, translate);
        mat4.scale(modelMatrix, modelMatrix, scale);
        // View transform
        // console.log(this.lightPos);
        // console.log(this.focalPoint);
        // console.log(this.lightUp);
        mat4.lookAt(viewMatrix, lightPos, focalPoint, lightUp);
        // Projection transform
        const aspect = 1;
        const size = 500.0;
        mat4.ortho(projectionMatrix,
            -aspect * size/2,
            aspect * size/2,
            -size/2,
            size/2,
            1.0, 300.0);
        
        mat4.multiply(lightMVP, projectionMatrix, viewMatrix);
        mat4.multiply(lightMVP, lightMVP, modelMatrix);

        return lightMVP;
    }
}
