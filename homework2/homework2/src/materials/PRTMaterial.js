class PRTMaterial extends Material {

    constructor(precomputeLr, precomputeLg, precomputeLb, vertexShader, fragmentShader) {
        super({
            'uPrecomputeLr': { type:'matrix3fv', value: precomputeLr },
            'uPrecomputeLg': { type:'matrix3fv', value: precomputeLg },
            'uPrecomputeLb': { type:'matrix3fv', value: precomputeLb },
        }, [
            'aPrecomputeLT'
        ], vertexShader, fragmentShader, null);
    }
}

async function buildPRTMaterial(precomputeLr, precomputeLg, precomputeLb, vertexPath, fragmentPath) {


    let vertexShader = await getShaderString(vertexPath);
    let fragmentShader = await getShaderString(fragmentPath);
    return new PRTMaterial(precomputeLr, precomputeLg, precomputeLb, vertexShader, fragmentShader);

}