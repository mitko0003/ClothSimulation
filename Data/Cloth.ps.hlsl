__import Shading;

struct ClothPSIn
{
    float4 pos : SV_POSITION;
    linear float3 posW : POS_W;
    linear float3 normal : NORMAL;
    nointerpolation float3 color : COLOR;
};

cbuffer PerFrameCB : register(b0)
{
    LightData gDirLight;
};

float3 AdjustNormalsToTranslucency(float3 light_dir, float3 normal)
{
    return normalize(normal + light_dir * 1.0f);
}

float4 main(ClothPSIn pIn, bool bFrontFace : SV_IsFrontFace) : SV_TARGET
{
    ShadingData sd = initShadingData();

    sd.posW = pIn.posW;
    sd.V = normalize(gCamera.posW - pIn.posW);
    sd.N = AdjustNormalsToTranslucency(-gDirLight.dirW, bFrontFace ? pIn.normal : -pIn.normal);
    sd.B = normalize(cross(sd.N, float3(1.0f, 0.0, 0.0f)));
    sd.T = normalize(cross(sd.B, sd.N));
    sd.uv = float2(0.0f, 0.0f);
    sd.NdotV = dot(sd.N, sd.V);

    sd.diffuse = pIn.color;
    sd.opacity = 1.0f;
    sd.specular = float3(0.0f, 0.0f, 0.0f);
    sd.linearRoughness = 1.0f;
    sd.roughness = 1.0f;
    sd.emissive = float3(0.0f, 0.0f, 0.0f);
    sd.IoR = 0.0f;
    sd.doubleSidedMaterial = true;

    float4 finalColor;
    finalColor.a = 1.0f;
    finalColor.rgb = evalMaterial(sd, gDirLight, 1).color.rgb;

    return finalColor;
}
