__import Shading;

Texture2D<float4> gBaseColorMap;
Texture2D<float> gRoughnessMap;
Texture2D<float3> gNormalMap;
SamplerState gSampler;

struct ClothPSIn
{
    float4 pos : SV_POSITION;
    linear float3 posW : POS_W;
    linear float3 normal : NORMAL;
    linear float2 uv : TEXCOORD;
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

    float3 normalMap = normalize(gNormalMap.Sample(gSampler, pIn.uv));

    sd.N = bFrontFace ? pIn.normal : -pIn.normal;
    normalMap = bFrontFace ? normalMap : -normalMap;

    sd.N = normalize(float3(normalMap.xy + sd.N.xy, normalMap.z * sd.N.z));
    sd.N = AdjustNormalsToTranslucency(-gDirLight.dirW, sd.N);

    sd.B = normalize(cross(sd.N, float3(1.0f, 0.0, 0.0f)));
    sd.T = normalize(cross(sd.B, sd.N));
    sd.uv = float2(0.0f, 0.0f);
    sd.NdotV = dot(sd.N, sd.V);

    float4 baseColorMap = gBaseColorMap.Sample(gSampler, pIn.uv);
    sd.diffuse = lerp(pIn.color, baseColorMap.xyz, baseColorMap.w);

    sd.opacity = 1.0f;
    sd.specular = float3(0.0f, 0.0f, 0.0f);
    sd.linearRoughness = gRoughnessMap.Sample(gSampler, pIn.uv);
    sd.roughness = sd.linearRoughness;
    sd.emissive = float3(0.0f, 0.0f, 0.0f);
    sd.IoR = 0.0f;
    sd.doubleSidedMaterial = true;

    float4 finalColor;
    finalColor.a = 1.0f;
    finalColor.rgb = evalMaterial(sd, gDirLight, 1).color.rgb;

    return finalColor;
}
