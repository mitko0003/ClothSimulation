__import ShaderCommon;

struct ClothVSIn
{
    float3 position : POSITION;
    float3 normal : NORMAL;
    float3 color : COLOR;
    float2 uv : TEXCOORD;
};

struct ClothVSOut
{
    float4 pos : SV_POSITION;
    linear float3 posW : POS_W;
    linear float3 normal : NORMAL;
    linear float2 uv : TEXCOORD;
    nointerpolation float3 color : COLOR;
};

//cbuffer DataCB : register(b0)
//{
//    float4x4 gWorldMat;
//}

ClothVSOut main(ClothVSIn vIn)
{
    ClothVSOut vOut;
    vOut.pos.xyz = vIn.position;
    vOut.pos.w = 1.0f;

    vOut.posW = vIn.position;
    vOut.pos = mul(vOut.pos, gCamera.viewProjMat);

    vOut.normal = vIn.normal;
    vOut.color = vIn.color;
    vOut.uv = vIn.uv;

    return vOut;
}
