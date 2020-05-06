// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#include "TriangleIntersection.h"

static bool intersects(glm::vec3 U[3], glm::vec3 V[3], glm::vec3 segment[2], glm::vec3 &response)
{
    // Compute the plane normal for triangle U.
    glm::vec3 edge1 = U[1] - U[0];
    glm::vec3 edge2 = U[2] - U[0];
    glm::vec3 normal = glm::normalize(glm::cross(edge1, edge2));

    // Test whether the edges of triangle V transversely intersect the
    // plane of triangle U.
    glm::float32_t d[3];
    glm::int32_t positive = 0, negative = 0, zero = 0;
    for (glm::int32_t i = 0; i < 3; ++i)
    {
        d[i] = glm::dot(normal, V[i] - U[0]);
        if (d[i] > 0.0f)
        {
            ++positive;
        }
        else if (d[i] < 0.0f)
        {
            ++negative;
        }
        else
        {
            ++zero;
        }
    }
    // positive + negative + zero == 3

    if (positive > 0 && negative > 0)
    {
        if (positive == 2)  // and negative == 1
        {
            if (d[0] < 0.0f)
            {
                segment[0] = (d[1] * V[0] - d[0] * V[1]) / (d[1] - d[0]);
                segment[1] = (d[2] * V[0] - d[0] * V[2]) / (d[2] - d[0]);
                response = normal * d[0];
            }
            else if (d[1] < 0.0f)
            {
                segment[0] = (d[0] * V[1] - d[1] * V[0]) / (d[0] - d[1]);
                segment[1] = (d[2] * V[1] - d[1] * V[2]) / (d[2] - d[1]);
                response = normal * d[1];
            }
            else  // d[2] < 0.0f
            {
                segment[0] = (d[0] * V[2] - d[2] * V[0]) / (d[0] - d[2]);
                segment[1] = (d[1] * V[2] - d[2] * V[1]) / (d[1] - d[2]);
                response = normal * d[2];
            }
        }
        else if (negative == 2)  // and positive == 1
        {
            if (d[0] > 0.0f)
            {
                segment[0] = (d[1] * V[0] - d[0] * V[1]) / (d[1] - d[0]);
                segment[1] = (d[2] * V[0] - d[0] * V[2]) / (d[2] - d[0]);
                response = normal * d[0];
            }
            else if (d[1] > 0.0f)
            {
                segment[0] = (d[0] * V[1] - d[1] * V[0]) / (d[0] - d[1]);
                segment[1] = (d[2] * V[1] - d[1] * V[2]) / (d[2] - d[1]);
                response = normal * d[1];
            }
            else  // d[2] > 0.0f
            {
                segment[0] = (d[0] * V[2] - d[2] * V[0]) / (d[0] - d[2]);
                segment[1] = (d[1] * V[2] - d[2] * V[1]) / (d[1] - d[2]);
                response = normal * d[2];
            }
        }
        else  // positive == 1, negative == 1, zero == 1
        {
            if (d[0] == 0.0f)
            {
                segment[0] = V[0];
                segment[1] = (d[2] * V[1] - d[1] * V[2]) / (d[2] - d[1]);
                if (d[1] < 0.0f)
                    response = normal * (-d[1] < d[2] ? d[1] : d[2]);
                else
                    response = normal * (d[1] < -d[2] ? d[1] : d[2]);

            }
            else if (d[1] == 0.0f)
            {
                segment[0] = V[1];
                segment[1] = (d[0] * V[2] - d[2] * V[0]) / (d[0] - d[2]);
                if (d[0] < 0.0f)
                    response = normal * (-d[0] < d[2] ? d[0] : d[2]);
                else
                    response = normal * (d[0] < -d[2] ? d[0] : d[2]);
            }
            else  // d[2] == 0.0f
            {
                segment[0] = V[2];
                segment[1] = (d[1] * V[0] - d[0] * V[1]) / (d[1] - d[0]);
                if (d[0] < 0.0f)
                    response = normal * (-d[0] < d[1] ? d[0] : d[1]);
                else
                    response = normal * (d[0] < -d[1] ? d[0] : d[1]);
            }
        }
        return true;
    }

    // Triangle V does not transversely intersect triangle U, although it is
    // possible a vertex or edge of V is just touching U.  In this case, we
    // do not call this an intersection.
    return false;
}

// TODO(doc): Eberly, David H - GPGPU Programming for Games and Science-CRC Press (2014).pdf
bool intersectTriangles(glm::vec3 U[3], glm::vec3 V[3], glm::vec3 &response)
{
    glm::vec3 S0[2], S1[2];
    glm::vec3 R0, R1;
    if (intersects(V, U, S0, R0) && intersects(U, V, S1, R1))
    {
        // Theoretically, the segments lie on the same line.  A direction D
        // of the line is the Cross(NormalOf(U),NormalOf(V)).  We choose the
        // average A of the segment endpoints as the line origin.
        glm::vec3 uNormal = glm::cross(U[1] - U[0], U[2] - U[0]);
        glm::vec3 vNormal = glm::cross(V[1] - V[0], V[2] - V[0]);
        glm::vec3 D = glm::normalize(glm::cross(uNormal, vNormal));
        glm::vec3 A = 0.25f*(S0[0] + S0[1] + S1[0] + S1[1]);

        // Each segment endpoint is of the form A + t*D.  Compute the
        // t-values to obtain I0 = [t0min,t0max] for S0 and I1 = [t1min,t1max]
        // for S1.  The segments intersect when I0 overlaps I1.  Although this
        // application acts as a "test intersection" query, in fact the
        // construction here is a "find intersection" query.
        auto t00 = glm::dot(D, S0[0] - A), t01 = glm::dot(D, S0[1] - A);
        auto t10 = glm::dot(D, S1[0] - A), t11 = glm::dot(D, S1[1] - A);
        auto I0 = std::minmax(t00, t01);
        auto I1 = std::minmax(t10, t11);

        response = glm::dot(R0, R0) < glm::dot(R1, R1) ? R0 : R1;
        return (I0.second > I1.first && I0.first < I1.second);
    }
    return false;
}

// Möller–Trumbore ray-triangle intersection algorithm
// TODO(doc): https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection
bool intersectTriangle(glm::vec3 vertices[3], glm::vec3 rayOrigin, glm::vec3 rayDirection, glm::vec3 &intersection)
{
    const auto edge1 = vertices[1] - vertices[0];
    const auto edge2 = vertices[2] - vertices[0];
    const auto DxE2 = cross(rayDirection, edge2);
    const auto det = dot(edge1, DxE2);

    // ray and triangle are parallel if det is close to 0
    if (fabs(det) < std::numeric_limits<float>::epsilon())
        return false;

    const auto T = rayOrigin - vertices[0];
    const auto u = dot(T, DxE2) / det;
    if (u < 0.0f || u > 1.0f)
        return false;

    const auto TxE1 = cross(T, edge1);
    const auto v = dot(rayDirection, TxE1) / det;
    if (v < 0.0f || u + v > 1.0f)
        return false;

    const auto t = dot(edge2, TxE1) / det;

    if (t >= 0.0f)
    {
        intersection = rayOrigin + t * rayDirection;
        return true;
    }
    else
    {
        return false;
    }
}

bool intersectTriangle(glm::vec3 V[3], glm::vec3 rayOrigin, glm::vec3 rayDirection)
{
    glm::vec3 dummyIntersection;
    return intersectTriangle(V, rayOrigin, rayDirection, dummyIntersection);
}
