// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2020
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2019.08.13

#pragma once

#include <algorithm>
#include "glm/glm.hpp"

bool intersectTriangles(const glm::vec3 U[3], glm::vec3 V[3], glm::vec3 &response);
bool intersectTriangle(const glm::vec3 V[3], glm::vec3 rayOrigin, glm::vec3 rayDirection);
bool intersectTriangle(const glm::vec3 V[3], glm::vec3 rayOrigin, glm::vec3 rayDirection, glm::vec3 &intersection);
glm::vec3 triangleClosestPoint(glm::vec3 V[3], const glm::vec3 &point);
