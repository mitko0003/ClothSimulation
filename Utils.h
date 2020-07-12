#pragma once

#include "glm/glm.hpp"

class NonCopyable
{
public:
	NonCopyable() = default;
	NonCopyable(const NonCopyable&) = delete;
	NonCopyable& operator=(const NonCopyable&) = delete;
};

vec3 getYawPitchRoll(vec3 direction)
{
	direction = normalize(direction);
	const auto phi = atan2(direction[0], direction[2]);
	const auto theta = asin(-direction[1]);
	return vec3(phi, theta + radians(90.0f), 0.0f);
}