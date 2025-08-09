#pragma once

#define GLM_ENABLE_EXPERIMENTAL

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtc/constants.hpp> // For glm::pi
#include <vector>
#include <random>
#include <iostream>

// Constants for Poisson Disk Sampling
// 建议在主程序中根据网格大小动态计算此值，以保证均匀性
const float MIN_DISTANCE_SQUARED = 2000.0f;
const int K = 30;

// Function prototypes
bool isValid(const glm::vec2& p, const std::vector<glm::vec2>& points, float min_dist_sq, int grid_width, int grid_height);
std::vector<glm::vec2> generatePoissonPoints(int grid_width, int grid_height, int num_points);

// Check if a new point 'p' is valid (i.e., not too close to existing points)
bool isValid(const glm::vec2& p, const std::vector<glm::vec2>& points, float min_dist_sq, int grid_width, int grid_height) {
    if (p.x < 0 || p.x >= grid_width || p.y < 0 || p.y >= grid_height) {
        return false;
    }

    for (const auto& existing_point : points) {
        if (glm::distance2(p, existing_point) < min_dist_sq) {
            return false;
        }
    }
    return true;
}

std::vector<glm::vec2> generatePoissonPoints(int grid_width, int grid_height, int num_points) {
    std::vector<glm::vec2> points;
    std::vector<glm::vec2> active_list;
    points.reserve(num_points);
    active_list.reserve(num_points);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<float> width_dist(0.0f, (float)grid_width);
    std::uniform_real_distribution<float> height_dist(0.0f, (float)grid_height);

    // 1. Pick a random point and add it to the lists
    glm::vec2 first_point(width_dist(gen), height_dist(gen));
    points.push_back(first_point);
    active_list.push_back(first_point);

    int max_attempts = 1000; // 设置一个最大尝试次数，防止无限循环
    while (!active_list.empty() && points.size() < num_points) {
        if (max_attempts-- <= 0) {
            std::cerr << "Error: Failed to generate all " << num_points << " Poisson points." << std::endl;
            return {};
        }

        std::uniform_int_distribution<size_t> active_dist(0, active_list.size() - 1);
        int randomIndex = active_dist(gen);
        glm::vec2 current_point = active_list[randomIndex];

        bool found_new = false;
        for (int i = 0; i < K; ++i) {
            std::uniform_real_distribution<float> angle_dist(0.0f, 2.0f * glm::pi<float>());
            std::uniform_real_distribution<float> radius_dist(sqrt(MIN_DISTANCE_SQUARED), 2.0f * sqrt(MIN_DISTANCE_SQUARED));

            float angle = angle_dist(gen);
            float radius = radius_dist(gen);

            glm::vec2 new_point;
            new_point.x = current_point.x + radius * glm::cos(angle);
            new_point.y = current_point.y + radius * glm::sin(angle);

            if (isValid(new_point, points, MIN_DISTANCE_SQUARED, grid_width, grid_height)) {
                points.push_back(new_point);
                active_list.push_back(new_point);
                found_new = true;
                break;
            }
        }

        if (!found_new) {
            active_list.erase(active_list.begin() + randomIndex);
        }
    }
    return points;
}