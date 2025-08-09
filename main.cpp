#define GLM_ENABLE_EXPERIMENTAL
#include "GLM/glm.hpp"
#include "GLM/vec3.hpp"
#include "GLM/gtc/matrix_transform.hpp"
#include "GLM/gtx/norm.hpp"
#include "PoissonDiskSampling.hpp"
#include "FastNoiseLite/FastNoiseLite.h"
#include "RayLib/include/raylib.h"

#include <vector>
#include <queue>
#include <iostream>
#include <limits>
#include <random>
#include <map>
#include <set>
#include <algorithm>
#include <string> // Needed for std::to_string
#include <omp.h>  // Include OpenMP for multi-threading

// =================================================================
// 全局定义 (Global Definitions)
// =================================================================
const int MAP_TILE_WIDTH = 1800;
const int MAP_TILE_HEIGHT = 1800;
const int PARTICLE_GRID_SIZE = 300;
const int PARTICLE_SPACING_KM = 6;
const bool VISUALIZE_FAULTS = true;

// Staged simulation control
enum SimulationState {
    TECTONICS,          // Stage 1: Crustal thickening and stress accumulation
    FRACTURING,         // Stage 2: One-shot fault generation
    ISOSTATIC_SETTLING, // Stage 3: Let the thickened crust rise to equilibrium
    DONE                // Simulation finished
};
// MODIFIED: Increased the duration of the tectonics phase to allow for more mature mountain formation.
const int TECTONICS_STEPS = 3000; // Was 1500
const int ISOSTATIC_SETTLING_STEPS = 500;

struct Particle {
    glm::vec3 position;
    glm::vec3 velocity;
    float mass; // Represents density now
    int plateID;
    float cohesion = 1.0f;
    glm::vec2 stress = glm::vec2(0.0f);
    float crustThickness = 0.0f;
};

struct NoiseOctave {
    FastNoiseLite noise;
    float frequency;
    float amplitude;
};

class TerrainGenerator {
private:
    int m_screenWidth = 900;
    int m_screenHeight = 900;

    std::vector<glm::vec2> m_plateSeeds;
    std::vector<std::vector<Particle>> m_particleGrid;
    std::vector<std::vector<float>> m_faultZoneMap;
    std::vector<NoiseOctave> m_mantleNoiseX;
    std::vector<NoiseOctave> m_mantleNoiseY;
    std::mt19937 m_rng;

    SimulationState m_currentState = TECTONICS;
    int m_stepCount = 0;

    Color lerp(Color a, Color b, float t) {
        t = std::max(0.0f, std::min(1.0f, t));
        Color result;
        result.r = static_cast<unsigned char>(a.r + (b.r - a.r) * t);
        result.g = static_cast<unsigned char>(a.g + (b.g - a.g) * t);
        result.b = static_cast<unsigned char>(a.b + (b.b - a.b) * t);
        result.a = 255;
        return result;
    }

    void drawTerrain() {
        const auto& particles = m_particleGrid;
        float particleWidth = (float)m_screenWidth / PARTICLE_GRID_SIZE;
        float particleHeight = (float)m_screenHeight / PARTICLE_GRID_SIZE;

        float minElevation = std::numeric_limits<float>::max();
        float maxElevation = std::numeric_limits<float>::lowest();
        for (const auto& row : particles) {
            for (const auto& p : row) {
                if (p.position.z < minElevation) minElevation = p.position.z;
                if (p.position.z > maxElevation) maxElevation = p.position.z;
            }
        }
        if (maxElevation <= minElevation) maxElevation = minElevation + 1.0f;

        const float seaLevelRatio = 0.15f;
        const float seaLevel = minElevation + (maxElevation - minElevation) * seaLevelRatio;

        const Color DEEP_OCEAN = { 10, 40, 80, 255 };
        const Color SHALLOW_OCEAN = { 80, 150, 200, 255 };
        const Color BEACH = { 210, 195, 150, 255 };
        const Color LOWLAND_GREEN = { 70, 110, 40, 255 };
        const Color UPLAND_ARID = { 160, 140, 90, 255 };
        const Color ROCKY_MOUNTAIN = { 130, 120, 110, 255 };
        const Color SNOW_CAP = { 255, 255, 255, 255 };

        for (int i = 0; i < PARTICLE_GRID_SIZE; ++i) {
            for (int j = 0; j < PARTICLE_GRID_SIZE; ++j) {
                const auto& p = particles[i][j];
                float elevation = p.position.z;
                Color color;

                if (elevation < seaLevel) {
                    float waterDepthRatio = 1.0f;
                    if (seaLevel > minElevation) {
                        waterDepthRatio = (elevation - minElevation) / (seaLevel - minElevation);
                    }
                    color = lerp(DEEP_OCEAN, SHALLOW_OCEAN, waterDepthRatio);
                }
                else {
                    float landHeightRatio = 0.0f;
                    if (maxElevation > seaLevel) {
                        landHeightRatio = (elevation - seaLevel) / (maxElevation - seaLevel);
                    }
                    if (landHeightRatio < 0.02f) {
                        color = BEACH;
                    }
                    else if (landHeightRatio < 0.25f) {
                        color = lerp(BEACH, LOWLAND_GREEN, (landHeightRatio - 0.02f) / 0.23f);
                    }
                    else if (landHeightRatio < 0.5f) {
                        color = lerp(LOWLAND_GREEN, UPLAND_ARID, (landHeightRatio - 0.25f) / 0.25f);
                    }
                    else if (landHeightRatio < 0.75f) {
                        color = lerp(UPLAND_ARID, ROCKY_MOUNTAIN, (landHeightRatio - 0.75f) / 0.25f);
                    }
                    else {
                        color = lerp(ROCKY_MOUNTAIN, SNOW_CAP, (landHeightRatio - 0.75f) / 0.25f);
                    }
                }

                DrawRectangle(
                    static_cast<int>((float)i * particleWidth), static_cast<int>((float)j * particleHeight),
                    static_cast<int>(particleWidth + 1.0f), static_cast<int>(particleHeight + 1.0f), color
                );

                if (VISUALIZE_FAULTS && m_faultZoneMap[i][j] > 0.0f && m_currentState != TECTONICS) {
                    float strength = m_faultZoneMap[i][j];
                    Color faultColor = { 255, (unsigned char)(255 * (1.0f - strength)), 0, (unsigned char)(strength * 150) };
                    DrawRectangle(
                        static_cast<int>((float)i * particleWidth), static_cast<int>((float)j * particleHeight),
                        static_cast<int>(particleWidth + 1.0f), static_cast<int>(particleHeight + 1.0f), faultColor
                    );
                }
            }
        }

        std::string statusText;
        int progress = 0;
        switch (m_currentState) {
        case TECTONICS:
            progress = (int)(((float)m_stepCount / TECTONICS_STEPS) * 100);
            statusText = "Stage 1/3: Tectonics (" + std::to_string(progress) + "%)";
            break;
        case FRACTURING:
            statusText = "Stage 2/3: Generating Faults...";
            break;
        case ISOSTATIC_SETTLING:
            progress = (int)(((float)(m_stepCount - TECTONICS_STEPS) / ISOSTATIC_SETTLING_STEPS) * 100);
            statusText = "Stage 3/3: Isostatic Settling (" + std::to_string(progress) + "%)";
            break;
        case DONE:
            statusText = "Generation Complete.";
            break;
        }
        DrawText(statusText.c_str(), 10, 10, 20, LIME);
    }

public:
    TerrainGenerator() : m_particleGrid(PARTICLE_GRID_SIZE, std::vector<Particle>(PARTICLE_GRID_SIZE)), m_rng(std::random_device{}()) {}

    void initialize() {
        initParticleGrid();
        generateVoronoiPlates();
        finalizeInitialization();
        initMantleNoise();
    }

    void runRenderer() {
        InitWindow(m_screenWidth, m_screenHeight, "Emergent Tectonics Simulation");
        SetTargetFPS(60);
        float dt = 0.02f;

        omp_set_num_threads(12);

        const int SIM_STEPS_PER_FRAME = 10;

        while (!WindowShouldClose()) {
            if (m_currentState != DONE) {
                for (int i = 0; i < SIM_STEPS_PER_FRAME; ++i) {
                    if (m_currentState == TECTONICS) {
                        simulateStep_Tectonics(dt);
                        if (m_stepCount >= TECTONICS_STEPS) {
                            m_currentState = FRACTURING;
                            break;
                        }
                    }
                    else if (m_currentState == FRACTURING) {
                        generateFaultsFromStressMap();
                        m_currentState = ISOSTATIC_SETTLING;
                        break;
                    }
                    else if (m_currentState == ISOSTATIC_SETTLING) {
                        simulateStep_Isostasy(dt);
                        if (m_stepCount >= TECTONICS_STEPS + ISOSTATIC_SETTLING_STEPS) {
                            m_currentState = DONE;
                            break;
                        }
                    }
                    m_stepCount++;
                }
            }

            BeginDrawing();
            ClearBackground(RAYWHITE);
            drawTerrain();
            EndDrawing();
        }
        CloseWindow();
    }

private:
    float getIsostaticEquilibriumZ(const Particle& p) {
        const float MANTLE_DENSITY = 3300.0f;
        float crust_density = p.mass;
        float equilibrium_z = p.crustThickness * (1.0f - (crust_density / MANTLE_DENSITY));
        return equilibrium_z * 100.0f;
    }

    void initParticleGrid() {
        m_faultZoneMap.assign(PARTICLE_GRID_SIZE, std::vector<float>(PARTICLE_GRID_SIZE, 0.0f));
        for (int i = 0; i < PARTICLE_GRID_SIZE; ++i) {
            for (int j = 0; j < PARTICLE_GRID_SIZE; ++j) {
                Particle& p = m_particleGrid[i][j];
                p.position = glm::vec3(static_cast<float>(i * PARTICLE_SPACING_KM), static_cast<float>(j * PARTICLE_SPACING_KM), 0.0f);
                p.velocity = glm::vec3(0.0f, 0.0f, 0.0f);
                p.plateID = -1;
                p.cohesion = 1.0f;
                p.stress = glm::vec2(0.0f);
                p.crustThickness = 0.0f;
            }
        }
    }

    void generateVoronoiPlates() {
        const int NUM_PLATE_SEEDS = 8;
        m_plateSeeds = generatePoissonPoints(PARTICLE_GRID_SIZE, PARTICLE_GRID_SIZE, NUM_PLATE_SEEDS);
        if (m_plateSeeds.empty()) { std::cerr << "Error: Poisson disk sampling failed." << std::endl; return; }

        FastNoiseLite domainWarp;
        domainWarp.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
        domainWarp.SetSeed(1337);
        domainWarp.SetFrequency(0.008f);
        domainWarp.SetFractalType(FastNoiseLite::FractalType_Ridged);
        domainWarp.SetFractalOctaves(4);
        domainWarp.SetFractalLacunarity(2.0f);
        domainWarp.SetFractalGain(0.5f);
        domainWarp.SetDomainWarpType(FastNoiseLite::DomainWarpType_OpenSimplex2);
        domainWarp.SetDomainWarpAmp(200.0f);

#pragma omp parallel for
        for (int i = 0; i < PARTICLE_GRID_SIZE; ++i) {
            for (int j = 0; j < PARTICLE_GRID_SIZE; ++j) {
                float warped_x = static_cast<float>(i);
                float warped_y = static_cast<float>(j);
                domainWarp.DomainWarp(warped_x, warped_y);
                glm::vec2 warped_pos(warped_x, warped_y);
                float min_dist_sq = std::numeric_limits<float>::max();
                int closest_seed_id = -1;
                for (size_t k = 0; k < m_plateSeeds.size(); ++k) {
                    float dist_sq = glm::distance2(warped_pos, m_plateSeeds[k]);
                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                        closest_seed_id = static_cast<int>(k);
                    }
                }
                m_particleGrid[i][j].plateID = closest_seed_id;
            }
        }
    }

    void finalizeInitialization() {
        const float MIN_DENSITY = 2700.0f;
        const float MAX_DENSITY = 3000.0f;
        const float MIN_THICKNESS = 7.0f;
        const float MAX_THICKNESS = 35.0f;

        FastNoiseLite property_noise;
        property_noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
        property_noise.SetFrequency(0.01f);
        property_noise.SetSeed(2026);

#pragma omp parallel for
        for (int i = 0; i < PARTICLE_GRID_SIZE; ++i) {
            for (int j = 0; j < PARTICLE_GRID_SIZE; ++j) {
                Particle& p = m_particleGrid[i][j];
                float noise_val = property_noise.GetNoise((float)i, (float)j);
                float t = (noise_val + 1.0f) * 0.5f;
                p.mass = MIN_DENSITY + (MAX_DENSITY - MIN_DENSITY) * t;
                p.crustThickness = MAX_THICKNESS - (MAX_THICKNESS - MIN_THICKNESS) * t;
                p.position.z = getIsostaticEquilibriumZ(p);
            }
        }
    }

    void generateFaultsFromStressMap() {
        std::vector<std::pair<float, glm::ivec2>> stressed_particles;
        const float MAX_BRITTLENESS = 1.5f;
        const float MIN_DENSITY = 2700.0f;
        const float MAX_DENSITY = 3000.0f;

        for (int i = 0; i < PARTICLE_GRID_SIZE; ++i) {
            for (int j = 0; j < PARTICLE_GRID_SIZE; ++j) {
                Particle& p = m_particleGrid[i][j];
                float stress_mag = glm::length(p.stress);
                float density_t = (p.mass - MIN_DENSITY) / (MAX_DENSITY - MIN_DENSITY);
                float brittleness = 1.0f + (1.0f - density_t) * (MAX_BRITTLENESS - 1.0f);
                float effective_stress = stress_mag * brittleness;
                if (effective_stress > 1.0f) {
                    stressed_particles.push_back({ effective_stress, {i, j} });
                }
            }
        }

        std::sort(stressed_particles.begin(), stressed_particles.end(),
            [](const auto& a, const auto& b) {
                return a.first > b.first;
            });

        const int num_main_faults = 15;
        const int main_fault_length = 150;
        const float main_fault_cohesion = 0.2f;
        for (int i = 0; i < num_main_faults && i < stressed_particles.size(); ++i) {
            glm::ivec2 start_pos = stressed_particles[i].second;
            glm::vec2 stress_dir = m_particleGrid[start_pos.x][start_pos.y].stress;
            createFaultPath(start_pos, stress_dir, main_fault_length, main_fault_cohesion);
        }

        const int num_secondary_faults = 40;
        const int secondary_fault_length = 60;
        const float secondary_fault_cohesion = 0.4f;
        int start_index = num_main_faults;
        for (int i = start_index; i < start_index + num_secondary_faults && i < stressed_particles.size(); ++i) {
            glm::ivec2 start_pos = stressed_particles[i].second;
            glm::vec2 stress_dir = m_particleGrid[start_pos.x][start_pos.y].stress;
            createFaultPath(start_pos, stress_dir, secondary_fault_length, secondary_fault_cohesion);
        }

        createFaultZoneMap();
    }

    void createFaultPath(glm::ivec2 start_pos, glm::vec2 stress_dir, int max_length, float cohesion) {
        if (glm::length(stress_dir) < 1e-6) return;
        stress_dir = glm::normalize(stress_dir);
        glm::vec2 perp_dir = glm::vec2(-stress_dir.y, stress_dir.x);
        std::uniform_real_distribution<float> dist_float(0.0f, 1.0f);
        glm::vec2 current_dir;
        if (dist_float(m_rng) < 0.7f) {
            current_dir = (dist_float(m_rng) < 0.5f) ? perp_dir : -perp_dir;
        }
        else {
            current_dir = (dist_float(m_rng) < 0.5f) ? stress_dir : -stress_dir;
        }
        glm::ivec2 current_pos = start_pos;
        for (int i = 0; i < max_length; ++i) {
            if (current_pos.x < 0 || current_pos.x >= PARTICLE_GRID_SIZE || current_pos.y < 0 || current_pos.y >= PARTICLE_GRID_SIZE) break;
            m_particleGrid[current_pos.x][current_pos.y].cohesion = std::min(m_particleGrid[current_pos.x][current_pos.y].cohesion, cohesion);
            glm::vec2 random_dir = glm::normalize(glm::vec2(dist_float(m_rng) - 0.5f, dist_float(m_rng) - 0.5f));
            current_dir = glm::normalize(current_dir * 0.8f + random_dir * 0.2f);
            current_pos.x += static_cast<int>(round(current_dir.x));
            current_pos.y += static_cast<int>(round(current_dir.y));
        }
    }

    void createFaultZoneMap() {
        const int influence_radius = 20;
        std::vector<glm::ivec2> fault_points;
        for (int i = 0; i < PARTICLE_GRID_SIZE; ++i) {
            for (int j = 0; j < PARTICLE_GRID_SIZE; ++j) {
                if (m_particleGrid[i][j].cohesion < 1.0f) {
                    fault_points.push_back({ i, j });
                }
            }
        }

        if (fault_points.empty()) return;

#pragma omp parallel for
        for (int i = 0; i < PARTICLE_GRID_SIZE; ++i) {
            for (int j = 0; j < PARTICLE_GRID_SIZE; ++j) {
                float min_dist_sq = std::numeric_limits<float>::max();
                for (const auto& fault_pos : fault_points) {
                    float dist_sq = glm::distance2(glm::vec2(i, j), glm::vec2(fault_pos));
                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                    }
                }
                if (min_dist_sq < influence_radius * influence_radius) {
                    float dist = sqrt(min_dist_sq);
                    m_faultZoneMap[i][j] = cos(dist / influence_radius * (3.14159f / 2.0f));
                }
            }
        }
    }

    void initMantleNoise() {
        m_mantleNoiseX.push_back({ FastNoiseLite(), 0.002f, 1.0f });
        m_mantleNoiseY.push_back({ FastNoiseLite(), 0.002f, 1.0f });
        m_mantleNoiseX.push_back({ FastNoiseLite(), 0.008f, 0.5f });
        m_mantleNoiseY.push_back({ FastNoiseLite(), 0.008f, 0.5f });
        m_mantleNoiseX.push_back({ FastNoiseLite(), 0.02f, 0.25f });
        m_mantleNoiseY.push_back({ FastNoiseLite(), 0.02f, 0.25f });
        int seed = 2025;
        for (auto& octave : m_mantleNoiseX) {
            octave.noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
            octave.noise.SetSeed(seed++);
            octave.noise.SetFrequency(octave.frequency);
        }
        for (auto& octave : m_mantleNoiseY) {
            octave.noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
            octave.noise.SetSeed(seed++);
            octave.noise.SetFrequency(octave.frequency);
        }
    }

    void applyUniversalForces(std::vector<std::vector<glm::vec3>>& forces, bool useIsostasy) {
        const float DAMPING_COEFFICIENT = 0.5f;
        const float MANTLE_DRAG_SCALE = 7000.0f;
        const float ISOSTATIC_REBOUND_SCALE = 0.05f;

#pragma omp parallel for
        for (int i = 0; i < PARTICLE_GRID_SIZE; ++i) {
            for (int j = 0; j < PARTICLE_GRID_SIZE; ++j) {
                Particle& p = m_particleGrid[i][j];
                if (p.plateID == -1) continue;

                float total_mantle_vel_x = 0.0f, total_mantle_vel_y = 0.0f;
                for (const auto& octave : m_mantleNoiseX) { total_mantle_vel_x += octave.noise.GetNoise(p.position.x, p.position.y) * octave.amplitude; }
                for (const auto& octave : m_mantleNoiseY) { total_mantle_vel_y += octave.noise.GetNoise(p.position.x, p.position.y) * octave.amplitude; }
                forces[i][j].x += total_mantle_vel_x * MANTLE_DRAG_SCALE;
                forces[i][j].y += total_mantle_vel_y * MANTLE_DRAG_SCALE;

                if (useIsostasy) {
                    float target_z = getIsostaticEquilibriumZ(p);
                    forces[i][j].z += (target_z - p.position.z) * ISOSTATIC_REBOUND_SCALE;
                }

                forces[i][j] -= p.velocity * DAMPING_COEFFICIENT;
            }
        }
    }

    void updatePositions(float dt, const std::vector<std::vector<glm::vec3>>& forces) {
#pragma omp parallel for
        for (int i = 0; i < PARTICLE_GRID_SIZE; ++i) {
            for (int j = 0; j < PARTICLE_GRID_SIZE; ++j) {
                Particle& p = m_particleGrid[i][j];
                if (p.plateID == -1) continue;
                glm::vec3 acceleration = forces[i][j] / (p.mass / 1000.0f);
                p.velocity += acceleration * dt;
                float max_vel = 50.0f;
                if (glm::length(p.velocity) > max_vel) { p.velocity = glm::normalize(p.velocity) * max_vel; }
                p.position += p.velocity * dt;
                if (p.position.x < 0) p.position.x += MAP_TILE_WIDTH;
                if (p.position.x >= MAP_TILE_WIDTH) p.position.x -= MAP_TILE_WIDTH;
                if (p.position.y < 0) p.position.y += MAP_TILE_HEIGHT;
                if (p.position.y >= MAP_TILE_HEIGHT) p.position.y -= MAP_TILE_HEIGHT;
            }
        }
    }

    void simulateStep_Tectonics(float dt) {
        std::vector<std::vector<glm::vec3>> forces(PARTICLE_GRID_SIZE, std::vector<glm::vec3>(PARTICLE_GRID_SIZE, glm::vec3(0.0f)));
        applyUniversalForces(forces, true);

        const float SPRING_CONSTANT = 90.0f;
        const float REST_LENGTH = PARTICLE_SPACING_KM;
        const float HORIZONTAL_PRESSURE_SCALE = 90.0f;
        const float THICKENING_RATE = 0.05f;
        const float DENSITY_DIFF_THRESHOLD = 1.05f;
        const float THICKNESS_RESISTANCE_SCALE = 120.0f;


#pragma omp parallel for
        for (int i = 0; i < PARTICLE_GRID_SIZE; ++i) {
            for (int j = 0; j < PARTICLE_GRID_SIZE; ++j) {
                Particle& p = m_particleGrid[i][j];
                if (p.plateID == -1) continue;
                glm::ivec2 directions[] = { {0, 1}, {1, 0}, {1, 1}, {-1, 1} };
                for (auto& dir : directions) {
                    int nx = i + dir.x, ny = j + dir.y;
                    if (nx >= 0 && nx < PARTICLE_GRID_SIZE && ny >= 0 && ny < PARTICLE_GRID_SIZE) {
                        Particle& neighbor = m_particleGrid[nx][ny];
                        if (neighbor.plateID == -1) continue;
                        glm::vec3 distance_vector = p.position - neighbor.position;
                        float dist = glm::length(distance_vector);
                        if (dist < 1e-6) continue;
                        glm::vec3 dir_vec = distance_vector / dist;
                        float displacement = dist - REST_LENGTH;
                        glm::vec3 spring_force = -dir_vec * displacement * SPRING_CONSTANT;

#pragma omp atomic
                        forces[i][j].x += spring_force.x;
#pragma omp atomic
                        forces[i][j].y += spring_force.y;
#pragma omp atomic
                        forces[nx][ny].x -= spring_force.x;
#pragma omp atomic
                        forces[nx][ny].y -= spring_force.y;

                        if (displacement < 0) {
                            float compression_depth = -displacement;
                            glm::vec3 pressure_force = -dir_vec * compression_depth * HORIZONTAL_PRESSURE_SCALE;

#pragma omp atomic
                            forces[i][j].x += pressure_force.x;
#pragma omp atomic
                            forces[i][j].y += pressure_force.y;
#pragma omp atomic
                            forces[nx][ny].x -= pressure_force.x;
#pragma omp atomic
                            forces[nx][ny].y -= pressure_force.y;

                            float base_thickness_to_add = compression_depth * THICKENING_RATE;

                            if (p.mass > neighbor.mass * DENSITY_DIFF_THRESHOLD) {
                                float resistance = neighbor.crustThickness / THICKNESS_RESISTANCE_SCALE;
                                float thickness_to_add = base_thickness_to_add * std::max(0.0f, 1.0f - resistance);
#pragma omp atomic
                                neighbor.crustThickness += thickness_to_add;
                            }
                            else if (neighbor.mass > p.mass * DENSITY_DIFF_THRESHOLD) {
                                float resistance = p.crustThickness / THICKNESS_RESISTANCE_SCALE;
                                float thickness_to_add = base_thickness_to_add * std::max(0.0f, 1.0f - resistance);
#pragma omp atomic
                                p.crustThickness += thickness_to_add;
                            }
                            else {
                                float p_resistance = p.crustThickness / THICKNESS_RESISTANCE_SCALE;
                                float n_resistance = neighbor.crustThickness / THICKNESS_RESISTANCE_SCALE;
                                float p_thickness_to_add = base_thickness_to_add * 0.5f * std::max(0.0f, 1.0f - p_resistance);
                                float n_thickness_to_add = base_thickness_to_add * 0.5f * std::max(0.0f, 1.0f - n_resistance);
#pragma omp atomic
                                p.crustThickness += p_thickness_to_add;
#pragma omp atomic
                                neighbor.crustThickness += n_thickness_to_add;
                            }
                        }
                    }
                }
            }
        }
        updatePositions(dt, forces);
    }

    void simulateStep_Isostasy(float dt) {
        std::vector<std::vector<glm::vec3>> forces(PARTICLE_GRID_SIZE, std::vector<glm::vec3>(PARTICLE_GRID_SIZE, glm::vec3(0.0f)));

        const float DAMPING_COEFFICIENT = 2.0f;
        const float ISOSTATIC_REBOUND_SCALE = 0.2f;

#pragma omp parallel for
        for (int i = 0; i < PARTICLE_GRID_SIZE; ++i) {
            for (int j = 0; j < PARTICLE_GRID_SIZE; ++j) {
                Particle& p = m_particleGrid[i][j];
                if (p.plateID == -1) continue;

                float target_z = getIsostaticEquilibriumZ(p);

                if (m_faultZoneMap[i][j] > 0.0f) {
                    float fault_influence = m_faultZoneMap[i][j];
                    target_z -= fault_influence * 20.0f * p.cohesion;
                }

                forces[i][j].z += (target_z - p.position.z) * ISOSTATIC_REBOUND_SCALE;
                forces[i][j] -= p.velocity * DAMPING_COEFFICIENT;
            }
        }
        updatePositions(dt, forces);
    }
};

int main() {
    TerrainGenerator generator;
    generator.initialize();
    generator.runRenderer();
    return 0;
}
