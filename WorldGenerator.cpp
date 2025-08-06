#include "WorldGenerator.h"
#include "World.h"
#include "FastNoiseLite.h"
#include <cmath>

#ifndef PI
#define PI 3.1415926535f
#endif

// File: WorldGenerator.cpp
void WorldGenerator::Generate(World& world, int seed, float scale, int octaves, float lacunarity, float persistence, float radius, float waterLevel) {
    FastNoiseLite noise(seed);
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise.SetFrequency(scale);

    noise.SetFractalType(FastNoiseLite::FractalType_FBm);
    noise.SetFractalOctaves(octaves);
    noise.SetFractalLacunarity(lacunarity);
    noise.SetFractalGain(persistence);

    // FIX 1: Introduce a radius to scale our coordinates, let's start with 4.0.
    // This "zooms in" on the 3D noise space.

    for (int y = 0; y < world.GetHeight(); ++y) {
        for (int x = 0; x < world.GetWidth(); ++x) {
            float lon = ((float)x / world.GetWidth()) * 2.0f * PI;
            float lat = (((float)y / world.GetHeight()) * PI) - (PI / 2.0f);

            // FIX 2: Correct the math for py.
            float px = std::cos(lon) * std::cos(lat);
            float py = std::sin(lon) * std::cos(lat); // Corrected from sin(lat) to cos(lat)
            float pz = std::sin(lat);

            // Apply the radius to our unit sphere coordinates before getting the noise.
            float noiseValue = noise.GetNoise(px * radius, py * radius, pz * radius);

            float normalizedValue = (noiseValue + 1.0f) / 2.0f;

            if (normalizedValue > waterLevel) {
                world.SetTile(x, y, ETileType::Land);
            }
            else {
                world.SetTile(x, y, ETileType::Ocean);
            }
        }
    }
}