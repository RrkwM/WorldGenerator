#pragma once

// Forward declaration to avoid circular include.
class World;

// File: WorldGenerator.h
class WorldGenerator {
public:
    // Update the function signature to accept all noise parameters.
    void Generate(World& world,
        int seed,
        float scale,
        int octaves,
        float lacunarity,
        float persistence, // We'll call it persistence, it maps to "gain" in FastNoiseLite
        float radius, // Resize the noise map
        float waterLevel);
};