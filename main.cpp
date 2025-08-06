#include "World.h"
#include "WorldGenerator.h"
#include <iostream>

// File: main.cpp
int main() {
    // 1. Define all our world generation parameters in one place.
    //    This makes tweaking and experimenting super easy!
    const int worldWidth = 80;
    const int worldHeight = 40;

    // General Settings
    const int seed = 13370;
    const float waterLevel = 0.7f; // 70% of the world is ocean

    // Noise Settings
    const float scale = 0.05f;   // Lower value = larger continents. Let's use this as frequency directly.
    const int octaves = 5;       // Number of noise layers to combine.
    const float lacunarity = 2.0f; // How much detail increases with each octave.
    const float persistence = 0.5f;// How much each octave contributes to the overall shape.
    const float radius = 4.0f;

    // 2. Create the world and generator objects.
    World myWorld(worldWidth, worldHeight);
    WorldGenerator generator;

    // 3. Generate the world using the updated function call.
    std::cout << "Generating world with new parameters..." << std::endl;
    generator.Generate(myWorld, seed, scale, octaves, lacunarity, persistence, radius, waterLevel);
    std::cout << "Generation complete." << std::endl;

    // 4. Print the result.
    myWorld.PrintToConsole();

    return 0;
}