#pragma once

#include "Tile.h"
#include <vector>
#include <iostream>

// File: World.h
// Holds the grid of tiles that represents our world.
class World {
public:
    World(int width, int height)
        : m_Width(width), m_Height(height), m_Grid(height, std::vector<ETileType>(width)) {
    }

    ETileType GetTile(int x, int y) const {
        if (x >= 0 && x < m_Width && y >= 0 && y < m_Height) {
            return m_Grid[y][x];
        }
        return ETileType::Ocean; // Return ocean for out-of-bounds access.
    }

    void SetTile(int x, int y, ETileType type) {
        if (x >= 0 && x < m_Width && y >= 0 && y < m_Height) {
            m_Grid[y][x] = type;
        }
    }

    int GetWidth() const { return m_Width; }
    int GetHeight() const { return m_Height; }

    // A simple function to print the world to the console for debugging.
    void PrintToConsole() const {
        for (int y = 0; y < m_Height; ++y) {
            for (int x = 0; x < m_Width; ++x) {
                // Use '#' for land and '.' for ocean.
                std::cout << (GetTile(x, y) == ETileType::Land ? '#' : '.');
            }
            std::cout << std::endl;
        }
    }

private:
    int m_Width;
    int m_Height;
    std::vector<std::vector<ETileType>> m_Grid;
};