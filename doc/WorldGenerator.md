1. World is based on Tiles.
2. World Generator uses perlin noise. It ensures:
    - Land : Water ~= 3:7
    - The generated 2D noise map is based on sphere to simulate actual earth plate
    - Water body is always continuous.
3. ImGui is used for visualization.