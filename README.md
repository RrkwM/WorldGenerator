# WorldGenerator TODO List: Procedural Planetary Genesis Engine

This list outlines the implementation plan for the procedural world generation project, from initial setup to a running tectonic simulation.

## Phase 0: Project Setup & Core Infrastructure

- [ ] **Initialize Project Environment**
    - [ ] Create a new C++ solution in Visual Studio.
    - [ ] Initialize a Git repository.
    - [ ] Create a comprehensive `.gitignore` file for Visual Studio C++ projects.
    - [ ] Push the initial project structure to the GitHub repository.

- [ ] **Dependency Management**
    - [ ] Decide on a dependency management strategy (vcpkg recommended, or a manual `vendor` folder).
    - [ ] Integrate required third-party libraries:
        - [ ] **FastNoiseLite**: For all noise generation.
        - [ ] **ImGui**: For the debug UI.
        - [ ] **GLFW**: For window and input management.
        - [ ] **Glad**: For loading OpenGL functions.

- [ ] **Core Data Structures**
    - [ ] Implement `Tile.h`: Define the `Tile` class/struct with fields like `m_Height`, `m_OwnerPlate`, etc.
    - [ ] Implement `Plate.h`: Define the `Plate` class with all identification, physical, and kinematic properties.
    * [ ] Implement `World.h`: Define the `World` class to hold the grid of `Tile` objects.

## Phase 1: Initial Planet Generation

- [ ] **Implement `WorldGenerator` class**
    - [ ] Create the `GenerateInitialPlanet()` function signature that accepts all necessary configuration parameters.

- [ ] **Implement Pangaea Generation Logic within `GenerateInitialPlanet()`**
    - [ ] **Spherical Mapping:** Convert 2D grid coordinates `(x, y)` to 3D spherical coordinates `(px, py, pz)`.
    - [ ] **Noise Sampling:** Integrate `FastNoiseLite` to sample 3D fractal noise using the spherical coordinates. Remember to use the `radius` parameter to scale the sampling domain.
    - [ ] **Continent Gradient:** Implement the radial bias logic to ensure a single continent forms. This involves calculating the distance from a center point and creating a falloff value.
    - [ ] **Value Combination:** Combine the `normalized_noise` and `continent_bias` values (e.g., via multiplication).
    - [ ] **Height Initialization:** Implement the function to map the `final_value` to a realistic initial `m_Height` for each `Tile`, based on the `waterLevel`.

- [ ] **Implement Post-processing for Water Connectivity (Optional but Recommended)**
    - [ ] Implement the `EnforceSingleOcean()` function.
    - [ ] This involves a flood-fill (BFS/DFS) algorithm to find all separate bodies of water.
    - [ ] Identify the largest body of water and convert all smaller ones to land by adjusting their `m_Height`.

## Phase 2: Tectonic Initialization

- [ ] **Implement the `InitializeTectonics()` function**
    - [ ] This function will take the `World` object and `plateCount` as input.

- [ ] **Plate Fracturing (Voronoi Diagram)**
    - [ ] Implement logic to generate `N` random seed points on the sphere's surface.
    - [ ] For each `Tile` in the `World`, calculate which seed point is closest.
    - [ ] Assign each `Tile` to its corresponding new `Plate` object.

- [ ] **Plate Property Assignment**
    - [ ] For each `Plate`, determine its `m_Type` (Continental/Oceanic) by analyzing the land-to-water ratio of its constituent tiles.
    - [ ] Assign physical properties (`m_Density`, `m_Thickness`, `m_Age`) based on the plate's type, using random values within realistic ranges.
    - [ ] Assign kinematic properties (`m_EulerPole`, `m_AngularSpeed`) with random values.

## Phase 3: Tectonic Simulation Loop

- [ ] **Implement the `SimulateTick(deltaTime)` function**

- [ ] **Data Management**
    - [ ] Implement the **Double Buffer** pattern for the height map to ensure stable updates. All reads should come from a `read_buffer` and all writes should go to a `write_buffer`.

- [ ] **Plate Kinematics**
    - [ ] For each plate, calculate its rotation for the current `tick` based on its Euler Pole and speed.
    - [ ] Apply this rotation to all tiles belonging to the plate.

- [ ] **Boundary Analysis & Geology**
    - [ ] Implement a method to detect boundary tiles and their neighbors on other plates.
    - [ ] Calculate the relative velocity vectors at these boundary points.
    - [ ] Implement the geological rules engine:
        - [ ] **Convergent Logic:** Calculate height changes for uplift (mountains) and subduction (trenches).
        - [ ] **Divergent Logic:** Calculate height changes for rifting and new crust formation.
        - [ ] **Transform Logic:** (Optional) Add minor crust deformation.
    - [ ] Ensure all geological changes are written to the `write_buffer`.

- [ ] **Erosion Simulation**
    - [ ] Implement a simple hydraulic erosion pass that moves height from higher tiles to lower neighbors.
    * [ ] Implement a simple thermal erosion pass that simulates material sliding down steep slopes.
    - [ ] Ensure erosion changes are also written to the `write_buffer`.

- [ ] **End of Tick**
    - [ ] Swap the `read_buffer` and `write_buffer`.

## Phase 4: UI & Visualization Harness

- [ ] **Implement `Application` class**
    - [ ] Set up the main application loop.
    - [ ] Handle window creation (GLFW) and destruction.
    - [ ] Initialize and shut down OpenGL and ImGui contexts.

- [ ] **Implement `UIManager` class**
    - [ ] Create ImGui widgets (sliders, buttons, etc.) for every single generation and simulation parameter.
    - [ ] Implement the logic for the "Generate" button to call `GenerateInitialPlanet()` and `InitializeTectonics()`.
    - [ ] Implement a texture buffer (`std::vector<unsigned char>`) for map visualization.
    - [ ] Create a function `UpdateTexture()` that reads the `World`'s height/tile data and converts it to colors in the texture buffer.
    - [ ] Use OpenGL to create and update a `GLuint` texture from the buffer.
    - [ ] Use `ImGui::Image()` to display the texture in a window.

## Future Enhancements (Idea Parking Lot)

- [ ] **Refactor to a common `IWorldDataSource` interface** to support both procedural generation and manual painting/loading from a file.
- [ ] Research and implement more advanced erosion models (e.g., flow-based hydraulic erosion).
- [ ] Performance optimizations for the simulation loop (e.g., focusing calculations only on plate boundaries).
- [ ] Add climate simulation based on height, latitude, and ocean currents.
