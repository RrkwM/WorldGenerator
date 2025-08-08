# WorldGenerator TODO List

## 1. Vendor Utility (The "Harness")

This section covers all the scaffolding, third-party libraries, and UI tools needed to run, visualize, and debug the core simulation.

- [ ] **Project & Environment Setup**
    - [ ] Initialize a Visual Studio C++ Solution and Git repository.
    - [ ] Create and configure the `.gitignore` file.
    - [ ] Set up a dependency management solution (vcpkg or a manual `vendor` folder).
        - [ ] Integrate **FastNoiseLite**.
        - [ ] Integrate **ImGui**.
        - [ ] Integrate **GLFW**.
        - [ ] Integrate **Glad**.

- [ ] **Application Layer (`Application` class)**
    - [ ] Implement the main application class.
    - [ ] Handle window creation and the main program loop.
    - [ ] Initialize and manage the lifecycles of OpenGL and ImGui.

- [ ] **UI & Visualization Layer (`UIManager` class)**
    - [ ] Implement the UI manager class.
    - [ ] Create ImGui widgets (sliders, buttons) to control all simulation parameters.
    - [ ] Link the "Generate" and "Tick" buttons to the Main Logic functions.
    - [ ] Implement the visualization pipeline:
        - [ ] Create a texture buffer (`std::vector<unsigned char>`).
        - [ ] Write a function to update this buffer with colors based on the `World`'s data (e.g., height or tile type).
        - [ ] Manage an OpenGL texture (`GLuint`) for displaying the map.
        - [ ] Render the texture within an ImGui window using `ImGui::Image()`.

---

## 2. Main Logic (The "Engine")

This section contains the pure, portable C++ code that constitutes the core planetary simulation. It has no knowledge of ImGui or any other vendor utility.

- [ ] **Core Data Structures**
    - [ ] Define `Tile.h` (with `m_Height`, `m_OwnerPlate`, etc.).
    - [ ] Define `Plate.h` (with all identification, physical, and kinematic properties).
    - [ ] Define `World.h` (to hold the grid of tiles and manage the double-buffered height map).

- [ ] **Initial State Generation (`GenerateInitialPlanet` function)**
    - [ ] Implement spherical mapping (2D grid to 3D coordinates).
    - [ ] Implement the Pangaea-shaping logic by combining 3D fractal noise with a radial gradient.
    - [ ] Implement the function to set the initial `m_Height` of every tile.

- [ ] **Tectonic Initialization (`InitializeTectonics` function)**
    - [ ] Implement the plate fracturing algorithm using a Spherical Voronoi Diagram.
    - [ ] Implement the property assignment logic for each new plate (Type, Density, Thickness, Age, Euler Pole, Speed).

- [ ] **Tectonic Simulation (`SimulateTick` function)**
    - [ ] Implement the double-buffering pattern for all height map updates.
    - [ ] **Plate Kinematics:** Move plates by applying a rotation based on their Euler Pole and speed.
    - [ ] **Boundary Analysis:** Detect boundary interactions and calculate relative velocity.
    - [ ] **Geological Rules Engine:**
        - [ ] Implement height changes for **convergent** boundaries (uplift, subduction).
        - [ ] Implement height changes for **divergent** boundaries (rifting).
    - [ ] **Erosion Engine:**
        - [ ] Implement the hydraulic erosion pass.
        - [ ] Implement the thermal erosion pass.

---

## 3. Enhancement (The "Future")

This section is the "idea parking lot" for future features and architectural improvements once the main logic is complete and stable.

- [ ] **Architectural Refinements**
    - [ ] Design and implement a common `IWorldDataSource` interface to abstract the world creation process (procedural vs. loaded from file).

- [ ] **Advanced Simulation Features**
    - [ ] Research and implement more complex, flow-based hydraulic erosion models.
    - [ ] Implement a basic climate simulation (e.g., temperature based on latitude and altitude, precipitation based on wind and mountains).

- [ ] **Performance & Art-Directability**
    - [ ] Profile the simulation and optimize bottlenecks (e.g., focus calculations only on active plate boundaries).
    - [ ] Add tools for artistic control, such as painting "weight maps" to influence geological activity (e.g., uplift or erosion strength).
