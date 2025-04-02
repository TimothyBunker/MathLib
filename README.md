# C Framebuffer Graphics Engine

## Project Description

This project is a low-level graphics rendering engine implemented in C. It provides a foundation for 2D and 3D graphics by combining a custom-built linear algebra library with direct framebuffer manipulation on Linux systems.  The engine bypasses higher-level graphics APIs, allowing for a direct understanding of the graphics pipeline from mathematical transformations to pixel-level rendering.

This project serves as a learning exercise and a demonstration of core graphics programming principles. It's designed to be lightweight, with minimal external dependencies, making it potentially suitable for educational purposes or resource-constrained environments.

## Features

*   **Custom Linear Algebra Library:**
    *   Implementation of vector (2D, 3D, 4D) and matrix (3x3, 4x4) types.
    *   Essential vector operations: addition, subtraction, dot product, cross product, magnitude, normalization.
    *   Key matrix operations: multiplication, identity matrix generation, transpose. *(Inverse matrix calculation may be included depending on implementation).*
    *   Functions to create transformation matrices (translation, rotation, scaling, and potentially projection matrices for 3D).
*   **Direct Framebuffer Access:**
    *   Opens and memory-maps the Linux framebuffer device (`/dev/fb*`).
    *   Queries framebuffer properties like resolution, color depth, and line stride.
    *   `draw_pixel()` function to directly write color values to the framebuffer memory, adapting to different pixel formats.
*   **Rendering Primitives:**
    *   Efficient line drawing using Bresenham's algorithm.

## Getting Started

### Prerequisites

*   **Linux Operating System:**  Requires a Linux distribution with framebuffer support enabled in the kernel.
*   **C Compiler:**  GCC or Clang (and standard build tools).
*   **Make:**  Build automation tool (usually comes with development tools).
*   **Development Libraries (Minimal):** Standard C library is the primary dependency. No external graphics libraries are needed.

### Building

1.  **Clone the Repository:**
    ```bash
    git clone <YOUR_GITHUB_REPO_LINK_HERE> # <<< Replace with your repo link
    cd <your-repository-directory>
    ```

2.  **Compile the Code:**
    Use the provided `Makefile` to compile the project.
    ```bash
    make
    ```
    This should create an executable file (e.g., `framebuffer_demo` or similar, check your Makefile).

### Running

1.  **Permissions:** Ensure you have permissions to access the framebuffer device (usually `/dev/fb0`). You might need to be part of the `video` group or use `sudo` for testing (be cautious with `sudo` and direct hardware access).

## Future Work

*   Implement solid triangle rasterization with Z-buffering for 3D rendering.
*   Add lighting models (ambient, diffuse, specular shading).
*   Explore texture mapping.
*   Implement more advanced clipping algorithms.
*   Optimize performance for real-time rendering.
*   Integrate basic input handling for interactivity.
*   Refine the API and add documentation.
*   Potentially explore shader-like functionality in software.

## License

This project is open-source and available under the [Apache 2.0](LICENSE). See the `LICENSE` file for details.

## Author

[Timothy Bunker] - [bunkertimothy0@gmail.com]

---